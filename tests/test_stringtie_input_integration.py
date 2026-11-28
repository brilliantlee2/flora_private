"""Actual release binaries and processRead state, opt in via STRINGTIE_INPUT_BUILD."""

import hashlib
import json
import os
from pathlib import Path
import random
import shlex
import subprocess
import tempfile
import time
import unittest

ROOT = Path(__file__).resolve().parents[1]
PACKAGE = ROOT / "packaging/stringtie_input"


def records(pysam, header, adversarial=False):
    rng = random.Random(222)
    for chrom in range(2):
        for locus in (1000, 20000, 40000):
            for i in range(180):
                r = pysam.AlignedSegment(header)
                r.query_name = f"pair{chrom}_{locus}_{i % 20}" if i < 40 else f"r{chrom}_{locus}_{i}"
                r.reference_id = chrom
                r.reference_start = locus + (5 if 20 <= i < 40 else 0 if i < 140 else i % 7)
                r.mapping_quality = 60
                left = 100 + (i % 9 if i < 60 else rng.randrange(5))
                # Same exon coordinates, but different deletion-junction signatures.
                deletion = i % 4 if i % 3 == 0 else 0
                r.cigartuples = ([(0, left - deletion), (2, deletion)] if deletion else [(0, left)]) + [(3, 200), (0, 140 + i % 5)]
                if i % 11 == 0:
                    r.cigartuples = [(4, 5)] + r.cigartuples + [(4, 3)]
                r.query_sequence = "ACGT" * (sum(n for op, n in r.cigartuples if op in (0, 1, 4)) // 4 + 1)
                r.query_sequence = r.query_sequence[:sum(n for op, n in r.cigartuples if op in (0, 1, 4))]
                r.query_qualities = pysam.qualitystring_to_array("I" * len(r.query_sequence))
                r.flag = (99 if i < 20 else 147) if i < 40 else (256 if i % 17 == 0 else 2048 if i % 19 == 0 else 0)
                r.next_reference_id = chrom
                r.next_reference_start = locus + (5 if i < 20 else 0)
                r.template_length = 440 if i % 2 == 0 else -440
                r.set_tag("NH", (1, 2, 4)[i % 3])
                r.set_tag("HI", 1)
                r.set_tag("NM", i % 6)
                r.set_tag("YC", 1 + i % 3)
                if i % 3 != 2:
                    r.set_tag("XS", "+" if i % 3 == 0 else "-", "A")
                r.set_tag("YK", (-2.0 if adversarial and i % 7 == 0 else 2.5 if i % 5 == 0 else 0.0), "f")
                r.set_tag("ZG", int(i % 2 == 0))
                r.set_tag("ZU", i % 2)
                if adversarial and i < 6:
                    # Three guide-adjusted nonmatching reads share a stored
                    # signature. A zero-deletion query must update the newest.
                    r.cigartuples = ([(0, 98), (2, 2)] if i < 3 else [(0, 100)]) + [(3, 200), (0, 140 if i < 4 else 141)]
                    r.flag = 0
                    for tag, value in (("NH", 1), ("NM", 0), ("YC", 1), ("ZG", 1), ("ZU", 1)):
                        r.set_tag(tag, value)
                    r.set_tag("XS", "+", "A")
                    r.set_tag("YK", -2.0 if i == 4 else 3.0 if i == 5 else 0.0, "f")
                    r.query_sequence = "A" * sum(n for op, n in r.cigartuples if op == 0)
                    r.query_qualities = pysam.qualitystring_to_array("I" * len(r.query_sequence))
                yield r


def make_bam(pysam, path, adversarial=False):
    header = pysam.AlignmentHeader.from_dict({"HD": {"VN": "1.6", "SO": "coordinate"},
                                            "SQ": [{"SN": c, "LN": 100000} for c in ("chr2", "chr22")]})
    data = list(records(pysam, header, adversarial))
    data.sort(key=lambda r: (r.reference_id, r.reference_start))
    with pysam.AlignmentFile(path, "wb", header=header) as out:
        for record in data:
            out.write(record)
    return len(data)


@unittest.skipUnless(os.environ.get("STRINGTIE_INPUT_BUILD"), "set STRINGTIE_INPUT_BUILD to a completed optional build")
class ReleaseParityTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        import pysam
        cls.pysam = pysam
        cls.build = Path(os.environ["STRINGTIE_INPUT_BUILD"]).resolve()
        cls.stock = cls.build / "stringtie-2.2.2-stock"
        cls.cached = cls.build / "stringtie-2.2.2-same-start"
        cls.metadata = json.loads((cls.build / "build-provenance.json").read_text())
        for binary in (cls.stock, cls.cached):
            if subprocess.check_output([str(binary), "--version"], text=True).strip() != "2.2.2":
                raise AssertionError(f"unexpected version: {binary}")

    def modes(self):
        return [(self.stock, "linear"), (self.cached, "linear"), (self.cached, "indexed")]

    def test_real_process_read_state_after_every_call(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            bam = root / "state.bam"
            count = make_bam(self.pysam, bam, adversarial=True)
            harnesses = []
            dependency = next((self.build / "dependency").glob("libdeflate-*/libdeflate.a"))
            cxx = shlex.split(os.environ.get("CXX", "c++"))
            flags = self.metadata.get("flags", {})
            for name in ("stock", "accelerated"):
                source = self.build / name / "stringtie-2.2.2"
                includes = ["-I" + str(source), "-I" + str(source / "gclib"), "-I" + str(source / "htslib")]
                main = root / (name + "-main.o")
                command = cxx + ["-std=c++11", "-O3", "-pthread", "-fsigned-char"] + includes + shlex.split(flags.get("cppflags", ""))
                subprocess.run(command + ["-Dmain=stringtie_original_main", "-c", str(source / "stringtie.cpp"), "-o", str(main)], check=True, capture_output=True)
                exe = root / (name + "-process-read")
                objects = list((source / "gclib").glob("*.o")) + [source / (n + ".o") for n in ("rlink", "tablemaker", "tmerge")]
                subprocess.run(command + [str(PACKAGE / "tests/process_read.cpp"), str(main)] + [str(p) for p in objects]
                               + [str(source / "htslib/libhts.a"), str(dependency), "-lz", "-lbz2", "-llzma", "-lm"]
                               + shlex.split(flags.get("ldflags", "")) + ["-o", str(exe)], check=True, capture_output=True)
                harnesses.append(exe)
            for mode in ("short", "long", "mixed", "unitig", "nomulti"):
                expected = None
                for exe, lookup in [(harnesses[0], "linear"), (harnesses[1], "linear"), (harnesses[1], "indexed")]:
                    result = subprocess.run([str(exe), str(bam), mode], env=dict(os.environ, STRINGTIE_INPUT_LOOKUP=lookup),
                                            check=True, capture_output=True, timeout=120)
                    self.assertEqual(result.stdout.count(b"CALL:"), count)
                    self.assertIn(b" P:", result.stdout, "fixture must exercise actual pair-count mutations")
                    if expected is None:
                        expected = result.stdout
                    else:
                        self.assertEqual(result.stdout, expected, f"processRead mutation mismatch: {mode}/{lookup}")
                print(f"processRead {mode}: {count} exact state snapshots, sha256={hashlib.sha256(expected).hexdigest()}", flush=True)
                if mode == "long":
                    fourth = next(line for line in expected.splitlines() if line.startswith(b"CALL:4 "))
                    self.assertEqual(fourth.count(b" R:"), 3, "guide deletion duplicates must remain separate")
                    self.assertIn(b"0x1p+1", fourth.split(b" R:2,")[1].split(b" E:")[0], "newest match must receive count 2")

    def test_synthetic_bam_assemblies_and_quantification(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            bam = root / "reads.bam"
            count = make_bam(self.pysam, bam)
            guide = root / "guide.gtf"
            guide.write_text("".join(
                f'{chrom}\tfixture\texon\t{start}\t{end}\t.\t+\t.\tgene_id "G{chrom}_{locus}"; transcript_id "T{chrom}_{locus}";\n'
                for chrom in ("chr2", "chr22") for locus in (1000, 20000, 40000)
                for start, end in ((locus + 1, locus + 100), (locus + 301, locus + 440))))
            scenarios = [("short", []), ("long", ["-L"]), ("guided", ["-L", "-G", str(guide)]),
                         ("quantify", ["-L", "-e", "-G", str(guide)]), ("nomulti", ["-L", "-u"]),
                         ("threads", ["-L", "-p", "2"]), ("mixed", ["--mix"])]
            for scenario, options in scenarios:
                expected = None
                for index, (binary, lookup) in enumerate(self.modes()):
                    out = root / f"{scenario}-{index}"
                    out.mkdir()
                    command = [str(binary)] + options + ["-o", str(out / "assembly.gtf"), "-A", str(out / "abundance.tsv")]
                    command += [str(bam), str(bam)] if scenario == "mixed" else [str(bam)]
                    start = time.perf_counter()
                    subprocess.run(command, env=dict(os.environ, STRINGTIE_INPUT_LOOKUP=lookup), check=True, capture_output=True, timeout=120)
                    elapsed = time.perf_counter() - start
                    # Only CLI comments differ (binary/output paths); retain every data row verbatim.
                    data = {p.name: [line for line in p.read_text().splitlines() if not line.startswith("#")]
                            for p in out.iterdir() if p.is_file()}
                    self.assertTrue(any("\ttranscript\t" in line for line in data["assembly.gtf"]), scenario)
                    if expected is None:
                        expected = data
                    else:
                        self.assertEqual(data, expected, f"assembly/abundance mismatch: {scenario}/{lookup}")
                    print(f"assembly {scenario}/{index} {lookup}: reads={count}, elapsed={elapsed:.4f}s", flush=True)

    def test_dense_same_start_and_merge(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            bam = root / "dense.bam"
            header = self.pysam.AlignmentHeader.from_dict({"HD": {"SO": "coordinate"}, "SQ": [{"SN": "chr2", "LN": 100000}]})
            distinct = 12000
            with self.pysam.AlignmentFile(bam, "wb", header=header) as out:
                for i in range(distinct * 2):
                    r = self.pysam.AlignedSegment(header)
                    r.query_name = f"dense{i}"
                    r.reference_id = 0
                    r.reference_start = 1000
                    r.mapping_quality = 60
                    r.cigartuples = [(0, 200 + i % distinct)]
                    # SEQ=* is legal BAM; this isolates geometry, not sequence compression.
                    r.set_tag("NH", 1)
                    out.write(r)
            expected = None
            assemblies = []
            for i, (binary, lookup) in enumerate(self.modes()):
                gtf = root / f"dense-{i}.gtf"
                start = time.perf_counter()
                subprocess.run([str(binary), "-L", "-p", "1", "-o", str(gtf), str(bam)],
                               env=dict(os.environ, STRINGTIE_INPUT_LOOKUP=lookup), check=True, capture_output=True, timeout=120)
                elapsed = time.perf_counter() - start
                data = [line for line in gtf.read_text().splitlines() if not line.startswith("#")]
                self.assertTrue(data)
                if expected is None:
                    expected = data
                else:
                    self.assertEqual(data, expected)
                assemblies.append(gtf)
                print(f"dense same-start {i}/{lookup}: distinct={distinct}, reads={distinct*2}, elapsed={elapsed:.4f}s", flush=True)
            merged = None
            for i, (binary, lookup) in enumerate(self.modes()):
                out = root / f"merge-{i}.gtf"
                subprocess.run([str(binary), "--merge", "-o", str(out), str(assemblies[0]), str(assemblies[0])],
                               env=dict(os.environ, STRINGTIE_INPUT_LOOKUP=lookup), check=True, capture_output=True, timeout=120)
                data = [line for line in out.read_text().splitlines() if not line.startswith("#")]
                self.assertTrue(data)
                if merged is None:
                    merged = data
                else:
                    self.assertEqual(data, merged)


if __name__ == "__main__":
    unittest.main()
