"""Bounded execution and semantic parity for the ONT performance backend."""

import contextlib
import argparse
import gzip
import importlib.util
import io
from pathlib import Path
import sqlite3
import subprocess
import sys
import tempfile
import threading
import time
import unittest
from unittest.mock import Mock, patch

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("ont_performance", ROOT / "scripts/ont_features.py")
MOD = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = MOD
SPEC.loader.exec_module(MOD)


class BatchTests(unittest.TestCase):
    def test_feature_shards_are_read_only_chromosome_specific_and_bounded(self):
        with tempfile.TemporaryDirectory() as temp:
            paths = {}
            for chrom in ["chr1", "chr2"]:
                paths[chrom] = Path(temp) / f"{chrom}.sqlite"
                with MOD.open_scratch_database(paths[chrom]) as db:
                    db.execute("INSERT INTO features VALUES (?,?,?,?,?)", (chrom, "same", chrom, "T", "candidate"))
                db.close()
            with MOD.ChromosomeFeatureLookup(paths) as lookup:
                self.assertEqual(lookup.lookup("chr1", ["same", "missing"]), {"same": ("chr1", "T", "candidate")})
                first = lookup.db
                self.assertEqual(lookup.lookup("chr1", ["same"])["same"][0], "chr1")
                self.assertIs(lookup.db, first)
                self.assertEqual(lookup.lookup("chr2", ["same"])["same"][0], "chr2")
                with self.assertRaises(sqlite3.ProgrammingError):
                    first.execute("SELECT 1")
                with self.assertRaises(sqlite3.OperationalError):
                    lookup.db.execute("DELETE FROM features")
                self.assertEqual(lookup.lookup("absent", ["same"]), {})
            with self.assertRaises(sqlite3.ProgrammingError):
                lookup.db.execute("SELECT 1")

    def test_genome_mapq_database_is_separate_from_retained_features(self):
        with tempfile.TemporaryDirectory() as temp:
            features, genome = Path(temp) / "features.sqlite", Path(temp) / "genome.sqlite"
            db = MOD.open_scratch_database(features, genome_path=genome)
            budget = abs(db.execute("PRAGMA main.cache_size").fetchone()[0]) + abs(db.execute("PRAGMA maps.cache_size").fetchone()[0])
            self.assertLessEqual(budget, 32768, "attached databases must share the original worker cache budget")
            db.execute("INSERT INTO genome VALUES ('r1', 60)")
            db.commit()
            self.assertEqual(MOD.lookup_batch(db, "genome", ["r1"]), {"r1": (60,)})
            self.assertIsNone(db.execute("SELECT name FROM main.sqlite_master WHERE name='genome'").fetchone())
            db.close()
            genome.unlink()
            with sqlite3.connect(features) as db:
                self.assertEqual(db.execute("SELECT count(*) FROM features").fetchone(), (0,))

    def test_both_workflows_expose_and_forward_performance_options(self):
        for name in ["run_all.sh", "run_all_mixed_species.sh"]:
            source = (ROOT / name).read_text()
            self.assertIn("[--ont-jobs 1]", source)
            self.assertIn("--ont-temp-dir", source)
            self.assertIn('ONT_JOBS=1', source)
            self.assertIn('--jobs "${ONT_JOBS}"', source)
            self.assertIn('ONT_ARGS+=(--temp-dir "${ONT_TEMP_DIR}")', source)
            self.assertIn('FLORA_TRIM_CACHE_DIR="${ONT_TEMP_DIR}" run_stage flora', source)
        runtime = (ROOT / "src/workflow_runtime.rs").read_text()
        self.assertIn('"--ont-temp-dir",', runtime)

    def test_job_count_is_capped_by_total_thread_budget(self):
        self.assertTrue(hasattr(MOD, "execution_layout"), "need total-budget job allocation")
        self.assertEqual(MOD.execution_layout(16, 2), (2, 8))
        self.assertEqual(MOD.execution_layout(16, 20), (2, 8))
        self.assertEqual(MOD.execution_layout(4, 2), (1, 4))
        for threads, jobs in [(0, 1), (1, 0), (-1, 1)]:
            with self.assertRaises(ValueError):
                MOD.execution_layout(threads, jobs)

    def test_lookup_is_one_join_and_handles_missing_and_other_chromosome(self):
        db = sqlite3.connect(":memory:")
        self.addCleanup(db.close)
        db.execute("CREATE TABLE features(chrom,read_id,gene,transcript,reason,PRIMARY KEY(chrom,read_id)) WITHOUT ROWID")
        db.executemany("INSERT INTO features VALUES (?,?,?,?,?)", [
            ("chr1", "same", "G1", "T1", "candidate"),
            ("chr2", "same", "G2", "T2", "candidate"),
        ])
        statements = []
        db.set_trace_callback(statements.append)
        self.assertTrue(hasattr(MOD, "lookup_batch"), "need bounded batch joins")
        result = MOD.lookup_batch(db, "features", ["same", "missing"], "chr1")
        self.assertEqual(result, {"same": ("G1", "T1", "candidate")})
        self.assertEqual(sum(s.lstrip().upper().startswith("SELECT") for s in statements), 1)

    def test_duplicate_validation_catches_cross_batch_and_same_batch(self):
        self.assertTrue(hasattr(MOD, "validate_read_ids"), "need batch duplicate validation")
        for batches in [(["a", "a"],), (["a", "b"], ["c", "a"])]:
            db = sqlite3.connect(":memory:")
            try:
                db.execute("CREATE TABLE counted_reads(read_id TEXT PRIMARY KEY) WITHOUT ROWID")
                with self.assertRaisesRegex(ValueError, "Duplicate primary read name: a"):
                    for batch in batches:
                        MOD.validate_read_ids(db, batch)
            finally:
                db.close()

    def test_disposable_database_retains_savepoint_rollback(self):
        with tempfile.TemporaryDirectory() as temp:
            db = MOD.open_scratch_database(Path(temp) / "scratch.sqlite")
            self.addCleanup(db.close)
            db.execute("CREATE TABLE counted_reads(read_id TEXT PRIMARY KEY) WITHOUT ROWID")
            MOD.validate_read_ids(db, ["z"])
            db.commit()
            with self.assertRaisesRegex(ValueError, "Duplicate primary read name: z"):
                MOD.validate_read_ids(db, ["a", "b", "z"])

    def test_parallel_map_is_ordered_and_bounded(self):
        self.assertTrue(hasattr(MOD, "ordered_parallel_map"), "need bounded chromosome scheduling")
        lock = threading.Lock()
        active = peak = 0
        def work(n):
            nonlocal active, peak
            with lock:
                active += 1
                peak = max(peak, active)
            time.sleep(0.02 if n == 0 else 0.005)
            with lock:
                active -= 1
            return n * 2
        self.assertEqual(list(MOD.ordered_parallel_map(work, range(6), 2)), [0, 2, 4, 6, 8, 10])
        self.assertEqual(peak, 2)

    def test_parallel_failure_propagates_and_joins_workers(self):
        self.assertTrue(hasattr(MOD, "ordered_parallel_map"))
        def fail(n):
            if n == 1:
                raise ValueError("worker failed")
            return n
        with self.assertRaisesRegex(ValueError, "worker failed"):
            list(MOD.ordered_parallel_map(fail, range(8), 2))


class SchedulerTests(unittest.TestCase):
    def consume_in_thread(self, iterator):
        results, errors = [], []
        finished = threading.Event()

        def consume():
            try:
                for result in iterator:
                    results.append(result)
            except BaseException as error:
                errors.append(error)
            finally:
                finished.set()

        thread = threading.Thread(target=consume, daemon=True)
        thread.start()
        return thread, results, errors, finished

    def test_replenishes_while_first_task_is_blocked_and_keeps_order(self):
        for jobs in (2, 3):
            with self.subTest(jobs=jobs):
                release_first = threading.Event()
                later_started = threading.Event()
                initial_workers = threading.Barrier(jobs)
                lock = threading.Lock()
                active = peak = 0

                def work(n):
                    nonlocal active, peak
                    with lock:
                        active += 1
                        peak = max(peak, active)
                    try:
                        if n < jobs:
                            initial_workers.wait(2)
                        if n == 0:
                            release_first.wait()
                        if n == jobs:
                            later_started.set()
                        return n * 2
                    finally:
                        with lock:
                            active -= 1

                thread, results, errors, _ = self.consume_in_thread(
                    MOD.ordered_parallel_map(work, range(12), jobs))
                try:
                    self.assertTrue(later_started.wait(2), "idle worker must start a replacement before task 0 finishes")
                    self.assertEqual(results, [], "exports must still wait for task 0")
                finally:
                    release_first.set()
                    thread.join(5)
                self.assertFalse(thread.is_alive())
                self.assertEqual(errors, [])
                self.assertEqual(results, list(range(0, 24, 2)))
                self.assertEqual(peak, jobs)
                self.assertEqual(active, 0)

    def test_bounds_outstanding_plus_buffered_results_and_input_consumption(self):
        for jobs in (2, 3):
            with self.subTest(jobs=jobs):
                release_first = threading.Event()
                window_started = threading.Event()
                consumed = []
                results = []
                limit = 2 * jobs

                def items():
                    for n in range(30):
                        self.assertLess(len(consumed) - len(results), limit,
                                        "running, queued and buffered tasks must share one bound")
                        consumed.append(n)
                        yield n

                def work(n):
                    if n == 0:
                        release_first.wait()
                    if n == limit - 1:
                        window_started.set()
                    return n

                iterator = MOD.ordered_parallel_map(work, items(), jobs)
                self.assertEqual(consumed, [], "constructing the iterator must be lazy")
                thread, results, errors, _ = self.consume_in_thread(iterator)
                try:
                    self.assertTrue(window_started.wait(2), "completion-driven work should fill the bounded window")
                    self.assertEqual(consumed, list(range(limit)))
                    self.assertEqual(results, [])
                finally:
                    release_first.set()
                    thread.join(5)
                self.assertFalse(thread.is_alive())
                self.assertEqual(errors, [])
                self.assertEqual(results, list(range(30)))

    def test_lazy_serial_and_parallel_maps_accept_none_and_empty_input(self):
        for jobs in (1, 2):
            with self.subTest(jobs=jobs):
                consumed, called = [], []
                values = [0, 1, None, 3, None, 5]

                def items():
                    for item in values:
                        consumed.append(item)
                        yield item

                def work(item):
                    called.append(item)
                    return item

                iterator = MOD.ordered_parallel_map(work, items(), jobs)
                self.assertEqual(consumed, [])
                self.assertEqual(called, [])
                self.assertEqual(next(iterator), 0)
                self.assertLessEqual(len(consumed), 1 if jobs == 1 else 2 * jobs)
                self.assertEqual(list(iterator), values[1:])
                self.assertCountEqual(called, values)
                self.assertEqual(list(MOD.ordered_parallel_map(work, [], jobs)), [])
                self.assertEqual(len(called), len(values))

    def test_worker_error_propagates_after_running_workers_are_joined(self):
        release_first = threading.Event()
        first_started = threading.Event()
        failed = threading.Event()
        active = set()
        lock = threading.Lock()
        error = ValueError("worker failed")

        def work(n):
            with lock:
                active.add(n)
            try:
                if n == 0:
                    first_started.set()
                    release_first.wait()
                elif n == 1:
                    self.assertTrue(first_started.wait(2))
                    failed.set()
                    raise error
                return n
            finally:
                with lock:
                    active.remove(n)

        thread, _, errors, finished = self.consume_in_thread(
            MOD.ordered_parallel_map(work, range(20), 2))
        try:
            self.assertTrue(failed.wait(2))
            self.assertFalse(finished.is_set(), "a failed map must join the still-running worker")
        finally:
            release_first.set()
            thread.join(5)
        self.assertFalse(thread.is_alive())
        self.assertEqual(errors, [error])
        self.assertEqual(active, set())

    def test_source_error_propagates_and_joins_initial_workers(self):
        started = threading.Event()
        release = threading.Event()
        source_failed = threading.Event()
        worker_finished = threading.Event()
        error = RuntimeError("input failed")

        def items():
            yield 0
            self.assertTrue(started.wait(2))
            source_failed.set()
            raise error

        def work(n):
            started.set()
            try:
                release.wait()
                return n
            finally:
                worker_finished.set()

        thread, _, errors, finished = self.consume_in_thread(
            MOD.ordered_parallel_map(work, items(), 2))
        try:
            self.assertTrue(source_failed.wait(2))
            self.assertFalse(finished.is_set())
        finally:
            release.set()
            thread.join(5)
        self.assertFalse(thread.is_alive())
        self.assertEqual(errors, [error])
        self.assertTrue(worker_finished.is_set())

    def test_closing_map_joins_running_workers_without_consuming_more_input(self):
        started = threading.Event()
        release = threading.Event()
        worker_finished = threading.Event()
        consumed = []

        def items():
            for n in range(20):
                consumed.append(n)
                yield n

        def work(n):
            if n == 1:
                started.set()
                try:
                    release.wait()
                finally:
                    worker_finished.set()
            return n

        iterator = MOD.ordered_parallel_map(work, items(), 2)
        thread = None
        errors = []

        def close():
            try:
                iterator.close()
            except BaseException as error:
                errors.append(error)

        try:
            self.assertEqual(next(iterator), 0)
            self.assertTrue(started.wait(2))
            before_close = list(consumed)
            thread = threading.Thread(target=close, daemon=True)
            thread.start()
        finally:
            release.set()
            if thread is not None:
                thread.join(5)
            else:
                iterator.close()
        self.assertFalse(thread.is_alive())
        self.assertEqual(errors, [])
        self.assertTrue(worker_finished.is_set())
        self.assertEqual(consumed, before_close)


class AssemblyLoggingTests(unittest.TestCase):
    def test_stringtie_start_logs_context_before_execution_without_changing_command(self):
        for options in ("-c 2", "-c 2 -A 'coverage values.tsv'"):
            with self.subTest(options=options), tempfile.TemporaryDirectory() as temp:
                directory = Path(temp)
                (directory / "stringtie.gff").write_text("")
                log = io.StringIO()
                view = Mock(stdout=io.BytesIO(), args=["samtools", "view", "-h", "reads.bam", "chr7"])
                view.wait.return_value = 0

                def run(command, **kwargs):
                    message = log.getvalue()
                    self.assertIn("Stage stringtie (chr7) starting:", message)
                    self.assertIn(f"options={options!r}", message)
                    self.assertIn(f"work_dir={directory}", message)
                    self.assertNotIn("Stage stringtie (chr7) elapsed:", message)
                    self.assertIs(kwargs["stdin"], view.stdout)
                    return subprocess.CompletedProcess(command, 0)

                with patch.object(MOD.subprocess, "Popen", return_value=view) as popen, \
                        patch.object(MOD.subprocess, "run", side_effect=run) as execute, \
                        contextlib.redirect_stderr(log):
                    MOD.assemble(Path("reads.bam"), "chr7", Path("genes.gtf"), Path("genome.fa"),
                                 directory, 8, options, "flora_chr7")
                extra = [] if options == "-c 2" else ["-A", "coverage values.tsv"]
                execute.assert_called_once_with(
                    ["stringtie", "-L", "-c", "2", *extra, "-p", "8", "-G", "genes.gtf",
                     "-l", "flora_chr7", "-o", str(directory / "stringtie.gff"), "-"], stdin=view.stdout)
                popen.assert_called_once_with(view.args, stdout=subprocess.PIPE)
                self.assertTrue(view.stdout.closed)
                view.wait.assert_called_once_with()
                self.assertIn("Stage stringtie (chr7) elapsed:", log.getvalue())


class BamParityTests(unittest.TestCase):
    def test_run_isolates_worker_databases_and_uses_requested_scratch(self):
        import pysam
        import subprocess
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            scratch = (root / "scratch").resolve()
            scratch.mkdir()
            genome = root / "genome.fa"
            genome.write_text(">chr1\n" + "A" * 1000 + "\n>chr2\n" + "C" * 1000 + "\n")
            pysam.faidx(str(genome))
            gtf = root / "genes.gtf"
            gtf.write_text(''.join(f'{c}\tt\ttranscript\t1\t100\t.\t+\t.\tgene_name "G"; transcript_id "T";\n' for c in ["chr1", "chr2"]))
            bam = root / "input.bam"
            header = {"HD": {"SO": "coordinate"}, "SQ": [{"SN": c, "LN": 1000} for c in ["chr1", "chr2"]]}
            with pysam.AlignmentFile(bam, "wb", header=header) as out:
                for tid in [0, 1]:
                    r = pysam.AlignedSegment(out.header)
                    r.query_name, r.query_sequence = f"r{tid}", "A" * 50
                    r.reference_id, r.reference_start, r.mapping_quality, r.cigarstring = tid, 0, 60, "50M"
                    out.write(r)
            pysam.index(str(bam))
            worker_paths = set()
            def assemble(bam, chrom, gtf, genome, directory, threads, options, label):
                self.assertTrue(directory.is_relative_to(scratch))
                (directory / "transcriptome.fa").write_text(f">{label}\n" + "A" * 100 + "\n")
                (directory / "gffcompare.annotated.gtf").write_text(gtf.read_text())
            def remap(db, bam, chrom, fasta, output, threads):
                path = Path(db.execute("PRAGMA database_list").fetchone()[2])
                worker_paths.add(path)
                self.assertEqual(threads, 8)
                with pysam.AlignmentFile(output, "wb", header={"HD": {"SO": "queryname"}, "SQ": [{"SN": "T", "LN": 100}]}) as out:
                    r = pysam.AlignedSegment(out.header)
                    r.query_name, r.query_sequence = ("r0" if chrom == "chr1" else "r1"), "A" * 50
                    r.reference_id, r.reference_start, r.mapping_quality, r.cigarstring = 0, 0, 60, "50M"
                    r.set_tag("AS", 50)
                    out.write(r)
                    db.execute("INSERT INTO genome VALUES (?,?)", (r.query_name, 60))
                    db.commit()
            args = argparse.Namespace(bam=bam, tagged_bam=bam, genome=genome, gtf=gtf,
                output_dir=root / "out", output_bam=root / "out.bam", gene_assigns=root / "g.tsv",
                transcript_assigns=root / "t.tsv", sample="sample", threads=16, mapq=30,
                stringtie_opts="-c 2", save_intermediate=True, temp_dir=scratch, jobs=2)
            def version(command, **kwargs):
                return subprocess.CompletedProcess(command, 0, "2.2.2" if command[0] == "stringtie" else "test", "")
            self.assertTrue(hasattr(MOD, "execution_layout"), "need isolated chromosome workers")
            with patch.object(MOD.shutil, "which", return_value="tool"), patch.object(MOD.subprocess, "run", side_effect=version), \
                    patch.object(MOD, "assemble", side_effect=assemble), patch.object(MOD, "remap_chromosome", side_effect=remap), \
                    patch.object(MOD, "load_candidates", return_value={"T": ("G", "T")}), contextlib.redirect_stderr(io.StringIO()):
                MOD.run(args)
            self.assertEqual(len(worker_paths), 2)
            self.assertFalse(list(scratch.iterdir()), "scratch must be removed after success")
            self.assertEqual((root / "g.tsv").read_text(), "r0\tAssigned\t60\tG\nr1\tAssigned\t60\tG\n")
            with gzip.open(args.output_dir / "sample.transcriptome.fa.gz", "rt") as handle:
                self.assertEqual(handle.read(), "".join(f">flora_chr{i}\n" + "A" * 100 + "\n" for i in range(2)))
            with gzip.open(args.output_dir / "sample.transcriptome.gtf.gz", "rt") as handle:
                self.assertEqual(handle.read(), gtf.read_text())
            self.assertEqual(len(list((args.output_dir / "intermediates").glob("ont-work-*/chrom_*/reference.gtf"))), 2)

    def test_native_fastq_matches_original_for_reverse_missing_and_hard_clips(self):
        import pysam
        self.assertTrue(hasattr(MOD, "fastq_bytes"), "need native quality conversion")
        header = pysam.AlignmentHeader.from_dict({"SQ": [{"SN": "chr1", "LN": 1000}]})
        for reverse in (False, True):
            for missing in (False, True):
                record = pysam.AlignedSegment(header)
                record.query_name = "read1"
                record.query_sequence = "ACGTN" * 10
                record.flag = 16 if reverse else 0
                record.reference_id = 0
                record.reference_start = 0
                record.cigarstring = "2H50M3H"
                if not missing:
                    record.query_qualities = pysam.qualitystring_to_array("!5I0B" * 10)
                seq = record.get_forward_sequence()
                qual = record.get_forward_qualities() if not missing else None
                expected = "".join(chr(q + 33) for q in qual) if qual is not None else "B" * len(seq)
                self.assertEqual(MOD.fastq_bytes(record), f"@read1\n{seq}\n+\n{expected}\n".encode())

    def test_batched_writer_preserves_record_order_status_and_tags(self):
        import pysam
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            bam = root / "input.bam"
            header = {"HD": {"SO": "coordinate"}, "SQ": [{"SN": "chr1", "LN": 1000}, {"SN": "chr2", "LN": 1000}]}
            db = sqlite3.connect(":memory:")
            self.addCleanup(db.close)
            db.execute("CREATE TABLE features(chrom,read_id,gene,transcript,reason,PRIMARY KEY(chrom,read_id)) WITHOUT ROWID")
            cases = [("r1", 60, "G", "T", "candidate", "Assigned", "Assigned"),
                     ("r2", 1, "-", "T", "candidate", "Unassigned_mapq", "Assigned"),
                     ("r3", 60, "G", "-", "candidate", "Assigned", "Unassigned_ambiguous"),
                     ("r4", 60, "-", "-", "ambiguous_or_zero_transcript_mapq", "Unassigned_ambiguous", "Unassigned_ambiguous"),
                     ("r5", 60, "-", "-", "no_reference_candidate", "Unassigned_no_features", "Unassigned_no_features")]
            with pysam.AlignmentFile(bam, "wb", header=header) as out:
                for i, (rid, mapq, gene, tx, reason, _, _) in enumerate(cases):
                    rec = pysam.AlignedSegment(out.header)
                    rec.query_name, rec.query_sequence = rid, "ACGT" * 10
                    rec.reference_id, rec.reference_start, rec.mapping_quality = 0, i * 50, mapq
                    rec.cigarstring = "40M"
                    rec.set_tags([("CB", "cell1"), ("UR", "ACGT"), ("GN", "old")])
                    out.write(rec)
                    db.execute("INSERT INTO features VALUES (?,?,?,?,?)", ("chr1", rid, gene, tx, reason))
                rec.flag, rec.reference_id, rec.reference_start = 2048, 1, 0
                out.write(rec)
            trace = []
            db.set_trace_callback(trace.append)
            log = io.StringIO()
            with contextlib.redirect_stderr(log):
                MOD.write_tagged_assignments(db, bam, root / "out.bam", root / "g.tsv", root / "t.tsv", 30, 2)
            expected_g = "".join(f"{rid}\t{gs}\t{q}\t{g if g != '-' else 'NA'}\n" for rid, q, g, t, r, gs, ts in cases)
            expected_t = "".join(f"{rid}\t{ts}\t{q}\t{g if g != '-' else 'NA'}\t{t if t != '-' else 'NA'}\n" for rid, q, g, t, r, gs, ts in cases)
            self.assertEqual((root / "g.tsv").read_text(), expected_g)
            self.assertEqual((root / "t.tsv").read_text(), expected_t)
            with pysam.AlignmentFile(root / "out.bam") as result:
                rows = list(result)
            self.assertEqual([r.query_name for r in rows], [x[0] for x in cases])
            self.assertEqual([r.get_tag("GN") for r in rows], [x[2] if x[2] != "-" else "NA" for x in cases])
            selects = [s for s in trace if s.lstrip().upper().startswith("SELECT")]
            self.assertLessEqual(len(selects), 1, "must not query features once per read")
            self.assertIn("Stage tagged_assignment_output elapsed:", log.getvalue())
            self.assertIn("Stage feature_bam_index elapsed:", log.getvalue())

            shards = {"chr1": root / "chr1.sqlite"}
            shard = MOD.open_scratch_database(shards["chr1"])
            shard.executemany("INSERT INTO features VALUES (?,?,?,?,?)", db.execute("SELECT * FROM features"))
            shard.commit()
            shard.close()
            db2 = MOD.open_scratch_database(root / "global.sqlite")
            self.addCleanup(db2.close)
            with MOD.ChromosomeFeatureLookup(shards) as lookup, contextlib.redirect_stderr(io.StringIO()):
                stats = MOD.write_tagged_assignments(db2, bam, root / "sharded.bam", root / "sg.tsv", root / "st.tsv", 30, 2, feature_lookup=lookup)
            self.assertEqual((root / "sg.tsv").read_bytes(), (root / "g.tsv").read_bytes())
            self.assertEqual((root / "st.tsv").read_bytes(), (root / "t.tsv").read_bytes())
            with pysam.AlignmentFile(root / "sharded.bam") as result:
                self.assertEqual([r.to_string() for r in result], [r.to_string() for r in rows])
            self.assertEqual(sum(v for k, v in stats.items() if k.startswith("gene:")), len(cases))

    def test_sharded_writer_matches_global_index_across_batch_and_chromosome_boundaries(self):
        import itertools
        import pysam
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            header = {"HD": {"SO": "coordinate"}, "SQ": [{"SN": c, "LN": 3000000} for c in ["chr1", "chr2"]]}
            bam = root / "input.bam"
            db = MOD.open_scratch_database(root / "global.sqlite")
            self.addCleanup(db.close)
            paths = {c: root / f"{c}.sqlite" for c in ["chr1", "chr2"]}
            for c, path in paths.items():
                shard = MOD.open_scratch_database(path)
                with pysam.AlignmentFile(bam if c == "chr1" else root / "part2.bam", "wb", header=header) as out:
                    for i in range(10003):
                        r = pysam.AlignedSegment(out.header)
                        r.query_name, r.query_sequence = f"{c}_r{i}", "ACGT" * 10
                        r.reference_id, r.reference_start = int(c[-1]) - 1, i * 50
                        r.mapping_quality, r.cigarstring = (1 if i % 7 == 0 else 60), "40M"
                        out.write(r)
                        if i % 11:
                            row = (c, r.query_name, "G" if i % 7 else "-", "T" if i % 3 else "-", "candidate")
                            db.execute("INSERT INTO features VALUES (?,?,?,?,?)", row)
                            shard.execute("INSERT INTO features VALUES (?,?,?,?,?)", row)
                shard.commit()
                shard.close()
            pysam.cat("-o", str(root / "combined.bam"), str(bam), str(root / "part2.bam"))
            db.commit()
            with contextlib.redirect_stderr(io.StringIO()):
                expected = MOD.write_tagged_assignments(db, root / "combined.bam", root / "old.bam", root / "g.tsv", root / "t.tsv", 30, 2)
            empty = MOD.open_scratch_database(root / "empty.sqlite")
            self.addCleanup(empty.close)
            with MOD.ChromosomeFeatureLookup(paths) as lookup, contextlib.redirect_stderr(io.StringIO()):
                actual = MOD.write_tagged_assignments(empty, root / "combined.bam", root / "new.bam", root / "sg.tsv", root / "st.tsv", 30, 2, feature_lookup=lookup)
            self.assertEqual(actual, expected)
            self.assertEqual((root / "g.tsv").read_bytes(), (root / "sg.tsv").read_bytes())
            self.assertEqual((root / "t.tsv").read_bytes(), (root / "st.tsv").read_bytes())
            with pysam.AlignmentFile(root / "old.bam") as old, pysam.AlignmentFile(root / "new.bam") as new:
                for a, b in itertools.zip_longest(old, new):
                    self.assertIsNotNone(a)
                    self.assertIsNotNone(b)
                    self.assertEqual(a.to_string(), b.to_string())
            # A repeated primary ID on another chromosome still fails globally.
            with self.assertRaisesRegex(ValueError, "Duplicate primary read name: chr1_r0"):
                MOD.validate_read_ids(empty, ["chr1_r0"])


if __name__ == "__main__":
    unittest.main()
