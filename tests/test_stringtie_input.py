"""Offline contracts and compiled candidate-index regressions; no runtime wiring."""

import importlib.util
import contextlib
import io
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
PACKAGE = ROOT / "packaging/stringtie_input"


def builder():
    spec = importlib.util.spec_from_file_location("stringtie_input_builder", PACKAGE / "build.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class BuilderTests(unittest.TestCase):
    def test_optional_builder_exists(self):
        self.assertTrue((PACKAGE / "build.py").is_file(), "missing opt-in source-pinned builder")

    def test_same_binary_linear_control_is_explicit(self):
        self.assertIn("STRINGTIE_INPUT_LOOKUP", (PACKAGE / "same_start_index.h").read_text())
        self.assertIn("StringTieSameStartIndex::enabled()", builder().NEW_LOOP)

    def test_build_root_cannot_be_inside_repository(self):
        with contextlib.redirect_stderr(io.StringIO()), self.assertRaises(SystemExit) as error:
            builder().main(["--build-root", str(ROOT / "build-stringtie-test-never-created"), "--prepare-only"])
        self.assertEqual(error.exception.code, 2)
        self.assertFalse((ROOT / "build-stringtie-test-never-created").exists())

    def test_build_flags_keep_environment_prefixes(self):
        mod = builder()
        flags = mod.build_flags({"CPPFLAGS": "-I/conda/include", "CFLAGS": "-O2 -fPIC",
                                 "CXXFLAGS": "-O3 -fPIC", "LDFLAGS": "-L/conda/lib -Wl,-rpath,/conda/lib"})
        self.assertIn("-I/conda/include", flags["cppflags"])
        self.assertIn("-O2 -fPIC", flags["cflags"])
        self.assertIn("-O3 -fPIC", flags["cxxflags"])
        self.assertIn("-L/conda/lib", flags["ldflags"])

    def test_bad_archive_is_rejected_before_extraction(self):
        mod = builder()
        with tempfile.TemporaryDirectory() as tmp:
            archive = Path(tmp) / "bad.tar.gz"
            archive.write_bytes(b"not an upstream archive")
            destination = Path(tmp) / "source"
            with self.assertRaisesRegex(ValueError, "SHA256"):
                mod.extract_verified(archive, destination, mod.STRINGTIE)
            self.assertFalse(destination.exists())

    def test_source_patch_is_exact_and_refuses_drift(self):
        archive = os.environ.get("STRINGTIE_INPUT_ARCHIVE")
        if not archive:
            self.skipTest("set STRINGTIE_INPUT_ARCHIVE for pinned-source checks")
        mod = builder()
        with tempfile.TemporaryDirectory() as tmp:
            source = mod.extract_verified(Path(archive), Path(tmp) / "source", mod.STRINGTIE)
            cpp = (source / "rlink.cpp").read_text()
            mod.apply_accelerator(source)
            patched = (source / "rlink.cpp").read_text()
            self.assertIn("bdata.same_start_candidates.find", patched)
            # Everything before and after the candidate enumeration stays byte-identical.
            self.assertEqual(cpp.split("\tif(!mergeMode) while(")[0], patched.split("\tif(!mergeMode) {")[0])
            self.assertEqual(cpp.split("\telse if((alndata.tinfo->cov")[1], patched.split("\telse if((alndata.tinfo->cov")[1])
            self.assertIn("match=exonmatch(readlist[n]->segs,brec.exons);", patched)
            self.assertIn("match=deljuncmatch(readlist[n],brec.juncsdel);", patched)
            header = (source / "rlink.h").read_text()
            self.assertEqual(header.count("same_start_candidates.clear();"), 2)
            with self.assertRaisesRegex(ValueError, "SHA256"):
                mod.apply_accelerator(source)

    def test_compiled_index_preserves_reverse_predicate_and_lifecycle(self):
        cxx = shutil.which(os.environ.get("CXX", "c++"))
        if not cxx:
            self.skipTest("C++11 compiler unavailable")
        with tempfile.TemporaryDirectory() as tmp:
            exe = Path(tmp) / "index-test"
            subprocess.run([cxx, "-std=c++11", "-O2", "-Wall", "-Wextra", "-Werror",
                            "-I", str(PACKAGE), str(PACKAGE / "tests/same_start.cpp"),
                            "-o", str(exe)], check=True)
            result = subprocess.run([str(exe)], check=True, text=True, capture_output=True)
            self.assertIn("candidate regressions passed", result.stdout)


if __name__ == "__main__":
    unittest.main()
