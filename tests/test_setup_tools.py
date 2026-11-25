import hashlib
import importlib.util
import io
import json
import os
from pathlib import Path
import tarfile
import tempfile
import unittest
from unittest.mock import patch
from urllib.error import URLError


ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("setup_tools", ROOT / "packaging/setup_tools.py")
SETUP = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(SETUP)


def executable(version):
    return f"#!/bin/sh\nprintf '%s\\n' '{version}'\n".encode()


class SetupToolsTests(unittest.TestCase):
    def test_directory_precedence(self):
        with patch.dict(os.environ, {"FLORA_TOOLS_DIR": "/chosen", "XDG_DATA_HOME": "/xdg"}, clear=True):
            self.assertEqual(SETUP.resolve_tools_dir(Path("/explicit")), Path("/explicit"))
            self.assertEqual(SETUP.resolve_tools_dir(None), Path("/chosen"))
            del os.environ["FLORA_TOOLS_DIR"]
            self.assertEqual(SETUP.resolve_tools_dir(None), Path("/xdg/flora/tools"))
            del os.environ["XDG_DATA_HOME"]
            with patch.object(Path, "home", return_value=Path("/user")):
                self.assertEqual(SETUP.resolve_tools_dir(None), Path("/user/.local/share/flora/tools"))

    def fixtures(self, root):
        asset_root = root / "assets"
        asset_root.mkdir()
        archives = root / "archives"
        archives.mkdir()
        config = json.loads((ROOT / "runtime-tools.json").read_text())
        bonefish = config["bonefish"]["releases"]["linux-x86_64"]
        archive = archives / bonefish["archive"]
        with tarfile.open(archive, "w:gz") as tar:
            data = executable("bonefish 0.5.0")
            item = tarfile.TarInfo(bonefish["member"])
            item.size = len(data)
            tar.addfile(item, io.BytesIO(data))
        bonefish["sha256"] = hashlib.sha256(archive.read_bytes()).hexdigest()
        for name in ("source", "libdeflate"):
            spec = config["stringtie"][name]
            data = name.encode()
            (archives / spec["archive"]).write_bytes(data)
            spec["sha256"] = hashlib.sha256(data).hexdigest()
        (asset_root / "runtime-tools.json").write_text(json.dumps(config))
        return asset_root, archives

    def test_offline_install_repeat_and_integrity_failure(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            assets, archives = self.fixtures(root)
            tools = root / "tools"
            tools.mkdir()
            (tools / "bin").mkdir()
            note = tools / "bin/independent"
            note.write_text("keep")

            def build(_assets, _source, _libdeflate, _build_dir, _jobs, output):
                output.write_bytes(executable("2.2.2") + b"# STRINGTIE_INPUT_LOOKUP\n")
                output.chmod(0o755)
                return {"patch": "same-start-candidate-index-v1",
                        "index_sha256": "5035220a2872d4155c9569345f713bdc65d176772c48592705c29b2a1dc187bb"}

            with patch.object(SETUP, "build_stringtie", side_effect=build) as builder:
                self.assertTrue(SETUP.install(tools, archives, 1, assets, "linux-x86_64"))
                self.assertFalse(SETUP.install(tools, None, 1, assets, "linux-x86_64"))
                self.assertEqual(builder.call_count, 1)
                self.assertTrue((tools / "bin/bonefish").is_file())
                self.assertTrue((tools / "bin/stringtie-2.2.2-same-start").is_file())
                self.assertEqual(note.read_text(), "keep")
                (tools / "bin/bonefish").write_text("modified")
                with self.assertRaisesRegex(ValueError, "integrity"):
                    SETUP.install(tools, archives, 1, assets, "linux-x86_64")
                self.assertEqual((tools / "bin/bonefish").read_text(), "modified")

    def test_unmanaged_executable_is_preserved(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            assets, archives = self.fixtures(root)
            tools = root / "tools"
            (tools / "bin").mkdir(parents=True)
            existing = tools / "bin/bonefish"
            existing.write_text("external")
            with self.assertRaisesRegex(ValueError, "already exists"):
                SETUP.install(tools, archives, 1, assets, "linux-x86_64")
            self.assertEqual(existing.read_text(), "external")
            self.assertFalse((tools / "installed-tools.json").exists())

    def test_online_download_retries_a_transport_failure(self):
        data = b"pinned archive"
        spec = {
            "archive": "test.tar.gz",
            "url": "https://example.invalid/test.tar.gz",
            "sha256": hashlib.sha256(data).hexdigest(),
        }
        with tempfile.TemporaryDirectory() as temp:
            with patch.object(SETUP, "build_opener") as opener, patch.object(SETUP.time, "sleep"):
                opener.return_value.open.side_effect = [URLError("temporary failure"), io.BytesIO(data)]
                path = SETUP.obtain_archive(spec, None, Path(temp))
                self.assertEqual(path.read_bytes(), data)
                self.assertEqual(opener.return_value.open.call_count, 2)

    def test_bad_archive_does_not_publish_outputs(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            assets, archives = self.fixtures(root)
            (archives / "libdeflate.tar.gz").write_bytes(b"changed")
            with self.assertRaisesRegex(ValueError, "SHA-256 mismatch"):
                SETUP.install(root / "tools", archives, 1, assets, "linux-x86_64")
            self.assertFalse((root / "tools/bin").exists())
            self.assertFalse((root / "tools/.setup-tools.lock").exists())

    def test_publish_removes_its_own_output_on_a_racing_external_collision(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            stage, bin_dir = root / "stage", root / "bin"
            stage.mkdir()
            names = (SETUP.BONEFISH, SETUP.STRINGTIE)
            for name in (*names, SETUP.RECEIPT):
                (stage / name).write_text("new")
            external = bin_dir / SETUP.STRINGTIE
            link = os.link

            def introduce_collision(source, target):
                if target == external:
                    external.write_text("external")
                link(source, target)

            receipt = root / SETUP.RECEIPT
            with patch.object(SETUP.os, "link", side_effect=introduce_collision):
                with self.assertRaises(FileExistsError):
                    SETUP.publish(stage, bin_dir, receipt, names)
            self.assertFalse((bin_dir / SETUP.BONEFISH).exists())
            self.assertEqual(external.read_text(), "external")
            self.assertFalse(receipt.exists())

    def test_archive_symlink_is_rejected(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            archive = root / "unsafe.tar.gz"
            member = "bonefish/bonefish"
            with tarfile.open(archive, "w:gz") as tar:
                item = tarfile.TarInfo(member)
                item.type = tarfile.SYMTYPE
                item.linkname = "/some/other/executable"
                tar.addfile(item)
            with self.assertRaisesRegex(ValueError, "regular"):
                SETUP.extract_bonefish(archive, member, root / "output")
            self.assertFalse((root / "output").exists())


if __name__ == "__main__":
    unittest.main()
