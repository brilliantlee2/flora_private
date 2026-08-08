import unittest
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]


class ReleaseLayoutTests(unittest.TestCase):
    def test_release_manifests_exist(self):
        for relative_path in [
            ".gitignore",
            "Cargo.toml",
            "Cargo.lock",
            "environment.yml",
            "environment.runtime.yml",
            "requirements.txt",
            "README.md",
            "README_zh-CN.md",
            "THIRD_PARTY_NOTICES.md",
            "licenses/Glycine-MIT.txt",
            "run_all.sh",
            "run_all_mixed_species.sh",
            "scripts/report_template.html",
            "scripts/plotly-2.26.0.min.js",
            "packaging/build_binary_release.sh",
            "packaging/compile_python_assets.py",
            "packaging/refresh_binary_release_metadata.sh",
            "packaging/prepare_public_repository.sh",
            "docs/repository-templates/public/README.md",
            "docs/repository-templates/public/README_zh-CN.md",
            "docs/repository-templates/private/README.md",
            "docs/repository-templates/private/README_zh-CN.md",
            "src/glycine/mod.rs",
            "runtime_manifest.json",
            "vendor/edlib_rs-0.1.2/Cargo.toml",
            "vendor/rust-htslib/Cargo.toml",
        ]:
            self.assertTrue((PROJECT_ROOT / relative_path).is_file(), relative_path)

    def test_binary_release_builds_and_stages_only_flora(self):
        script = (PROJECT_ROOT / "packaging/build_binary_release.sh").read_text(
            encoding="utf-8"
        )

        self.assertIn("--target x86_64-unknown-linux-gnu --bin flora", script)
        self.assertIn('"${STAGE_DIR}/flora"', script)
        self.assertIn('"${STAGE_DIR}/flora" run --help', script)
        self.assertIn('"${STAGE_DIR}/flora" run-mixed --help', script)
        self.assertIn("'^  flora run([[:space:]]|$)'", script)
        self.assertIn("'^  flora run-mixed([[:space:]]|$)'", script)
        self.assertIn("'^Usage: flora run([[:space:]]|$)'", script)
        self.assertIn("'^Usage: flora run-mixed([[:space:]]|$)'", script)
        self.assertNotIn("RUST_BINS=", script)
        self.assertNotIn('install -m 0755 "${ROOT_DIR}/run_all.sh"', script)
        self.assertNotIn('install -m 0755 "${ROOT_DIR}/run_all_mixed_species.sh"', script)
        self.assertNotIn('"${STAGE_DIR}/target/release"', script)

    def test_binary_release_uses_allowlisted_python_assets(self):
        script = (PROJECT_ROOT / "packaging/build_binary_release.sh").read_text(
            encoding="utf-8"
        )

        self.assertIn("PYTHON_RUNTIME_ASSETS=(", script)
        self.assertNotIn('for source_path in "${ROOT_DIR}"/scripts/*.py', script)
        self.assertIn("runtime_manifest.json", script)
        self.assertIn("compile_python_assets.py", script)

    def test_runners_require_report_static_assets(self):
        for runner_name in ["run_all.sh", "run_all_mixed_species.sh"]:
            runner = (PROJECT_ROOT / runner_name).read_text(encoding="utf-8")
            self.assertIn('report_template.html', runner)
            self.assertIn('plotly-2.26.0.min.js', runner)

    def test_glycine_is_integrated_into_flora(self):
        for runner_name in ["run_all.sh", "run_all_mixed_species.sh"]:
            runner = (PROJECT_ROOT / runner_name).read_text(encoding="utf-8")
            self.assertIn('"${FLORA_BIN}" glycine', runner)
            self.assertNotIn("--glycine-bin-dir", runner)
            self.assertNotIn("require_cmd glycine", runner)

    def test_readmes_have_language_switches_and_core_commands(self):
        public_docs = PROJECT_ROOT / "docs" / "repository-templates" / "public"
        english = (public_docs / "README.md").read_text(encoding="utf-8")
        chinese = (public_docs / "README_zh-CN.md").read_text(encoding="utf-8")

        for readme in [english, chinese]:
            self.assertIn("README.md", readme)
            self.assertIn("README_zh-CN.md", readme)
            self.assertIn("conda env create -f environment.yml", readme)
            self.assertIn("bash run_all.sh", readme)
            self.assertIn("bash run_all_mixed_species.sh", readme)

        for private_name in ["README.md", "README_zh-CN.md"]:
            private_readme = (PROJECT_ROOT / private_name).read_text(encoding="utf-8")
            self.assertIn("README.md", private_readme)
            self.assertIn("README_zh-CN.md", private_readme)
            self.assertIn("cargo build --release", private_readme)

    def test_environment_uses_python_311_and_scanpy_dependencies(self):
        environment = (PROJECT_ROOT / "environment.yml").read_text(encoding="utf-8")
        runtime_environment = (PROJECT_ROOT / "environment.runtime.yml").read_text(encoding="utf-8")
        requirements = (PROJECT_ROOT / "requirements.txt").read_text(encoding="utf-8")
        english = (PROJECT_ROOT / "README.md").read_text(encoding="utf-8")
        chinese = (PROJECT_ROOT / "README_zh-CN.md").read_text(encoding="utf-8")

        self.assertIn("python>=3.11,<3.12", environment)
        self.assertNotIn("python-gil", environment)
        self.assertIn("- nodefaults", environment)
        for dependency in [
            "scanpy>=1.11,<1.12",
            "numba>=0.66,<0.67",
            "python-igraph",
            "leidenalg",
        ]:
            self.assertIn(f"- {dependency}", environment)
        for dependency in ["scanpy>=1.11,<1.12", "numba>=0.66,<0.67", "igraph", "leidenalg"]:
            self.assertIn(dependency, requirements)
        for readme in [english, chinese]:
            self.assertIn("Python 3.11", readme)
            self.assertNotIn("3.14t", readme)
        self.assertRegex(environment, r"(?m)^\s*- rust(?:[=<>]|\s*$)")
        self.assertNotRegex(environment, r"(?m)^\s*- cargo(?:[=<>]|\s*$)")
        for dependency in ["samtools", "minimap2", "bedtools", "clang", "libclang", "pip"]:
            self.assertIn(f"- {dependency}", environment)
        for dependency in ["samtools", "minimap2", "bedtools", "scanpy>=1.11,<1.12", "pip"]:
            self.assertIn(f"- {dependency}", runtime_environment)
        for build_dependency in ["rust", "cmake", "c-compiler", "cxx-compiler", "clang", "libclang"]:
            self.assertNotRegex(runtime_environment, rf"(?m)^\s*- {build_dependency}(?:[=<>]|\s*$)")

    def test_large_runtime_data_are_ignored(self):
        gitignore = (PROJECT_ROOT / ".gitignore").read_text(encoding="utf-8")

        for pattern in ["target/", "*.fastq.gz", "*.bam", "report_new/", "report_new_2/", "vendor.zip"]:
            self.assertIn(pattern, gitignore)
        self.assertIn("/glycine/", gitignore)
        self.assertNotIn("\nglycine/", gitignore)

    def test_integrated_glycine_source_is_complete(self):
        for source_name in [
            "args.rs",
            "file_system.rs",
            "identifier.rs",
            "mod.rs",
            "qc.rs",
            "reader.rs",
            "utils.rs",
            "writer.rs",
        ]:
            self.assertTrue(
                (PROJECT_ROOT / "src" / "glycine" / source_name).is_file(),
                source_name,
            )

    def test_flora_package_metadata_are_consistent(self):
        cargo_toml = (PROJECT_ROOT / "Cargo.toml").read_text(encoding="utf-8")
        cargo_lock = (PROJECT_ROOT / "Cargo.lock").read_text(encoding="utf-8")
        runners = "\n".join(
            (PROJECT_ROOT / name).read_text(encoding="utf-8")
            for name in ["run_all.sh", "run_all_mixed_species.sh"]
        )

        self.assertRegex(cargo_toml, r'(?m)^name = "flora"$')
        self.assertRegex(cargo_toml, r'(?m)^version = "0\.1\.0"$')
        self.assertRegex(cargo_toml, r'(?ms)^\[profile\.release\].*?^strip = "symbols"$')
        self.assertIn('name = "flora"\nversion = "0.1.0"', cargo_lock)
        self.assertNotIn("StrintRust", runners)


if __name__ == "__main__":
    unittest.main()
