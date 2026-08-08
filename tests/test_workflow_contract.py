import json
import math
import unittest
from pathlib import Path


FIXTURE_DIR = Path(__file__).parent / "fixtures" / "cli"
COMMANDS = ("single", "mixed", "glycine", "analyze")
OPTION_FIELDS = {
    "name",
    "aliases",
    "type",
    "required",
    "required_when",
    "default",
    "conflicts",
    "description",
}
INVOCATION_FIELDS = {"argv", "reason"}
ALLOWED_TYPES = {"bool", "float", "int", "path", "string", "string_list"}


def load_contract(name):
    with (FIXTURE_DIR / f"{name}.json").open(encoding="utf-8") as handle:
        return json.load(handle)


def option_map(contract):
    return {option["name"]: option for option in contract["options"]}


class WorkflowContractTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.contracts = {name: load_contract(name) for name in COMMANDS}

    def test_fixture_schema_is_complete(self):
        for fixture_name, contract in self.contracts.items():
            with self.subTest(fixture=fixture_name):
                self.assertEqual(1, contract["schema_version"])
                self.assertEqual(fixture_name, contract["id"])
                self.assertIsInstance(contract["command"], list)
                self.assertTrue(contract["command"])
                self.assertIsInstance(contract["description"], str)
                self.assertTrue(contract["description"])
                self.assertIsInstance(contract["options"], list)
                self.assertTrue(contract["options"])
                self.assertEqual(
                    {"accepted", "rejected"}, set(contract["invocations"])
                )
                for option in contract["options"]:
                    self.assertTrue(OPTION_FIELDS.issubset(option))
                    self.assertTrue(option["name"].startswith("--"))
                    self.assertIsInstance(option["aliases"], list)
                    self.assertIn(option["type"], ALLOWED_TYPES)
                    self.assertIsInstance(option["required"], bool)
                    self.assertIsInstance(option["required_when"], list)
                    self.assertIsInstance(option["conflicts"], list)
                    self.assertTrue(option["description"])
                for outcome in ("accepted", "rejected"):
                    self.assertTrue(contract["invocations"][outcome])
                    for invocation in contract["invocations"][outcome]:
                        self.assertTrue(INVOCATION_FIELDS.issubset(invocation))
                        self.assertIsInstance(invocation["argv"], list)
                        self.assertTrue(invocation["reason"])

    def test_option_names_and_aliases_are_unique_per_command(self):
        for fixture_name, contract in self.contracts.items():
            seen = {}
            for option in contract["options"]:
                for spelling in (option["name"], *option["aliases"]):
                    with self.subTest(fixture=fixture_name, spelling=spelling):
                        self.assertNotIn(spelling, seen)
                    seen[spelling] = option["name"]

    def test_aliases_defaults_conditions_and_conflicts_reference_known_options(self):
        for fixture_name, contract in self.contracts.items():
            options = option_map(contract)
            all_spellings = {
                spelling
                for option in contract["options"]
                for spelling in (option["name"], *option["aliases"])
            }
            for option in contract["options"]:
                with self.subTest(fixture=fixture_name, option=option["name"]):
                    for alias in option["aliases"]:
                        self.assertTrue(alias.startswith("-"))
                    for conflict in option["conflicts"]:
                        self.assertIn(conflict, options)
                        self.assertIn(option["name"], options[conflict]["conflicts"])
                    for condition in option["required_when"]:
                        self.assertIn(condition["option"], all_spellings)
                        self.assertIn(condition["operator"], {"equals", "not_equals"})
                        self.assertIn("value", condition)
                    if option["required"]:
                        self.assertIsNone(option["default"])

    def test_legacy_aliases_and_representative_defaults_are_frozen(self):
        workflow_aliases = {
            "--tso-seq": ["--tso_seq"],
            "--rtp-seq": ["--rtp_seq"],
            "--barcode-list-10bp": ["--barcode_list_10bp"],
            "--ref-dir": ["--ref_dir"],
            "--gene-fasta": ["--genome-fa"],
            "--out-dir": ["--outdir"],
            "--sample-id": ["--sample"],
            "--threads": ["--thread"],
            "--barcode-extract-mode": ["--barcode_extract_mode"],
            "--glycine-err": ["--err"],
            "--glycine-shift": ["--shift"],
            "--min-len": ["--min_len"],
            "--umi-len": ["--umi_len"],
            "--help": ["-h"],
            "--version": ["-V"],
        }
        for fixture_name in ("single", "mixed"):
            options = option_map(self.contracts[fixture_name])
            aliases = {
                name: option["aliases"]
                for name, option in options.items()
                if option["aliases"]
            }
            with self.subTest(fixture=fixture_name):
                self.assertEqual(workflow_aliases, aliases)
                self.assertEqual(96, options["--threads"]["default"])
                self.assertTrue(options["--light-output"]["default"])
                self.assertFalse(options["--full-output"]["default"])
                self.assertIsNone(options["--python"]["default"])

        glycine = option_map(self.contracts["glycine"])
        self.assertEqual("0.25,0.25", glycine["--err"]["default"])
        self.assertEqual("100,100", glycine["--shift"]["default"])
        self.assertEqual(100, glycine["--min_len"]["default"])
        self.assertEqual(4, glycine["--thread"]["default"])

        analyze = option_map(self.contracts["analyze"])
        self.assertEqual("Flora", analyze["--out_dir"]["default"])
        self.assertEqual("fixed_seq", analyze["--barcode_extract_mode"]["default"])
        self.assertEqual(20, analyze["--min_reads_per_cell"]["default"])
        self.assertTrue(analyze["--absorb_unassigned_paired"]["default"])

    def test_legacy_toggle_pairs_are_ordered_overrides_not_conflicts(self):
        expected_rules = [
            "full-pipeline and upstream-only are accepted together; the last occurrence wins",
            "full-output and light-output are accepted together; the last occurrence wins",
        ]
        for fixture_name in ("single", "mixed"):
            contract = self.contracts[fixture_name]
            with self.subTest(fixture=fixture_name):
                self.assertEqual(expected_rules, contract["resolution_rules"][:2])
                self.assertFalse(
                    any(option["conflicts"] for option in contract["options"])
                )

    def test_each_contract_has_representative_acceptance_and_rejection(self):
        for fixture_name, contract in self.contracts.items():
            accepted = contract["invocations"]["accepted"]
            rejected = contract["invocations"]["rejected"]
            with self.subTest(fixture=fixture_name):
                self.assertGreaterEqual(len(accepted), 2)
                self.assertGreaterEqual(len(rejected), 2)
                self.assertTrue(all(item["argv"][: len(contract["command"])] == contract["command"] for item in accepted + rejected))

    def test_full_workflows_freeze_python_precedence(self):
        for fixture_name in ("single", "mixed"):
            contract = self.contracts[fixture_name]
            options = option_map(contract)
            with self.subTest(fixture=fixture_name):
                self.assertIn("--python", options)
                self.assertEqual("path", options["--python"]["type"])
                self.assertEqual(
                    ["--python", "FLORA_PYTHON", "python3"],
                    contract["interpreter"]["precedence"],
                )
                self.assertEqual("FLORA_PYTHON", contract["interpreter"]["environment"])
                self.assertEqual("python3", contract["interpreter"]["fallback"])
                self.assertEqual(
                    "all_python_stages",
                    contract["interpreter"]["scope"],
                )

    def test_glycine_bin_dir_is_an_explicit_full_workflow_rejection(self):
        for fixture_name in ("single", "mixed"):
            contract = self.contracts[fixture_name]
            rejected = contract["compatibility_exceptions"]["rejected_options"]
            invocations = contract["invocations"]["rejected"]
            with self.subTest(fixture=fixture_name):
                self.assertEqual(
                    "must_error_never_ignore", rejected["--glycine-bin-dir"]
                )
                self.assertNotIn("--glycine-bin-dir", option_map(contract))
                self.assertTrue(
                    any("--glycine-bin-dir" in item["argv"] for item in invocations)
                )

    def test_mixed_singlet_threshold_preserves_python_float_semantics(self):
        contract = self.contracts["mixed"]
        option = option_map(contract)["--singlet-threshold"]
        semantics = contract["singlet_threshold_semantics"]
        self.assertEqual("float", option["type"])
        self.assertEqual(0.9, option["default"])
        self.assertEqual("python_float", semantics["parser"])
        self.assertEqual(
            ["barnyard_qc", "parameters_tsv"], semantics["propagates_to"]
        )
        self.assertEqual([], semantics["range_constraints"])

        for value in semantics["accepted_values"]:
            with self.subTest(value=value):
                parsed = float(value)
                self.assertTrue(math.isfinite(parsed) or math.isnan(parsed) or math.isinf(parsed))
        for value in semantics["rejected_values"]:
            with self.subTest(value=value):
                with self.assertRaises(ValueError):
                    float(value)

        single = self.contracts["single"]
        self.assertNotIn("--singlet-threshold", option_map(single))

    def test_command_boundaries_are_isolated(self):
        glycine = self.contracts["glycine"]["execution"]
        self.assertEqual(["glycine"], glycine["stages"])
        self.assertEqual(
            ["reference", "alignment", "python"], glycine["does_not_require"]
        )

        analyze = self.contracts["analyze"]["execution"]
        self.assertEqual(
            ["barcode_assignment", "umi_assignment", "cell_assignment"],
            analyze["stages"],
        )
        self.assertEqual(
            ["glycine", "reference", "alignment", "python"],
            analyze["does_not_require"],
        )

        self.assertNotIn("--python", option_map(self.contracts["glycine"]))
        self.assertNotIn("--python", option_map(self.contracts["analyze"]))


if __name__ == "__main__":
    unittest.main()
