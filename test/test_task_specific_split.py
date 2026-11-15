import hashlib
import json
from pathlib import Path
import random
import shutil
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import SplitTaskSpecific as split


@unittest.skipUnless(shutil.which("ltlfilt"), "Spot's ltlfilt is required")
class SemanticSplitTests(unittest.TestCase):
    def fixture(self, path):
        formulas = ["F(blue)", "F((X(purple) U blue))", "F(blue)",
                    "G(yellow)", "!F(!yellow)", "F(green)",
                    "X(purple)", "brown", "G(blue)"]
        records = []
        for i, formula in enumerate(formulas, 1):
            record = {"id": i, "LTL": formula, "metadata": {"original": [i, None, True]}}
            for field in split.ACTIVE_NL_FIELDS + ["back_translation", "natural_paraphrase",
                    "paraphrase_gpt5.4-mini", "paraphrase_gemini-2.5-flash", "paraphrase_deepseek"]:
                record[field] = f"Original NL text: record {i}, field {field}."
            records.append(record)
        path.write_text(json.dumps(records))
        return records

    def test_exact_semantic_groups_uniform_assignment_and_preserved_records(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / "source.json"
            records = self.fixture(source)
            source_hash = hashlib.sha256(source.read_bytes()).hexdigest()
            output = root / "split"
            result = split.split_dataset(source, output, seed=20261027)
            self.assertEqual(result["semantic_class_count"], 6)
            self.assertEqual(result["unique_formula_strings"], 8)
            self.assertEqual(result["splits"]["calibration"]["semantic_class_count"], 4)
            self.assertEqual(result["splits"]["test"]["semantic_class_count"], 2)
            groups = json.loads((output / "semantic_groups.json").read_text())
            by_record = {rid: group for group in groups for rid in group["record_ids"]}
            self.assertEqual(by_record[1]["class_id"], by_record[2]["class_id"])
            self.assertEqual(by_record[1]["class_id"], by_record[3]["class_id"])
            self.assertEqual(by_record[4]["class_id"], by_record[5]["class_id"])
            self.assertNotEqual(by_record[1]["class_id"], by_record[9]["class_id"])
            self.assertEqual(by_record[1]["split"], by_record[2]["split"])
            self.assertEqual(by_record[1]["split"], by_record[3]["split"])
            expected_shuffle = [group["class_id"] for group in groups]
            random.Random(20261027).shuffle(expected_shuffle)
            self.assertEqual(result["shuffled_class_ids"], expected_shuffle)
            self.assertEqual(set(result["splits"]["calibration"]["class_ids"]), set(expected_shuffle[:4]))
            calibration = json.loads((output / "calibration.json").read_text())
            test = json.loads((output / "test.json").read_text())
            self.assertEqual(sorted(calibration + test, key=lambda r: r["id"]), records)
            self.assertFalse({r["id"] for r in calibration} & {r["id"] for r in test})
            self.assertFalse({r["LTL"] for r in calibration} & {r["LTL"] for r in test})
            self.assertEqual(sum(s["active_nl_ltl_pair_count"] for s in result["splits"].values()), 36)
            for name, digest in result["file_sha256"].items():
                self.assertEqual(hashlib.sha256((output / name).read_bytes()).hexdigest(), digest)
            self.assertEqual(hashlib.sha256(source.read_bytes()).hexdigest(), source_hash)
            # Same seed produces the identical datasets and group assignment.
            other = root / "same-seed"
            repeated = split.split_dataset(source, other, seed=20261027, workers=1)
            self.assertEqual(repeated["shuffled_class_ids"], result["shuffled_class_ids"])
            for name in ("calibration.json", "test.json", "semantic_groups.json"):
                self.assertEqual((output / name).read_bytes(), (other / name).read_bytes())
            before = (output / "manifest.json").read_bytes()
            with self.assertRaises(FileExistsError):
                split.split_dataset(source, output, seed=1)
            self.assertEqual((output / "manifest.json").read_bytes(), before)

    def test_invalid_benchmark_never_creates_a_split(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / "source.json"
            records = self.fixture(source)
            records[1]["id"] = records[0]["id"]
            source.write_text(json.dumps(records))
            with self.assertRaises(ValueError):
                split.split_dataset(source, root / "duplicate")
            self.assertFalse((root / "duplicate").exists())
            records = self.fixture(source)
            records[0]["LTL"] = "F("
            source.write_text(json.dumps(records))
            with self.assertRaises(RuntimeError):
                split.split_dataset(source, root / "invalid")
            self.assertFalse((root / "invalid").exists())


if __name__ == "__main__":
    unittest.main()
