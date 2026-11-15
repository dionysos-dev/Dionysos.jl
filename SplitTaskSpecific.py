#!/usr/bin/env python3
"""Uniform calibration/test split of complete LTL semantic-equivalence classes.

Requires only Python's standard library and Spot's ltlfilt. No model/API access.
All record fields are retained; the four rewritten NL fields define the active
NL-LTL pair counts used by the full-benchmark experiments.
"""
import argparse
import collections
import concurrent.futures
import hashlib
import json
import math
import os
from pathlib import Path
import random
import shutil
import subprocess
import tempfile

ACTIVE_NL_FIELDS = [
    "natural_paraphrase_gpt54mini_rewritten",
    "paraphrase_gpt5.4-mini_gpt54mini_rewritten",
    "paraphrase_gemini-2.5-flash_gpt54mini_rewritten",
    "paraphrase_deepseek_gpt54mini_rewritten",
]
ATOMS = ["blue", "brown", "green", "purple", "yellow"]
PROBE_SEED = 20261026


def encoded(data):
    return (json.dumps(data, ensure_ascii=False, indent=2, allow_nan=False) + "\n").encode("utf-8")


def sha(data):
    return hashlib.sha256(data).hexdigest()


def spot_matches(spot, options, formulas):
    result = subprocess.run(
        [spot] + options + ["--format=%L"],
        input="\n".join(formulas) + "\n", text=True, capture_output=True, timeout=120,
    )
    if result.returncode not in (0, 1) or result.stderr.strip():
        raise RuntimeError("Spot failed; no split will be written: " + result.stderr.strip())
    indices = [int(line) - 1 for line in result.stdout.splitlines() if line.strip()]
    if any(i < 0 or i >= len(formulas) for i in indices):
        raise RuntimeError("Invalid Spot output indices")
    return indices


def probe_words():
    def valuation(mask):
        return " & ".join(ap if mask & (1 << i) else "!" + ap for i, ap in enumerate(ATOMS))
    words = ["cycle{" + valuation(m) + "}" for m in range(32)]
    for i in range(5):
        words += [valuation(0) + ";cycle{" + valuation(1 << i) + "}",
                  valuation(1 << i) + ";cycle{" + valuation(0) + "}"]
    rng = random.Random(PROBE_SEED)
    for _ in range(64):
        prefix = [valuation(rng.randrange(32)) for _ in range(rng.randrange(0, 5))]
        cycle = [valuation(rng.randrange(32)) for _ in range(rng.randrange(2, 8))]
        words.append(";".join(prefix + ["cycle{" + ";".join(cycle) + "}"]))
    return words


def semantic_groups(records, spot, workers=4):
    formulas = sorted({r["LTL"] for r in records})
    words = probe_words()
    signatures = [0] * len(formulas)
    def probe(word):
        return spot_matches(spot, ["--accept-word=" + word], formulas)
    with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as pool:
        for bit, matches in enumerate(pool.map(probe, words)):
            for index in matches:
                signatures[index] |= 1 << bit
    buckets = collections.defaultdict(list)
    for index, signature in enumerate(signatures):
        buckets[signature].append(index)

    # Different word-acceptance signatures prove inequivalence. Equal signatures
    # do NOT prove equivalence: resolve each remaining bucket with exact Spot tests.
    def resolve(bucket):
        remaining = list(bucket)
        groups, checks = [], 0
        while remaining:
            if len(remaining) == 1:
                groups.append(remaining)
                break
            local = spot_matches(spot, ["--equivalent-to=" + formulas[remaining[0]]],
                                 [formulas[i] for i in remaining])
            checks += len(remaining)
            if 0 not in local:
                raise RuntimeError("Spot did not match a formula with itself")
            matches = set(local)
            groups.append([i for j, i in enumerate(remaining) if j in matches])
            remaining = [i for j, i in enumerate(remaining) if j not in matches]
        return groups, checks

    classes, comparisons = [], 0
    with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as pool:
        for groups, checks in pool.map(resolve, list(buckets.values())):
            classes.extend(groups)
            comparisons += checks
    assert sorted(i for group in classes for i in group) == list(range(len(formulas)))
    ids = collections.defaultdict(list)
    for record in records:
        ids[record["LTL"]].append(record["id"])
    groups = []
    for group in classes:
        members = sorted(formulas[i] for i in group)
        groups.append({
            "class_id": "ltl-" + sha(encoded(members))[:16],
            "representative": members[0], "formulas": members,
            "record_ids": [rid for formula in members for rid in ids[formula]],
        })
    groups.sort(key=lambda g: tuple(g["formulas"]))
    assert len({g["class_id"] for g in groups}) == len(groups)
    return groups, {"lasso_probes": len(words), "probe_seed": PROBE_SEED,
                    "exact_spot_comparisons": comparisons}


def split_dataset(dataset_path, output_directory, seed=20261027,
                  calibration_fraction=0.7, workers=4):
    dataset_path, output_directory = Path(dataset_path).resolve(), Path(output_directory).absolute()
    if not isinstance(seed, int) or isinstance(seed, bool):
        raise ValueError("seed must be an integer")
    if not math.isfinite(calibration_fraction) or not 0 < calibration_fraction < 1:
        raise ValueError("calibration_fraction must lie strictly between 0 and 1")
    if workers < 1:
        raise ValueError("workers must be positive")
    if output_directory.exists():
        raise FileExistsError("Split directory already exists; refusing to overwrite: " + str(output_directory))
    source = dataset_path.read_bytes()
    records = json.loads(source)
    if not isinstance(records, list) or not records:
        raise ValueError("Benchmark must be a nonempty record array")
    seen = set()
    for record in records:
        rid = record.get("id")
        if isinstance(rid, bool) or not isinstance(rid, (str, int)) or rid in seen:
            raise ValueError("Missing, invalid, or duplicate record id: " + str(rid))
        seen.add(rid)
        formula = record.get("LTL")
        if not isinstance(formula, str) or not formula.strip() or "\n" in formula or "\r" in formula:
            raise ValueError("Missing or multiline reference LTL at record " + str(rid))
        for field in ACTIVE_NL_FIELDS:
            text = record.get(field)
            if not isinstance(text, str) or not text.strip():
                raise ValueError("Missing rewritten NL field at record " + str(rid) + ": " + field)
    spot = shutil.which("ltlfilt")
    if not spot:
        raise RuntimeError("Spot's ltlfilt must be available on PATH")
    groups, grouping_audit = semantic_groups(records, spot, workers)
    n_calibration = math.floor(calibration_fraction * len(groups) + 0.5)
    if not 0 < n_calibration < len(groups):
        raise ValueError("Not enough semantic classes to form both nonempty sets")
    shuffled = [g["class_id"] for g in groups]
    # One uniform shuffle; no stratification, balancing, or seed selection.
    random.Random(seed).shuffle(shuffled)
    calibration_ids = set(shuffled[:n_calibration])
    assignment = {}
    for group in groups:
        group["split"] = "calibration" if group["class_id"] in calibration_ids else "test"
        for rid in group["record_ids"]:
            assert rid not in assignment
            assignment[rid] = group["split"]
    assert set(assignment) == seen
    split_records = {name: [r for r in records if assignment[r["id"]] == name]
                     for name in ("calibration", "test")}
    nl_fields = [key for key in records[0]
                 if "paraphrase" in key or key == "back_translation"]
    manifest = {
        "schema_version": 1, "dataset_path": str(dataset_path), "dataset_sha256": sha(source),
        "reference_field": "LTL", "seed": seed,
        "calibration_fraction_requested": calibration_fraction,
        "test_fraction_requested": 1 - calibration_fraction,
        "rounding": "nearest integer number of classes; halves round upward",
        "randomization_unit": "semantic_equivalence_class",
        "sampling": "one uniform random shuffle without replacement; first n classes assigned to calibration",
        "grouping_semantics": "standard infinite-trace LTL, with atomic proposition names fixed",
        "grouping_method": "word-acceptance partition followed by exact Spot equivalence within each bucket",
        "spot_executable": spot,
        "spot_version": subprocess.check_output([spot, "--version"], text=True).strip(),
        "grouping_audit": grouping_audit,
        "record_count": len(records), "unique_formula_strings": len({r["LTL"] for r in records}),
        "semantic_class_count": len(groups), "active_nl_fields": ACTIVE_NL_FIELDS,
        "all_nl_fields_preserved": nl_fields,
        "records_preserved": "all original fields, values, and NL variants; original order within each set",
        "shuffled_class_ids": shuffled,
        "splits": {},
    }
    for name, subset in split_records.items():
        manifest["splits"][name] = {
            "file": name + ".json", "semantic_class_count": sum(g["split"] == name for g in groups),
            "record_count": len(subset), "active_nl_ltl_pair_count": len(subset) * len(ACTIVE_NL_FIELDS),
            "record_ids": [r["id"] for r in subset],
            "class_ids": [g["class_id"] for g in groups if g["split"] == name],
        }
    cal_ids = set(manifest["splits"]["calibration"]["record_ids"])
    test_ids = set(manifest["splits"]["test"]["record_ids"])
    assert not cal_ids & test_ids and cal_ids | test_ids == seen
    # The complete partition is exact: each equivalent formula and all records
    # containing it share one class, hence cannot cross the split boundary.
    manifest["validation"] = {
        "all_records_in_exactly_one_set": True, "semantic_class_overlap": 0,
        "formula_string_overlap": 0, "all_record_fields_preserved": True,
    }
    # Report known prior development usage without excluding or reassigning any
    # class: doing that would change the requested uniform full-benchmark split.
    development_path = dataset_path.parent.parent / "Results/task-specific-temperature-prompt-v2/manifest.json"
    if development_path.is_file():
        dev_bytes = development_path.read_bytes()
        dev = json.loads(dev_bytes)
        development_counts = collections.Counter(assignment.get(p["record_id"], "outside_benchmark")
                                                 for p in dev["prompts"])
        manifest["prior_development_overlap"] = {
            "manifest_path": str(development_path), "manifest_sha256": sha(dev_bytes),
            "prompt_count": len(dev["prompts"]),
            "prompt_counts_by_new_split": dict(development_counts),
            "assignment_changed_for_development_membership": False,
        }
    files = {name + ".json": encoded(subset) for name, subset in split_records.items()}
    files["semantic_groups.json"] = encoded(groups)
    manifest["file_sha256"] = {name: sha(content) for name, content in files.items()}
    files["manifest.json"] = encoded(manifest)
    cal, test = (manifest["splits"][name] for name in ("calibration", "test"))
    readme = f"""# Task-specific calibration/test split

Uniform random assignment of complete semantic-equivalence classes; seed {seed}.
Calibration receives the nearest integer to {calibration_fraction:.0%} of classes.
Equivalent formulas and every NL variant remain together. Proposition names are
fixed; equivalence is standard LTL over unconstrained infinite traces.

| Set | Semantic classes | Benchmark records | Active rewritten NL-LTL pairs |
|---|---:|---:|---:|
| Calibration | {cal['semantic_class_count']} | {cal['record_count']} | {cal['active_nl_ltl_pair_count']} |
| Test | {test['semantic_class_count']} | {test['record_count']} | {test['active_nl_ltl_pair_count']} |

The two dataset JSON files preserve complete original records, including all
original and rewritten NL fields and metadata. The pair counts above use the
four rewritten fields from the full-benchmark evaluation. Do not deduplicate
records or distribute individual paraphrases independently.

`semantic_groups.json` records equivalent formula strings, record IDs, and their
assigned set. `manifest.json` records the source hash, seed, full shuffled class
order, file hashes, counts, and checks. Record order inside each dataset follows
the source dataset. Different group sizes mean record/pair proportions need not
be exactly 70/30. This partitions the complete existing benchmark, including
examples previously used in development and full-benchmark evaluations.

Reproduce from the Dionysos repository root, using a new output directory:

```sh
python3 SplitTaskSpecific.py --seed {seed} --output /tmp/reproduced-task-specific-split
```

No model/API calls or controller synthesis are involved. Existing LLM outputs
can be assigned using the record IDs in the manifest.
"""
    files["README.md"] = readme.encode("utf-8")
    if "prior_development_overlap" in manifest:
        overlap = manifest["prior_development_overlap"]
        files["README.md"] += ("\nKnown prior development overlap: " +
            json.dumps(overlap["prompt_counts_by_new_split"], sort_keys=True) +
            " prompts. No class was reassigned based on this overlap.\n").encode("utf-8")
    output_directory.parent.mkdir(parents=True, exist_ok=True)
    staging = Path(tempfile.mkdtemp(prefix=".semantic-split-", dir=output_directory.parent))
    try:
        for name, content in files.items():
            (staging / name).write_bytes(content)
        os.rename(staging, output_directory)
    finally:
        if staging.exists():
            shutil.rmtree(staging)
    return manifest


def main():
    root = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dataset", type=Path, default=root / "Benchmark/Task_Spescific_Dataset.json")
    parser.add_argument("--output", type=Path,
                        default=root / "Benchmark/TaskSpecificSplits/semantic-70-30-seed-20261027")
    parser.add_argument("--seed", type=int, default=20261027)
    parser.add_argument("--calibration-fraction", type=float, default=0.7)
    parser.add_argument("--workers", type=int, default=4)
    args = parser.parse_args()
    result = split_dataset(args.dataset, args.output, args.seed, args.calibration_fraction, args.workers)
    print(json.dumps({"output": str(args.output.absolute()), "semantic_classes": result["semantic_class_count"],
                      "seed": result["seed"], "splits": {k: {n: v[n] for n in (
                          "semantic_class_count", "record_count", "active_nl_ltl_pair_count")}
                          for k, v in result["splits"].items()}}, indent=2))


if __name__ == "__main__":
    main()
