#!/usr/bin/env python3
"""Frozen DeepSeek development ABCD experiment. Python standard library only.

Julia runs in the destination Dionysos environment. No model/API calls exist in
this runner. A persistent, isolated Julia process is killed on a job timeout.
"""
import argparse
import collections
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import signal
import subprocess
import sys
import time

OUTCOMES = ("controller_issued", "logically_inconsistent", "plant_unrealizable",
            "computational_failure", "no_valid_candidate")
LABELS = ("Controller issued", "Logically inconsistent", "Plant-unrealizable",
          "Timeout / computational failure", "No valid candidate")


def read(path):
    return json.loads(Path(path).read_text())


def save(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, ensure_ascii=False) + "\n")
    temporary.replace(path)


def digest(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def signature(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, ensure_ascii=False).encode()).hexdigest()


def freeze(source, root, output, timeout):
    manifest = read(source / "manifest.json")
    assert manifest["model"] == "deepseek-v4-pro"
    assert manifest["prompt_count"] == len(manifest["prompts"]) == 40
    assert manifest["repetitions"] == 10 and 0.7 in manifest["temperatures"]
    assert manifest["thinking_type"] == "disabled"
    assert manifest["prompt_version"] == 2
    assert set(manifest["atomic_propositions"]) == {"blue", "brown", "green", "purple", "yellow"}
    analysis = read(source / "conjunction-curves-sat-filtered" / "summary.json")
    report = next(r for r in analysis["temperature_results"] if r["temperature"] == 0.7)
    assert report["candidate_filter"] == "parse_then_allowed_APs_then_SAT_before_clustering"
    assert report["missing_policy"] == "saturate"
    assert len(report["prompt_results"]) == 40 and report["max_k"] >= 5
    cases, raw_hashes = [], {}
    for prompt, result in zip(manifest["prompts"], report["prompt_results"]):
        assert result["prompt_index"] == prompt["prompt_index"]
        assert result["record_id"] == prompt["record_id"]
        assert result["reference_ltl"] == prompt["ground_truth_ltl"]
        relative = f'details/prompt_{prompt["prompt_index"]:03d}_T_0.7.json'
        detail = read(source / relative)
        assert digest(source / relative) == result["source_sha256"]
        assert detail["manifest_sha256"] == report["manifest_sha256"]
        assert detail["prompt"] == prompt and detail["temperature"] == 0.7
        assert [s["repetition"] for s in detail["samples"]] == list(range(1, 11))
        raw_hashes[relative] = digest(source / relative)
        clusters = result["ranked_clusters"]
        assert clusters == sorted(clusters, key=lambda c: (-c["frequency"], min(c["sample_indices"])))
        assert all(c["satisfiable"] for c in clusters)
        assert sum(c["frequency"] for c in clusters) == result["satisfiable_samples"]
        for k in range(1, 6):
            included = clusters[:k]
            formula = " & ".join(f'({c["representative"]})' for c in included) if included else None
            if included:
                assert formula == result["conjunctions"][min(k, len(clusters)) - 1]["phi_k"]
            cases.append(dict(prompt_index=prompt["prompt_index"], record_id=prompt["record_id"],
                input_field=prompt["input_field"], source_text=prompt["source_text"],
                reference_ltl=prompt["ground_truth_ltl"], k=k, formula=formula,
                plant_id="interactive_vehicle_shared", available_candidates=result["satisfiable_samples"],
                available_clusters=len(clusters), included_clusters=len(included)))
    code = {}
    for p in [Path(__file__), Path(__file__).with_name("worker.jl"),
              root / "evaluate_claude_buchi_synthesis.jl", root / "System_helper.jl", root / "InteractiveWithFeedBack.jl",
              root / "Project.toml", root / "Manifest.toml", *sorted((root / "src").rglob("*.jl"))]:
        key = ("runner/" + p.name) if p.parent == Path(__file__).parent else str(p.relative_to(root))
        code[key] = digest(p)
    config = dict(schema_version=1, purpose="development_only", primary_model="deepseek-v4-pro",
        temperature=0.7, repetitions=10, prompt_count=40, budgets=list(range(1, 6)),
        ranking="semantic_cluster_frequency_descending_then_earliest_sample",
        candidate_filter="parse_then_allowed_APs_then_individual_SAT_before_clustering",
        missing_policy="saturate", issuance_rule="initial_state_winning_and_universal_infinite_verification_and_export",
        abstention_outcomes=list(OUTCOMES[1:]), timeout_seconds=timeout,
        timeout_scope="per_unique_formula_SAT_translation_game_verification_export; startup_and_JIT_warmup_excluded",
        plant=dict(id="interactive_vehicle_shared", assignment="same plant for all 40 prompts",
            source="fresh abstraction matching InteractiveWithFeedBack.jl configuration",
            state_box=[[0, 0, -math.pi - 0.4], [10, 10, math.pi + 0.4]],
            grid_origin=[0, 0, math.pi / 2], state_spacing=[0.2, 0.2, 0.3],
            input_spacing=[0.3, 0.3], time_step=0.3, obstacle_margin=0.12,
            initial_state=[4.8, 5.0, math.pi / 2], semantics="robust",
            labelling="existing LTL2Automata.labeler_and_inner_sets; finite labelled abstraction"),
        reserved_cross_model=dict(model="mistral-large-4", temperature=0.7, max_k=4, run=False),
        source_manifest_sha256=digest(source / "manifest.json"), raw_sample_sha256=raw_hashes,
        conjunction_report_sha256=digest(source / "conjunction-curves-sat-filtered" / "summary.json"),
        code_sha256=code, cases=cases)
    path = output / "frozen-development.json"
    if path.exists() and read(path) != config:
        raise RuntimeError("Frozen experiment changed; choose a new output directory")
    save(path, config)
    return config


def summarize(config, results, output):
    rows = []
    for k in range(1, 6):
        selected = [r for r in results if r["k"] == k]
        counts = {outcome: sum(r["outcome"] == outcome for r in selected) for outcome in OUTCOMES}
        assert sum(counts.values()) == len(selected)
        issued = counts["controller_issued"]
        rows.append(dict(k=k, completed=len(selected), planned=40, counts=counts,
            rates={name: count / 40 for name, count in counts.items()},
            issuance_rate=issued / 40, abstention_rate=(len(selected) - issued) / 40,
            pending=40-len(selected),
            failure_subreasons=dict(collections.Counter(r["subreason"] for r in selected
                if r["outcome"] == "computational_failure"))))
    summary = dict(schema_version=1, model="deepseek-v4-pro", temperature=0.7,
        frozen_config_sha256=signature(config), complete=len(results) == 200,
        denominator=40, caveat="plant-unrealizable means losing on this frozen labelled abstraction and initial state; not a proof of physical impossibility on the continuous plant",
        summary=rows, results=results)
    save(output / "summary.json", summary)
    with (output / "summary.csv").open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["k", "completed", "pending", *OUTCOMES, "issuance_rate", "abstention_rate"])
        for r in rows:
            writer.writerow([r["k"], r["completed"], r["pending"], *[r["counts"][o] for o in OUTCOMES],
                             r["issuance_rate"], r["abstention_rate"]])
    text = ["# Frozen DeepSeek development ABCD", "",
        "40 development prompts; T=0.7; R=10; k=1..5. Existing model outputs only; no API calls.", "",
        "All rates use 40 prompts per k. Only verified, exported controllers count as issuance. Every other completed case abstains; pending cases are shown separately.", "",
        "| k | Complete | Issued | Inconsistent | Plant-unrealizable | Timeout / computational failure | No candidate | Issuance | Abstention |",
        "|---|---:|---:|---:|---:|---:|---:|---:|---:|"]
    for r in rows:
        text.append(f'| {r["k"]} | {r["completed"]}/40 | ' + " | ".join(str(r["counts"][o]) for o in OUTCOMES)
                    + f' | {100*r["issuance_rate"]:.1f}% | {100*r["abstention_rate"]:.1f}% |')
    text += ["", f'Complete: {summary["complete"]}. Timeout: {config["timeout_seconds"]} seconds per distinct formula, including SAT, translation, game, verification, and export; backend initialization and warm-up excluded.', "",
        "The same fresh interactive vehicle abstraction is used for all prompts. Plant-unrealizability denotes loss from its configured initial state on this abstraction; abstraction conservatism can prevent synthesis even when a continuous controller exists.", "",
        "Unsupported specification automata, resource errors, timeouts, and verification failures belong to computational failure, with detailed subreasons in summary.json. No inference of plant-unrealizability is made from these failures.", "",
        "Saturated duplicate conjunctions reuse the same plant/formula job and exported controller; results remain prompt-weighted. Controller .jls artifacts store the product policy, automaton memory/update, winning set, initial states, and verification certificate; plant.jls stores adjacency, labels, state centers, and physical input vectors. Guarantees refer to this labelled abstraction.", "",
        "Mistral T=0.7, k_max=4 is reserved and was not run."]
    (output / "comparison.md").write_text("\n".join(text) + "\n")
    return summary


class Worker:
    def __init__(self, root, output, startup_timeout=300):
        self.root, self.output = root, output
        self.log = (output / "worker.log").open("a")
        ready = output / "worker-ready.json"
        ready.unlink(missing_ok=True)
        self.process = subprocess.Popen(["julia", f"--project={root}", "--threads=4",
            str(Path(__file__).with_name("worker.jl")), str(root), str(output)],
            stdin=subprocess.PIPE, stdout=self.log, stderr=subprocess.STDOUT,
            text=True, start_new_session=True)
        end = time.monotonic() + startup_timeout
        while not ready.exists():
            if self.process.poll() is not None:
                self.stop()
                raise RuntimeError("Julia backend startup failed; inspect worker.log")
            if time.monotonic() > end:
                self.stop()
                raise RuntimeError("Julia backend startup exceeded 300 seconds; inspect worker.log")
            time.sleep(0.1)

    def stop(self):
        if self.process.poll() is None:
            os.killpg(self.process.pid, signal.SIGTERM)
            try:
                self.process.wait(timeout=5)
            except subprocess.TimeoutExpired:
                os.killpg(self.process.pid, signal.SIGKILL)
                self.process.wait()
        if self.process.stdin is not None and not self.process.stdin.closed:
            try:
                self.process.stdin.close()
            except BrokenPipeError:
                pass
        self.log.close()

    def solve(self, job, job_path, result_path, timeout):
        self.process.stdin.write(json.dumps(dict(job_file=str(job_path), result_file=str(result_path))) + "\n")
        self.process.stdin.flush()
        end = time.monotonic() + timeout
        while not result_path.exists():
            if self.process.poll() is not None:
                return dict(job_id=job["job_id"], formula=job["formula"], outcome="computational_failure",
                    issued=False, subreason="worker_process_failure", controller_file=None,
                    error=f"Julia exited with code {self.process.returncode}")
            if time.monotonic() > end:
                stage = read(self.output / "worker-stage.json") if (self.output / "worker-stage.json").exists() else {}
                self.stop()
                return dict(job_id=job["job_id"], formula=job["formula"], outcome="computational_failure",
                    issued=False, subreason="timeout", controller_file=None, elapsed_seconds=timeout,
                    last_stage=stage.get("stage"), error=f"Exceeded frozen {timeout}-second wall-clock limit")
            time.sleep(0.1)
        return read(result_path)


def run(root, source, output, timeout):
    output.mkdir(parents=True, exist_ok=True)
    import fcntl
    with (output / "experiment.lock").open("w") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        config = freeze(source, root, output, timeout)
        results, worker = [], None
        try:
            summarize(config, results, output)
            print("Frozen: 40 prompts, 200 cases; building/loading the shared plant.", flush=True)
            worker = Worker(root, output)
            plant_hash = digest(output / "plant.jls")
            plant_identity = dict(sha256=plant_hash, backend=read(output / "worker-ready.json"))
            plant_identity["backend"].pop("pid", None)
            plant_path = output / "plant-identity.json"
            if plant_path.exists() and read(plant_path) != plant_identity:
                raise RuntimeError("Frozen plant identity changed")
            save(plant_path, plant_identity)
            for case in config["cases"]:
                if case["formula"] is None:
                    outcome = dict(outcome="no_valid_candidate", issued=False, subreason="empty_filtered_pool",
                        controller_file=None, error=None)
                else:
                    job_id = signature(dict(formula=case["formula"],plant_sha256=plant_hash,
                        frozen_config_sha256=signature(config)))
                    job = dict(job_id=job_id, formula=case["formula"], plant_sha256=plant_hash)
                    job_path, result_path = output / "jobs" / f"{job_id}.json", output / "job-results" / f"{job_id}.json"
                    save(job_path, job)
                    if result_path.exists():
                        outcome = read(result_path)
                        assert outcome["formula"] == case["formula"] and outcome["job_id"] == job_id
                    else:
                        if worker.process.poll() is not None:
                            worker = Worker(root, output)
                            assert digest(output / "plant.jls") == plant_hash
                        outcome = worker.solve(job, job_path, result_path, timeout)
                        save(result_path, outcome)
                    assert outcome["outcome"] in OUTCOMES
                    assert outcome["issued"] == (outcome["outcome"] == "controller_issued")
                    if outcome["issued"]:
                        assert digest(output / outcome["controller_file"]) == outcome["controller_sha256"]
                row = dict(case, **{k: v for k, v in outcome.items() if k != "formula"})
                results.append(row)
                summarize(config, results, output)
                print(f'{len(results):03d}/200 prompt={case["prompt_index"]:02d} k={case["k"]} {row["outcome"]} ({row["subreason"]})', flush=True)
        finally:
            if worker is not None:
                worker.stop()
            summarize(config, results, output)
    return summarize(config, results, output)


def install(root, source):
    destination = root / "DevelopmentABCD"
    destination.mkdir(exist_ok=True)
    for name in ("run.py", "worker.jl", "README.md", "test_runner.py", "test_worker.jl"):
        src = Path(__file__).with_name(name)
        if src.exists() and src.resolve() != (destination / name).resolve():
            shutil.copy2(src, destination / name)
    target = root / "Results" / source.name
    if source.resolve() != target.resolve():
        for src in source.rglob("*"):
            if not src.is_file():
                continue
            dst = target / src.relative_to(source)
            if dst.exists():
                if digest(dst) != digest(src):
                    raise RuntimeError(f"Refusing to replace differing experiment data: {dst}")
            else:
                dst.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(src, dst)
    return target


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dionysos-root", type=Path, required=True)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--timeout", type=float, default=120)
    parser.add_argument("--install", action="store_true")
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    if not math.isfinite(args.timeout) or args.timeout <= 0:
        parser.error("timeout must be finite and positive")
    root, source = args.dionysos_root.resolve(), args.source.resolve()
    if args.install:
        source = install(root, source)
    if args.self_test:
        subprocess.run([sys.executable,str(Path(__file__).with_name("test_runner.py"))],check=True)
        with __import__("tempfile").TemporaryDirectory() as folder:
            subprocess.run(["julia",f"--project={root}","--threads=4",
                str(Path(__file__).with_name("test_worker.jl")),str(root),folder],check=True)
    output = args.output.resolve() if args.output else root / "Results" / "deepseek-development-abcd"
    run(root, source, output, args.timeout)


if __name__ == "__main__":
    main()
