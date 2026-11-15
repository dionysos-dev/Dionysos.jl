"""Validate the completed run and install its independent audit documentation."""
import collections
import importlib.util
import json
from pathlib import Path
import shutil
import sys

root=Path(sys.argv[1]).resolve()
output=Path(sys.argv[2]).resolve()
spec=importlib.util.spec_from_file_location("frozen_abcd",root/"DevelopmentABCD/run.py")
runner=importlib.util.module_from_spec(spec)
spec.loader.exec_module(runner)
config=runner.read(output/"frozen-development.json")
summary=runner.read(output/"summary.json")
assert runner.freeze(root/"Results/task-specific-temperature-prompt-v2",root,output,120)==config
assert summary["complete"] and len(summary["results"])==200
assert len({(r["prompt_index"],r["k"]) for r in summary["results"]})==200
jobs=list((output/"job-results").glob("*.json"))
assert len(jobs)==83
assert all(sum(row["counts"].values())==40 and row["pending"]==0 for row in summary["summary"])
assert all(row["issuance_rate"]+row["abstention_rate"]==1 for row in summary["summary"])
for a in summary["results"]:
    for b in summary["results"]:
        if a["prompt_index"]==b["prompt_index"] and a["k"]<b["k"] and a["outcome"] in ("plant_unrealizable","logically_inconsistent"):
            assert not b["issued"]
issued=[r for r in summary["results"] if r["issued"]]
files=sorted({r["controller_file"] for r in issued})
assert len(files)==24 and len(issued)==76
for row in issued:
    assert runner.digest(output/row["controller_file"])==row["controller_sha256"]
audit_log=Path("/tmp/deepseek-controller-audit.log")
assert "448    448" in audit_log.read_text()
assert "Independently re-verified 24 controller artifacts for 76 issued prompt/budget cases." in audit_log.read_text()
shutil.copy2(audit_log,output/"controller-audit.log")
for name in ("README.md","audit_controllers.jl","finalize.py"):
    shutil.copy2(Path(__file__).with_name(name),root/"DevelopmentABCD"/name)
validation=dict(complete=True,prompt_budget_cases=200,unique_synthesis_jobs=83,
    issued_cases=76,independently_reverified_controller_artifacts=24,
    independent_artifact_checks_passed=448,source_and_frozen_config_verified=True,
    timeout_cases=sum(r["subreason"]=="timeout" for r in summary["results"]),
    computational_failure_subreasons=dict(collections.Counter(r["subreason"] for r in summary["results"] if r["outcome"]=="computational_failure")),
    no_additional_model_api_calls=True,mistral_experiment_run=False,
    plant_sha256=runner.digest(output/"plant.jls"),
    frozen_config_file_sha256=runner.digest(output/"frozen-development.json"))
runner.save(output/"validation.json",validation)
print(json.dumps(validation,indent=2))
