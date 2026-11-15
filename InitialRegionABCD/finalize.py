"""Audit a completed all-initial-region run and write prompt-level diagnostics."""
import argparse
import collections
import csv
import hashlib
import json
from pathlib import Path
import shutil
import subprocess


def read(path):return json.loads(path.read_text())
def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dionysos-root",type=Path,required=True)
    parser.add_argument("--experiment",type=Path,required=True)
    parser.add_argument("--semantic-audit",type=Path,required=True)
    args=parser.parse_args()
    root=args.dionysos_root.resolve();output=args.experiment.resolve()
    config=read(output/"frozen-development.json");summary=read(output/"summary.json")
    assert summary["complete"] and len(summary["results"])==200
    for key,expected in config["code_sha256"].items():
        relative="DevelopmentABCD/worker.jl" if key=="plant_builder" else "InitialRegionABCD/"+key[7:] if key.startswith("runner/") else key
        assert digest(root/relative)==expected,(relative,"code changed")
    semantic=read(args.semantic_audit)
    assert semantic["complete"] and semantic["reference_checks"]==200 and semantic["prefix_checks"]==160
    assert not semantic["reference_mismatches"] and not semantic["prefix_violations"]
    shutil.copy2(args.semantic_audit,output/"semantic-audit.json")
    # Add the reusable review utilities without changing any frozen solver file.
    for name in ("finalize.py","verify_logic.py"):
        source=Path(__file__).with_name(name);target=root/"InitialRegionABCD"/name
        if source.resolve()!=target.resolve():shutil.copy2(source,target)
    with (output/"region-audit.log").open("w") as stream:
        subprocess.run(["julia",f"--project={root}","--threads=4",str(root/"InitialRegionABCD/audit_regions.jl"),str(root),str(output)],
                       stdout=stream,stderr=subprocess.STDOUT,check=True)
    audit=summary["monotonicity"]
    assert audit["adjacent_pairs_checked"]==160
    assert not audit["winning_region_inclusion_violations"]
    assert not audit["saturation_violations"]
    assert not audit["issuance_criterion_violations"]
    cases=summary["results"]
    index={(r["prompt_index"],r["k"]):r for r in cases}
    prompt_rows=[]
    for i in range(1,41):
        rs=[index[i,k] for k in range(1,6)]
        pairs=[p for p in audit["pairs"] if p["prompt_index"]==i]
        row=dict(prompt_index=i,record_id=rs[0]["record_id"],
            loss_nonincreasing=all(p["loss_nonincreasing"] for p in pairs),
            region_pairs_checked=sum(p["regions_comparable"] for p in pairs),
            region_pairs_unknown=sum(not p["regions_comparable"] for p in pairs),
            region_subset_violations=sum(p["winning_region_inclusion"] is False for p in pairs),
            saturated_pairs=sum(p["saturated"] for p in pairs))
        for r in rs:
            for field in ("g","loss","x_win_size","outcome","abstention_reason"):
                row[f"{field}_k{r['k']}"]=r.get(field)
        prompt_rows.append(row)
    with (output/"prompt-monotonicity.csv").open("w",newline="") as stream:
        writer=csv.DictWriter(stream,fieldnames=list(prompt_rows[0]));writer.writeheader();writer.writerows(prompt_rows)
    known=[r for r in cases if r["region_known"]]
    # Ensure a deadline did not leave an exact region unreported in the table.
    unpublished_regions=[r["job_id"] for r in cases if not r["region_known"] and
                         (output/"regions"/(r["job_id"]+".json")).exists()]
    assert not unpublished_regions,unpublished_regions
    for r in cases:
        assert r["g"]==int(r["issued"])
        assert r["loss"]==int(r["issued"] and not r["implies_reference"])
        if r["region_known"]:
            assert digest(output/r["region_file"])==r["region_sha256"]
            data=read(output/r["region_file"])
            assert len(set(data["x_win"]))==data["size"]==r["x_win_size"]
        else:assert r["x_win_size"] is None and not r["issued"]
        if r["issued"]:
            assert digest(output/r["controller_file"])==r["controller_sha256"]
    best=min(s["crc_corrected_risk"] for s in summary["summary"])
    result=dict(complete=True,cases=200,admissible_initial_states=read(output/"admissible-initial-states.json")["size"],
        original_fixed_start_results_preserved=True,known_region_cases=len(known),unknown_region_cases=200-len(known),
        known_region_size_min=min(r["x_win_size"] for r in known),known_region_size_max=max(r["x_win_size"] for r in known),
        loss_adjacent_pairs=160,loss_violations=audit["loss_monotonicity_violations"],
        region_pairs_checked=audit["region_pairs_comparable"],region_pairs_unknown=audit["region_pairs_unknown"],
        region_inclusion_violations=[],saturated_pairs=audit["saturated_pairs"],
        saturated_pairs_with_known_regions=audit["saturated_pairs_with_known_regions"],saturation_violations=[],
        minimum_crc_corrected_risk=best,
        best_k=[s["k"] for s in summary["summary"] if s["crc_corrected_risk"]==best],
        abstention_causes=dict(collections.Counter(r["abstention_reason"] for r in cases if not r["issued"])))
    (output/"final-audit.json").write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps(result,indent=2))


if __name__=="__main__":main()
