"""Apply cumulative issuance to saved synthesis outputs; no solver or model calls."""
import argparse
import collections
import csv
import hashlib
import json
from pathlib import Path

from policy import apply_prefix_issuance, audit_prefix_issuance


def read(path):return json.loads(path.read_text())


def digest(path):
    checksum=hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda:stream.read(1024*1024),b""):checksum.update(chunk)
    return checksum.hexdigest()


def save(path,value):
    temporary=path.with_name(path.name+".tmp")
    temporary.write_text(json.dumps(value,indent=2,ensure_ascii=False)+"\n")
    temporary.replace(path)


def analyze(source,output):
    source,output=Path(source).resolve(),Path(output).resolve()
    if source==output:raise ValueError("Use a new output directory to preserve raw synthesis records")
    base=read(source/"summary.json")
    frozen=read(source/"frozen-development.json")
    semantic=read(source/"semantic-audit.json")
    assert base["complete"] and base["schema_version"]==3
    assert frozen["initial_domain"]=="all_admissible_plant_states"
    n=frozen["prompt_count"];budgets=frozen["budgets"]
    assert n==40 and budgets==list(range(1,6)) and len(base["results"])==200
    assert semantic["complete"] and not semantic["reference_mismatches"] and not semantic["prefix_violations"]
    semantics={(r["prompt_index"],r["k"]):r["implies_reference"] for r in semantic["reference_results"]}
    # Verify each artifact once, including controllers later withheld by policy.
    verified=set()
    for row in base["results"]:
        assert row["implies_reference"]==semantics[row["prompt_index"],row["k"]]
        assert row["g"]==int(row["issued"])
        for file_key,hash_key in (("controller_file","controller_sha256"),("region_file","region_sha256")):
            relative=row.get(file_key)
            if relative and relative not in verified:
                assert digest(source/relative)==row[hash_key],relative
                verified.add(relative)
    result=apply_prefix_issuance(base["results"])
    audit=audit_prefix_issuance(result)
    assert audit["adjacent_pairs"]==160
    assert not audit["issuance_violations"] and not audit["implication_violations"] and not audit["loss_violations"]
    # The new policy can withhold a controller with a nonempty winning region.
    # X_win is deliberately left unchanged; it belongs to raw synthesis.
    raw_region_mismatches=[]
    for row in result:
        if row["region_known"] and bool(row["x_win_size"])!=bool(row["tilde_g"]):
            raw_region_mismatches.append(dict(prompt_index=row["prompt_index"],k=row["k"]))
    assert not raw_region_mismatches
    index={(r["prompt_index"],r["k"]):r for r in result}
    saturation=[]
    for p in sorted({r["prompt_index"] for r in result}):
        for k in budgets[:-1]:
            a,b=index[p,k],index[p,k+1]
            if a["formula"]==b["formula"]:
                equal=all(a[key]==b[key] for key in ("job_id","tilde_g","g","loss","outcome","region_file","x_win_size"))
                saturation.append(dict(prompt_index=p,lower_k=k,upper_k=k+1,identical=equal))
    assert all(r["identical"] for r in saturation)
    audit.update(saturated_pairs=len(saturation),saturation_violations=[r for r in saturation if not r["identical"]],
        synthesis_region_audit=base["monotonicity"],raw_success_region_violations=raw_region_mismatches)
    stats=[]
    for k in budgets:
        selected=[r for r in result if r["k"]==k]
        assert len(selected)==n
        issued=sum(r["g"] for r in selected);loss=sum(r["loss"] for r in selected)
        stats.append(dict(k=k,raw_synthesis_successes=sum(r["tilde_g"] for r in selected),
            issued=issued,prefix_blocked=sum(r["prefix_blocked"] for r in selected),
            issued_and_guaranteed=issued-loss,issued_but_not_guaranteed=loss,abstained=n-issued,
            joint_risk_rate=loss/n,crc_corrected_risk=(loss+1)/(n+1),
            conditional_diagnostic=loss/issued if issued else None,
            abstention_causes=dict(collections.Counter(r["abstention_reason"] for r in selected if not r["issued"])),
            raw_synthesis_outcomes=dict(collections.Counter(r["synthesis_outcome"] for r in selected))))
    output.mkdir(parents=True,exist_ok=True)
    config=dict(schema_version=4,purpose="development_only",model=frozen["primary_model"],
        temperature=frozen["temperature"],repetitions=frozen["repetitions"],prompt_count=n,budgets=budgets,
        source_experiment=str(source),artifact_base_directory=str(source),
        source_summary_sha256=digest(source/"summary.json"),
        source_frozen_config_sha256=digest(source/"frozen-development.json"),
        source_semantic_audit_sha256=digest(source/"semantic-audit.json"),
        initial_domain=frozen["initial_domain"],raw_synthesis_rule=frozen["issuance_rule"],
        issuance_rule="g_k = product_{j=1..k} tilde_g_j",
        loss_rule="loss_k = g_k * indicator(Phi_k does not imply reference)",
        no_new_synthesis=True,no_model_calls=True,
        code_sha256={name:digest(Path(__file__).with_name(name)) for name in ("policy.py","analyze.py")})
    config_path=output/"analysis-config.json"
    if config_path.exists() and read(config_path)!=config:
        raise ValueError("Analysis configuration changed; choose a new output directory")
    save(config_path,config)
    blocked=[dict(prompt_index=r["prompt_index"],record_id=r["record_id"],k=r["k"],
        first_failure_k=r["prefix_failure_k"],first_failure_reason=r["prefix_failure_reason"],
        implies_reference=r["implies_reference"],x_win_size=r["x_win_size"],
        synthesized_controller_file=r["synthesized_controller_file"]) for r in result if r["prefix_blocked"]]
    summary=dict(schema_version=4,complete=True,issuance_rule=config["issuance_rule"],
        raw_synthesis_indicator="tilde_g",source_experiment=str(source),artifact_base_directory=str(source),
        denominator=n,crc_denominator=n+1,summary=stats,results=result,
        monotonicity=audit,newly_withheld_cases=blocked,
        minimum_crc_corrected_risk=min(r["crc_corrected_risk"] for r in stats))
    save(output/"summary.json",summary)
    save(output/"merged-cases.json",result)
    save(output/"monotonicity.json",audit)
    save(output/"newly-withheld-cases.json",blocked)
    fields=["prompt_index","record_id","k","tilde_g","g","issued","loss","implies_reference",
        "prefix_blocked","prefix_failure_k","prefix_failure_outcome","prefix_failure_reason",
        "outcome","abstention_reason","synthesis_outcome","synthesis_abstention_reason",
        "x_win_size","region_known","issued_region_size","formula","reference_ltl",
        "controller_file","synthesized_controller_file","region_file","job_id"]
    with (output/"merged-cases.csv").open("w",newline="") as stream:
        writer=csv.DictWriter(stream,fieldnames=fields,extrasaction="ignore");writer.writeheader();writer.writerows(result)
    fields=["k","raw_synthesis_successes","issued","prefix_blocked","issued_and_guaranteed",
        "issued_but_not_guaranteed","abstained","joint_risk_rate","crc_corrected_risk","conditional_diagnostic"]
    with (output/"joint-risk.csv").open("w",newline="") as stream:
        writer=csv.DictWriter(stream,fieldnames=fields,extrasaction="ignore");writer.writeheader();writer.writerows(stats)
    text=["# Cumulative issuance: saved DeepSeek development results","",
        "g_k = product_{j=1..k} tilde_g_j. Raw synthesis success and its winning region are retained separately. All artifact paths are relative to the source experiment:",str(source),"",
        "| k | Issued and guaranteed | Issued but not guaranteed | Abstained | Joint-risk rate | CRC-corrected risk | Conditional diagnostic |",
        "|---|---:|---:|---:|---:|---:|---:|"]
    for r in stats:
        conditional="undefined" if r["conditional_diagnostic"] is None else f'{100*r["conditional_diagnostic"]:.2f}%'
        text.append(f'| {r["k"]} | {r["issued_and_guaranteed"]} | {r["issued_but_not_guaranteed"]} | {r["abstained"]} | {100*r["joint_risk_rate"]:.2f}% | {100*r["crc_corrected_risk"]:.2f}% | {conditional} |')
    text += ["",f'Newly withheld successful syntheses: {len(blocked)}. Issuance and loss increases: {len(audit["issuance_violations"])} and {len(audit["loss_violations"])} over {audit["adjacent_pairs"]} adjacent pairs.',
        "The mathematical X_win regions and raw synthesis outcomes are unchanged. A blocked row can have a verified controller and nonempty X_win while g=0. Such an artifact is recorded as synthesized_controller_file and is withheld from controller_file.",
        "Joint risk divides by 40; the requested correction is (losses+1)/41; the conditional diagnostic divides by cumulative issuance. These are development-set diagnostics. No synthesis or model calls were made."]
    (output/"joint-risk.md").write_text("\n".join(text)+"\n")
    print(json.dumps(dict(rows=len(result),newly_withheld=blocked,summary=stats,
        issuance_violations=len(audit["issuance_violations"]),loss_violations=len(audit["loss_violations"])),indent=2))
    return summary


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source",type=Path,required=True)
    parser.add_argument("--output",type=Path,required=True)
    args=parser.parse_args()
    analyze(args.source,args.output)


if __name__=="__main__":main()
