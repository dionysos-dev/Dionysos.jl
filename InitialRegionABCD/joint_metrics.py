"""Merge issuance with logical implication and audit exact initial regions."""
import collections
import csv
import hashlib
import json
from pathlib import Path

OUTCOMES=("controller_issued","logically_inconsistent","plant_unrealizable","computational_failure","no_valid_candidate")
_REGIONS={}
_INCLUSIONS={}


def save(path,value):
    path=Path(path);temporary=path.with_name(path.name+".tmp")
    temporary.write_text(json.dumps(value,indent=2,ensure_ascii=False)+"\n")
    temporary.replace(path)


def exact_region(row,output):
    relative=row.get("region_file")
    if not relative:
        return None
    path=Path(output)/relative
    if path not in _REGIONS:
        raw=path.read_bytes()
        assert hashlib.sha256(raw).hexdigest()==row["region_sha256"]
        data=json.loads(raw)
        assert data["job_id"]==row["job_id"] and data["formula"]==row["formula"]
        states=frozenset(data["x_win"])
        assert len(states)==data["size"]==row["x_win_size"]
        _REGIONS[path]=states
    return _REGIONS[path]


def monotonicity(results,output):
    index={(r["prompt_index"],r["k"]):r for r in results}
    loss_increases=[];region_increases=[];saturation_violations=[];pairs=[]
    known=0;unknown=0;saturated=0;saturated_known=0
    criterion_violations=[]
    for r in results:
        region=exact_region(r,output)
        if region is not None and bool(region)!=bool(r["issued"]):
            criterion_violations.append(dict(prompt_index=r["prompt_index"],k=r["k"],
                outcome=r["outcome"],subreason=r["subreason"],x_win_size=len(region)))
    for i in range(1,41):
        for k in range(1,5):
            if (i,k) not in index or (i,k+1) not in index:
                continue
            low,high=index[i,k],index[i,k+1]
            details=dict(prompt_index=i,record_id=low["record_id"],lower_k=k,upper_k=k+1,
                lower_loss=low["loss"],upper_loss=high["loss"],
                lower_outcome=low["outcome"],upper_outcome=high["outcome"],
                lower_reason=low.get("subreason"),upper_reason=high.get("subreason"))
            if high["loss"]>low["loss"]:
                loss_increases.append(details)
            a,b=exact_region(low,output),exact_region(high,output)
            check=dict(details,loss_nonincreasing=high["loss"]<=low["loss"],
                saturated=low["formula"]==high["formula"],regions_comparable=a is not None and b is not None,
                winning_region_inclusion=None)
            if a is not None and b is not None:
                known+=1
                key=(str(Path(output).resolve()),low["job_id"],high["job_id"])
                if key not in _INCLUSIONS:
                    _INCLUSIONS[key]=(b.issubset(a),sorted(b-a))
                inclusion,extra=_INCLUSIONS[key]
                check["winning_region_inclusion"]=inclusion
                if not inclusion:
                    region_increases.append(dict(details,added_states=extra,added_count=len(extra)))
            else:
                unknown+=1
            if check["saturated"]:
                saturated+=1
                same=low["job_id"]==high["job_id"] and low["outcome"]==high["outcome"] and low["subreason"]==high["subreason"]
                if a is not None and b is not None:
                    saturated_known+=1;same=same and a==b
                if not same:
                    saturation_violations.append(details)
            pairs.append(check)
    return dict(adjacent_pairs_checked=len(pairs),expected_pairs=160,
        loss_monotonicity_violations=loss_increases,winning_region_inclusion_violations=region_increases,
        region_pairs_comparable=known,region_pairs_unknown=unknown,
        saturated_pairs=saturated,saturated_pairs_with_known_regions=saturated_known,
        saturation_violations=saturation_violations,issuance_criterion_violations=criterion_violations,
        pairs=pairs)


def summarize(config,results,output):
    output=Path(output)
    rows=[]
    for k in range(1,6):
        selected=[r for r in results if r["k"]==k]
        issued=sum(r["issued"] for r in selected)
        loss=sum(r["loss"] for r in selected)
        counts={o:sum(r["outcome"]==o for r in selected) for o in OUTCOMES}
        assert sum(counts.values())==len(selected)
        rows.append(dict(k=k,completed=len(selected),planned=40,pending=40-len(selected),counts=counts,
            issued_and_guaranteed=issued-loss,issued_but_not_guaranteed=loss,abstained=len(selected)-issued,
            issuance_rate=issued/40,abstention_rate=(len(selected)-issued)/40,
            joint_risk_rate=loss/40,crc_corrected_risk=(loss+1)/41,
            conditional_diagnostic=loss/issued if issued else None,
            issued_by_backend=dict(collections.Counter(r.get("backend") for r in selected if r["issued"])),
            failure_subreasons=dict(collections.Counter(r["subreason"] for r in selected if r["outcome"]=="computational_failure")),
            initial_region_sizes=[dict(prompt_index=r["prompt_index"],size=r.get("x_win_size"),
                region_known=r.get("region_known",False)) for r in selected]))
    audit=monotonicity(results,output)
    summary=dict(schema_version=3,complete=len(results)==200,model="deepseek-v4-pro",temperature=.7,
        initial_domain=config["initial_domain"],issuance_rule=config["issuance_rule"],
        frozen_config_sha256=hashlib.sha256(json.dumps(config,sort_keys=True,ensure_ascii=False).encode()).hexdigest(),
        denominator=40,crc_denominator=41,conditional_diagnostic_definition="issued_but_not_guaranteed / issued",
        unknown_region_policy="null for computational abstention; never replace unknown by the empty set",
        summary=rows,monotonicity=audit,results=results)
    save(output/"summary.json",summary)
    save(output/"merged-cases.json",results)
    save(output/"monotonicity.json",audit)
    fields=["prompt_index","record_id","k","formula","reference_ltl","implies_reference",
            "conjunction_satisfiable","g","loss","outcome","abstention_reason","backend",
            "x_win_size","region_known","region_file","job_id"]
    with (output/"merged-cases.csv").open("w",newline="") as stream:
        writer=csv.DictWriter(stream,fieldnames=fields,extrasaction="ignore");writer.writeheader();writer.writerows(results)
    columns=["k","completed","issued_and_guaranteed","issued_but_not_guaranteed","abstained",
             "joint_risk_rate","crc_corrected_risk","conditional_diagnostic"]
    with (output/"joint-risk.csv").open("w",newline="") as stream:
        writer=csv.DictWriter(stream,fieldnames=columns,extrasaction="ignore");writer.writeheader();writer.writerows(rows)
    text=["# All-admissible-initial-state development synthesis", "",
        "DeepSeek V4 Pro; T=0.7; R=10; 40 development prompts; existing model outputs only.","",
        "The cached plant's fixed start is replaced by all admissible plant states. A controller is issued only after a nonempty label-consistent initial winning region is computed, universally verified, and exported. Unknown regions remain null on computational failure.","",
        "| k | Issued and guaranteed | Issued but not guaranteed | Abstained | Joint-risk rate | CRC-corrected risk | Conditional diagnostic |",
        "|---|---:|---:|---:|---:|---:|---:|"]
    for r in rows:
        conditional="undefined" if r["conditional_diagnostic"] is None else f'{100*r["conditional_diagnostic"]:.2f}%'
        text.append(f'| {r["k"]} | {r["issued_and_guaranteed"]} | {r["issued_but_not_guaranteed"]} | {r["abstained"]} | {100*r["joint_risk_rate"]:.2f}% | {100*r["crc_corrected_risk"]:.2f}% | {conditional} |')
    text += ["",f'Complete: {summary["complete"]}. Conditional diagnostic = losses / issued. Joint risk = losses / 40. CRC-corrected risk = (losses+1) / 41.',"",
        "## Synthesis outcomes","","| k | Issued | Inconsistent | Plant-unrealizable | Computational failure | No candidate |","|---|---:|---:|---:|---:|---:|"]
    for r in rows:text.append(f'| {r["k"]} | '+" | ".join(str(r["counts"][o]) for o in OUTCOMES)+" |")
    text += ["","## Monotonicity","",
        f'Loss increases: {len(audit["loss_monotonicity_violations"])} / {audit["adjacent_pairs_checked"]} adjacent pairs.',
        f'Winning-region inclusion violations: {len(audit["winning_region_inclusion_violations"])} among {audit["region_pairs_comparable"]} comparable pairs; {audit["region_pairs_unknown"]} pairs have unknown regions.',
        f'Saturated pairs: {audit["saturated_pairs"]}; mismatched outcomes/regions: {len(audit["saturation_violations"])}.',
        f'Known nonempty-region / issuance mismatches: {len(audit["issuance_criterion_violations"])}.',"",
        "The numerical CRC correction is reported as requested on development data. This is not a held-out calibration guarantee; operational abstention from timeouts/resource limits can also break loss monotonicity. See monotonicity.json for every prompt/pair."]
    (output/"joint-risk.md").write_text("\n".join(text)+"\n")
    (output/"comparison.md").write_text("\n".join(text)+"\n")
    return summary
