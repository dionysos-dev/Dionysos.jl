"""Recheck cached implication labels and conjunction-prefix implication with Spot."""
import argparse
import json
from pathlib import Path
import subprocess


def implies(formula,reference):
    result=subprocess.run(["ltlfilt","-f",formula,"--imply="+reference,"--count"],
                          capture_output=True,text=True,timeout=120)
    if result.returncode not in (0,1):raise RuntimeError(result.stderr)
    return result.returncode==0


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("experiment",type=Path)
    parser.add_argument("--output",type=Path,required=True)
    args=parser.parse_args()
    config=json.loads((args.experiment/"frozen-development.json").read_text())
    cases=config["cases"]
    cache={};mismatches=[];reference=[];prefixes=[]
    def check(a,b):
        key=a,b
        if key not in cache:cache[key]=implies(a,b)
        return cache[key]
    for row in cases:
        if row["formula"] is None:continue
        actual=check(row["formula"],row["reference_ltl"])
        result=dict(prompt_index=row["prompt_index"],k=row["k"],implies_reference=actual)
        reference.append(result)
        if actual!=row["implies_reference"]:mismatches.append(result)
    index={(r["prompt_index"],r["k"]):r for r in cases}
    for i in range(1,41):
        for k in range(1,5):
            low,high=index[i,k],index[i,k+1]
            if low["formula"] is None or high["formula"] is None:continue
            actual=check(high["formula"],low["formula"])
            prefixes.append(dict(prompt_index=i,lower_k=k,upper_k=k+1,higher_implies_lower=actual))
    result=dict(complete=True,reference_checks=len(reference),prefix_checks=len(prefixes),
        unique_spot_queries=len(cache),reference_mismatches=mismatches,
        prefix_violations=[r for r in prefixes if not r["higher_implies_lower"]],
        reference_results=reference,prefix_results=prefixes)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2)+"\n")
    print({k:v for k,v in result.items() if k not in ("reference_results","prefix_results")})
    assert not mismatches and not result["prefix_violations"]


if __name__=="__main__":main()
