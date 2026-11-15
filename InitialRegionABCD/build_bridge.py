#!/usr/bin/env python3
"""Build the thin bridge against an installed Spot C++ library."""
import argparse
from pathlib import Path
import os
import shutil
import subprocess

def build(output=None,prefix=None):
    source=Path(__file__).with_name("spot_parity_bridge.cpp")
    output=Path(output) if output else source.parent/"build"/"spot_parity_bridge"
    exe=shutil.which("ltl2tgba")
    if not exe:
        raise RuntimeError("Spot's ltl2tgba is required on PATH")
    prefix=Path(prefix or os.environ.get("SPOT_PREFIX") or Path(exe).resolve().parents[1])
    if not (prefix/"include"/"spot"/"twaalgos"/"game.hh").exists():
        raise RuntimeError("Spot development headers missing; supply --prefix or SPOT_PREFIX")
    output.parent.mkdir(parents=True,exist_ok=True)
    command=[os.environ.get("CXX","c++"),"-std=c++17","-O2",f"-I{prefix/'include'}",
        str(source),f"-L{prefix/'lib'}",f"-Wl,-rpath,{prefix/'lib'}","-lspot","-lbddx","-o",str(output)]
    subprocess.run(command,check=True)
    return output.resolve()

if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output",type=Path)
    parser.add_argument("--prefix",type=Path)
    args=parser.parse_args()
    print(build(args.output,args.prefix))
