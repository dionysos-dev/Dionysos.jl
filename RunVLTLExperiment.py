#!/usr/bin/env python3
"""Collect VLTL, analyze saved responses, and produce the two-benchmark figure."""
import fcntl
import json
import os
from pathlib import Path
import subprocess
import sys

root=Path(__file__).resolve().parent
output=root/'Results/vltl-experiment'
output.mkdir(parents=True,exist_ok=True)


def status(stage,**values):
    path=output/'progress.json'
    temporary=path.with_suffix('.json.tmp')
    temporary.write_text(json.dumps(dict(stage=stage,**values),indent=2)+'\n')
    temporary.replace(path)


def main():
    with (output/'experiment.lock').open('w') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        env=os.environ.copy()
        env['PATH']='/Users/abayat/.juliaup/bin:/opt/homebrew/opt/spot/bin:'+env.get('PATH','')
        try:
            status('collecting_and_analyzing',benchmark='VLTL',planned_queries=420)
            subprocess.run(['/Users/abayat/.juliaup/bin/julia','--startup-file=no',
                '--project='+str(root),'--threads=1',str(root/'VLTLEvaluate.jl'),'--run'],
                cwd=root,env=env,check=True)
            status('plotting',new_controller_synthesis_jobs=0)
            python=Path(os.environ.get('VLTL_PLOT_PYTHON','/tmp/vltl-plot-20261109/bin/python'))
            subprocess.run([str(python if python.exists() else sys.executable),
                str(root/'PlotSemanticInclusionComparison.py')],cwd=root,env=env,check=True)
            status('complete',complete=True,figure_directory=str(root/'Figures/semantic-inclusion-r10-with-VLTL'))
        except BaseException as error:
            status('stopped',complete=False,error=type(error).__name__+': '+str(error))
            raise


if __name__=='__main__':
    main()
