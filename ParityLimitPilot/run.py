#!/usr/bin/env python3
"""Bounded higher-node-limit parity pilot; no algorithm changes or LLM calls."""
import argparse
import collections
import copy
import hashlib
import json
import os
from pathlib import Path
import re
import signal
import subprocess
import time


def read(path):return json.loads(Path(path).read_text())


def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(1048576),b''):h.update(block)
    return h.hexdigest()


def save(path,data):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    tmp=path.with_name(path.name+'.tmp')
    tmp.write_text(json.dumps(data,indent=2,ensure_ascii=False)+'\n');tmp.replace(path)


def stop(process):
    if process.poll() is None:
        os.killpg(process.pid,signal.SIGTERM)
        try:process.wait(timeout=5)
        except subprocess.TimeoutExpired:
            os.killpg(process.pid,signal.SIGKILL);process.wait()


def group_rss_kib(pgid):
    out=subprocess.run(['ps','-axo','pgid=,rss='],capture_output=True,text=True,check=True).stdout
    total=0
    for line in out.splitlines():
        fields=line.split()
        if len(fields)==2 and int(fields[0])==pgid:total+=int(fields[1])
    return total


def memory_free_percent():
    process=subprocess.run(['memory_pressure','-Q'],capture_output=True,text=True,timeout=10)
    match=re.search(r'System-wide memory free percentage:\s*(\d+)%',process.stdout)
    if process.returncode or not match:raise RuntimeError('Cannot monitor system memory pressure')
    return int(match[1])


def run_bounded(command,folder,remaining,max_rss_kib,min_free):
    started=time.monotonic();peak=0;reason=None;last_pressure=0
    if remaining<=0:return dict(exit_code=None,guard='timeout',seconds=0,peak_rss_kib=0)
    with (folder/'run.log').open('a') as log:
        process=subprocess.Popen(command,stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
        save(folder/'active-process.json',dict(pid=process.pid,command=command))
        try:
            while process.poll() is None:
                elapsed=time.monotonic()-started
                peak=max(peak,group_rss_kib(process.pid))
                if peak>max_rss_kib:reason='process_memory_guard';break
                if elapsed>remaining:reason='timeout';break
                if time.monotonic()-last_pressure>2:
                    if memory_free_percent()<min_free:reason='system_memory_guard';break
                    last_pressure=time.monotonic()
                time.sleep(0.15)
        finally:
            stop(process)
            save(folder/'active-process.json',dict(pid=process.pid,stopped=True))
    return dict(exit_code=process.returncode,guard=reason,seconds=time.monotonic()-started,peak_rss_kib=peak)


def classify_bridge_failure(log):
    if 'maximum nodes exceeded' in log:return 'node_limit'
    if 'maximum edges exceeded' in log:return 'edge_limit'
    return 'native_bridge_failure'


def select_cases(rows,count):
    cases=[]
    for row in rows:
        result=row['synthesis']
        if result.get('subreason')!='resource_limit':continue
        hoa=Path(result['artifact_base_directory'])/'parity-jobs'/result['job_id']/'specification.hoa'
        match=re.search(r'^States:\s*(\d+)',hoa.read_text(),re.M)
        if not match:raise ValueError('Missing DPA state count')
        cases.append(dict(class_id=row['group']['class_id'],formula=result['formula'],
            source_job_id=result['job_id'],hoa_path=str(hoa),hoa_sha256=sha(hoa),
            automaton_states=int(match[1]),record_ids=row['group']['record_ids'],
            plant_binary=str(Path(result['artifact_base_directory'])/'plant.parity.bin')))
    cases.sort(key=lambda r:(r['automaton_states'],len(r['formula']),r['formula'],r['class_id']))
    return cases[:count]


def run(root,output,count=6,max_nodes=4_000_000,max_edges=40_000_000,timeout=600,
        memory_gib=3,min_free=15):
    import fcntl
    root,output=Path(root).resolve(),Path(output).resolve()
    output.mkdir(parents=True,exist_ok=True)
    with (output/'experiment.lock').open('w') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        source=root/'Results/ground-truth-controller-audit-timeout600-v2'
        summary=read(source/'summary.json')
        if not summary['complete']:raise ValueError('Source audit is incomplete')
        old=read(source/'frozen-config.json')['config']
        for path,checksum in old['solver']['code_sha256'].items():
            if sha(root/path)!=checksum:raise ValueError('Production algorithm changed: '+path)
        if sha(old['plant_cache'])!=old['plant_sha256']:raise ValueError('Plant changed')
        cases=select_cases(read(source/'reference-results.json'),count)
        if len(cases)!=count:raise ValueError('Not enough eligible pilot cases')
        bridge=root/'InitialRegionABCD/build/spot_parity_bridge'
        plant_hashes={sha(c['plant_binary']) for c in cases}
        if len(plant_hashes)!=1:raise ValueError('Native plant caches differ')
        for c in cases:
            c.update(plant_cache=old['plant_cache'],plant_sha256=old['plant_sha256'])
        config=dict(source_reference_results_sha256=sha(source/'reference-results.json'),
            production_solver=old['solver'],native_bridge_sha256=sha(bridge),
            pilot_runner_sha256=sha(__file__),pilot_verifier_sha256=sha(root/'ParityLimitPilot/verify.jl'),
            max_nodes=max_nodes,max_edges=max_edges,per_case_timeout_seconds=timeout,
            process_group_rss_guard_gib=memory_gib,system_memory_free_floor_percent=min_free,
            workers=1,verification_threads=1,max_simulation_steps=1200,
            verification_horizon='infinite',algorithm_change=False,
            memory_guard='polled process-group resident memory and system pressure; not an OS hard allocation quota',
            selection='six smallest existing failed DPAs, then shortest formula; exploratory not a representative sample',
            native_plant_sha256=next(iter(plant_hashes)),cases=cases)
        frozen=output/'frozen-config.json'
        if frozen.exists() and read(frozen)!=config:raise ValueError('Frozen pilot changed; use another directory')
        save(frozen,config)
        results=[]
        report=dict(complete=False,case_count=count,new_api_calls=0,algorithm_change=False,config=str(frozen))
        save(output/'summary.json',report)
        try:
            for index,case in enumerate(cases,1):
                folder=output/'cases'/case['class_id'];folder.mkdir(parents=True,exist_ok=True)
                result_path=folder/'result.json'
                if result_path.exists():result=read(result_path)
                else:
                    if memory_free_percent()<min_free:raise RuntimeError('System memory guard: pilot not started')
                    save(folder/'case.json',case)
                    save(output/'progress.json',dict(stage='native_parity_game',case_index=index,total_cases=count,
                        completed_cases=len(results),formula=case['formula']))
                    print(f"Case {index}/{count}: {case['formula']} (DPA states={case['automaton_states']})",flush=True)
                    start=time.monotonic()
                    command=[str(bridge),case['plant_binary'],case['hoa_path'],str(folder/'spot-result.json'),
                        str(folder/'spot-policy.bin'),'1200',str(max_nodes),str(max_edges)]
                    native=run_bounded(command,folder,timeout,int(memory_gib*1024*1024),min_free)
                    result=dict(class_id=case['class_id'],formula=case['formula'],issued=False,
                        outcome='computational_failure',native_execution=native)
                    if native['guard']:result['subreason']=native['guard']
                    elif native['exit_code']!=0:result['subreason']=classify_bridge_failure((folder/'run.log').read_text())
                    else:
                        info=read(folder/'spot-result.json')
                        result['native_statistics']={k:v for k,v in info.items() if not isinstance(v,list)}
                        if not info['winning']:
                            if info['winning_initial_state_count']!=0:raise ValueError('Inconsistent empty region')
                            result.update(outcome='plant_unrealizable',subreason='empty_admissible_initial_winning_region',x_win_size=0)
                        else:
                            if not info['verified']:raise ValueError('Native universal verification did not pass')
                            save(output/'progress.json',dict(stage='independent_parity_verification',case_index=index,
                                total_cases=count,completed_cases=len(results),formula=case['formula']))
                            command=['/Users/abayat/.juliaup/bin/julia','--startup-file=no','--project='+str(root),'--threads=1',
                                str(root/'ParityLimitPilot/verify.jl'),str(root),str(folder),'1200',str(folder/'case.json')]
                            verification=run_bounded(command,folder,timeout-(time.monotonic()-start),int(memory_gib*1024*1024),min_free)
                            result['verification_execution']=verification
                            if verification['guard']:result['subreason']=verification['guard']
                            elif verification['exit_code']!=0:result['subreason']='independent_verification_or_export_failure'
                            else:result.update(read(folder/'verified-result.json'))
                    result['elapsed_seconds']=time.monotonic()-start
                    save(result_path,result)
                results.append(result)
                report['results']=results
                report['outcomes']=dict(collections.Counter(r['outcome'] for r in results))
                report['failure_reasons']=dict(collections.Counter(r.get('subreason') for r in results if r['outcome']=='computational_failure'))
                save(output/'summary.json',report)
                save(output/'progress.json',dict(stage='case_complete',completed_cases=len(results),total_cases=count,outcomes=report['outcomes']))
                print('COMPLETED '+json.dumps(result),flush=True)
            report['complete']=True
            lines=['# Higher-node-limit parity pilot','',
                f'Same native Spot game construction and solver; {max_nodes:,} nodes (previously 2,000,000), {max_edges:,} edges, {timeout} seconds per case.',
                f'Sequential execution; {memory_gib} GiB monitored RSS guard; system-free floor {min_free}%.',
                'Only independently verified and exported strategies count as controller issuance.','',
                '| Formula | Outcome | Reason | Peak process memory (GiB) | Time (s) |',
                '|---|---|---|---:|---:|']
            for r in results:
                peak=max(r['native_execution']['peak_rss_kib'],r.get('verification_execution',{}).get('peak_rss_kib',0))/(1024*1024)
                lines.append(f"| {r['formula']} | {r['outcome']} | {r.get('subreason','verified controller')} | {peak:.2f} | {r['elapsed_seconds']:.1f} |")
            lines+=['','This deliberately selected small-automaton pilot does not estimate a recovery percentage for all 165 parity failures. Original benchmark reports are unchanged.','']
            (output/'comparison.md').write_text('\n'.join(lines))
        except BaseException as error:
            report['interruption_or_error']=type(error).__name__+': '+str(error)
            raise
        finally:save(output/'summary.json',report)
        save(output/'progress.json',dict(stage='complete',completed_cases=len(results),total_cases=count,outcomes=report['outcomes']))
        print('FINAL '+json.dumps(report),flush=True)
        return report


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root',type=Path,required=True)
    parser.add_argument('--output',type=Path)
    args=parser.parse_args()
    run(args.root,args.output or args.root/'Results/parity-limit-pilot-4m-v1')
