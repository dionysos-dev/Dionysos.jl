#!/usr/bin/env python3
"""Continue remaining parity-limit cases sequentially, with explicit rest breaks.

Unchanged production bridge and verification; 4M nodes, 60M edges, 600 seconds.
Reuses the six completed pilot decisions and saves every new class separately.
"""
import argparse
import collections
import copy
import hashlib
import importlib.util
import json
from pathlib import Path
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


def validate_result(result):
    if bool(result['issued'])!=(result['outcome']=='controller_issued'):
        raise ValueError('Inconsistent issuance flag')
    if result['issued']:
        if not result.get('region_known') or result.get('x_win_size',0)<=0 or not result.get('certificate',{}).get('satisfied'):
            raise ValueError('Issued strategy lacks verified admissible winning states')
        for file_key,hash_key in [('controller_file','controller_sha256'),('region_file','region_sha256')]:
            if sha(result[file_key])!=result[hash_key]:raise ValueError('Artifact hash mismatch')
    elif result['outcome']=='plant_unrealizable':
        info=result['native_statistics']
        if info['winning'] or info['winning_initial_state_count']!=0:
            raise ValueError('Unrealizable outcome does not have an empty native winning region')


def remaining_cases(cases,pilot_results):
    known={r['class_id']:r for r in pilot_results}
    if len(known)!=len(pilot_results):raise ValueError('Duplicate pilot classes')
    for r in pilot_results:
        if r['outcome'] not in ('controller_issued','plant_unrealizable'):
            raise ValueError('Only resolved pilot decisions may be reused')
    ids={c['class_id'] for c in cases}
    if not set(known)<=ids:raise ValueError('Pilot class is not in the source parity failures')
    return [c for c in cases if c['class_id'] not in known]


def aggregate(source_rows,overrides):
    counts,records=collections.Counter(),collections.Counter()
    for row in source_rows:
        result=overrides.get(row['group']['class_id'],row['synthesis'])
        counts[result['outcome']]+=1
        records[result['outcome']]+=len(row['group']['record_ids'])
    return dict(semantic_classes=len(source_rows),benchmark_records=sum(records.values()),
        class_counts=dict(counts),record_counts=dict(records),
        class_rates={k:v/len(source_rows) for k,v in counts.items()},
        record_rates={k:v/sum(records.values()) for k,v in records.items()})


def cooldown(output,seconds,completed,total,reason):
    deadline=time.monotonic()+seconds
    print(f'COOLDOWN {seconds:g} seconds ({reason}); no synthesis process is active',flush=True)
    while True:
        left=deadline-time.monotonic()
        if left<=0:break
        save(output/'progress.json',dict(stage='cooldown',reason=reason,completed_cases=completed,
            total_cases=total,cooldown_seconds_remaining=round(left),synthesis_process_active=False))
        time.sleep(min(2,left))


def run(root,output,cooldown_seconds=60,long_break_seconds=300,long_break_after=1800):
    import fcntl
    root,output=Path(root).resolve(),Path(output).resolve()
    output.mkdir(parents=True,exist_ok=True)
    with (output/'experiment.lock').open('w') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        helper_path=root/'ParityLimitPilot/run_edge60m.py'
        spec=importlib.util.spec_from_file_location('bounded_parity_pilot',helper_path)
        helper=importlib.util.module_from_spec(spec);spec.loader.exec_module(helper)
        source=root/'Results/ground-truth-controller-audit-timeout600-v2'
        original=read(source/'summary.json')
        if not original['complete']:raise ValueError('Source audit is incomplete')
        source_rows=read(source/'reference-results.json')
        old=read(source/'frozen-config.json')['config']
        for name,checksum in old['solver']['code_sha256'].items():
            if sha(root/name)!=checksum:raise ValueError('Production algorithm changed: '+name)
        if sha(old['plant_cache'])!=old['plant_sha256']:raise ValueError('Plant changed')
        pilot_path=root/'Results/parity-limit-pilot-4m-60m-v2/combined-pilot-summary.json'
        pilot=read(pilot_path)
        if not pilot['complete']:raise ValueError('Pilot is incomplete')
        pilot_results=pilot['results']
        for result in pilot_results:validate_result(result)
        all_cases=helper.select_cases(source_rows,165)
        cases=remaining_cases(all_cases,pilot_results)
        if len(all_cases)!=165 or len(cases)!=159:raise ValueError('Unexpected parity audit coverage')
        plant_hashes={p:sha(p) for p in {c['plant_binary'] for c in cases}}
        if len(set(plant_hashes.values()))!=1:raise ValueError('Native plant caches differ')
        for c in cases:c.update(plant_cache=old['plant_cache'],plant_sha256=old['plant_sha256'])
        bridge=root/'InitialRegionABCD/build/spot_parity_bridge'
        config=dict(source_reference_sha256=sha(source/'reference-results.json'),pilot_sha256=sha(pilot_path),
            production_solver=old['solver'],runner_sha256=sha(__file__),helper_sha256=sha(helper_path),
            verifier_sha256=sha(root/'ParityLimitPilot/verify.jl'),native_plant_hashes=plant_hashes,
            cases=cases,max_nodes=4000000,max_edges=60000000,timeout_seconds=600,
            workers=1,verification_threads=1,process_group_rss_guard_gib=3,
            system_memory_free_floor_percent=15,simulation_steps=1200,verification_horizon='infinite',
            cooldown_seconds=cooldown_seconds,long_break_seconds=long_break_seconds,
            long_break_after_processing_seconds=long_break_after,temperature_measured=False,
            algorithm_change=False,pilot_cases_reused=len(pilot_results),new_cases=len(cases))
        frozen=output/'frozen-config.json'
        if frozen.exists() and read(frozen)!=config:raise ValueError('Frozen continuation changed')
        save(frozen,config)
        if (output/'summary.json').exists() and read(output/'summary.json').get('complete'):
            print('Already complete',flush=True);return read(output/'summary.json')
        overrides={r['class_id']:r for r in pilot_results}
        completed=[]
        for c in cases:
            path=output/'cases'/c['class_id']/'result.json'
            if path.exists():
                result=read(path);validate_result(result)
                if result['class_id']!=c['class_id'] or result['formula']!=c['formula']:
                    raise ValueError('Continuation checkpoint identity mismatch')
                completed.append(result);overrides[c['class_id']]=result
        report=dict(complete=False,new_api_calls=0,algorithm_change=False,pilot_cases_reused=len(pilot_results),
            new_cases_total=len(cases),new_cases_completed=len(completed),
            new_case_outcomes=dict(collections.Counter(r['outcome'] for r in completed)),
            known_results=aggregate(source_rows,overrides),cooldown_seconds=cooldown_seconds)
        save(output/'summary.json',report)
        save(output/'progress.json',dict(stage='initializing',completed_cases=len(completed),total_cases=len(cases)))
        (output/'README.md').write_text('''# Remaining parity-limit cases, with cooldown breaks

159 remaining ground-truth parity cases, plus six resolved pilot decisions
reused with provenance. Same production Spot bridge and independent verifier.
Four-million-node and sixty-million-edge caps; 600 seconds per new case,
including independent verification and export. Cases run one at a time.
Julia verification uses one thread. No LLM/model requests.

Monitored process-group RSS guard: 3 GiB. System-wide free-memory floor: 15%.
These are memory guards, not temperature sensors or OS allocation quotas.
Every new completed case is followed by 60 seconds with no synthesis process.
After roughly 30 minutes of case-processing time, use a 5-minute break instead.
The final case needs no cooldown because all processing has finished.
Cool-downs reduce sustained load; CPU temperature is not measured.

Pause by sending SIGINT to this runner. Its active child process group is
stopped in finally; interrupted cases are not saved as failed decisions.
Resume with the same command; completed case/result.json files are reused.
progress.json distinguishes active synthesis, verification, cooldown and
waiting for memory. summary.json complete=true marks final report completion.
Original audits and pilot files remain unchanged. Controller artifacts retain
the frozen plant-cache path and hash.
''')
        active_since_break=0
        try:
            for index,case in enumerate(cases,1):
                folder=output/'cases'/case['class_id'];folder.mkdir(parents=True,exist_ok=True)
                result_path=folder/'result.json'
                if result_path.exists():continue
                while helper.memory_free_percent()<15:
                    cooldown(output,60,len(completed),len(cases),'waiting for system memory')
                save(folder/'case.json',case)
                save(output/'progress.json',dict(stage='native_parity_game',case_index=index,
                    completed_cases=len(completed),total_cases=len(cases),formula=case['formula'],synthesis_process_active=True))
                print(f"CASE {index}/{len(cases)}: {case['formula']} (DPA states={case['automaton_states']})",flush=True)
                start=time.monotonic()
                command=[str(bridge),case['plant_binary'],case['hoa_path'],str(folder/'spot-result.json'),
                    str(folder/'spot-policy.bin'),'1200','4000000','60000000']
                native=helper.run_bounded(command,folder,600,3*1024*1024,15)
                result=dict(class_id=case['class_id'],formula=case['formula'],issued=False,
                    outcome='computational_failure',native_execution=native,artifact_base_directory=str(folder))
                if native['guard']:result['subreason']=native['guard']
                elif native['exit_code']!=0:result['subreason']=helper.classify_bridge_failure((folder/'run.log').read_text())
                else:
                    info=read(folder/'spot-result.json')
                    result['native_statistics']={k:v for k,v in info.items() if not isinstance(v,list)}
                    if not info['winning']:
                        if info['winning_initial_state_count']!=0:raise ValueError('Inconsistent empty region')
                        result.update(outcome='plant_unrealizable',subreason='empty_admissible_initial_winning_region',x_win_size=0)
                    else:
                        if not info['verified']:raise ValueError('Native strategy is unverified')
                        save(output/'progress.json',dict(stage='independent_parity_verification',case_index=index,
                            completed_cases=len(completed),total_cases=len(cases),formula=case['formula'],synthesis_process_active=True))
                        command=['/Users/abayat/.juliaup/bin/julia','--startup-file=no','--project='+str(root),'--threads=1',
                            str(root/'ParityLimitPilot/verify.jl'),str(root),str(folder),'1200',str(folder/'case.json')]
                        verification=helper.run_bounded(command,folder,600-(time.monotonic()-start),3*1024*1024,15)
                        result['verification_execution']=verification
                        if verification['guard']:result['subreason']=verification['guard']
                        elif verification['exit_code']!=0:result['subreason']='independent_verification_or_export_failure'
                        else:result.update(read(folder/'verified-result.json'))
                result['elapsed_seconds']=time.monotonic()-start
                validate_result(result);save(result_path,result)
                completed.append(result);overrides[case['class_id']]=result
                active_since_break+=result['elapsed_seconds']
                report.update(new_cases_completed=len(completed),
                    new_case_outcomes=dict(collections.Counter(r['outcome'] for r in completed)),
                    failure_reasons=dict(collections.Counter(r.get('subreason') for r in completed if r['outcome']=='computational_failure')),
                    known_results=aggregate(source_rows,overrides))
                save(output/'summary.json',report)
                save(output/'progress.json',dict(stage='case_complete',completed_cases=len(completed),total_cases=len(cases),outcomes=report['new_case_outcomes'],synthesis_process_active=False))
                print('COMPLETED '+json.dumps(result),flush=True)
                if len(completed)<len(cases):
                    if active_since_break>=long_break_after:
                        cooldown(output,long_break_seconds,len(completed),len(cases),'periodic longer break')
                        active_since_break=0
                    else:cooldown(output,cooldown_seconds,len(completed),len(cases),'between cases')
            merged=copy.deepcopy(source_rows)
            for row in merged:
                replacement=overrides.get(row['group']['class_id'])
                if replacement is not None:
                    row['previous_synthesis']=row['synthesis']
                    row['synthesis']=replacement
            save(output/'reference-results.json',merged)
            report.update(complete=True,known_results=aggregate(source_rows,overrides))
            lines=['# Ground-truth audit after larger parity limits','',
                'Production algorithms unchanged. Six resolved pilot cases reused; all 159 remaining parity cases attempted sequentially with cooldown breaks.','',
                '| Outcome | Semantic classes | % of classes | Records | % of records |',
                '|---|---:|---:|---:|---:|']
            stats=report['known_results']
            for outcome in ('controller_issued','plant_unrealizable','logically_inconsistent','computational_failure'):
                n=stats['class_counts'].get(outcome,0);nr=stats['record_counts'].get(outcome,0)
                lines.append(f"| {outcome} | {n} | {100*n/stats['semantic_classes']:.2f}% | {nr} | {100*nr/stats['benchmark_records']:.2f}% |")
            lines+=['','Computational failures remain unknown controller existence. All successful parity controllers passed both native and separate Julia infinite-horizon verification. Original audit reports are preserved.','']
            (output/'comparison.md').write_text('\n'.join(lines))
            save(output/'progress.json',dict(stage='complete',completed_cases=len(completed),total_cases=len(cases),outcomes=report['new_case_outcomes'],synthesis_process_active=False))
        except BaseException as error:
            report['interruption_or_error']=type(error).__name__+': '+str(error)
            save(output/'progress.json',dict(stage='stopped',completed_cases=len(completed),total_cases=len(cases),synthesis_process_active=False))
            raise
        finally:save(output/'summary.json',report)
        print('FINAL '+json.dumps(report),flush=True)
        return report


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root',type=Path,required=True)
    parser.add_argument('--output',type=Path)
    args=parser.parse_args()
    run(args.root,args.output or args.root/'Results/parity-limit-continuation-4m-60m-v1')
