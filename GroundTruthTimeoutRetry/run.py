#!/usr/bin/env python3
"""Retry ground-truth timeouts with the SAME fixed-point solver and more time.

No worklist installation, no parity changes, no LLM calls. Completed source
decisions and hard node-limit outcomes remain intact and retain provenance.
"""
import argparse
import collections
import copy
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import subprocess


def read(path):
    return json.loads(Path(path).read_text())


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def save(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + '.tmp')
    tmp.write_text(json.dumps(value, indent=2, ensure_ascii=False) + '\n')
    tmp.replace(path)


def merge_references(rows, replacements):
    merged = copy.deepcopy(rows)
    for row in merged:
        previous = row['synthesis']
        if previous.get('subreason') == 'timeout':
            new = replacements[row['group']['class_id']]
            if new['formula'] != previous['formula']:
                raise ValueError('Retry formula differs from the reference job')
            if bool(new['issued']) != (new['outcome'] == 'controller_issued'):
                raise ValueError('Inconsistent issuance flag')
            if new['issued'] and not (new.get('region_known') and new.get('x_win_size',0)>0
                    and new.get('certificate',{}).get('satisfied')):
                raise ValueError('Recovered controller lacks verified winning initial states')
            row['previous_synthesis'] = previous
            row['synthesis'] = copy.deepcopy(new)
    return merged


def summarize(rows):
    classes, records = collections.Counter(), collections.Counter()
    for row in rows:
        outcome = row['synthesis']['outcome']
        classes[outcome] += 1
        records[outcome] += len(row['group']['record_ids'])
    n, nr = len(rows), sum(records.values())
    return dict(semantic_classes=n, benchmark_records=nr, class_counts=dict(classes),
        record_counts=dict(records), class_rates={k:v/n for k,v in classes.items()},
        record_rates={k:v/nr for k,v in records.items()},
        controller_existence_lower_bound=classes['controller_issued']/n,
        controller_existence_upper_bound=(classes['controller_issued']+classes['computational_failure'])/n,
        computational_failures_are_unknown=True)


def run(root, output, timeout=600.0, workers=1):
    import fcntl
    if not math.isfinite(timeout) or timeout <= 120 or workers < 1:
        raise ValueError('Require timeout >120 seconds and positive workers')
    root, output = Path(root).resolve(), Path(output).resolve()
    output.mkdir(parents=True, exist_ok=True)
    with (output/'experiment.lock').open('w') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        source = root/'Results/ground-truth-controller-audit-threadsafe-v1'
        original = read(source/'summary.json')
        if not original['complete']:
            raise ValueError('Source ground-truth audit is incomplete')
        frozen_old = read(source/'frozen-config.json')
        source_rows = read(source/'reference-results.json')
        failures = [r for r in source_rows if r['synthesis'].get('subreason') == 'timeout']
        failures.sort(key=lambda r:(len(r['synthesis']['formula']),r['group']['class_id']))
        old_config = frozen_old['config']
        for name, sha in old_config['solver']['code_sha256'].items():
            if digest(root/name) != sha:
                raise ValueError('Production solver changed; this audit changes ONLY computation time: '+name)
        if digest(old_config['plant_cache']) != old_config['plant_sha256']:
            raise ValueError('Plant changed')
        spot_version = subprocess.check_output(['ltl2tgba','--version'],text=True).splitlines()[0]
        if spot_version != old_config['solver']['spot_version']:
            raise ValueError('Spot version changed')
        engine = root/'ControllerCalibration/run.py'
        spec = importlib.util.spec_from_file_location('ground_truth_timeout_backend',engine)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        config = copy.deepcopy(old_config)
        config['timeout_seconds'] = timeout
        config['workers'] = workers
        config['solver']['timeout_seconds'] = timeout
        config['solver_identity'] = module.signature(config['solver'])
        config['code_sha256'] = digest(__file__)
        frozen = dict(config=config,source_directory=str(source),
            source_summary_sha256=digest(source/'summary.json'),
            source_reference_results_sha256=digest(source/'reference-results.json'),
            source_frozen_sha256=digest(source/'frozen-config.json'),
            backend_runner_sha256=digest(engine), retry_count=len(failures),
            algorithm_change=False, simulation_horizon_change=False,
            policy='retry only timeouts with more wall-clock time; preserve completed decisions and node-limit failures')
        path = output/'frozen-config.json'
        if path.exists() and read(path) != frozen:
            raise ValueError('Frozen experiment changed; choose a fresh output directory')
        save(path,frozen)
        if (output/'summary.json').exists() and read(output/'summary.json').get('complete'):
            print('Already complete',flush=True)
            return read(output/'summary.json')
        save(output/'retry-plan.json',[dict(class_id=r['group']['class_id'],formula=r['synthesis']['formula'],
            source_job_id=r['synthesis']['job_id'],old_stage=r['synthesis'].get('last_stage')) for r in failures])
        report = dict(complete=False,new_api_calls=0,retry_count=len(failures),timeout_seconds=timeout,
            workers=workers,algorithm='unchanged production nested Buchi fixed-point solver with existing parity fallback',
            simulation_steps=config['max_steps'],verification_horizon='infinite',
            original_class_counts=original['class_counts'])
        save(output/'summary.json',report)
        save(output/'progress.json',dict(stage='initializing',completed_jobs=0,total_jobs=len(failures),outcomes={}))
        (output/'README.md').write_text(f'''# Ground-truth timeout retry: unchanged solver, {timeout:g} seconds

Production solver and plant hashes must match the completed ground-truth audit.
Only the wall-clock budget and worker count change: {timeout:g} seconds per
formula, {workers} concurrent Julia worker(s), four Julia threads per worker.
The 1,200-step illustrative simulation and infinite-horizon verification remain
unchanged. No worklist or parity action reduction is enabled. No model requests.

Only the {len(failures)} old timeout classes are rerun. Old positive and negative
completed decisions remain valid under the same solver semantics. Hard parity
node-limit failures retain the same graph construction and same node/edge caps;
more computation time cannot remove those errors, so they are retained with
provenance. All original experiment files remain untouched.

Check progress.json or run.log. Each job-results JSON is a resumable checkpoint.
summary.json becomes complete=true after all retries, verification and merged
report generation finish. comparison.md then gives updated class/record rates.
Use SIGINT on the runner to pause gracefully; repeat the same command to resume.
''')
        backend = module.Backend(root,output,config,workers)
        backend.old.clear()
        backend.stage = 'ground-truth-timeout-retry'
        try:
            outcomes = backend.solve_many([r['synthesis']['formula'] for r in failures])
            replacements = {r['group']['class_id']:outcomes[r['synthesis']['formula']] for r in failures}
            merged = merge_references(source_rows,replacements)
            save(output/'reference-results.json',merged)
            report.update(summarize(merged))
            report['retry_outcomes'] = dict(collections.Counter(r['outcome'] for r in replacements.values()))
            report['remaining_failure_reasons'] = dict(collections.Counter(
                r['synthesis'].get('subreason') for r in merged if r['synthesis']['outcome']=='computational_failure'))
            lines = ['# Ground-truth controllers: longer computation budget', '',
                f'Unchanged production fixed-point solver; {timeout:g}-second computation limit, {workers} worker(s).',
                'Simulation remains 1,200 steps; synthesis and verification are infinite-horizon.', '',
                '| Outcome | Semantic classes | % of classes | Records | % of records |',
                '|---|---:|---:|---:|---:|']
            for outcome in ('controller_issued','plant_unrealizable','logically_inconsistent','computational_failure'):
                n = report['class_counts'].get(outcome,0)
                nr = report['record_counts'].get(outcome,0)
                lines.append(f"| {outcome} | {n} | {100*n/report['semantic_classes']:.2f}% | {nr} | {100*nr/report['benchmark_records']:.2f}% |")
            lines += ['',f'Retried {len(failures)} old timeout classes. Other source decisions retained with provenance.',
                'Computational failures remain unknown controller existence, not proven unrealizability.', '']
            (output/'comparison.md').write_text('\n'.join(lines))
            report['complete'] = True
            save(output/'progress.json',dict(stage='complete',completed_jobs=len(failures),total_jobs=len(failures),outcomes=report['retry_outcomes']))
        except BaseException as error:
            report['interruption_or_error'] = type(error).__name__+': '+str(error)
            raise
        finally:
            backend.close()
            save(output/'summary.json',report)
        print('FINAL '+json.dumps(report),flush=True)
        return report


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root',type=Path,required=True)
    parser.add_argument('--output',type=Path)
    parser.add_argument('--timeout',type=float,default=600)
    parser.add_argument('--workers',type=int,default=1)
    args=parser.parse_args()
    label=f'{args.timeout:g}'.replace('.','p')
    run(args.root,args.output or args.root/f'Results/ground-truth-controller-audit-timeout{label}-v2',args.timeout,args.workers)
