#!/usr/bin/env python3
"""Finalize completed timeout retries, checking normalized LTL equivalence.

No synthesis, Julia workers, or LLM calls. Keeps the original frozen runner
and config untouched; records a separate provenance trail for this repair.
"""
import argparse
import collections
import copy
import hashlib
import importlib.util
import json
from pathlib import Path
import subprocess


def read(path):
    return json.loads(Path(path).read_text())


def digest(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(1048576),b''):
            h.update(block)
    return h.hexdigest()


def save(path,data):
    path=Path(path)
    temporary=path.with_name(path.name+'.tmp')
    temporary.write_text(json.dumps(data,indent=2,ensure_ascii=False)+'\n')
    temporary.replace(path)


def merge_verified_references(rows,replacements,equivalent):
    merged=copy.deepcopy(rows)
    changes=[]
    for row in merged:
        old=row['synthesis']
        if old.get('subreason')!='timeout':
            continue
        result=replacements[row['group']['class_id']]
        changed=result['formula']!=old['formula']
        if changed:
            if not equivalent(old['formula'],result['formula']):
                raise ValueError('Retry is not semantically equivalent to the source specification')
            changes.append(dict(class_id=row['group']['class_id'],source_formula=old['formula'],
                retry_formula=result['formula'],spot_equivalent=True))
        if bool(result['issued'])!=(result['outcome']=='controller_issued'):
            raise ValueError('Inconsistent issuance flag')
        if result['issued'] and not (result.get('region_known') and result.get('x_win_size',0)>0
                and result.get('certificate',{}).get('satisfied')):
            raise ValueError('Recovered controller lacks a verified initial winning region')
        row['previous_synthesis']=old
        row['synthesis']=copy.deepcopy(result)
        row['retry_formula_equivalence_verified']=True
        row['retry_formula_spelling_changed']=changed
    return merged,changes


def validate_artifact(result):
    outcome=result['outcome']
    if outcome not in ('controller_issued','computational_failure','plant_unrealizable','logically_inconsistent'):
        raise ValueError('Unexpected synthesis outcome')
    if result.get('region_known'):
        base=Path(result['artifact_base_directory'])
        path=base/result['region_file']
        if digest(path)!=result['region_sha256']:
            raise ValueError('Winning-region hash mismatch')
        region=read(path)
        if region['size']!=len(set(region['x_win'])) or region['size']!=result['x_win_size']:
            raise ValueError('Winning-region size mismatch')
        if region['job_id']!=result.get('artifact_job_id',result['job_id']):
            raise ValueError('Winning-region job mismatch')
        if region['formula']!=result.get('artifact_formula',result['formula']):
            raise ValueError('Winning-region formula mismatch')
    if result['issued']:
        path=Path(result['artifact_base_directory'])/result['controller_file']
        if digest(path)!=result['controller_sha256']:
            raise ValueError('Controller hash mismatch')


def finalize(root,output):
    import fcntl
    root,output=Path(root).resolve(),Path(output).resolve()
    with (output/'experiment.lock').open('w') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        report=read(output/'summary.json')
        if report.get('complete'):
            return report
        frozen=read(output/'frozen-config.json')
        source=Path(frozen['source_directory'])
        for name,key in [('summary.json','source_summary_sha256'),
                         ('reference-results.json','source_reference_results_sha256'),
                         ('frozen-config.json','source_frozen_sha256')]:
            if digest(source/name)!=frozen[key]:
                raise ValueError('Source audit changed: '+name)
        if digest(root/'GroundTruthTimeoutRetry/run.py')!=frozen['config']['code_sha256']:
            raise ValueError('Frozen retry runner changed')
        version=subprocess.check_output(['ltl2tgba','--version'],text=True).splitlines()[0]
        if version!=frozen['config']['solver']['spot_version']:
            raise ValueError('Spot version changed')
        plan=read(output/'retry-plan.json')
        if len(plan)!=frozen['retry_count']:
            raise ValueError('Retry plan size mismatch')
        normalized=subprocess.run(['ltlfilt','-r1','--format=%f'],
            input='\n'.join(row['formula'].replace('\n',' ').replace('\r',' ') for row in plan)+'\n',
            text=True,capture_output=True,timeout=60)
        lines=normalized.stdout.splitlines()
        if normalized.returncode or normalized.stderr.strip() or len(lines)!=len(plan):
            raise ValueError('Cannot reproduce batch formula normalization')
        replacements={}
        for row,formula in zip(plan,lines):
            identity=dict(formula=formula,solver_identity=frozen['config']['solver_identity'])
            jid=hashlib.sha256(json.dumps(identity,sort_keys=True).encode()).hexdigest()
            result=read(output/'job-results'/(jid+'.json'))
            if result['job_id']!=jid or result['formula']!=formula:
                raise ValueError('Retry checkpoint identity mismatch')
            validate_artifact(result)
            replacements[row['class_id']]=result
        def equivalent(a,b):
            check=subprocess.run(['ltlfilt','-f',a,'--equivalent-to='+b,'--count'],
                capture_output=True,text=True,timeout=60)
            if check.returncode not in (0,1) or check.stderr.strip():
                raise ValueError('Spot equivalence check failed')
            return check.returncode==0 and check.stdout.strip()=='1'
        source_rows=read(source/'reference-results.json')
        if {r['group']['class_id'] for r in source_rows if r['synthesis'].get('subreason')=='timeout'}!=set(replacements):
            raise ValueError('Retry classes do not exactly match old timeouts')
        merged,changes=merge_verified_references(source_rows,replacements,equivalent)
        spec=importlib.util.spec_from_file_location('frozen_retry_statistics',root/'GroundTruthTimeoutRetry/run.py')
        module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
        stats=module.summarize(merged)
        if stats['semantic_classes']!=1126 or stats['benchmark_records']!=1184:
            raise ValueError('Merged benchmark coverage mismatch')
        save(output/'reference-results.json',merged)
        if not (output/'summary.before-finalization-fix.json').exists():
            (output/'summary.before-finalization-fix.json').write_bytes((output/'summary.json').read_bytes())
        provenance=dict(issue='Strict string equality rejected reordered but equivalent formulas',
            remedy='Reproduce saved normalization and verify all changed spellings with exact Spot equivalence',
            changed_spellings=len(changes),equivalence_checks=changes,
            finalizer_sha256=digest(__file__),source_reference_results_sha256=frozen['source_reference_results_sha256'],
            job_results_verified=len(replacements),new_synthesis_jobs=0,new_api_calls=0)
        save(output/'finalization.json',provenance)
        report.update(stats)
        report['retry_outcomes']=dict(collections.Counter(r['outcome'] for r in replacements.values()))
        report['remaining_failure_reasons']=dict(collections.Counter(
            r['synthesis'].get('subreason') for r in merged if r['synthesis']['outcome']=='computational_failure'))
        report['finalization_repair']=dict(changed_formula_spellings=len(changes),all_spot_equivalent=True,
            provenance_file='finalization.json',previous_error=report.pop('interruption_or_error',None))
        lines=['# Ground-truth controllers: longer computation budget','',
            f"Unchanged production fixed-point solver; {report['timeout_seconds']:g}-second computation limit; {report['workers']} worker.",
            'Simulation remains 1,200 steps; synthesis and verification are infinite-horizon.','',
            '| Outcome | Semantic classes | % of classes | Records | % of records |',
            '|---|---:|---:|---:|---:|']
        for outcome in ('controller_issued','plant_unrealizable','logically_inconsistent','computational_failure'):
            n=stats['class_counts'].get(outcome,0);nr=stats['record_counts'].get(outcome,0)
            lines.append(f"| {outcome} | {n} | {100*n/1126:.2f}% | {nr} | {100*nr/1184:.2f}% |")
        lines+=['',f"Retried {len(replacements)} timeouts: {report['retry_outcomes']}. Other source decisions are retained with provenance.",
            f"Remaining computational failures: {report['remaining_failure_reasons']}.",
            'Computational failures mean unknown controller existence, not proven unrealizability.',
            f'{len(changes)} normalized spellings were checked with exact Spot equivalence before merging. No synthesis rerun was needed to repair finalization.','']
        (output/'comparison.md').write_text('\n'.join(lines))
        report['complete']=True
        save(output/'summary.json',report)
        save(output/'progress.json',dict(stage='complete',completed_jobs=len(replacements),
            total_jobs=len(plan),outcomes=report['retry_outcomes']))
        with (output/'run.log').open('a') as log:
            log.write('\nFINALIZATION REPAIRED: verified normalized formula equivalence; no new synthesis\n')
            log.write('FINAL '+json.dumps(report)+'\n')
        return report


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root',type=Path,required=True)
    parser.add_argument('--output',type=Path)
    args=parser.parse_args()
    report=finalize(args.root,args.output or args.root/'Results/ground-truth-controller-audit-timeout600-v2')
    print(json.dumps(report,indent=2))
