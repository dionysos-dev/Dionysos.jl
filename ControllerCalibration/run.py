#!/usr/bin/env python3
"""Calibrate the original cumulative-controller-issuance loss; no model calls.

Only logical-zero losses and already-failed prefixes skip synthesis. Their g or
tilde_g is recorded as unknown unless it is determined; no success is fabricated.
"""
import argparse
import collections
import concurrent.futures
from fractions import Fraction
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import threading
import time


def read(path):
    return json.loads(Path(path).read_text())


def save(path, data):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + '.tmp')
    temporary.write_text(json.dumps(data, indent=2, ensure_ascii=False) + '\n')
    temporary.replace(path)


def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1048576), b''):
            h.update(block)
    return h.hexdigest()


def signature(data):
    return hashlib.sha256(json.dumps(data, sort_keys=True).encode()).hexdigest()


def evaluate_level(prompts, k, previous, solve_many, class_by_record):
    rows, needed = [], []
    for p in prompts:
        values = p['conjunctions']
        c = values[min(k, len(values))-1] if values else None
        prev = previous.get(p['prompt_index'])
        row = dict(prompt_index=p['prompt_index'], record_id=p['record_id'],
            input_field=p['input_field'], class_id=class_by_record[p['record_id']], k=k,
            formula=c['phi_k'] if c else None,
            implies_reference=c['implies_reference'] if c else None,
            conjunction_satisfiable=c['satisfiable'] if c else None,
            tilde_g=None, g=None, loss=0, controller_job_id=None)
        if c is None:
            row.update(g=0, tilde_g=0, evaluation='no_valid_candidate')
        elif prev is not None and prev['g'] == 0:
            row.update(g=0, evaluation='prefix_blocked')
        elif not c['satisfiable']:
            assert c['implies_reference'], 'Unsatisfiable formulas imply the reference vacuously'
            row.update(g=0, tilde_g=0, evaluation='logically_inconsistent')
        elif c['implies_reference']:
            row.update(evaluation='loss_zero_by_implication')
        else:
            if k > 1 and (prev is None or prev['g'] != 1):
                raise ValueError('A non-implying specification needs its full successful synthesis prefix')
            needed.append(row)
        rows.append(row)
    outcomes = solve_many([r['formula'] for r in needed])
    for row in needed:
        outcome = outcomes[row['formula']]
        raw = int(outcome['issued'])
        assert bool(raw) == (outcome['outcome'] == 'controller_issued')
        row.update(tilde_g=raw, g=raw, loss=raw, evaluation='synthesis_evaluated',
            synthesis_outcome=outcome['outcome'], synthesis_subreason=outcome.get('subreason'),
            controller_job_id=outcome['job_id'], artifact_base_directory=outcome.get('artifact_base_directory'),
            controller_file=outcome.get('controller_file'), region_file=outcome.get('region_file'),
            x_win_size=outcome.get('x_win_size'), region_known=outcome.get('region_known', False))
    for row in rows:
        if row['g'] is not None and row['implies_reference'] is not None:
            assert row['loss'] == int(bool(row['g']) and not row['implies_reference'])
        prev = previous.get(row['prompt_index'])
        if prev is not None:
            assert row['loss'] <= prev['loss'], 'Original prefix loss increased'
        upper = int(bool(row['conjunction_satisfiable']) and row['implies_reference'] is False)
        assert row['loss'] <= upper
    return rows


def statistics(rows):
    groups = collections.defaultdict(list)
    for row in rows:
        groups[row['class_id']].append(row['loss'])
    n = len(groups)
    total = sum((Fraction(sum(v), len(v)) for v in groups.values()), Fraction())
    pair_total = sum(row['loss'] for row in rows)
    return dict(n_classes=n, n_prompts=len(rows), class_loss_sum=str(total),
        class_mean_empirical_risk=float(total/n), class_mean_corrected_risk=float((total+1)/(n+1)),
        class_mean_corrected_exact=str((total+1)/(n+1)), pair_losses=pair_total,
        pair_empirical_risk=pair_total/len(rows), pair_corrected_risk=(pair_total+1)/(len(rows)+1),
        class_worst_losses=sum(max(v) for v in groups.values()),
        known_synthesis_outcomes=dict(collections.Counter(r.get('synthesis_outcome') for r in rows if r['evaluation']=='synthesis_evaluated')),
        evaluations=dict(collections.Counter(r['evaluation'] for r in rows)))


class Backend:
    def __init__(self, root, output, config, workers):
        self.root, self.output, self.config = root, output, config
        self.pool = concurrent.futures.ThreadPoolExecutor(max_workers=workers)
        self.local, self.mutex = threading.local(), threading.Lock()
        self.workers, self.norm, self.verified = [], {}, set()
        self.abort = threading.Event()
        self.stage = 'initializing'
        sys.path.insert(0, str(root/'InitialRegionABCD'))
        spec = importlib.util.spec_from_file_location('original_region_worker_runner', root/'InitialRegionABCD/run.py')
        runner = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(runner)
        self.Worker = runner.Worker
        self.old = {}
        old_dir = root/'Results/deepseek-development-abcd-all-initial-v3'
        old_config = read(old_dir/'frozen-development.json')
        compatible = all(digest(self.code_path(root, key)) == value
                         for key, value in old_config['code_sha256'].items())
        if compatible and digest(old_dir/'plant.jls') == config['plant_sha256']:
            old_jobs = {r['job_id']: r for r in read(old_dir/'summary.json')['results'] if r.get('job_id')}
            candidates = [r for r in old_jobs.values() if r['outcome'] in ('controller_issued', 'plant_unrealizable')]
            self.normalize([r['formula'] for r in candidates])
            for r in sorted(candidates, key=lambda r:r['job_id']):
                f = self.norm[r['formula']]
                self.old.setdefault(f, dict(r, artifact_base_directory=str(old_dir),
                    artifact_job_id=r['job_id'], artifact_formula=r['formula']))

    @staticmethod
    def code_path(root, key):
        if key.startswith('runner/'):
            return root/'InitialRegionABCD'/key[len('runner/'):]
        if key == 'plant_builder':
            return root/'DevelopmentABCD/worker.jl'
        return root/key

    def normalize(self, formulas):
        missing = list(dict.fromkeys(f for f in formulas if f not in self.norm))
        if missing:
            lines=[f.replace('\n',' ').replace('\r',' ') for f in missing]
            p = subprocess.run(['ltlfilt', '-r1', '--format=%f'], input='\n'.join(lines)+'\n',
                               capture_output=True, text=True, timeout=60)
            if p.returncode != 0 or p.stderr.strip() or len(p.stdout.splitlines()) != len(missing):
                raise RuntimeError('Basic Spot normalization failed: '+p.stderr)
            self.norm.update(zip(missing, p.stdout.splitlines()))

    def validate(self, result, path):
        if path in self.verified:
            return
        assert result['outcome'] in ('controller_issued', 'plant_unrealizable', 'computational_failure', 'logically_inconsistent')
        assert bool(result['issued']) == (result['outcome']=='controller_issued')
        if result.get('region_known') and result.get('region_file'):
            base = Path(result['artifact_base_directory'])
            assert digest(base/result['region_file']) == result['region_sha256']
            region = read(base/result['region_file'])
            assert region['size'] == len(set(region['x_win'])) == result['x_win_size']
            assert region['job_id'] == result.get('artifact_job_id', result['job_id'])
            assert region['formula'] == result.get('artifact_formula', result['formula'])
        if result['issued']:
            assert result['region_known'] and result['x_win_size'] > 0
            assert result['certificate']['satisfied']
            base = Path(result['artifact_base_directory'])
            assert digest(base/result['controller_file']) == result['controller_sha256']
        elif result['outcome'] in ('plant_unrealizable','logically_inconsistent'):
            assert result['region_known'] and result['x_win_size'] == 0
        self.verified.add(path)

    def get_worker(self):
        worker = getattr(self.local, 'worker', None)
        if worker is None or worker.process.poll() is not None:
            if worker is None:
                with self.mutex:
                    number = len(self.workers)
                    self.workers.append(None)
                folder = self.output/'workers'/f'worker-{number:02d}'
                folder.mkdir(parents=True, exist_ok=True)
                shutil.copy2(self.config['plant_cache'], folder/'plant.jls')
                self.local.number, self.local.folder = number, folder
            worker = self.Worker(self.root, self.local.folder, max_steps=self.config['max_steps'])
            self.local.worker = worker
            self.workers[self.local.number] = worker
            if self.abort.is_set():
                worker.stop()
                raise InterruptedError('Backend initialization interrupted')
            assert digest(self.local.folder/'plant.jls') == self.config['plant_sha256']
            ready = read(self.local.folder/'worker-ready.json')
            assert ready['initial_domain']=='all_admissible_plant_states'
            assert ready['admissible_initial_state_count']==44643
        return worker

    def execute(self, formula):
        jid = signature(dict(formula=formula, solver_identity=self.config['solver_identity']))
        path = self.output/'job-results'/f'{jid}.json'
        if path.exists():
            result = read(path)
            if 'artifact_base_directory' not in result:
                job = read(self.output/'jobs'/f'{jid}.json')
                result['artifact_base_directory'] = job['worker_output']
        elif formula in self.old:
            result = dict(self.old[formula], job_id=jid, formula=formula, reused_verified_development_artifact=True)
            save(path, result)
        else:
            if self.abort.is_set():
                raise InterruptedError('Synthesis interrupted')
            worker = self.get_worker()
            job = dict(job_id=jid, formula=formula, plant_sha256=self.config['plant_sha256'], worker_output=str(worker.output))
            job_path = self.output/'jobs'/f'{jid}.json'
            save(job_path, job)
            result = worker.solve(job, job_path, path, self.config['timeout_seconds'])
            if self.abort.is_set():
                raise InterruptedError('Synthesis interrupted; incomplete jobs are not scored')
            result['artifact_base_directory'] = str(worker.output)
            save(path, result)
        assert result['job_id']==jid and result['formula']==formula
        self.validate(result, path)
        return result

    def solve_many(self, formulas):
        self.normalize(formulas)
        unique = list(dict.fromkeys(self.norm[f] for f in formulas))
        if not unique:
            return {}
        print(f'{self.stage}: {len(formulas)} affected prompts, {len(unique)} distinct normalized controller jobs', flush=True)
        futures = {self.pool.submit(self.execute, f): f for f in unique}
        outcomes, last = {}, time.monotonic()
        try:
            for future in concurrent.futures.as_completed(futures):
                outcomes[futures[future]] = future.result()
                if len(outcomes)%10==0 or time.monotonic()-last>30:
                    progress=dict(stage=self.stage, completed_jobs=len(outcomes), total_jobs=len(unique),
                        outcomes=dict(collections.Counter(r['outcome'] for r in outcomes.values())))
                    save(self.output/'progress.json', progress)
                    print(json.dumps(progress), flush=True)
                    last=time.monotonic()
        except BaseException:
            self.abort.set()
            for worker in self.workers:
                if worker is not None:
                    worker.stop()
            raise
        return {f: outcomes[self.norm[f]] for f in formulas}

    def close(self):
        self.abort.set()
        for worker in self.workers:
            if worker is not None:
                worker.stop()
        self.pool.shutdown(wait=True, cancel_futures=True)


def run(root, output, alpha=Fraction(1,5), workers=2, timeout=120):
    root, output = Path(root).resolve(), Path(output).resolve()
    if not 0<alpha<1 or workers<1 or timeout<=0:
        raise ValueError('Require 0<alpha<1, positive workers, and positive timeout')
    output.mkdir(parents=True, exist_ok=True)
    import fcntl
    with (output/'experiment.lock').open('w') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX|fcntl.LOCK_NB)
        split=root/'Benchmark/TaskSpecificSplits/semantic-70-30-seed-20261027'
        semantics=root/'Results/deepseek-v4-pro-task-specific-full-T0.7-K10-conjunctions/summary.json'
        manifest, groups, analysis = read(split/'manifest.json'), read(split/'semantic_groups.json'), read(semantics)
        assert analysis['complete'] and analysis['repetitions']==10 and analysis['prompt_count']==4736
        assert analysis['model']=='deepseek-v4-pro' and analysis['temperature']==0.7
        raw_path=Path(analysis['source_path'])
        assert digest(raw_path)==analysis['analysis_identity']['source_sha256']
        raw=read(raw_path)
        assert raw['config']['dataset_sha256']==manifest['dataset_sha256'] and raw['summary']['complete']
        for filename, checksum in manifest['file_sha256'].items():
            assert digest(split/filename)==checksum
        class_by_record={rid:g['class_id'] for g in groups for rid in g['record_ids']}
        sets={name:set(s['record_ids']) for name,s in manifest['splits'].items()}
        prompts={name:[p for p in analysis['prompt_results'] if p['record_id'] in ids] for name,ids in sets.items()}
        originals={r['id']:r for name in sets for r in read(split/(name+'.json'))}
        keys=set()
        for p in analysis['prompt_results']:
            key=(p['record_id'],p['input_field'])
            assert key not in keys
            keys.add(key)
            assert p['reference_ltl']==originals[p['record_id']]['LTL']
            assert p['source_text']==originals[p['record_id']][p['input_field']]
        assert len(keys)==4736
        for name, values in prompts.items():
            assert len(values)==manifest['splits'][name]['active_nl_ltl_pair_count']
        plant=root/'Results/deepseek-development-abcd-all-initial-v3/plant.jls'
        frozen_old=read(plant.parent/'frozen-development.json')
        assert digest(plant)==read(plant.parent/'plant-identity.json')['sha256']
        code={str(p.relative_to(root)):digest(p) for p in sorted((root/'src').rglob('*.jl'))}
        for name in ['InitialRegionABCD/worker.jl','InitialRegionABCD/ParityBackend.jl',
                     'InitialRegionABCD/spot_parity_bridge.cpp','InitialRegionABCD/build/spot_parity_bridge',
                     'InitialRegionABCD/run.py','DevelopmentABCD/worker.jl','evaluate_claude_buchi_synthesis.jl',
                     'System_helper.jl','Project.toml','Manifest.toml']:
            code[name]=digest(root/name)
        solver=dict(policy='deterministic_buchi_first_then_verified_global_spot_parity',
            initial_domain='all_admissible_plant_states', semantics='robust',
            projection='(x,delta(q_init,label(x))) in product winning set', timeout_seconds=timeout,
            max_steps=1200, simulation_is_proof=False, verification_horizon='infinite',
            normalization='Spot ltlfilt -r1; semantics-preserving basic normalization',
            plant_sha256=digest(plant), code_sha256=code,
            spot_version=subprocess.check_output(['ltl2tgba','--version'],text=True).splitlines()[0])
        config=dict(schema_version=1, model=analysis['model'], temperature=analysis['temperature'], repetitions=10,
            alpha=str(alpha), max_k=10, workers=workers, timeout_seconds=timeout, max_steps=1200,
            plant_cache=str(plant), plant_sha256=solver['plant_sha256'], plant=frozen_old['plant'],
            solver=solver, solver_identity=signature(solver),
            loss='g_k * indicator(Phi_k does not imply reference)', issuance='g_k = product_{j<=k} tilde_g_j',
            risk_unit='uniform semantic class; mean original binary prompt losses within class',
            logical_zero_shortcut='loss=0 if implication true, without inventing g',
            prefix_shortcut='g=0 after first failed synthesis; unattempted tilde_g stays null',
            source_semantics_sha256=digest(semantics), split_manifest_sha256=digest(split/'manifest.json'),
            code_sha256=digest(__file__))
        config_path=output/'frozen-config.json'
        if config_path.exists() and read(config_path)!=config:
            raise RuntimeError('Frozen experiment changed; use a fresh output directory')
        save(config_path, config)
        backend=Backend(root, output, config, workers)
        report=dict(complete=False, config=config, calibration=[], selected_k=None, test=None)
        try:
            previous={}
            for k in range(1,11):
                backend.stage=f'calibration-k{k}'
                rows=evaluate_level(prompts['calibration'], k, previous, backend.solve_many, class_by_record)
                stats=statistics(rows)
                assert stats['n_classes']==788 and stats['n_prompts']==3288
                stats['k']=k
                report['calibration'].append(stats)
                save(output/f'cases-calibration-k{k:02d}.json',rows)
                save(output/'summary.json', report)
                previous={r['prompt_index']:r for r in rows}
                print('CALIBRATION COMPLETE '+json.dumps(stats),flush=True)
                if Fraction(stats['class_mean_corrected_exact'])<=alpha:
                    report['selected_k']=k
                    selection=dict(k=k,alpha=str(alpha),calibration=stats,selection_data='calibration_only',
                        frozen_config_sha256=digest(config_path))
                    selected=output/'selected-threshold.json'
                    if selected.exists() and read(selected)!=selection:
                        raise RuntimeError('Previously frozen threshold changed')
                    save(selected,selection)
                    break
            if report['selected_k'] is not None:
                previous={}
                for k in range(1,report['selected_k']+1):
                    backend.stage=f'test-k{k}'
                    rows=evaluate_level(prompts['test'], k, previous, backend.solve_many, class_by_record)
                    save(output/f'cases-test-k{k:02d}.json',rows)
                    previous={r['prompt_index']:r for r in rows}
                report['test']=statistics(rows)
                assert report['test']['n_classes']==338 and report['test']['n_prompts']==1448
                report['test']['k']=report['selected_k']
                report['test']['empirically_below_alpha']=report['test']['class_mean_empirical_risk']<=float(alpha)
            report['complete']=True
        finally:
            backend.close()
            save(output/'summary.json',report)
        lines=['# Original controller-aware CRC calibration','',
            'Loss = cumulative issuance times non-implication. Complete: '+str(report['complete']),
            'Calibration unit: semantic LTL class, averaging binary prompt losses within each class.','',
            '| k | Calibration corrected class-mean risk | Pair-weighted corrected diagnostic |',
            '|---|---:|---:|']
        for s in report['calibration']:
            lines.append(f"| {s['k']} | {100*s['class_mean_corrected_risk']:.4f}% | {100*s['pair_corrected_risk']:.4f}% |")
        lines+=['',f"Selected k: {report['selected_k']}. Selection uses calibration only.",'',
            'Test result: '+json.dumps(report['test'],indent=2),'',
            'Synthesis skips logical-zero losses and already-failed prefixes. Such skips do not fabricate issuance outcomes; null g/tilde_g means not measured. All computed losses equal the original loss. The report does not estimate overall issuance coverage for the skipped logical-zero cases.']
        (output/'comparison.md').write_text('\n'.join(lines)+'\n')
        print('FINAL '+json.dumps({k:report[k] for k in ['complete','calibration','selected_k','test']}),flush=True)
        return report


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root',type=Path,default=Path(__file__).resolve().parent.parent)
    parser.add_argument('--output',type=Path)
    parser.add_argument('--alpha',default='0.20')
    parser.add_argument('--workers',type=int,default=2)
    parser.add_argument('--timeout',type=float,default=120.0)
    args=parser.parse_args()
    run(args.root,args.output or args.root/'Results/deepseek-r10-crc-original-alpha0.20',
        Fraction(args.alpha),args.workers,args.timeout)
