#!/usr/bin/env python3
"""Multi-alpha CRC table, existing DeepSeek samples, cumulative issuance.

No LLM requests. Calibration alone selects k; complete test issuance is then
measured. Every completed distinct formula is checkpointed. SIGINT pauses safely.
"""
import argparse
import collections
import copy
import fcntl
from fractions import Fraction
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import re
import subprocess
import time

ALPHAS = tuple(Fraction(x) for x in ('0.05', '0.10', '0.15', '0.20', '0.25'))


def read(path):
    return json.loads(Path(path).read_text())


def save(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + '.tmp')
    tmp.write_text(json.dumps(value, indent=2, ensure_ascii=False) + '\n')
    tmp.replace(path)


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1048576), b''):
            h.update(block)
    return h.hexdigest()


def key(value):
    return hashlib.sha256(value.encode()).hexdigest()


def load_module(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def conjunction(prompt, k):
    values = prompt['conjunctions']
    return values[min(k, len(values))-1] if values else None


def validate_ranking(prompts):
    seen = set()
    for p in prompts:
        assert p['prompt_index'] not in seen
        seen.add(p['prompt_index'])
        assert p['samples'] == 10
        clusters = p['ranked_clusters']
        assert clusters == sorted(clusters, key=lambda c: (-c['frequency'], min(c['sample_indices'])))
        assert sum(c['frequency'] for c in clusters) == p['satisfiable_samples']
        indices = []
        assert len(p['conjunctions']) == len(clusters)
        for rank, c in enumerate(clusters, 1):
            assert c['rank'] == rank and c['satisfiable']
            assert len(c['sample_indices']) == c['frequency']
            indices.extend(c['sample_indices'])
            expected = ' & '.join('('+r['representative']+')' for r in clusters[:rank])
            assert p['conjunctions'][rank-1]['phi_k'] == expected
            assert p['conjunctions'][rank-1]['cluster_ranks'] == list(range(1, rank+1))
        assert len(set(indices)) == len(indices)
        assert all(1 <= i <= 10 for i in indices)
    return dict(prompts_checked=len(seen), ranking='frequency_descending_then_earliest_query',
                conjunction='top min(k,m) semantic-class representatives', repetitions=10)


def evaluate(prompts, k, previous, solve, class_by_record, complete_issuance=False):
    rows = []
    for p in prompts:
        c = conjunction(p, k)
        prev = previous.get(p['prompt_index'])
        row = dict(prompt_index=p['prompt_index'], record_id=p['record_id'],
                   class_id=class_by_record[p['record_id']], input_field=p['input_field'], k=k,
                   formula=c['phi_k'] if c else None,
                   implies_reference=c['implies_reference'] if c else None,
                   conjunction_satisfiable=c['satisfiable'] if c else None,
                   g=None, tilde_g=None, loss=0, guaranteed_issuance=None)
        if c is None:
            row.update(g=0, tilde_g=0, evaluation='no_valid_candidate')
        elif prev is not None and prev['g'] == 0:
            row.update(g=0, evaluation='prefix_blocked')
        elif not c['satisfiable']:
            assert c['implies_reference']
            row.update(g=0, tilde_g=0, evaluation='logically_inconsistent')
        elif c['implies_reference'] and not complete_issuance:
            row.update(evaluation='loss_zero_by_implication')
        else:
            if k > 1 and (prev is None or prev['g'] != 1):
                raise ValueError('An active issuance decision requires its successful prefix')
            result = solve(c['phi_k'], p, k)
            issued = int(result['issued'])
            row.update(g=issued, tilde_g=issued, evaluation='synthesis_evaluated',
                       synthesis_outcome=result['outcome'], synthesis_subreason=result.get('subreason'),
                       controller_job_id=result['job_id'], x_win_size=result.get('x_win_size'),
                       result_file=result['checkpoint_file'])
            row['loss'] = int(issued and not c['implies_reference'])
        if row['g'] is not None:
            row['guaranteed_issuance'] = int(row['g'] and bool(row['implies_reference']))
            assert row['loss'] == int(row['g'] and not row['implies_reference'])
        if prev is not None:
            assert row['loss'] <= prev['loss'], 'Nonmonotone prefix loss'
            if complete_issuance:
                assert row['g'] <= prev['g']
        rows.append(row)
    return rows


def statistics(rows):
    groups = collections.defaultdict(list)
    for r in rows:
        groups[r['class_id']].append(r)
    n = len(groups)
    def total(field):
        return sum((Fraction(sum(r[field] for r in rs), len(rs)) for rs in groups.values()), Fraction())
    loss = total('loss')
    out = dict(n_classes=n, n_prompts=len(rows), class_loss_sum=str(loss),
               empirical_risk=float(loss/n), corrected_risk=float((loss+1)/(n+1)),
               corrected_exact=str((loss+1)/(n+1)),
               pair_empirical_risk=sum(r['loss'] for r in rows)/len(rows))
    if all(r['g'] is not None for r in rows):
        issuance, guaranteed = total('g')/n, total('guaranteed_issuance')/n
        assert issuance == loss/n + guaranteed
        out.update(issuance=float(issuance), guaranteed_issuance=float(guaranteed),
                   conditional_error=float(loss/n/issuance) if issuance else None,
                   pair_issuance=sum(r['g'] for r in rows)/len(rows),
                   pair_guaranteed_issuance=sum(r['guaranteed_issuance'] for r in rows)/len(rows))
    return out


def select(calibration, alpha, max_k=10):
    for s in calibration:
        if Fraction(s['corrected_exact']) <= alpha:
            return s['k']
    return 'infinity' if len(calibration) == max_k else None


def absolute_result(result, base):
    result = copy.deepcopy(result)
    base = Path(result.get('artifact_base_directory', base))
    for field in ('controller_file', 'region_file'):
        if result.get(field):
            result[field] = str(base/Path(result[field]))
    return result


def validate_result(r):
    assert bool(r['issued']) == (r['outcome'] == 'controller_issued')
    if r['issued']:
        assert r.get('region_known') and r.get('x_win_size', 0) > 0
        assert r.get('certificate', {}).get('satisfied')
        for f, h in [('controller_file', 'controller_sha256'), ('region_file', 'region_sha256')]:
            assert sha(r[f]) == r[h], 'Saved proof artifact changed'
        region = read(r['region_file'])
        assert len(set(region['x_win'])) == region['size'] == r['x_win_size']
    elif r['outcome'] == 'plant_unrealizable':
        assert r['x_win_size'] == 0


class Engine:
    def __init__(self, root, output, config, formulas, sources, reference_results):
        self.root, self.output, self.config = root, output, config
        self.helper = load_module(root/'ParityLimitPilot/run_edge60m.py', 'crc_bounded_helper')
        # One normalization batch fixes Spot's AP ordering for the entire run.
        p = subprocess.run(['ltlfilt', '-r1', '--format=%f'],
                           input='\n'.join(f.replace('\n', ' ').replace('\r', ' ') for f in formulas)+'\n',
                           text=True, capture_output=True, timeout=120, check=True)
        normalized = p.stdout.splitlines()
        if p.stderr.strip() or len(normalized) != len(formulas):
            raise RuntimeError('Spot normalization failed')
        self.norm = dict(zip(formulas, normalized))
        self.cache, self.refs, self.validated = {}, {}, set()
        for result, provenance in sources:
            n = self.norm[result['formula']]
            old = self.cache.get(n)
            # A mathematical decision takes precedence over a resource failure.
            if old is None or (old[0]['outcome'] == 'computational_failure' and result['outcome'] != 'computational_failure'):
                self.cache[n] = (result, provenance)
        for rid, (result, provenance) in reference_results.items():
            self.refs[rid] = (result, provenance)
        self.new_count = 0
        self.active_seconds = 0
        self.stage = 'preparing'
        self.context = {}
        save(output/'normalization.json', self.norm)

    def jid(self, formula):
        return key(self.norm[formula])

    def wait(self, seconds, reason):
        end = time.monotonic()+seconds
        while time.monotonic() < end:
            save(self.output/'progress.json', dict(self.context, stage='cooldown',
                 reason=reason, seconds_remaining=round(end-time.monotonic()), synthesis_process_active=False))
            time.sleep(min(2, max(0, end-time.monotonic())))

    def solve(self, formula, prompt, k):
        self.context.update(prompt_index=prompt['prompt_index'], record_id=prompt['record_id'])
        jid = self.jid(formula)
        folder = self.output/'cases'/jid
        folder.mkdir(parents=True, exist_ok=True)
        path = folder/'result.json'
        if path.exists():
            result = read(path)
            assert result['normalized_formula'] == self.norm[formula]
            if jid not in self.validated:
                validate_result(result)
                self.validated.add(jid)
            return result
        reuse = self.cache.get(self.norm[formula])
        # The saved semantic-equivalence flag is an exact Spot check, not an
        # implication shortcut. It permits reuse when syntax normalization differs.
        if reuse is None and len(prompt['ranked_clusters']) and (
                k == 1 or len(prompt['ranked_clusters']) == 1) and prompt['ranked_clusters'][0]['equivalent_to_reference']:
            reuse = self.refs[prompt['record_id']]
        if reuse is not None:
            result, provenance = reuse
            result = copy.deepcopy(result)
            validate_result(result)
            result.update(reused_from=provenance, reused_artifact_formula=result['formula'])
        else:
            while self.helper.memory_free_percent() < 15:
                self.wait(60, 'waiting for system memory')
            result = self.compute(formula, jid, folder)
            self.new_count += 1
        result.update(job_id=jid, requested_formula=formula, normalized_formula=self.norm[formula],
                      checkpoint_file=str(path))
        validate_result(result)
        save(path, result)
        self.validated.add(jid)
        # Cache exact language aliases, including resource failures, within this run.
        self.cache[self.norm[formula]] = (result, str(path))
        print('DECISION '+json.dumps(dict(self.context, formula=formula, outcome=result['outcome'],
              reused=bool(result.get('reused_from')), new_synthesis_jobs=self.new_count)), flush=True)
        save(self.output/'progress.json', dict(self.context, stage='case_complete',
             formula=formula, outcome=result['outcome'], new_synthesis_jobs=self.new_count,
             synthesis_process_active=False))
        if reuse is None:
            rest = 300 if self.active_seconds >= 1800 else self.config['cooldown_seconds']
            self.wait(rest, 'scheduled rest after new synthesis')
            if rest == 300:
                self.active_seconds = 0
        return result

    def compute(self, formula, jid, folder):
        start = time.monotonic()
        result = dict(formula=formula, issued=False, outcome='computational_failure',
                      artifact_base_directory=str(folder), subreason=None)
        def bounded(command):
            return self.helper.run_bounded(command, folder, 600-(time.monotonic()-start), 3*1024*1024, 15)
        save(self.output/'progress.json', dict(self.context, stage='translation', formula=formula,
             synthesis_process_active=True))
        hoa = folder/'buchi.hoa'
        execution = bounded(['ltl2tgba', '-D', '-B', '-S', '--hoa', '-f', formula, '-o', str(hoa)])
        if execution['guard'] or execution['exit_code']:
            result['subreason'] = execution['guard'] or 'translation_failure'
        else:
            check = subprocess.run(['autfilt', '--is-deterministic', '--count', str(hoa)],
                                   capture_output=True, text=True, timeout=30)
            if check.returncode not in (0, 1):
                raise RuntimeError('Spot determinism check failed: '+check.stderr)
            deterministic = check.returncode == 0
            job = dict(job_id=jid, formula=formula, plant_cache=self.config['plant_cache'],
                       plant_sha256=self.config['plant_sha256'], plant_binary=self.config['plant_binary'])
            save(folder/'case.json', job)
            if deterministic:
                save(self.output/'progress.json', dict(self.context, stage='buchi_synthesis_and_verification',
                     formula=formula, synthesis_process_active=True))
                # The production artifact refers to plant.jls in its own folder.
                plant = folder/'plant.jls'
                if not plant.exists():
                    plant.symlink_to(self.config['plant_cache'])
                execution = bounded(['/Users/abayat/.juliaup/bin/julia', '--startup-file=no',
                    '--project='+str(self.root), '--threads=4', str(self.root/'CalibrationTable/buchi_case.jl'),
                    str(self.root), str(folder), '1200', str(folder/'case.json')])
                if execution['guard'] or execution['exit_code']:
                    result['subreason'] = execution['guard'] or 'buchi_worker_failure'
                else:
                    result.update(absolute_result(read(folder/'verified-result.json'), folder))
                    assert result['backend'] == 'buchi'
            else:
                result['backend'] = 'parity'
                dpa = folder/'specification.hoa'
                execution = bounded(['ltl2tgba', '-D', '-C', '--parity=max odd', '-S', '--hoa',
                                     '-f', formula, '-o', str(dpa)])
                if execution['guard'] or execution['exit_code']:
                    result['subreason'] = execution['guard'] or 'parity_translation_failure'
                else:
                    job.update(hoa_path=str(dpa))
                    save(folder/'case.json', job)
                    save(self.output/'progress.json', dict(self.context, stage='native_parity_game',
                         formula=formula, synthesis_process_active=True))
                    execution = bounded([str(self.root/'InitialRegionABCD/build/spot_parity_bridge'),
                        self.config['plant_binary'], str(dpa), str(folder/'spot-result.json'),
                        str(folder/'spot-policy.bin'), '1200', '4000000', '60000000'])
                    if execution['guard']:
                        result['subreason'] = execution['guard']
                    elif execution['exit_code']:
                        result['subreason'] = self.helper.classify_bridge_failure((folder/'run.log').read_text())
                    else:
                        info = read(folder/'spot-result.json')
                        result['native_statistics'] = {k:v for k,v in info.items() if not isinstance(v, list)}
                        if not info['winning']:
                            assert info['winning_initial_state_count'] == 0
                            result.update(outcome='plant_unrealizable', x_win_size=0, region_known=True,
                                          subreason='empty_admissible_initial_winning_region')
                        else:
                            assert info['verified']
                            save(self.output/'progress.json', dict(self.context, stage='independent_parity_verification',
                                 formula=formula, synthesis_process_active=True))
                            execution = bounded(['/Users/abayat/.juliaup/bin/julia', '--startup-file=no',
                                '--project='+str(self.root), '--threads=1', str(self.root/'ParityLimitPilot/verify.jl'),
                                str(self.root), str(folder), '1200', str(folder/'case.json')])
                            if execution['guard'] or execution['exit_code']:
                                result['subreason'] = execution['guard'] or 'independent_verification_or_export_failure'
                            else:
                                result.update(read(folder/'verified-result.json'))
        result['last_execution'] = execution
        result['elapsed_seconds'] = time.monotonic()-start
        self.active_seconds += result['elapsed_seconds']
        return result


def prepare(root, output, cooldown):
    split = root/'Benchmark/TaskSpecificSplits/semantic-cal300-seed-20261028'
    manifest = read(split/'manifest.json')
    for name, checksum in manifest['file_sha256'].items():
        assert sha(split/name) == checksum
    groups = read(split/'semantic_groups.json')
    classes = {rid:g['class_id'] for g in groups for rid in g['record_ids']}
    source = root/'Results/deepseek-v4-pro-task-specific-full-T0.7-K10-conjunctions/summary.json'
    analysis = read(source)
    assert analysis['complete'] and analysis['repetitions'] == 10
    assert analysis['model'] == 'deepseek-v4-pro' and analysis['temperature'] == 0.7
    assert sha(analysis['source_path']) == analysis['analysis_identity']['source_sha256']
    check = validate_ranking(analysis['prompt_results'])
    assert check['prompts_checked'] == 4736
    prompts = {name:[p for p in analysis['prompt_results'] if p['record_id'] in set(s['record_ids'])]
               for name,s in manifest['splits'].items()}
    assert len(prompts['calibration']) == 1256 and len(prompts['test']) == 3480
    old = read(root/'Results/ground-truth-controller-audit-timeout600-v2/frozen-config.json')['config']
    for name, checksum in old['solver']['code_sha256'].items():
        assert sha(root/name) == checksum, 'Production solver changed: '+name
    assert sha(old['plant_cache']) == old['plant_sha256']
    refs_path = root/'Results/parity-limit-continuation-4m-60m-v1/reference-results.json'
    assert read(refs_path.parent/'summary.json')['complete']
    ref_rows = read(refs_path)
    sources = []
    refs = {}
    source_hashes = {str(refs_path):sha(refs_path)}
    for row in ref_rows:
        result = absolute_result(row['synthesis'], refs_path.parent)
        sources.append((result, str(refs_path)))
        for rid in row['group']['record_ids']:
            refs[rid] = (result, str(refs_path))
    # Only verified positive results may be imported from the pre-race-fix run.
    # An old empty winning-set result is not a sound negative after that fix.
    for name, positives_only in [('deepseek-r10-crc-original-alpha0.20', True),
                                 ('deepseek-r10-computational-failures-threadsafe-v1', False)]:
        directory = root/'Results'/name
        for path in sorted((directory/'job-results').glob('*.json')):
            result = read(path)
            if result['outcome'] == 'computational_failure':
                continue  # old 120s/2M limits do not define this experiment
            if positives_only and not result['issued']:
                continue
            if not result.get('artifact_base_directory'):
                job = read(directory/'jobs'/(path.stem+'.json'))
                result['artifact_base_directory'] = job['worker_output']
            result = absolute_result(result, directory)
            sources.append((result, str(path)))
            source_hashes[str(path)] = sha(path)
    formulas = list(dict.fromkeys([c['phi_k'] for p in analysis['prompt_results'] for c in p['conjunctions']]
                                 +[r['formula'] for r,_ in sources]))
    native_plants = sorted((root/'Results/ground-truth-controller-audit-threadsafe-v1/workers').glob('*/plant.parity.bin'))
    assert native_plants
    config = dict(model='deepseek-v4-pro', temperature=0.7, repetitions=10, max_k=10,
                  alphas=[str(a) for a in ALPHAS], plant_cache=old['plant_cache'], plant_sha256=old['plant_sha256'],
                  plant_binary=str(native_plants[0]), plant_binary_sha256=sha(native_plants[0]),
                  code_sha256={**old['solver']['code_sha256'],
                      'CalibrationTable/run.py':sha(root/'CalibrationTable/run.py'),
                      'CalibrationTable/buchi_case.jl':sha(root/'CalibrationTable/buchi_case.jl'),
                      'ParityLimitPilot/run_edge60m.py':sha(root/'ParityLimitPilot/run_edge60m.py'),
                      'ParityLimitPilot/verify.jl':sha(root/'ParityLimitPilot/verify.jl')},
                  spot_version=subprocess.check_output(['ltl2tgba','--version'],text=True).splitlines()[0],
                  source_sha256=sha(source), split_manifest_sha256=sha(split/'manifest.json'),
                  reused_source_sha256=source_hashes, timeout_seconds=600, max_nodes=4000000,
                  max_edges=60000000, workers=1, buchi_threads=4, parity_verification_threads=1,
                  memory_guard_gib=3, system_free_floor_percent=15, cooldown_seconds=cooldown,
                  timeout_scope='new cases: translation, process startup, solving, verification and export',
                  reuse='resolved decisions with compatible plant and verified proof; latest matching-budget reference failures; no old negative Buchi results',
                  risk_unit='uniform semantic class then uniform rewritten NL prompt within class',
                  loss='product_{j<=k}(tilde_g_j) * indicator(Phi_k does not imply reference)',
                  threshold_selection='calibration only; first corrected risk <= alpha; infinity if none through k=10',
                  infinity_policy='always abstain, risk=issuance=guaranteed_issuance=0',
                  algorithm_change=False, new_api_calls=0)
    frozen = output/'frozen-config.json'
    if frozen.exists():
        assert read(frozen) == config, 'Frozen experiment changed; use a new output directory'
    save(frozen, config)
    save(output/'ranking-validation.json', check)
    engine = Engine(root, output, config, formulas, sources, refs)
    inventory = dict(prompts=4736, calibration_prompts=1256, test_prompts=3480,
                     distinct_normalized_conjunctions=len({engine.norm[c['phi_k']] for p in analysis['prompt_results'] for c in p['conjunctions']}),
                     existing_normalized_decisions=len(engine.cache), new_api_calls=0)
    for name, values in prompts.items():
        top = [p for p in values if conjunction(p,1) is not None]
        unmatched = {engine.norm[conjunction(p,1)['phi_k']] for p in top
                     if engine.norm[conjunction(p,1)['phi_k']] not in engine.cache
                     and not p['ranked_clusters'][0]['equivalent_to_reference']}
        inventory[name+'_k1_missing_distinct_formulas'] = len(unmatched)
    save(output/'inventory.json', inventory)
    return prompts, classes, engine, inventory


def table(report, output):
    entries = []
    for alpha in ALPHAS:
        selected = select(report['calibration'], alpha)
        row = dict(alpha=float(alpha), selected_k=selected, test_risk=None, issuance=None, guaranteed_issuance=None)
        if selected == 'infinity':
            row.update(test_risk=0, issuance=0, guaranteed_issuance=0, status='complete_always_abstain')
        elif selected is not None and str(selected) in report['test']:
            s = report['test'][str(selected)]
            row.update(test_risk=s['empirical_risk'], issuance=s['issuance'],
                       guaranteed_issuance=s['guaranteed_issuance'], status='complete',
                       empirically_below_alpha=s['empirical_risk'] <= float(alpha))
        else:
            row['status'] = 'pending_calibration' if selected is None else 'pending_test_synthesis'
        entries.append(row)
    report['table'] = entries
    save(output/'summary.json', report)
    lines = ['# DeepSeek controller-aware risk and issuance', '',
             'All percentages use uniform semantic-class means, matching calibration.',
             '300 calibration classes; 826 test classes. No new LLM calls.', '',
             '| Target alpha | Selected k | Test risk (%) | Issuance (%) | Guaranteed issuance (%) | Status |',
             '|---:|---:|---:|---:|---:|---|']
    for r in entries:
        k = 'infinity' if r['selected_k'] == 'infinity' else r['selected_k']
        cells = ['pending' if r[f] is None else f'{100*r[f]:.2f}' for f in ('test_risk','issuance','guaranteed_issuance')]
        lines.append(f"| {r['alpha']:.2f} | {k if k is not None else 'pending'} | "+' | '.join(cells)+f" | {r['status']} |")
    lines += ['', 'Issuance uses the successful prefix product, not merely synthesis at selected k.',
              'Test risk + guaranteed issuance = issuance. Computational failure causes abstention.',
              'The historical 16.32% used the earlier solver/resources and is not inserted into this new evaluation.',
              'This resplit has prior benchmark evaluation history; it is not a pristine held-out test.', '']
    (output/'table.md').write_text('\n'.join(lines))
    latex = []
    for r in entries:
        k = '$\\infty$' if r['selected_k'] == 'infinity' else str(r['selected_k'] or '--')
        cells = ['--' if r[f] is None else f'{100*r[f]:.2f}' for f in ('test_risk','issuance','guaranteed_issuance')]
        latex.append(f"{r['alpha']:.2f} & {k} & "+' & '.join(cells)+r' \\')
    (output/'table.tex').write_text('\n'.join(latex)+'\n')


def run(root, output, cooldown=60, prepare_only=False):
    root, output = Path(root).resolve(), Path(output).resolve()
    output.mkdir(parents=True, exist_ok=True)
    with (output/'experiment.lock').open('w') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX|fcntl.LOCK_NB)
        prompts, classes, engine, inventory = prepare(root, output, cooldown)
        print('INVENTORY '+json.dumps(inventory), flush=True)
        if prepare_only:
            return inventory
        report = dict(complete=False, new_api_calls=0, calibration=[], test={}, inventory=inventory)
        try:
            previous = {}
            for k in range(1,11):
                engine.context = dict(phase='calibration', k=k)
                rows = evaluate(prompts['calibration'], k, previous, engine.solve, classes)
                stats = dict(statistics(rows), k=k)
                assert stats['n_classes'] == 300
                save(output/f'cases-calibration-k{k:02d}.json', rows)
                report['calibration'].append(stats)
                previous = {r['prompt_index']:r for r in rows}
                table(report, output)
                print('CALIBRATION '+json.dumps(stats), flush=True)
                if all(select(report['calibration'], a) is not None for a in ALPHAS):
                    break
            selections = {str(a):select(report['calibration'], a) for a in ALPHAS}
            assert all(v is not None for v in selections.values())
            selection = dict(selections=selections, selection_data='calibration_only',
                             frozen_config_sha256=sha(output/'frozen-config.json'))
            if (output/'selected-thresholds.json').exists():
                assert read(output/'selected-thresholds.json') == selection
            save(output/'selected-thresholds.json', selection)
            finite = [v for v in selections.values() if isinstance(v,int)]
            max_k = max(finite, default=0)
            previous = {}
            for k in range(1, max_k+1):
                engine.context = dict(phase='test', k=k)
                rows = evaluate(prompts['test'], k, previous, engine.solve, classes, complete_issuance=True)
                stats = dict(statistics(rows), k=k)
                assert stats['n_classes'] == 826
                save(output/f'cases-test-k{k:02d}.json', rows)
                report['test'][str(k)] = stats
                previous = {r['prompt_index']:r for r in rows}
                table(report, output)
                print('TEST '+json.dumps(stats), flush=True)
            report.update(complete=True, new_synthesis_jobs_total=sum(
                not read(p).get('reused_from') for p in (output/'cases').glob('*/result.json')))
            table(report, output)
            save(output/'progress.json', dict(stage='complete', synthesis_process_active=False))
            print('FINAL '+json.dumps(report['table']), flush=True)
        except BaseException as error:
            report['interruption_or_error'] = type(error).__name__+': '+str(error)
            table(report, output)
            save(output/'progress.json', dict(engine.context, stage='stopped', synthesis_process_active=False,
                 interruption_or_error=report['interruption_or_error']))
            raise
        return report


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, default=Path(__file__).resolve().parent.parent)
    parser.add_argument('--output', type=Path)
    parser.add_argument('--cooldown', type=float, default=60)
    parser.add_argument('--prepare-only', action='store_true')
    args = parser.parse_args()
    if args.cooldown < 0:
        parser.error('Cooldown must be nonnegative')
    run(args.root, args.output or args.root/'Results/deepseek-r10-crc-alpha-grid-cal300-current-v1',
        args.cooldown, args.prepare_only)
