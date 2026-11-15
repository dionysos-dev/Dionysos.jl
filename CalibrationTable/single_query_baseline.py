#!/usr/bin/env python3
"""Evaluate query 1 on the frozen test split, without making LLM calls."""
import argparse
import collections
import fcntl
from fractions import Fraction
import json
from pathlib import Path
import subprocess

import run as core


def normalize(formulas, level):
    formulas = list(dict.fromkeys(formulas))
    result = subprocess.run(['ltlfilt', f'-r{level}', '--format=%f'],
                            input='\n'.join(f.replace('\n', ' ') for f in formulas)+'\n',
                            text=True, capture_output=True, check=True, timeout=300)
    values = result.stdout.splitlines()
    assert len(values) == len(formulas) and not result.stderr.strip()
    return dict(zip(formulas, values))


def select_first(analysis, raw, test_rows):
    raw = {p['prompt_index']: p for p in raw['prompts']}
    test = {r['prompt_index']: r for r in test_rows}
    rows, prompts = [], {}
    for p in analysis['prompt_results']:
        index = p['prompt_index']
        if index not in test:
            continue
        source = raw[index]
        samples = [s for s in source['samples'] if s['query_index'] == 1]
        assert len(samples) == 1
        sample = samples[0]
        cluster = next((c for c in p['ranked_clusters'] if 1 in c['sample_indices']), None)
        excluded = [e for e in p['excluded_candidates'] if e['repetition'] == 1]
        assert bool(cluster) != bool(excluded)
        row = dict(prompt_index=index, record_id=p['record_id'],
                   class_id=test[index]['class_id'], input_field=p['input_field'],
                   query_index=1, formula=sample.get('formula'), reference_ltl=p['reference_ltl'],
                   equivalent_to_reference=None, implies_reference=None,
                   issued=None, equivalence_loss=None, inclusion_loss=None,
                   status='pending_controller')
        if cluster:
            assert cluster['representative'] == sample['formula']
            assert not cluster['equivalent_to_reference'] or cluster['implies_reference']
            row.update(equivalent_to_reference=cluster['equivalent_to_reference'],
                       implies_reference=cluster['implies_reference'])
            # Engine's reference-equivalence shortcut must refer to query 1,
            # rather than the frequency-ranked first cluster.
            prompts[index] = dict(p, ranked_clusters=[cluster])
        else:
            row.update(issued=0, equivalence_loss=0, inclusion_loss=0,
                       status='abstain_invalid_or_unsatisfiable', exclusions=excluded)
        rows.append(row)
    assert len(rows) == len(test_rows)
    return rows, prompts


def assign(row, result, provenance, reuse_kind):
    issued = int(result['issued'])
    row.update(issued=issued,
               equivalence_loss=int(issued and not row['equivalent_to_reference']),
               inclusion_loss=int(issued and not row['implies_reference']),
               status='complete', outcome=result['outcome'], subreason=result.get('subreason'),
               result_file=provenance, reuse_kind=reuse_kind,
               artifact_formula=result['formula'])


def statistics(rows):
    groups = collections.Counter(r['class_id'] for r in rows)
    bounds = {name: [Fraction(), Fraction()] for name in ('equivalence_risk', 'inclusion_risk', 'issuance')}
    for r in rows:
        weight = Fraction(1, len(groups)*groups[r['class_id']])
        for name, field, can_err in (
                ('equivalence_risk', 'equivalence_loss', not r['equivalent_to_reference']),
                ('inclusion_risk', 'inclusion_loss', not r['implies_reference']),
                ('issuance', 'issued', True)):
            value = r[field]
            if value is None:
                bounds[name][1] += weight*int(can_err)
            else:
                bounds[name][0] += weight*value
                bounds[name][1] += weight*value
    return dict(n_classes=len(groups), n_prompts=len(rows),
                pending_prompts=sum(r['issued'] is None for r in rows),
                **{name: dict(lower=float(lo), upper=float(hi), lower_exact=str(lo), upper_exact=str(hi))
                   for name, (lo, hi) in bounds.items()})


def write_report(output, rows, framework, inventory):
    stats = statistics(rows)
    complete = stats['pending_prompts'] == 0
    report = dict(complete=complete, new_api_calls=0, sample_selection='query_index == 1, no replacement',
                  risk_unit='uniform semantic class then uniform test prompt within class',
                  invalid_policy='abstain; do not substitute a later query',
                  statistics=stats, inventory=inventory)
    core.save(output/'summary.json', report)
    core.save(output/'test-rows.json', rows)
    if not complete:
        return
    eq, inc, issuance = (100*stats[n]['lower'] for n in ('equivalence_risk', 'inclusion_risk', 'issuance'))
    table = [dict(method='Single LLM query (equivalence)', guarantee='No', alpha='--', k='--', risk=eq, issuance=issuance),
             dict(method='Single LLM query (language inclusion)', guarantee='No', alpha='--', k='--', risk=inc, issuance=issuance)]
    for r in framework['table']:
        if r['alpha'] == 0.05:
            continue
        table.append(dict(method='NL2RCtrl', guarantee='Yes', alpha=f"{r['alpha']:.2f}",
                          k=r['selected_k'], risk=100*r['test_risk'], issuance=100*r['issuance']))
    report['table'] = table
    core.save(output/'summary.json', report)
    markdown = ['All percentages use uniform semantic-class means; 300 calibration classes and 826 test classes.',
                'The baseline uses query 1 exactly. No new LLM calls.', '',
                '| Method | Statistical guarantee | α | k | Test risk (%) | Issuance (%) |',
                '|---|---|---:|---:|---:|---:|']
    for r in table:
        markdown.append(f"| {r['method']} | {r['guarantee']} | {r['alpha']} | {r['k']} | {r['risk']:.2f} | {r['issuance']:.2f} |")
    (output/'table.md').write_text('\n'.join(markdown)+'\n')
    latex = r'''\begin{table}[t]
\centering
\caption{Empirical test risk and controller issuance
for the single LLM query baseline and the proposed
NL2RCtrl framework at different risk tolerances $\alpha$.}
\label{tab:crc_results}
\small
\setlength{\tabcolsep}{3pt}
\begin{tabular}{lccccc}
\toprule
Method & Statistical guarantee & $\alpha$ & $\widehat{k}$
& Test risk (\%) & Issuance (\%) \\
\midrule
'''
    latex += f'Single LLM query (equivalence) & No & -- & -- & {eq:.2f} & {issuance:.2f} '+r'\\'+'\n'
    latex += f'Single LLM query (language inclusion) & No & -- & -- & {inc:.2f} & {issuance:.2f} '+r'\\'+'\n'
    latex += '\\midrule\n'
    for i, r in enumerate(table[2:]):
        alpha, k = r['alpha'], r['k']
        k = r'$\infty$' if k == 'infinity' else str(k)
        risk, rate = f"{r['risk']:.2f}", f"{r['issuance']:.2f}"
        if alpha == '0.15':
            alpha, k, risk, rate = [r'\textbf{'+v+'}' for v in (alpha, k, risk, rate)]
        latex += f"{'NL2RCtrl' if i == 0 else ''} & {'Yes' if i == 0 else ''} & {alpha} & {k} & {risk} & {rate} "+r'\\'+'\n'
    latex += '\\bottomrule\n\\end{tabular}\n\\end{table}\n'
    (output/'table.tex').write_text(latex)


def run(root, output, cooldown, prepare_only=False):
    output.mkdir(parents=True, exist_ok=True)
    with (output/'experiment.lock').open('w') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        source = root/'Results/deepseek-r10-crc-alpha-grid-cal300-current-v1'
        framework = core.read(source/'summary.json')
        assert framework['complete']
        frozen = core.read(source/'frozen-config.json')
        for name, digest in frozen['code_sha256'].items():
            assert core.sha(root/name) == digest, name
        for name, hash_name in (('plant_cache', 'plant_sha256'), ('plant_binary', 'plant_binary_sha256')):
            assert core.sha(frozen[name]) == frozen[hash_name]
        analysis_path = root/'Results/deepseek-v4-pro-task-specific-full-T0.7-K10-conjunctions/summary.json'
        analysis = core.read(analysis_path)
        assert core.sha(analysis_path) == frozen['source_sha256']
        assert core.sha(analysis['source_path']) == analysis['analysis_identity']['source_sha256']
        test_rows = core.read(source/'cases-test-k01.json')
        rows, prompts = select_first(analysis, core.read(analysis['source_path']), test_rows)
        assert len(rows) == 3480 and len({r['class_id'] for r in rows}) == 826
        sources, refs, hashes = [], {}, {}
        for path in sorted((source/'cases').glob('*/result.json')):
            result = core.read(path)
            # A prior exact-language alias is also a reusable cache key.
            for formula in dict.fromkeys((result['formula'], result['requested_formula'])):
                alias = dict(result, formula=formula)
                sources.append((alias, str(path)))
            hashes[str(path)] = core.sha(path)
        ref_path = root/'Results/parity-limit-continuation-4m-60m-v1/reference-results.json'
        for row in core.read(ref_path):
            result = core.absolute_result(row['synthesis'], ref_path.parent)
            sources.append((result, str(ref_path)))
            for rid in row['group']['record_ids']:
                refs[rid] = (result, str(ref_path))
        # Match the completed experiment's accepted historical source policy.
        for name, positive_only in [('deepseek-r10-crc-original-alpha0.20', True),
                                    ('deepseek-r10-computational-failures-threadsafe-v1', False)]:
            directory = root/'Results'/name
            for path in sorted((directory/'job-results').glob('*.json')):
                result = core.read(path)
                if result['outcome'] == 'computational_failure' or (positive_only and not result['issued']):
                    continue
                if not result.get('artifact_base_directory'):
                    result['artifact_base_directory'] = core.read(directory/'jobs'/path.name)['worker_output']
                sources.append((core.absolute_result(result, directory), str(path)))
        for path, digest in frozen['reused_source_sha256'].items():
            assert core.sha(path) == digest, path
        reusable_formulas = [r['formula'] for r in rows if r['issued'] is None] + [r['formula'] for r, _ in sources]
        formulas = list(core.read(source/'normalization.json')) + reusable_formulas
        norm1, norm3 = normalize(formulas, 1), normalize(reusable_formulas, 3)
        lookups = []
        for norm in (norm1, norm3):
            lookup = {}
            for result, path in sources:
                key = norm[result['formula']]
                old = lookup.get(key)
                if old is None or (old[0]['outcome'] == 'computational_failure' and result['outcome'] != 'computational_failure'):
                    lookup[key] = (result, path)
            lookups.append(lookup)
        config = dict(solver_config=dict(frozen, cooldown_seconds=cooldown),
                      framework_config_sha256=core.sha(source/'frozen-config.json'),
                      framework_summary_sha256=core.sha(source/'summary.json'),
                      test_rows_sha256=core.sha(source/'cases-test-k01.json'),
                      analysis_sha256=core.sha(analysis_path), runner_sha256=core.sha(Path(__file__)),
                      cache_sha256=hashes, query_index=1, invalid_policy='abstain_no_replacement',
                      equivalence_cache='Spot r3 identical simplifications; exact language equivalence',
                      new_api_calls=0)
        if (output/'frozen-config.json').exists():
            assert core.read(output/'frozen-config.json') == config, 'Baseline configuration changed'
        core.save(output/'frozen-config.json', config)
        engine = core.Engine(root, output, config['solver_config'], formulas, sources, refs)
        validated, pending, reuse_counts = set(), [], collections.Counter()
        for row in rows:
            if row['issued'] is not None:
                continue
            formula = row['formula']
            own = output/'cases'/engine.jid(formula)/'result.json'
            match = None
            if own.exists():
                match, kind = (core.read(own), str(own)), 'baseline_checkpoint'
            elif norm1[formula] in lookups[0]:
                match, kind = lookups[0][norm1[formula]], 'normalized_formula'
            elif row['equivalent_to_reference']:
                match, kind = refs[row['record_id']], 'reference_equivalence'
            elif norm3[formula] in lookups[1]:
                match, kind = lookups[1][norm3[formula]], 'spot_r3_equivalence'
            if match:
                result, path = match
                signature = (path, result.get('job_id'), result.get('controller_sha256'))
                if signature not in validated:
                    core.validate_result(result)
                    validated.add(signature)
                assign(row, result, path, kind)
                reuse_counts[kind] += 1
            else:
                pending.append(row)
        inventory = dict(test_prompts=len(rows), test_classes=826, reuse_counts=dict(reuse_counts),
                         invalid_or_unsatisfiable=sum(r['status']=='abstain_invalid_or_unsatisfiable' for r in rows),
                         pending_prompts=len(pending), missing_distinct_formulas=len({norm1[r['formula']] for r in pending}))
        core.save(output/'inventory.json', inventory)
        write_report(output, rows, framework, inventory)
        print('INVENTORY '+json.dumps(inventory), flush=True)
        if prepare_only:
            return
        try:
            for i, row in enumerate(pending, 1):
                engine.context = dict(phase='single_query_test', pending_index=i, pending_total=len(pending))
                result = engine.solve(row['formula'], prompts[row['prompt_index']], 1)
                assign(row, result, result['checkpoint_file'], 'baseline_synthesis')
                write_report(output, rows, framework, inventory)
            core.save(output/'progress.json', dict(stage='complete', synthesis_process_active=False))
            print('FINAL '+json.dumps(statistics(rows)), flush=True)
        except BaseException as error:
            write_report(output, rows, framework, inventory)
            core.save(output/'progress.json', dict(engine.context, stage='stopped', synthesis_process_active=False,
                                                   interruption_or_error=type(error).__name__+': '+str(error)))
            raise


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, default=Path(__file__).resolve().parent.parent)
    parser.add_argument('--output', type=Path)
    parser.add_argument('--cooldown', type=float, default=30)
    parser.add_argument('--prepare-only', action='store_true')
    args = parser.parse_args()
    if not 0 <= args.cooldown < float('inf'):
        parser.error('Cooldown must be finite and nonnegative')
    root = args.root.resolve()
    run(root, (args.output or root/'Results/deepseek-single-query-baseline-cal300-v1').resolve(), args.cooldown, args.prepare_only)
