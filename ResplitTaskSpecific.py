#!/usr/bin/env python3
"""Resplit audited LTL classes and reassign saved controller-aware k=1 losses.

Uses no model requests, Spot calls, or controller synthesis. Original datasets,
results, controller files, and winning regions remain in place.
"""
import argparse
import collections
import copy
from fractions import Fraction
import hashlib
import importlib.util
import json
from pathlib import Path
import random
import tempfile


def read(path):
    return json.loads(Path(path).read_text())


def encoded(value):
    return (json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False) + '\n').encode()


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def allocation(groups, count, seed):
    if isinstance(count, bool) or not isinstance(count, int) or not 0 < count < len(groups):
        raise ValueError('Calibration count must leave both sets nonempty')
    ordered = sorted(copy.deepcopy(groups), key=lambda g: tuple(g['formulas']))
    shuffled = [g['class_id'] for g in ordered]
    if len(set(shuffled)) != len(shuffled):
        raise ValueError('Duplicate semantic class ID')
    random.Random(seed).shuffle(shuffled)
    calibration = set(shuffled[:count])
    for group in ordered:
        group['split'] = 'calibration' if group['class_id'] in calibration else 'test'
    return ordered, shuffled


def reassign(rows, groups, records, fields):
    by_record = {rid: g for g in groups for rid in g['record_ids']}
    if len(by_record) != sum(len(g['record_ids']) for g in groups):
        raise ValueError('A record appears in multiple semantic classes')
    expected = {(r['id'], field) for r in records for field in fields}
    observed = set()
    indices = set()
    result = {'calibration': [], 'test': []}
    for row in sorted(rows, key=lambda r: r['prompt_index']):
        key = (row['record_id'], row['input_field'])
        if key not in expected or key in observed or row['prompt_index'] in indices:
            raise ValueError('Unexpected or duplicate controller result')
        group = by_record[row['record_id']]
        if row['k'] != 1 or row['class_id'] != group['class_id'] or row['loss'] not in (0, 1):
            raise ValueError('Controller case identity or loss mismatch')
        if row['g'] is not None and row['implies_reference'] is not None:
            if row['loss'] != int(bool(row['g']) and not row['implies_reference']):
                raise ValueError('Saved case does not use the original issuance-weighted loss')
        elif row['g'] is None and not (row['implies_reference'] is True and row['loss'] == 0):
            raise ValueError('Unmeasured synthesis does not determine a zero loss')
        result[group['split']].append(copy.deepcopy(row))
        observed.add(key)
        indices.add(row['prompt_index'])
    if observed != expected:
        raise ValueError('Saved cases do not cover every active NL-LTL pair')
    return result


def write_new_directory(path, files):
    path = Path(path)
    if path.exists():
        raise FileExistsError('Refusing to overwrite ' + str(path))
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix='.resplit-', dir=path.parent) as staging:
        staging = Path(staging)
        for name, content in files.items():
            (staging / name).write_bytes(content)
        staging.rename(path)


def run(root, count=300, seed=20261028, split_output=None, results_output=None):
    root = Path(root).resolve()
    old_split = root / 'Benchmark/TaskSpecificSplits/semantic-70-30-seed-20261027'
    old_results = root / 'Results/deepseek-r10-crc-original-alpha0.20'
    split_output = Path(split_output or root / f'Benchmark/TaskSpecificSplits/semantic-cal{count}-seed-{seed}').resolve()
    results_output = Path(results_output or root / f'Results/deepseek-r10-crc-original-alpha0.20-cal{count}-seed-{seed}').resolve()
    for path in (split_output, results_output):
        if path.exists():
            raise FileExistsError('Refusing to overwrite ' + str(path))
    old_manifest = read(old_split / 'manifest.json')
    for name, checksum in old_manifest['file_sha256'].items():
        if digest(old_split / name) != checksum:
            raise ValueError('Source split hash mismatch: ' + name)
    dataset = root / 'Benchmark/Task_Spescific_Dataset.json'
    if digest(dataset) != old_manifest['dataset_sha256']:
        raise ValueError('Benchmark changed since the semantic equivalence audit')
    records = read(dataset)
    old_groups = read(old_split / 'semantic_groups.json')
    groups, shuffled = allocation(old_groups, count, seed)
    by_record = {rid: g for g in groups for rid in g['record_ids']}
    if len(by_record) != len(records) or set(by_record) != {r['id'] for r in records}:
        raise ValueError('Semantic grouping does not partition the benchmark records')
    for group in groups:
        if group['class_id'] != 'ltl-' + hashlib.sha256(encoded(sorted(group['formulas']))).hexdigest()[:16]:
            raise ValueError('Semantic class identity changed')
    for record in records:
        if record['LTL'] not in by_record[record['id']]['formulas']:
            raise ValueError('Reference differs from audited semantic grouping')
    fields = old_manifest['active_nl_fields']
    subsets = {name: [r for r in records if by_record[r['id']]['split'] == name]
               for name in ('calibration', 'test')}
    manifest = copy.deepcopy(old_manifest)
    manifest.update(seed=seed, calibration_class_count_requested=count,
        calibration_fraction_requested=None, test_fraction_requested=None,
        calibration_class_fraction=count / len(groups),
        rounding='none; exact requested class count', shuffled_class_ids=shuffled,
        grouping_reused_from=str(old_split / 'semantic_groups.json'),
        grouping_source_sha256=digest(old_split / 'semantic_groups.json'),
        source_manifest_sha256=digest(old_split / 'manifest.json'),
        allocation_uses_controller_results=False,
        evaluation_history='A new random resplit of previously evaluated benchmark data; not an untouched test set',
        splits={})
    for name, subset in subsets.items():
        members = [g for g in groups if g['split'] == name]
        manifest['splits'][name] = dict(file=name + '.json', semantic_class_count=len(members),
            unique_formula_strings=len({r['LTL'] for r in subset}), record_count=len(subset),
            active_nl_ltl_pair_count=len(subset) * len(fields),
            record_ids=[r['id'] for r in subset], class_ids=[g['class_id'] for g in members])
    if 'prior_development_overlap' in manifest:
        overlap = manifest['prior_development_overlap']
        dev_path = Path(overlap['manifest_path'])
        if digest(dev_path) != overlap['manifest_sha256']:
            raise ValueError('Prior development manifest changed')
        overlap['prompt_counts_by_new_split'] = dict(collections.Counter(
            by_record[p['record_id']]['split'] for p in read(dev_path)['prompts']))
    files = {name + '.json': encoded(values) for name, values in subsets.items()}
    files['semantic_groups.json'] = encoded(groups)
    manifest['file_sha256'] = {name: hashlib.sha256(content).hexdigest() for name, content in files.items()}
    files['manifest.json'] = encoded(manifest)

    # Validate the completed source experiment before importing any outcomes.
    source = read(old_results / 'summary.json')
    if not source['complete'] or source['selected_k'] != 1:
        raise ValueError('Reassignment requires complete full-benchmark k=1 cases')
    if digest(old_split / 'manifest.json') != source['config']['split_manifest_sha256']:
        raise ValueError('Controller results use a different source split')
    engine = root / 'ControllerCalibration/run.py'
    if digest(engine) != source['config']['code_sha256']:
        raise ValueError('Source experiment code changed')
    semantics_path = root / 'Results/deepseek-v4-pro-task-specific-full-T0.7-K10-conjunctions/summary.json'
    if digest(semantics_path) != source['config']['source_semantics_sha256']:
        raise ValueError('Source semantic analysis changed')
    old_case_paths = [old_results / f'cases-{name}-k01.json' for name in ('calibration', 'test')]
    rows = [row for path in old_case_paths for row in read(path)]
    cases = reassign(rows, groups, records, fields)
    semantics = {p['prompt_index']: p for p in read(semantics_path)['prompt_results']}
    for row in rows:
        prompt = semantics[row['prompt_index']]
        if (row['record_id'], row['input_field']) != (prompt['record_id'], prompt['input_field']):
            raise ValueError('Saved case differs from the source LLM prompt')
        c = prompt['conjunctions'][0] if prompt['conjunctions'] else None
        if (row['formula'], row['implies_reference'], row['conjunction_satisfiable']) != (
                c['phi_k'] if c else None, c['implies_reference'] if c else None, c['satisfiable'] if c else None):
            raise ValueError('Saved case differs from the source conjunction')
    spec = importlib.util.spec_from_file_location('saved_crc_statistics', engine)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    alpha = Fraction(source['config']['alpha'])
    calibration = dict(module.statistics(cases['calibration']), k=1)
    selected = 1 if Fraction(calibration['class_mean_corrected_exact']) <= alpha else None
    # Selection is determined from calibration only; test cannot change it.
    test = dict(module.statistics(cases['test']), k=1) if selected else None
    if test is not None:
        test['empirically_below_alpha'] = test['class_mean_empirical_risk'] <= float(alpha)
    config = copy.deepcopy(source['config'])
    config['split_manifest_sha256'] = hashlib.sha256(files['manifest.json']).hexdigest()
    report = dict(complete=selected is not None, selected_k=selected, alpha=str(alpha),
        config=config, split_manifest=str(split_output / 'manifest.json'),
        calibration=[calibration], test=test,
        source_results_directory=str(old_results),
        source_summary_sha256=digest(old_results / 'summary.json'),
        source_case_sha256={str(p): digest(p) for p in old_case_paths},
        resplit_code_sha256=digest(__file__),
        new_api_calls=0, new_synthesis_jobs=0,
        controller_artifacts='Existing controller and winning-region paths retained; source results must remain available',
        interpretation='Uniform semantic-class means of original prompt losses; prior evaluation history preserved',
        pending='none' if selected else 'Need higher-k synthesis to select a threshold; k=1 data are fully reassigned')
    counts = ['| Set | Semantic LTL classes | Records | Active NL–LTL pairs |', '|---|---:|---:|---:|']
    for name, stats in manifest['splits'].items():
        counts.append(f"| {name.title()} | {stats['semantic_class_count']} | {stats['record_count']} | {stats['active_nl_ltl_pair_count']} |")
    common = '\n'.join(counts)
    files['README.md'] = (f'# Task-specific random split: {count} calibration classes\n\n'
        f'Uniform shuffle without replacement, seed {seed}; first {count} semantic classes assigned to calibration. '
        'No outcome-based balancing or seed search. Exact Spot equivalence groups reused from the hash-verified prior audit. '
        'Every equivalent formula and all NL representations stay together; all original record fields are preserved.\n\n'
        + common + '\n\nA new random resplit of previously evaluated data, not an untouched test set. '
        'The old split and experiment are preserved. Saved controller results are reassigned by record ID and input field.\n\n'
        f'Reassigned losses and risks: `{results_output}`.\n').encode()
    result_files = {f'cases-{name}-k01.json': encoded(values) for name, values in cases.items()}
    result_files['summary.json'] = encoded(report)
    if selected:
        result_files['selected-threshold.json'] = encoded(dict(k=selected, alpha=str(alpha),
            calibration=calibration, selection_data='new calibration only',
            split_manifest_sha256=config['split_manifest_sha256']))
    result_files['comparison.md'] = (f'# Reassigned controller-aware CRC results\n\nSeed: {seed}. '
        'No additional synthesis or LLM requests. All saved k=1 cases reassigned, with controller paths retained.\n\n'
        + common + '\n\n'
        f"Calibration CRC-corrected class-mean risk at k=1: {100 * calibration['class_mean_corrected_risk']:.4f}%.\n\n"
        f'Selected threshold: {selected}. Target alpha: {float(alpha):.0%}.\n\n'
        + (f"Test empirical class-mean risk: {100 * test['class_mean_empirical_risk']:.4f}%.\n\n" if test else '')
        + 'Loss: cumulative controller issuance multiplied by non-implication. Synthesis skipped for '
        'logically guaranteed formulas remains unmeasured (null g); their original losses are exactly zero. '
        'This report does not measure overall issuance coverage.\n\n'
        'This is a new random resplit of previously evaluated data; it is not an untouched test set. '
        'Source experiment outcomes, including timeout and resource-limit abstentions, are reused unchanged.\n').encode()
    write_new_directory(split_output, files)
    write_new_directory(results_output, result_files)
    return dict(split_output=str(split_output), results_output=str(results_output), seed=seed,
        splits={name: {key: value for key, value in stats.items() if key not in ('record_ids', 'class_ids')}
                for name, stats in manifest['splits'].items()},
        selected_k=selected, calibration=calibration, test=test,
        new_api_calls=0, new_synthesis_jobs=0)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, default=Path(__file__).resolve().parent)
    parser.add_argument('--calibration-count', type=int, default=300)
    parser.add_argument('--seed', type=int, default=20261028)
    parser.add_argument('--split-output', type=Path)
    parser.add_argument('--results-output', type=Path)
    args = parser.parse_args()
    print(json.dumps(run(args.root, args.calibration_count, args.seed,
                         args.split_output, args.results_output), indent=2))
