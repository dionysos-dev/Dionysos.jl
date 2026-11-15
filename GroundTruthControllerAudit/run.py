#!/usr/bin/env python3
"""Audit controller existence for reference LTLs, never generated candidates.

One representative per exact Spot-audited semantic class. Same frozen plant,
production solver and limits as the patched failure retries. Computational
failures remain unknown, not unrealizability. Checkpointed and resumable.
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
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def save(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + '.tmp')
    temporary.write_text(json.dumps(value, indent=2, ensure_ascii=False) + '\n')
    temporary.replace(path)


def outcome_counts(classes, results):
    by_class, by_record = collections.Counter(), collections.Counter()
    for group in classes:
        result = results.get(group['class_id'])
        outcome = result['outcome'] if result else 'pending'
        by_class[outcome] += 1
        by_record[outcome] += len(group['record_ids'])
    total_classes, total_records = len(classes), sum(len(g['record_ids']) for g in classes)
    return dict(semantic_classes=total_classes, benchmark_records=total_records,
        class_counts=dict(by_class), record_counts=dict(by_record),
        class_rates={key: count / total_classes for key, count in by_class.items()},
        record_rates={key: count / total_records for key, count in by_record.items()},
        controller_existence_lower_bound=by_class['controller_issued'] / total_classes,
        controller_existence_upper_bound=(by_class['controller_issued'] + by_class['computational_failure'] + by_class['pending']) / total_classes,
        computational_failures_are_unknown=True)


def run(root, output, workers=2):
    import fcntl
    root, output = Path(root).resolve(), Path(output).resolve()
    if workers < 1:
        raise ValueError('workers must be positive')
    output.mkdir(parents=True, exist_ok=True)
    with (output / 'experiment.lock').open('w') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        original_dir = root / 'Results/deepseek-r10-crc-original-alpha0.20'
        original = read(original_dir / 'summary.json')
        engine = root / 'ControllerCalibration/run.py'
        if digest(engine) != original['config']['code_sha256']:
            raise ValueError('Synthesis runner changed')
        spec = importlib.util.spec_from_file_location('ground_truth_crc_backend', engine)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        split = root / 'Benchmark/TaskSpecificSplits/semantic-cal300-seed-20261028'
        manifest = read(split / 'manifest.json')
        for name, checksum in manifest['file_sha256'].items():
            if digest(split / name) != checksum:
                raise ValueError('Audited semantic split changed: ' + name)
        dataset = root / 'Benchmark/Task_Spescific_Dataset.json'
        if digest(dataset) != manifest['dataset_sha256']:
            raise ValueError('Ground-truth benchmark changed')
        groups = read(split / 'semantic_groups.json')
        records = read(dataset)
        by_record = {r['id']: r for r in records}
        if {rid for g in groups for rid in g['record_ids']} != set(by_record):
            raise ValueError('Reference groups do not cover the benchmark')
        for group in groups:
            if group['representative'] not in group['formulas']:
                raise ValueError('Representative is not a reference formula')
            if any(by_record[rid]['LTL'] not in group['formulas'] for rid in group['record_ids']):
                raise ValueError('Group contains a mismatched ground truth')
        solver = copy.deepcopy(original['config']['solver'])
        solver['code_sha256'] = {name: digest(root / name) for name in solver['code_sha256']}
        fix = read(root / 'SolverDiagnostics/2026-10-28/production-fix.json')
        if solver['code_sha256']['src/LTL/LTL2Automata.jl'] != fix['after_sha256']['src/LTL/LTL2Automata.jl']:
            raise ValueError('Expected the tested production thread-safety fix')
        version = subprocess.check_output(['ltl2tgba', '--version'], text=True).splitlines()[0]
        if version != solver['spot_version']:
            raise ValueError('Spot version differs from the prior experiment')
        config = copy.deepcopy(original['config'])
        config.update(workers=workers, solver=solver, solver_identity=module.signature(solver),
            code_sha256=digest(__file__))
        if digest(config['plant_cache']) != config['plant_sha256']:
            raise ValueError('Plant cache changed')
        frozen = dict(schema_version=1, config=config,
            dataset_sha256=digest(dataset), semantic_groups_sha256=digest(split / 'semantic_groups.json'),
            semantic_class_count=len(groups), benchmark_record_count=len(records),
            specification_source='Ground-truth LTL field only, one audited semantic representative per class',
            success='verified nonempty admissible initial-state winning region',
            computational_failures='unknown controller existence',
            reuse='Only independently verified positive controllers for the identical or Spot-equivalent reference; no cached negative or timeout outcomes')
        frozen_path = output / 'frozen-config.json'
        if frozen_path.exists() and read(frozen_path) != frozen:
            raise ValueError('Frozen audit changed; use a fresh output directory')
        save(frozen_path, frozen)
        backend = module.Backend(root, output, config, workers)
        backend.old.clear()
        backend.stage = 'ground-truth-reference-controller-audit'
        report = dict(complete=False, new_api_calls=0, specification_source=frozen['specification_source'],
            solver_identity=config['solver_identity'], **outcome_counts(groups, {}))
        save(output / 'summary.json', report)
        save(output / 'progress.json', dict(stage='initializing', completed_jobs=0, total_jobs=len(groups), outcomes={}))
        (output / 'README.md').write_text('''# Ground-truth controller audit

Evaluates the reference LTLs in Benchmark/Task_Spescific_Dataset.json, one exact
Spot semantic class at a time. Equivalent formulas share the result. Reports
both class-weighted and record-weighted proportions. NL paraphrases do not
create new control problems. Same 44,643-state plant and robust semantics as
the existing experiment; every admissible plant state is a possible start.

Controller found means a verified nonempty projected initial winning region
and an exported controller. Computational failures are UNKNOWN existence,
not proofs of plant unrealizability. Infinite-horizon acceptance is verified;
only illustrative simulation is capped at 1,200 steps. The job timeout is
120 seconds; parity limits remain 2 million nodes / 40 million edges.

Uses the production thread-safety fix. Worklist and parity product reductions
are not yet enabled. No model API calls. Reuses only positive verified
controllers for exact or Spot-proved equivalent ground-truth formulas.
All new job outcomes are checkpointed; rerunning the same command resumes.
summary.json has complete=true only after the audit finishes. comparison.md
then contains the measured percentages and unknown-outcome bounds.
''')
        try:
            formulas = list(dict.fromkeys(f for g in groups for f in g['formulas']))
            backend.normalize(formulas)
            cache = {}
            source_positive_count = 0
            for directory in [original_dir, root / 'Results/deepseek-r10-computational-failures-threadsafe-v1']:
                for path in sorted((directory / 'job-results').glob('*.json')):
                    result = read(path)
                    if result.get('issued') and result.get('region_known') and result.get('x_win_size', 0) > 0 and result.get('certificate', {}).get('satisfied'):
                        cache.setdefault(result['formula'], (result, path))
            for group in groups:
                formula = backend.norm[group['representative']]
                jid = module.signature(dict(formula=formula, solver_identity=config['solver_identity']))
                target = output / 'job-results' / (jid + '.json')
                if target.exists():
                    continue
                hits = [cache[backend.norm[f]] for f in group['formulas'] if backend.norm[f] in cache]
                if not hits:
                    continue
                old, source_path = hits[0]
                # Verify exact equivalence even if a different audited member
                # or representation matched the cached positive controller.
                check = subprocess.run(['ltlfilt', '-f', formula, '--equivalent-to=' + old['formula'], '--count'],
                    capture_output=True, text=True, timeout=30)
                if check.returncode != 0 or check.stdout.strip() != '1':
                    raise ValueError('Positive artifact is not equivalent to the reference')
                result = dict(old, job_id=jid, formula=formula,
                    artifact_job_id=old.get('artifact_job_id', old['job_id']),
                    artifact_formula=old.get('artifact_formula', old['formula']),
                    reused_verified_reference_controller=True,
                    reuse_source_job=str(source_path), reuse_source_sha256=digest(source_path),
                    reference_equivalence_verified=True)
                backend.validate(result, target)
                save(target, result)
                source_positive_count += 1
            plan = [dict(class_id=g['class_id'], representative=g['representative'],
                normalized_formula=backend.norm[g['representative']], record_ids=g['record_ids']) for g in groups]
            save(output / 'reference-plan.json', plan)
            print(f'Ground truths: {len(groups)} semantic classes, {len(records)} records; {source_positive_count} positive controllers reused', flush=True)
            outcomes = backend.solve_many([g['representative'] for g in groups])
            results = {g['class_id']: outcomes[g['representative']] for g in groups}
            save(output / 'reference-results.json', [dict(group=g, synthesis=results[g['class_id']]) for g in groups])
            report.update(outcome_counts(groups, results))
            lines = ['# Ground-truth controller existence audit', '',
                'All benchmark reference formulas evaluated on the same frozen plant.', '',
                '| Outcome | Semantic classes | % of classes | Records | % of records |',
                '|---|---:|---:|---:|---:|']
            for outcome in ('controller_issued', 'plant_unrealizable', 'logically_inconsistent', 'computational_failure'):
                count = report['class_counts'].get(outcome, 0)
                nr = report['record_counts'].get(outcome, 0)
                lines.append(f"| {outcome} | {count} | {100*count/len(groups):.2f}% | {nr} | {100*nr/len(records):.2f}% |")
            lines += ['', 'Computational failures leave existence unknown. The verified controller rate is a lower bound on true controller existence.', '',
                f"Existence bounds over semantic classes: {100*report['controller_existence_lower_bound']:.2f}% to {100*report['controller_existence_upper_bound']:.2f}%.", '']
            (output / 'comparison.md').write_text('\n'.join(lines))
            report['complete'] = True
            save(output / 'progress.json', dict(stage='complete',completed_jobs=len(groups),total_jobs=len(groups),outcomes=report['class_counts']))
        except BaseException as error:
            report['interruption_or_error'] = type(error).__name__ + ': ' + str(error)
            raise
        finally:
            backend.close()
            save(output / 'summary.json', report)
        print('FINAL ' + json.dumps(report), flush=True)
        return report


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--output', type=Path)
    parser.add_argument('--workers', type=int, default=2)
    args = parser.parse_args()
    run(args.root, args.output or args.root / 'Results/ground-truth-controller-audit-threadsafe-v1', args.workers)
