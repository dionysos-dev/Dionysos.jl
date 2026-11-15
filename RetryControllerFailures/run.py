#!/usr/bin/env python3
"""Retry saved computational failures under the production thread-safety fix.

No LLM requests. Existing experiments are read-only. Checkpoints support safe
pause/resume; all resulting risk figures are labelled partial-rerun diagnostics.
"""
import argparse
import collections
import copy
from fractions import Fraction
import hashlib
import importlib.util
import json
from pathlib import Path
import subprocess


def read(path):
    return json.loads(Path(path).read_text())


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def save(path, data):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + '.tmp')
    tmp.write_text(json.dumps(data, indent=2, ensure_ascii=False) + '\n')
    tmp.replace(path)


def overlay_rows(rows, replacements, solver_identity):
    merged = []
    for row in rows:
        row = copy.deepcopy(row)
        if row.get('synthesis_outcome') == 'computational_failure':
            old_id = row['controller_job_id']
            result = replacements[old_id]
            if row['implies_reference'] is not False or row['k'] != 1:
                raise ValueError('Expected a non-implying k=1 source case')
            if bool(result['issued']) != (result['outcome'] == 'controller_issued'):
                raise ValueError('Invalid retry issuance flag')
            if result['issued'] and not (
                    result.get('region_known') and result.get('x_win_size', 0) > 0
                    and result.get('certificate', {}).get('satisfied')):
                raise ValueError('Recovered controller lacks a verified nonempty initial region')
            raw = int(result['issued'])
            row.update(tilde_g=raw, g=raw, loss=raw,
                synthesis_outcome=result['outcome'], synthesis_subreason=result.get('subreason'),
                controller_job_id=result['job_id'],
                previous_controller_job_id=old_id, retry_solver_identity=solver_identity,
                artifact_base_directory=result.get('artifact_base_directory'),
                controller_file=result.get('controller_file'), region_file=result.get('region_file'),
                x_win_size=result.get('x_win_size'), region_known=result.get('region_known', False))
        merged.append(row)
    return merged


def run(root, output, workers=2):
    import fcntl
    root, output = Path(root).resolve(), Path(output).resolve()
    if workers < 1:
        raise ValueError('workers must be positive')
    output.mkdir(parents=True, exist_ok=True)
    with (output / 'experiment.lock').open('w') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        source = root / 'Results/deepseek-r10-crc-original-alpha0.20'
        original = read(source / 'summary.json')
        if not original['complete'] or original['selected_k'] != 1:
            raise ValueError('Need the completed source k=1 experiment')
        engine = root / 'ControllerCalibration/run.py'
        if digest(engine) != original['config']['code_sha256']:
            raise ValueError('Source calibration engine changed')
        spec = importlib.util.spec_from_file_location('retry_crc_engine', engine)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        old_jobs = {p.stem: read(p) for p in (source / 'job-results').glob('*.json')}
        failures = [r for r in old_jobs.values() if r['outcome'] == 'computational_failure']
        failures.sort(key=lambda r: (r['subreason'] != 'resource_limit', len(r['formula']), r['job_id']))
        case_paths = [source / f'cases-{name}-k01.json' for name in ('calibration', 'test')]
        rows = sorted([r for p in case_paths for r in read(p)], key=lambda r: r['prompt_index'])
        if len(rows) != 4736 or len({r['prompt_index'] for r in rows}) != 4736:
            raise ValueError('Source cases do not cover the full benchmark')
        failed_rows = [r for r in rows if r.get('synthesis_outcome') == 'computational_failure']
        if {r['controller_job_id'] for r in failed_rows} != {r['job_id'] for r in failures}:
            raise ValueError('Prompt cases and failed jobs disagree')
        solver = copy.deepcopy(original['config']['solver'])
        code = {name: digest(root / name) for name in solver['code_sha256']}
        changes = [name for name, sha in code.items() if sha != solver['code_sha256'][name]]
        if changes != ['src/LTL/LTL2Automata.jl']:
            raise ValueError('Expected only the production race fix; changed files: ' + repr(changes))
        fix = read(root / 'SolverDiagnostics/2026-10-28/production-fix.json')
        if code[changes[0]] != fix['after_sha256'][changes[0]]:
            raise ValueError('Production source differs from the tested fix')
        solver['code_sha256'] = code
        version = subprocess.check_output(['ltl2tgba', '--version'], text=True).splitlines()[0]
        if version != solver['spot_version']:
            raise ValueError('Spot version changed')
        config = copy.deepcopy(original['config'])
        config.update(workers=workers, solver=solver, solver_identity=module.signature(solver),
            code_sha256=digest(__file__))
        if digest(config['plant_cache']) != config['plant_sha256']:
            raise ValueError('Plant cache changed')
        new_split = root / 'Benchmark/TaskSpecificSplits/semantic-cal300-seed-20261028'
        manifest = read(new_split / 'manifest.json')
        for name, sha in manifest['file_sha256'].items():
            if digest(new_split / name) != sha:
                raise ValueError('Resplit file changed: ' + name)
        frozen = dict(config=config, source_summary_sha256=digest(source / 'summary.json'),
            source_case_sha256={str(p): digest(p) for p in case_paths},
            source_failure_sha256={r['job_id']: digest(source / 'job-results' / (r['job_id'] + '.json')) for r in failures},
            split_manifest_sha256=digest(new_split / 'manifest.json'),
            unique_failures=len(failures), affected_prompts=len(failed_rows),
            policy='retry computational failures only; preserve original limits and all nonfailed source cases',
            interpretation='partial solver rerun diagnostics, not a fresh solver-wide CRC calibration')
        frozen_path = output / 'frozen-config.json'
        if frozen_path.exists() and read(frozen_path) != frozen:
            raise ValueError('Frozen retry experiment changed; use a fresh output directory')
        save(frozen_path, frozen)
        save(output / 'retry-plan.json', [dict(source_job_id=r['job_id'], formula=r['formula'],
            old_subreason=r['subreason'], old_last_stage=r.get('last_stage'),
            affected_prompt_indices=[p['prompt_index'] for p in failed_rows if p['controller_job_id'] == r['job_id']])
            for r in failures])
        report = dict(complete=False, unique_failures=len(failures), affected_prompts=len(failed_rows),
            new_api_calls=0, solver_identity=config['solver_identity'],
            interpretation=frozen['interpretation'])
        save(output / 'summary.json', report)
        save(output / 'progress.json', dict(stage='initializing', completed_jobs=0,
            total_jobs=len(failures), outcomes={}))
        (output / 'README.md').write_text('''# Retry computationally failed controller cases

Uses the tested production Büchi thread-safety fix. Same 120-second job timeout,
2-million-node / 40-million-edge parity caps, 1,200-step illustrative simulation,
robust plant semantics, and all-admissible-initial-state projection as before.
Acceptance verification remains infinite-horizon. No LLM/API calls.

`job-results/` checkpoints every completed distinct job. `progress.json` reports
progress. `summary.json` has `complete: true` only after all retries and report
generation finish. `comparison.md` then reports transitions from old outcomes.
Stop the runner with SIGINT to pause; rerun the same command to resume.

These are failure-only retry diagnostics. Nonfailed original cases are retained
and labelled by provenance. In particular, old negative winning-set decisions
have not been rerun, and the new risk numbers are not a fresh calibration of
the entire patched solver. Old controller artifacts and experiments are kept.
''')
        backend = module.Backend(root, output, config, workers)
        backend.old.clear()  # This experiment must genuinely rerun every failed job.
        backend.stage = 'retry-computational-failures-k1'
        try:
            outcomes = backend.solve_many([r['formula'] for r in failures])
            replacements = {r['job_id']: outcomes[r['formula']] for r in failures}
            merged = overlay_rows(rows, replacements, config['solver_identity'])
            unique_transitions = collections.Counter((r['subreason'], replacements[r['job_id']]['outcome'],
                replacements[r['job_id']].get('subreason')) for r in failures)
            prompt_transitions = collections.Counter((r['synthesis_subreason'], replacements[r['controller_job_id']]['outcome'],
                replacements[r['controller_job_id']].get('subreason')) for r in failed_rows)
            report.update(unique_outcomes=dict(collections.Counter(r['outcome'] for r in replacements.values())),
                prompt_outcomes=dict(collections.Counter(replacements[r['controller_job_id']]['outcome'] for r in failed_rows)),
                unique_transitions=[dict(old_subreason=a, new_outcome=b, new_subreason=c, count=n) for (a,b,c),n in unique_transitions.items()],
                prompt_transitions=[dict(old_subreason=a, new_outcome=b, new_subreason=c, count=n) for (a,b,c),n in prompt_transitions.items()],
                partial_rerun_diagnostics={})
            save(output / 'cases-full-k01.json', merged)
            for name in ('calibration', 'test'):
                ids = set(manifest['splits'][name]['record_ids'])
                subset = [r for r in merged if r['record_id'] in ids]
                save(output / f'cases-{name}-k01.json', subset)
                stats = module.statistics(subset)
                stats['k1_meets_alpha'] = Fraction(stats['class_mean_corrected_exact']) <= Fraction(config['alpha']) if name == 'calibration' else stats['class_mean_empirical_risk'] <= float(Fraction(config['alpha']))
                report['partial_rerun_diagnostics'][name] = stats
            lines = ['# Failure-case reevaluation', '',
                'Only old computational failures were rerun with the production race fix; original limits retained.', '',
                '| Previous reason | New outcome | New reason | Distinct jobs | Prompt cases |',
                '|---|---|---|---:|---:|']
            for (a,b,c),n in sorted(unique_transitions.items(),key=lambda item:repr(item[0])):
                lines.append(f'| {a} | {b} | {c} | {n} | {prompt_transitions[(a,b,c)]} |')
            lines += ['', 'No fresh threshold is selected by this failure-only diagnostic. Old nonfailed cases are retained; old negative winning regions require separate reassessment for a full patched-solver evaluation.', '',
                'Reassigned 300-class calibration / 826-class test diagnostics:', '', '```json',
                json.dumps(report['partial_rerun_diagnostics'],indent=2), '```', '']
            (output / 'comparison.md').write_text('\n'.join(lines))
            report['complete'] = True
            save(output / 'progress.json', dict(stage='complete',completed_jobs=len(failures),total_jobs=len(failures),outcomes=report['unique_outcomes']))
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
    run(args.root, args.output or args.root / 'Results/deepseek-r10-computational-failures-threadsafe-v1', args.workers)
