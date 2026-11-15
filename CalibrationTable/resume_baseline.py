#!/usr/bin/env python3
"""Resume the frozen single-query baseline with shorter operational rests."""
import argparse
from datetime import datetime, timezone
import math
from pathlib import Path

import single_query_baseline as baseline


def engine_with_rests(cooldown, long_rest):
    class RuntimeRestEngine(baseline.core.Engine):
        def __init__(self, root, output, config, *args, **kwargs):
            super().__init__(root, output, config, *args, **kwargs)
            policy = dict(
                started_at=datetime.now(timezone.utc).isoformat(),
                frozen_config_sha256=baseline.core.sha(output/'frozen-config.json'),
                launcher_sha256=baseline.core.sha(Path(__file__)),
                frozen_cooldown_seconds=config['cooldown_seconds'],
                effective_cooldown_seconds=cooldown,
                long_rest_seconds=long_rest,
                long_rest_after_active_seconds=1800,
                memory_waits='unchanged',
            )
            baseline.core.save(output/'runtime-policy.json', policy)
            print('RUNTIME_POLICY '+baseline.json.dumps(policy), flush=True)

        def wait(self, seconds, reason):
            if reason == 'scheduled rest after new synthesis':
                seconds = long_rest if seconds == 300 else cooldown
            return super().wait(seconds, reason)

    return RuntimeRestEngine


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, default=Path(__file__).resolve().parent.parent)
    parser.add_argument('--output', type=Path)
    parser.add_argument('--cooldown', type=float, default=5)
    parser.add_argument('--long-rest', type=float, default=30)
    args = parser.parse_args()
    if any(not math.isfinite(v) or v < 0 for v in (args.cooldown, args.long_rest)):
        parser.error('Rest durations must be finite and nonnegative')
    root = args.root.resolve()
    output = (args.output or root/'Results/deepseek-single-query-baseline-cal300-v1').resolve()
    frozen = baseline.core.read(output/'frozen-config.json')
    baseline.core.Engine = engine_with_rests(args.cooldown, args.long_rest)
    baseline.run(root, output, frozen['solver_config']['cooldown_seconds'])


if __name__ == '__main__':
    main()
