#!/usr/bin/env python3
"""Resume a frozen experiment with an audited operational cooldown override.

Keep the frozen runner, solver budgets, and checkpoint identities intact.
Only the ordinary post-synthesis rest changes; long rests and memory waits
remain controlled by the original runner.
"""
import argparse
from datetime import datetime, timezone
import math
from pathlib import Path

import run as runner


def engine_with_cooldown(seconds):
    class RuntimeCooldownEngine(runner.Engine):
        def __init__(self, root, output, config, *args, **kwargs):
            super().__init__(root, output, config, *args, **kwargs)
            # Do not mutate the frozen configuration object.
            self.config = dict(config, cooldown_seconds=seconds)
            policy = dict(
                started_at=datetime.now(timezone.utc).isoformat(),
                frozen_config_sha256=runner.sha(output / 'frozen-config.json'),
                launcher_sha256=runner.sha(Path(__file__)),
                frozen_cooldown_seconds=config['cooldown_seconds'],
                effective_cooldown_seconds=seconds,
                long_rest_seconds=300,
                long_rest_after_active_seconds=1800,
            )
            runner.save(output / 'runtime-policy.json', policy)
            print('RUNTIME_POLICY ' + runner.json.dumps(policy), flush=True)

    return RuntimeCooldownEngine


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, default=Path(__file__).resolve().parent.parent)
    parser.add_argument('--output', type=Path)
    parser.add_argument('--cooldown', type=float, required=True)
    args = parser.parse_args()
    if not math.isfinite(args.cooldown) or args.cooldown < 0:
        parser.error('Cooldown must be finite and nonnegative')
    root = args.root.resolve()
    output = (args.output or root / 'Results/deepseek-r10-crc-alpha-grid-cal300-current-v1').resolve()
    frozen = runner.read(output / 'frozen-config.json')
    runner.Engine = engine_with_cooldown(args.cooldown)
    runner.run(root, output, cooldown=frozen['cooldown_seconds'])


if __name__ == '__main__':
    main()
