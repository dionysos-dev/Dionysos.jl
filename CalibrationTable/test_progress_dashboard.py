"""The original missing-work denominator must survive evaluation resumes."""
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import progress_dashboard as dashboard


class ProgressDashboardTests(unittest.TestCase):
    def test_original_progress_includes_work_completed_before_resume(self):
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory)
            original = dict(test_prompts=100, pending_prompts=10, missing_distinct_formulas=10)
            resumed = dict(test_prompts=100, pending_prompts=6, missing_distinct_formulas=6)
            (output/'run.log').write_text('INVENTORY '+json.dumps(original)+'\nRESUME\nINVENTORY '+json.dumps(resumed)+'\n')
            (output/'inventory.json').write_text(json.dumps(resumed))
            (output/'summary.json').write_text(json.dumps(dict(complete=False,
                statistics=dict(n_prompts=100, pending_prompts=4))))
            with patch.object(dashboard, 'runner_alive', return_value=True):
                status = dashboard.ProgressReader(output).snapshot()
            self.assertEqual(status['original_missing_total'], 10)
            self.assertEqual(status['original_missing_done'], 6)
            self.assertEqual(status['completed_before_batch'], 4)
            self.assertEqual(status['batch_done'], 2)

            # Another resume keeps the original denominator while its batch shrinks.
            resumed['pending_prompts'] = 4
            resumed['missing_distinct_formulas'] = 4
            (output/'inventory.json').write_text(json.dumps(resumed))
            with (output/'run.log').open('a') as log:
                log.write('RESUME\nINVENTORY '+json.dumps(resumed)+'\n')
            (output/'summary.json').write_text(json.dumps(dict(complete=False,
                statistics=dict(n_prompts=100, pending_prompts=3))))
            with patch.object(dashboard, 'runner_alive', return_value=True):
                status = dashboard.ProgressReader(output).snapshot()
            self.assertEqual(status['original_missing_total'], 10)
            self.assertEqual(status['original_missing_done'], 7)
            self.assertEqual(status['completed_before_batch'], 6)
            self.assertEqual(status['batch_done'], 1)

    def test_complete_counts_all_original_jobs(self):
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory)
            (output/'run.log').write_text('INVENTORY '+json.dumps(dict(pending_prompts=10))+'\n')
            (output/'inventory.json').write_text(json.dumps(dict(pending_prompts=6)))
            (output/'summary.json').write_text(json.dumps(dict(complete=True,
                statistics=dict(n_prompts=100, pending_prompts=0))))
            status = dashboard.ProgressReader(output).snapshot()
            self.assertEqual(status['state'], 'complete')
            self.assertEqual(status['original_missing_done'], 10)
            self.assertEqual(status['original_missing_total'], 10)


if __name__ == '__main__':
    unittest.main()
