import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch
from run import run_bounded


class PilotTests(unittest.TestCase):
    def test_memory_guard_stops_the_child_instead_of_waiting(self):
        with tempfile.TemporaryDirectory() as tmp, patch('run.group_rss_kib',return_value=2048):
            result=run_bounded([sys.executable,'-c','import time; time.sleep(30)'],Path(tmp),
                remaining=10,max_rss_kib=1024,min_free=15)
            self.assertEqual(result['guard'],'process_memory_guard')
            self.assertIsNotNone(result['exit_code'])
            self.assertLess(result['seconds'],10)
            self.assertTrue(json.loads((Path(tmp)/'active-process.json').read_text())['stopped'])

    def test_timeout_stops_child_and_remains_a_computational_failure_reason(self):
        with tempfile.TemporaryDirectory() as tmp, patch('run.group_rss_kib',return_value=1), patch('run.memory_free_percent',return_value=100):
            result=run_bounded([sys.executable,'-c','import time; time.sleep(30)'],Path(tmp),
                remaining=0.1,max_rss_kib=1024,min_free=15)
            self.assertEqual(result['guard'],'timeout')
            self.assertIsNotNone(result['exit_code'])
            self.assertLess(result['seconds'],10)


if __name__=='__main__':unittest.main()
