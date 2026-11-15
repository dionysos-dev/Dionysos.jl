import importlib.util
import tempfile
import unittest
import subprocess
import sys
from pathlib import Path

spec = importlib.util.spec_from_file_location("abcd_runner", Path(__file__).with_name("run.py"))
runner = importlib.util.module_from_spec(spec)
spec.loader.exec_module(runner)


class OutcomeTests(unittest.TestCase):
    def test_timeout_kills_the_worker_and_abstains(self):
        with tempfile.TemporaryDirectory() as folder:
            output=Path(folder)
            worker=runner.Worker.__new__(runner.Worker)
            worker.output=output
            worker.log=(output/"log").open("w")
            worker.process=subprocess.Popen([sys.executable,"-c","import time; time.sleep(60)"],
                stdin=subprocess.PIPE,stdout=worker.log,stderr=worker.log,text=True,start_new_session=True)
            result=worker.solve({"job_id":"fixture","formula":"F(blue)"},
                output/"job.json",output/"result.json",0.1)
            self.assertEqual(result["outcome"],"computational_failure")
            self.assertEqual(result["subreason"],"timeout")
            self.assertFalse(result["issued"])
            self.assertIsNotNone(worker.process.poll())

    def test_five_outcomes_and_pending_use_all_40_prompts(self):
        with tempfile.TemporaryDirectory() as folder:
            results = [dict(k=1, outcome=o, issued=o=="controller_issued", subreason="fixture")
                       for o in runner.OUTCOMES]
            result = runner.summarize({"timeout_seconds":120},results,Path(folder))
            row = result["summary"][0]
            self.assertEqual(row["counts"],dict.fromkeys(runner.OUTCOMES,1))
            self.assertEqual(row["pending"],35)
            self.assertEqual(row["issuance_rate"],1/40)
            self.assertEqual(row["abstention_rate"],4/40)
            self.assertFalse(result["complete"])

    def test_complete_table_partitions_all_cases(self):
        with tempfile.TemporaryDirectory() as folder:
            results = [dict(k=k,outcome=runner.OUTCOMES[i%5],
                issued=i%5==0,subreason="timeout" if i%5==3 else "fixture")
                for k in range(1,6) for i in range(40)]
            result=runner.summarize({"timeout_seconds":120},results,Path(folder))
            self.assertTrue(result["complete"])
            for row in result["summary"]:
                self.assertEqual(sum(row["counts"].values()),40)
                self.assertEqual(row["issuance_rate"]+row["abstention_rate"],1)
                self.assertEqual(row["failure_subreasons"],{"timeout":8})

    def test_actual_input_freeze_is_immutable_and_saturates(self):
        root=Path(__file__).resolve().parents[2]/"Dio"/"Dionysos.jl_1"
        source=Path(__file__).resolve().parents[1]/"Results"/"task-specific-temperature-prompt-v2"
        if not root.is_dir() or not source.is_dir():
            self.skipTest("Local Dionysos and DeepSeek fixtures are unavailable")
        with tempfile.TemporaryDirectory() as folder:
            output=Path(folder)
            config=runner.freeze(source,root,output,120)
            self.assertEqual(len(config["cases"]),200)
            self.assertEqual(config["reserved_cross_model"]["max_k"],4)
            self.assertFalse(config["reserved_cross_model"]["run"])
            for case in config["cases"]:
                self.assertEqual(case["included_clusters"],min(case["k"],case["available_clusters"]))
            with self.assertRaises(RuntimeError):
                runner.freeze(source,root,output,121)


if __name__=="__main__":
    unittest.main()
