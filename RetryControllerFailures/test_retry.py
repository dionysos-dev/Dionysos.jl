import copy
import unittest
import run as retry


class RetryTests(unittest.TestCase):
    def fixture(self):
        old = dict(prompt_index=1, k=1, implies_reference=False, g=0, tilde_g=0,
                   loss=0, synthesis_outcome='computational_failure', controller_job_id='old')
        result = dict(job_id='new', issued=True, outcome='controller_issued',
                      region_known=True, x_win_size=5, certificate=dict(satisfied=True),
                      controller_file='controllers/new.jls',region_file='regions/new.json',
                      artifact_base_directory='/saved/new')
        return old,result

    def test_recovered_controller_changes_loss_and_keeps_provenance(self):
        old,result=self.fixture()
        before=copy.deepcopy(old)
        row=retry.overlay_rows([old],{'old':result},'solver2')[0]
        self.assertEqual((row['g'],row['tilde_g'],row['loss']),(1,1,1))
        self.assertEqual(row['previous_controller_job_id'],'old')
        self.assertEqual(row['controller_file'],'controllers/new.jls')
        self.assertEqual(old,before)

    def test_other_cases_unchanged_and_failed_retry_stays_zero(self):
        old,result=self.fixture()
        result.update(issued=False,outcome='computational_failure')
        untouched=dict(prompt_index=2, synthesis_outcome='plant_unrealizable',g=0,loss=0)
        rows=retry.overlay_rows([old,untouched],{'old':result},'solver2')
        self.assertEqual(rows[0]['loss'],0)
        self.assertEqual(rows[1],untouched)

    def test_success_without_verified_initial_region_rejected(self):
        old,result=self.fixture()
        result['x_win_size']=0
        with self.assertRaises(ValueError):retry.overlay_rows([old],{'old':result},'solver2')


if __name__=='__main__':unittest.main()
