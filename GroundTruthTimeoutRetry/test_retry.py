import unittest
from run import merge_references,summarize


class RetryTests(unittest.TestCase):
    def test_only_timeout_changes_and_preserves_reference_provenance(self):
        timeout=dict(group=dict(class_id='a',record_ids=[1,2]),synthesis=dict(
            formula='XFblue',outcome='computational_failure',subreason='timeout',issued=False))
        cap=dict(group=dict(class_id='b',record_ids=[3]),synthesis=dict(
            formula='FGblue',outcome='computational_failure',subreason='resource_limit',issued=False))
        positive=dict(formula='XFblue',outcome='controller_issued',issued=True,
            region_known=True,x_win_size=5,certificate=dict(satisfied=True))
        merged=merge_references([timeout,cap],{'a':positive})
        self.assertEqual(merged[0]['previous_synthesis'],timeout['synthesis'])
        self.assertEqual(merged[0]['synthesis'],positive)
        self.assertEqual(merged[1],cap)
        self.assertEqual(timeout['synthesis']['subreason'],'timeout')
        summary=summarize(merged)
        self.assertEqual(summary['class_counts'],dict(controller_issued=1,computational_failure=1))
        self.assertEqual(summary['record_counts']['controller_issued'],2)
        self.assertEqual(summary['controller_existence_upper_bound'],1)

    def test_different_formula_is_rejected(self):
        row=dict(group=dict(class_id='a'),synthesis=dict(formula='XFblue',subreason='timeout'))
        with self.assertRaises(ValueError):
            merge_references([row],{'a':dict(formula='Fblue')})

    def test_no_issue_without_verified_initial_region(self):
        row=dict(group=dict(class_id='a'),synthesis=dict(formula='XFblue',subreason='timeout'))
        result=dict(formula='XFblue',outcome='controller_issued',issued=True,region_known=True,x_win_size=0)
        with self.assertRaises(ValueError):merge_references([row],{'a':result})


if __name__=='__main__':unittest.main()
