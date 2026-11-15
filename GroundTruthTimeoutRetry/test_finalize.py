import unittest
from finalize import merge_verified_references


class FinalizationTests(unittest.TestCase):
    def fixture(self):
        row=dict(group=dict(class_id='a'),synthesis=dict(formula='Ggreen | Fbrown',subreason='timeout'))
        result=dict(formula='Fbrown | Ggreen',outcome='controller_issued',issued=True,
            region_known=True,x_win_size=5,certificate=dict(satisfied=True))
        return row,result

    def test_equivalent_reordering_is_accepted_with_provenance(self):
        row,result=self.fixture()
        merged,changes=merge_verified_references([row],{'a':result},lambda a,b:True)
        self.assertEqual(merged[0]['synthesis']['formula'],'Fbrown | Ggreen')
        self.assertEqual(merged[0]['previous_synthesis']['formula'],'Ggreen | Fbrown')
        self.assertEqual(len(changes),1)
        self.assertNotIn('previous_synthesis',row)

    def test_non_equivalent_retry_is_rejected(self):
        row,result=self.fixture()
        with self.assertRaises(ValueError):merge_verified_references([row],{'a':result},lambda a,b:False)

    def test_other_source_outcomes_unchanged(self):
        row=dict(group=dict(class_id='b'),synthesis=dict(formula='FGblue',subreason='resource_limit'))
        merged,changes=merge_verified_references([row],{},lambda a,b:False)
        self.assertEqual(merged,[row]);self.assertEqual(changes,[])


if __name__=='__main__':unittest.main()
