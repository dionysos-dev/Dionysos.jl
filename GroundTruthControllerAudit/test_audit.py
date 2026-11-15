import unittest
from run import outcome_counts


class AuditTests(unittest.TestCase):
    def test_unknown_is_not_counted_as_unrealizable(self):
        groups=[dict(class_id='a',record_ids=[1,2]),dict(class_id='b',record_ids=[3]),dict(class_id='c',record_ids=[4])]
        results={'a':dict(outcome='controller_issued'),'b':dict(outcome='computational_failure')}
        stats=outcome_counts(groups,results)
        self.assertEqual(stats['class_counts'],dict(controller_issued=1,computational_failure=1,pending=1))
        self.assertEqual(stats['record_counts']['controller_issued'],2)
        self.assertEqual(stats['controller_existence_lower_bound'],1/3)
        self.assertEqual(stats['controller_existence_upper_bound'],1)

    def test_proved_negatives_reduce_the_upper_bound(self):
        groups=[dict(class_id='a',record_ids=[1]),dict(class_id='b',record_ids=[2])]
        stats=outcome_counts(groups,{'a':dict(outcome='plant_unrealizable'),'b':dict(outcome='logically_inconsistent')})
        self.assertEqual(stats['controller_existence_lower_bound'],0)
        self.assertEqual(stats['controller_existence_upper_bound'],0)


if __name__=='__main__':unittest.main()
