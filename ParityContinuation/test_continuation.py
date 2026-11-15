import unittest
from run import remaining_cases,aggregate,validate_result


class ContinuationTests(unittest.TestCase):
    def test_resolved_pilot_is_not_resubmitted_and_input_is_preserved(self):
        cases=[dict(class_id='a'),dict(class_id='b'),dict(class_id='c')]
        pilot=[dict(class_id='b',outcome='plant_unrealizable')]
        self.assertEqual([c['class_id'] for c in remaining_cases(cases,pilot)],['a','c'])
        self.assertEqual(len(cases),3)

    def test_unknown_pilot_cannot_be_skipped(self):
        with self.assertRaises(ValueError):
            remaining_cases([dict(class_id='a')],[dict(class_id='a',outcome='computational_failure')])

    def test_aggregation_weights_all_equivalent_records_together(self):
        source=[dict(group=dict(class_id='a',record_ids=[1,2]),synthesis=dict(outcome='computational_failure')),
                dict(group=dict(class_id='b',record_ids=[3]),synthesis=dict(outcome='plant_unrealizable'))]
        stats=aggregate(source,{'a':dict(outcome='controller_issued')})
        self.assertEqual(stats['class_counts'],dict(controller_issued=1,plant_unrealizable=1))
        self.assertEqual(stats['record_counts'],dict(controller_issued=2,plant_unrealizable=1))

    def test_negative_decision_requires_empty_winning_initial_region(self):
        result=dict(issued=False,outcome='plant_unrealizable',native_statistics=dict(winning=True,winning_initial_state_count=1))
        with self.assertRaises(ValueError):validate_result(result)


if __name__=='__main__':unittest.main()
