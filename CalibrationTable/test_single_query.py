import tempfile
from pathlib import Path
import unittest

import single_query_baseline as baseline


class SingleQueryTests(unittest.TestCase):
    def test_first_query_is_not_frequency_winner(self):
        first = dict(representative='F(a)', sample_indices=[1],
                     equivalent_to_reference=False, implies_reference=True)
        winner = dict(representative='G(a)', sample_indices=[2, 3],
                      equivalent_to_reference=True, implies_reference=True)
        analysis = dict(prompt_results=[dict(prompt_index=7, record_id=3, input_field='nl',
                        reference_ltl='G(a)', ranked_clusters=[winner, first], excluded_candidates=[])])
        raw = dict(prompts=[dict(prompt_index=7, samples=[dict(query_index=2, formula='G(a)'),
                                                       dict(query_index=1, formula='F(a)')])])
        rows, prompts = baseline.select_first(analysis, raw, [dict(prompt_index=7, class_id='class')])
        self.assertEqual(rows[0]['formula'], 'F(a)')
        self.assertEqual(prompts[7]['ranked_clusters'], [first])
        baseline.assign(rows[0], dict(issued=True, outcome='controller_issued', formula='F(a)'), 'proof', 'cache')
        self.assertEqual((rows[0]['equivalence_loss'], rows[0]['inclusion_loss'], rows[0]['issued']), (1, 0, 1))

    def test_invalid_first_query_abstains_despite_valid_second(self):
        analysis = dict(prompt_results=[dict(prompt_index=1, record_id=1, input_field='nl', reference_ltl='a',
            ranked_clusters=[dict(representative='a', sample_indices=[2])],
            excluded_candidates=[dict(repetition=1, reason='invalid_ltl')])])
        raw = dict(prompts=[dict(prompt_index=1, samples=[dict(query_index=1, formula='oops'),
                                                       dict(query_index=2, formula='a')])])
        rows, prompts = baseline.select_first(analysis, raw, [dict(prompt_index=1, class_id='class')])
        self.assertFalse(prompts)
        self.assertEqual((rows[0]['issued'], rows[0]['equivalence_loss'], rows[0]['inclusion_loss']), (0, 0, 0))

    def test_class_means_and_unknown_bounds(self):
        rows = [dict(class_id='a', issued=1, equivalence_loss=1, inclusion_loss=0,
                     equivalent_to_reference=False, implies_reference=True),
                dict(class_id='a', issued=None, equivalence_loss=None, inclusion_loss=None,
                     equivalent_to_reference=False, implies_reference=False),
                dict(class_id='b', issued=0, equivalence_loss=0, inclusion_loss=0,
                     equivalent_to_reference=None, implies_reference=None)]
        stats = baseline.statistics(rows)
        self.assertEqual(stats['equivalence_risk']['lower_exact'], '1/4')
        self.assertEqual(stats['equivalence_risk']['upper_exact'], '1/2')
        self.assertEqual(stats['inclusion_risk']['upper_exact'], '1/4')
        self.assertEqual(stats['issuance']['upper_exact'], '1/2')
        self.assertEqual(stats['pending_prompts'], 1)

    def test_finished_table_uses_same_issuance_for_both_baselines(self):
        rows = [dict(class_id='a', issued=1, equivalence_loss=1, inclusion_loss=0,
                     equivalent_to_reference=False, implies_reference=True)]
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory)
            baseline.write_report(output, rows, {'table': []}, {})
            table = (output/'table.tex').read_text()
            self.assertIn('equivalence) & No & -- & -- & 100.00 & 100.00', table)
            self.assertIn('language inclusion) & No & -- & -- & 0.00 & 100.00', table)


if __name__ == '__main__':
    unittest.main()
