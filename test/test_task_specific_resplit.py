import copy
import random
import unittest

import ResplitTaskSpecific as split


class ResplitTests(unittest.TestCase):
    def test_uniform_exact_size_and_equivalent_records_stay_together(self):
        original = [dict(class_id=f'c{i}', formulas=[f'formula{i}'],
                         record_ids=[i, i + 10] if i == 0 else [i], split='test')
                    for i in range(5)]
        snapshot = copy.deepcopy(original)
        groups, order = split.allocation(original, 2, 1234)
        expected = [g['class_id'] for g in sorted(original, key=lambda g: tuple(g['formulas']))]
        random.Random(1234).shuffle(expected)
        self.assertEqual(order, expected)
        self.assertEqual({g['class_id'] for g in groups if g['split'] == 'calibration'}, set(expected[:2]))
        self.assertEqual(original, snapshot)
        self.assertEqual(next(g for g in groups if g['class_id'] == 'c0')['record_ids'], [0, 10])

    def test_reassignment_preserves_outcomes_and_artifact_paths(self):
        groups = [dict(class_id='a', record_ids=[1, 2], split='calibration'),
                  dict(class_id='b', record_ids=[3], split='test')]
        records = [dict(id=i) for i in range(1, 4)]
        rows = [dict(prompt_index=i, record_id=i, input_field='nl', class_id='a' if i < 3 else 'b',
                     k=1, loss=1 if i == 1 else 0, g=1 if i == 1 else None,
                     implies_reference=i != 1, controller_file='original/controllers/x.jls')
                for i in range(1, 4)]
        before = copy.deepcopy(rows)
        output = split.reassign(rows, groups, records, ['nl'])
        self.assertEqual(output['calibration'], rows[:2])
        self.assertEqual(output['test'], rows[2:])
        self.assertEqual(rows, before)
        self.assertEqual(split.reassign(rows, groups, records, ['nl']), output)

    def test_missing_duplicate_and_unknown_nonzero_losses_are_rejected(self):
        group = [dict(class_id='a', record_ids=[1], split='test')]
        record = [dict(id=1)]
        row = dict(prompt_index=1, record_id=1, input_field='nl', class_id='a',
                   k=1, loss=1, g=1, implies_reference=False)
        for rows in ([], [row, row], [dict(row, g=None)]):
            with self.assertRaises(ValueError):
                split.reassign(rows, group, record, ['nl'])


if __name__ == '__main__':
    unittest.main()
