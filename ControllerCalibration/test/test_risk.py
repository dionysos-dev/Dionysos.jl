import importlib.util
from pathlib import Path
from fractions import Fraction
import unittest

spec=importlib.util.spec_from_file_location('controller_crc',Path(__file__).resolve().parents[1]/'run.py')
crc=importlib.util.module_from_spec(spec)
spec.loader.exec_module(crc)


def prompt(index, record, specifications):
    return dict(prompt_index=index,record_id=record,input_field='fixture',conjunctions=[
        dict(phi_k=f,satisfiable=True,implies_reference=implies) for f,implies in specifications])


class OriginalLossTests(unittest.TestCase):
    def test_prefix_product_and_logical_zero_do_not_invent_success(self):
        data=[prompt(1,1,[('bad_A',False),('bad_A_and_more',False)]),
              prompt(2,1,[('good_A',True),('good_A_and_more',True)]),
              prompt(3,2,[('bad_B',False),('bad_B_and_more',False)]),
              prompt(4,3,[('bad_C',False),('good_C',True)])]
        classes={1:'A',2:'B',3:'C'}
        calls=[]
        def solve(formulas):
            calls.extend(formulas)
            return {f:dict(issued=f!='bad_B',outcome='controller_issued' if f!='bad_B' else 'computational_failure',job_id=f) for f in formulas}
        first=crc.evaluate_level(data,1,{},solve,classes)
        self.assertEqual([r['loss'] for r in first],[1,0,0,1])
        self.assertIsNone(first[1]['g'])
        self.assertEqual(first[2]['g'],0)
        self.assertEqual(first[2]['tilde_g'],0)
        stats=crc.statistics(first)
        self.assertEqual(stats['n_classes'],3)
        self.assertEqual(stats['class_loss_sum'],'3/2')
        self.assertEqual(stats['class_mean_corrected_exact'],'5/8')
        self.assertEqual(stats['pair_losses'],2)
        second=crc.evaluate_level(data,2,{r['prompt_index']:r for r in first},solve,classes)
        self.assertEqual([r['loss'] for r in second],[1,0,0,0])
        self.assertEqual(second[2]['evaluation'],'prefix_blocked')
        self.assertIsNone(second[2]['tilde_g'])
        self.assertEqual(second[2]['g'],0)
        self.assertEqual(calls,['bad_A','bad_B','bad_C','bad_A_and_more'])
        stats=crc.statistics(second)
        self.assertEqual(stats['class_loss_sum'],'1/2')
        self.assertEqual(stats['class_mean_corrected_exact'],'3/8')
        alpha=Fraction(1,2)
        self.assertGreater(Fraction(crc.statistics(first)['class_mean_corrected_exact']),alpha)
        self.assertLessEqual(Fraction(stats['class_mean_corrected_exact']),alpha)

    def test_implication_cannot_disappear_after_a_logical_zero_skip(self):
        data=[prompt(1,1,[('good',True),('bad',False)])]
        first=crc.evaluate_level(data,1,{},lambda _: {},{1:'A'})
        with self.assertRaises(ValueError):
            crc.evaluate_level(data,2,{1:first[0]},lambda _: {},{1:'A'})

    def test_no_candidate_and_inconsistency_have_original_loss_zero(self):
        data=[prompt(1,1,[]),dict(prompt_index=2,record_id=2,input_field='fixture',
              conjunctions=[dict(phi_k='false',satisfiable=False,implies_reference=True)])]
        result=crc.evaluate_level(data,1,{},lambda _: {},{1:'A',2:'B'})
        self.assertEqual([r['g'] for r in result],[0,0])
        self.assertEqual([r['loss'] for r in result],[0,0])
        self.assertEqual([r['evaluation'] for r in result],['no_valid_candidate','logically_inconsistent'])


if __name__=='__main__':unittest.main()
