import importlib.util
from pathlib import Path
from fractions import Fraction
import unittest

spec = importlib.util.spec_from_file_location('table_runner', Path(__file__).with_name('run.py'))
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)


def prompt(index=1, implying=(False, False), sat=(True, True)):
    return dict(prompt_index=index, record_id=index, input_field='NL',
        conjunctions=[dict(phi_k=f'f{j}', implies_reference=imp, satisfiable=s)
                      for j,(imp,s) in enumerate(zip(implying,sat),1)])


class TableTests(unittest.TestCase):
    def result(self, issued):
        return dict(issued=issued, outcome='controller_issued' if issued else 'computational_failure',
                    job_id='job', checkpoint_file='result.json')

    def test_failed_prefix_is_absorbing_even_when_next_synthesis_would_succeed(self):
        calls=[]
        def solve(*args):
            calls.append(args)
            return self.result(False)
        p=prompt()
        first=m.evaluate([p],1,{},solve,{1:'a'},True)
        second=m.evaluate([p],2,{1:first[0]},lambda *a:self.fail('Blocked prefix must skip synthesis'),{1:'a'},True)
        self.assertEqual((first[0]['g'],second[0]['g'],second[0]['loss']),(0,0,0))
        self.assertIsNone(second[0]['tilde_g'])
        self.assertEqual(len(calls),1)

    def test_guaranteed_formula_requires_synthesis_for_issuance(self):
        p=prompt(implying=(True,))
        shortcut=m.evaluate([p],1,{},lambda *a:self.fail('Logical zero needs no calibration solve'),{1:'a'})
        self.assertIsNone(shortcut[0]['g'])
        measured=m.evaluate([p],1,{},lambda *a:self.result(True),{1:'a'},True)
        self.assertEqual((measured[0]['g'],measured[0]['guaranteed_issuance'],measured[0]['loss']),(1,1,0))

    def test_saturation_and_empty_pool(self):
        p=prompt(implying=(False,))
        self.assertIs(m.conjunction(p,10),p['conjunctions'][0])
        p['conjunctions']=[]
        rows=m.evaluate([p],1,{},lambda *a:self.fail('No candidate must abstain'),{1:'a'},True)
        self.assertEqual((rows[0]['g'],rows[0]['loss'],rows[0]['guaranteed_issuance']),(0,0,0))

    def test_class_weighting_and_rate_identity(self):
        rows=[dict(class_id='a',loss=1,g=1,guaranteed_issuance=0),
              dict(class_id='a',loss=0,g=1,guaranteed_issuance=1),
              dict(class_id='b',loss=0,g=0,guaranteed_issuance=0)]
        s=m.statistics(rows)
        self.assertEqual(s['empirical_risk'],.25)
        self.assertEqual(s['issuance'],.5)
        self.assertEqual(s['guaranteed_issuance'],.25)
        self.assertEqual(Fraction(s['corrected_exact']),Fraction(1,2))
        self.assertEqual(s['pair_empirical_risk'],1/3)

    def test_first_threshold_exact_boundary_and_infinity(self):
        levels=[dict(k=1,corrected_exact='1/5'),dict(k=2,corrected_exact='1/10')]
        self.assertEqual(m.select(levels,Fraction(1,5)),1)
        self.assertEqual(m.select(levels,Fraction(1,10)),2)
        self.assertIsNone(m.select(levels,Fraction(1,20)))
        self.assertEqual(m.select(levels,Fraction(1,20),max_k=2),'infinity')

    def test_unsat_is_not_issued_despite_vacuous_implication(self):
        p=prompt(implying=(True,),sat=(False,))
        rows=m.evaluate([p],1,{},lambda *a:self.fail('UNSAT must not synthesize'),{1:'a'},True)
        self.assertEqual(rows[0]['evaluation'],'logically_inconsistent')
        self.assertEqual(rows[0]['guaranteed_issuance'],0)


if __name__=='__main__':
    unittest.main()
