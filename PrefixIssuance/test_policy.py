import itertools
import unittest
from policy import apply_prefix_issuance, audit_prefix_issuance


def fixture(raw, implication=None, prompt=1):
    implication=implication or [False]*len(raw)
    return [dict(prompt_index=prompt,k=k,issued=bool(success),
        outcome="controller_issued" if success else "computational_failure",
        subreason="verified_nonempty_initial_region" if success else "timeout",
        implies_reference=implies,x_win_size=12 if success else None,
        controller_file=f"controllers/{prompt}-{k}.jls" if success else None,
        controller_sha256=f"hash-{prompt}-{k}" if success else None)
        for k,(success,implies) in enumerate(zip(raw,implication),1)]


class PrefixTests(unittest.TestCase):
    def test_failure_withholds_all_later_controllers(self):
        original=fixture([1,0,1,1,1])
        result=apply_prefix_issuance(original)
        self.assertEqual([r["tilde_g"] for r in result],[1,0,1,1,1])
        self.assertEqual([r["g"] for r in result],[1,0,0,0,0])
        self.assertEqual([r["loss"] for r in result],[1,0,0,0,0])
        self.assertEqual(result[2]["prefix_failure_k"],2)
        self.assertEqual(result[2]["prefix_failure_reason"],"timeout")
        self.assertIsNone(result[2]["controller_file"])
        self.assertEqual(result[2]["synthesized_controller_file"],original[2]["controller_file"])
        self.assertEqual(result[2]["x_win_size"],12)
        self.assertEqual(result[2]["issued_region_size"],0)
        self.assertTrue(original[2]["issued"])

    def test_exhaustive_prefixes_and_monotone_implication(self):
        comparisons=0
        for raw in itertools.product((0,1),repeat=5):
            for first_correct in range(6):
                implication=[k>=first_correct for k in range(5)]
                result=apply_prefix_issuance(fixture(raw,implication))
                audit=audit_prefix_issuance(result)
                self.assertFalse(audit["issuance_violations"])
                self.assertFalse(audit["loss_violations"])
                self.assertFalse(audit["implication_violations"])
                comparisons+=1
        self.assertEqual(comparisons,192)

    def test_groups_do_not_share_failure_state_and_order_is_irrelevant(self):
        rows=fixture([0,1],prompt=1)+fixture([1,1],prompt=2)
        rows.reverse()
        result=apply_prefix_issuance(rows)
        self.assertEqual([r["g"] for r in result],[1,1,0,0])

    def test_reapplication_preserves_raw_success_and_all_fields(self):
        once=apply_prefix_issuance(fixture([0,1,1]))
        twice=apply_prefix_issuance(once)
        self.assertEqual(once,twice)

    def test_incomplete_prefix_and_duplicates_are_rejected(self):
        with self.assertRaises(ValueError):apply_prefix_issuance(fixture([1,1])[1:])
        row=fixture([1])[0]
        with self.assertRaises(ValueError):apply_prefix_issuance([row,row])

    def test_bad_implication_labels_are_reported(self):
        audit=audit_prefix_issuance(apply_prefix_issuance(fixture([1,1],[True,False])))
        self.assertEqual(len(audit["implication_violations"]),1)
        self.assertEqual(len(audit["loss_violations"]),1)


if __name__=="__main__":unittest.main()
