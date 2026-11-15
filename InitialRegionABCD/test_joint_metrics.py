import hashlib
import json
from pathlib import Path
import tempfile
import unittest
from joint_metrics import summarize

CONFIG=dict(initial_domain="all_admissible_plant_states",issuance_rule="nonempty projected region")


def row(folder,k,states,issued,implies,job=None,formula=None):
    job=job or f"k{k}"
    r=dict(prompt_index=1,record_id=1,k=k,formula=formula or f"formula{k}",reference_ltl="reference",
        job_id=job,issued=issued,g=int(issued),loss=int(issued and not implies),implies_reference=implies,
        outcome="controller_issued" if issued else "computational_failure" if states is None else "plant_unrealizable",
        subreason="verified" if issued else "resource_limit" if states is None else "empty_region",backend="parity",
        abstention_reason=None if issued else "fixture",region_known=states is not None,x_win_size=None if states is None else len(states))
    if states is not None:
        relative=f"{job}.json";path=folder/relative
        raw=json.dumps(dict(job_id=job,formula=r["formula"],x_win=states,size=len(states))).encode()
        path.write_bytes(raw);r.update(region_file=relative,region_sha256=hashlib.sha256(raw).hexdigest())
    return r


class MetricsTests(unittest.TestCase):
    def test_risk_denominators_and_conditional(self):
        with tempfile.TemporaryDirectory() as f:
            p=Path(f)
            result=summarize(CONFIG,[row(p,1,[1],True,False)],p)
            s=result["summary"][0]
            self.assertEqual(s["joint_risk_rate"],1/40)
            self.assertEqual(s["crc_corrected_risk"],2/41)
            self.assertEqual(s["conditional_diagnostic"],1)

    def test_inclusion_violation_is_detected_not_smoothed(self):
        with tempfile.TemporaryDirectory() as f:
            p=Path(f)
            s=summarize(CONFIG,[row(p,1,[1],True,False),row(p,2,[1,2],True,False)],p)
            self.assertEqual(len(s["monotonicity"]["winning_region_inclusion_violations"]),1)

    def test_resource_abstention_can_break_loss_monotonicity(self):
        with tempfile.TemporaryDirectory() as f:
            p=Path(f)
            s=summarize(CONFIG,[row(p,1,None,False,False),row(p,2,[1],True,False)],p)
            self.assertEqual(len(s["monotonicity"]["loss_monotonicity_violations"]),1)
            self.assertEqual(s["monotonicity"]["region_pairs_unknown"],1)

    def test_saturation_reuses_identical_known_regions(self):
        with tempfile.TemporaryDirectory() as f:
            p=Path(f)
            s=summarize(CONFIG,[row(p,1,[1],True,True,job="same",formula="same"),
                row(p,2,[1],True,True,job="same",formula="same")],p)
            self.assertEqual(s["monotonicity"]["saturated_pairs"],1)
            self.assertFalse(s["monotonicity"]["saturation_violations"])
            self.assertEqual(s["monotonicity"]["region_pairs_comparable"],1)


if __name__=="__main__":unittest.main()
