"""Cumulative issuance from verified per-level synthesis success.

The solver's success flag is tilde_g. Only g=prod_{j<=k} tilde_g_j is used
for controller issuance and the loss. The mathematical winning region remains
the result of the level-k game, including when issuance is withheld.
"""
from collections import defaultdict
from copy import deepcopy


def _binary(value, name):
    if value not in (False, True, 0, 1):
        raise ValueError(f"{name} must be binary")
    return int(value)


def apply_prefix_issuance(results):
    """Return new rows; require complete prefixes, independent of input order.

    Supports partial experiments only when each observed prompt has levels
    1..K with no gaps. Reapplying the policy reads the saved tilde_g, so it
    never mistakes a policy-withheld controller for a synthesis failure.
    """
    groups = defaultdict(list)
    seen = set()
    for original in results:
        row = deepcopy(original)
        key = row["prompt_index"], row["k"]
        if key in seen:
            raise ValueError(f"Duplicate prompt/level: {key}")
        if not isinstance(row["k"], int) or row["k"] < 1:
            raise ValueError("Levels must be positive integers")
        seen.add(key)
        groups[row["prompt_index"]].append(row)
    processed = {}
    for prompt, rows in groups.items():
        rows.sort(key=lambda r: r["k"])
        if [r["k"] for r in rows] != list(range(1, len(rows)+1)):
            raise ValueError(f"Missing prefix level for prompt {prompt}")
        issued_prefix = 1
        first_failure = None
        for row in rows:
            raw = _binary(row["tilde_g"] if "tilde_g" in row else row["issued"], "tilde_g")
            synthesis_outcome = row.get("synthesis_outcome", row["outcome"])
            synthesis_reason = row.get("synthesis_subreason", row.get("subreason"))
            raw_controller = row.get("synthesized_controller_file", row.get("controller_file"))
            raw_controller_hash = row.get("synthesized_controller_sha256", row.get("controller_sha256"))
            if bool(raw) != (synthesis_outcome == "controller_issued"):
                raise ValueError("Raw success must match the controller-issued solver outcome")
            if raw and not raw_controller:
                raise ValueError("Successful synthesis requires an exported controller")
            if raw and not isinstance(row.get("implies_reference"), bool):
                raise ValueError("Successful synthesis requires a Boolean implication label")
            if not raw and first_failure is None:
                first_failure = dict(k=row["k"], outcome=synthesis_outcome, reason=synthesis_reason)
            issued_prefix *= raw
            blocked = bool(raw and not issued_prefix)
            row.update(
                tilde_g=raw,
                synthesis_issued=bool(raw),
                synthesis_outcome=synthesis_outcome,
                synthesis_subreason=synthesis_reason,
                synthesis_abstention_reason=None if raw else synthesis_reason,
                synthesized_controller_file=raw_controller,
                synthesized_controller_sha256=raw_controller_hash,
                g=issued_prefix,
                issued=bool(issued_prefix),
                loss=int(bool(issued_prefix) and not row["implies_reference"]),
                prefix_blocked=blocked,
                prefix_failure_k=None if first_failure is None else first_failure["k"],
                prefix_failure_outcome=None if first_failure is None else first_failure["outcome"],
                prefix_failure_reason=None if first_failure is None else first_failure["reason"],
                issuance_rule="product_of_raw_synthesis_success_over_levels_1_through_k",
                # This is an issuance region, not a replacement for X_win.
                issued_region_size=row.get("x_win_size") if issued_prefix else 0,
                controller_file=raw_controller if issued_prefix else None,
                controller_sha256=raw_controller_hash if issued_prefix else None,
                outcome="prefix_blocked" if blocked else synthesis_outcome,
                subreason="earlier_prefix_synthesis_failed" if blocked else synthesis_reason,
                abstention_reason=None if issued_prefix else "earlier_prefix_synthesis_failed" if blocked else synthesis_reason,
            )
            processed[prompt, row["k"]] = row
    return [processed[r["prompt_index"], r["k"]] for r in results]


def audit_prefix_issuance(rows):
    groups = defaultdict(list)
    for row in rows:
        groups[row["prompt_index"]].append(row)
    checks = []
    for prompt, values in groups.items():
        values.sort(key=lambda r:r["k"])
        prefix = 1
        for row in values:
            prefix *= row["tilde_g"]
            if row["g"] != prefix or bool(row["g"]) != row["issued"]:
                raise AssertionError("Issuance does not equal the prefix product")
            if row["loss"] != int(bool(prefix) and not row["implies_reference"]):
                raise AssertionError("Loss does not use cumulative issuance")
        for low, high in zip(values,values[1:]):
            checks.append(dict(prompt_index=prompt, lower_k=low["k"], upper_k=high["k"],
                raw_success_low=low["tilde_g"], raw_success_high=high["tilde_g"],
                g_low=low["g"], g_high=high["g"], loss_low=low["loss"], loss_high=high["loss"],
                issuance_nonincreasing=high["g"]<=low["g"],
                implication_persists=not low["implies_reference"] or bool(high["implies_reference"]),
                loss_nonincreasing=high["loss"]<=low["loss"]))
    return dict(adjacent_pairs=len(checks),
        issuance_violations=[r for r in checks if not r["issuance_nonincreasing"]],
        implication_violations=[r for r in checks if not r["implication_persists"]],
        loss_violations=[r for r in checks if not r["loss_nonincreasing"]], pairs=checks)
