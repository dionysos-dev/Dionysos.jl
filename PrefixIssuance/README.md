# Cumulative controller issuance

`tilde_g_k` records successful synthesis, verification, and export at level k.
The issuance policy is now

```
g_k = product(tilde_g_j for j in 1..k)
loss_k = g_k * indicator(Phi_k does not imply reference)
```

The first failure makes every subsequent issuance zero for that prompt. Raw
successes, controllers, failure causes, and mathematical winning regions are
preserved. A later successful controller is withheld after an earlier failure;
its nonempty X_win is not changed to an empty region.

Fields in the new analysis:

- `tilde_g`, `synthesis_issued`, `synthesis_outcome`: original per-level success.
- `g`, `issued`, `loss`: cumulative issuance and the updated loss.
- `prefix_failure_k`, `prefix_failure_reason`: the first failing prefix level.
- `prefix_blocked`: a successful current synthesis withheld by an earlier failure.
- `x_win_size`: original mathematical winning-region size (null if unknown).
- `issued_region_size`: zero when issuance is withheld.
- `synthesized_controller_file`: original artifact, including withheld successes.
- `controller_file`: provided only when cumulative issuance is allowed.

All relative artifact paths refer to `artifact_base_directory`, the unchanged
source synthesis experiment. No controllers, regions or source results are
overwritten, and no synthesis or model calls are made.

```sh
python3 PrefixIssuance/test_policy.py
python3 PrefixIssuance/analyze.py \
  --source Results/deepseek-development-abcd-all-initial-v3 \
  --output Results/deepseek-development-abcd-prefix-issuance-v4
```

The policy is a reusable Python function, `apply_prefix_issuance(results)`, for
future evaluation pipelines. It requires every earlier level for each observed
prompt, sorts levels before applying the product, and preserves the raw success
flag on repeated application.

Since g_(k+1) = g_k * tilde_g_(k+1), issuance is nonincreasing even when raw
synthesis success changes from failure to success because of computational
limits. The loss is also nonincreasing when Phi_(k+1) implies Phi_k, as holds
for the cumulative conjunctions. Both properties are checked prompt by prompt.
