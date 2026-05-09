# Final scaling validation diagnostics

Monte Carlo stress test is not a proof. LMI feasibility is not enough if the ellipsoid leaves the local linearization box. A configuration is marked formally defensible only when certification succeeds and post-hoc local checks pass.

## Code state
- No adaptive boxes are used in the certifier core.
- The only refinement option is `state_scaling`.
- `state_scaling = nothing` calls the historical `UT.transition_backward` path.

## Identity scaling check
passed

## Algebraic round-trip checks
passed

## Ellipsoid convention
passed

## Stress test diagnosis
- Samples are generated in E0 and initial membership is checked in `final_validation_stress_debug.csv`.
- Low success rates are mainly chain exits, not target failures or simulation errors, in the current runs.
- Number of debug failures by step <= 3: 498.
- This strongly suggests the certified source ellipsoids are not always contained in the local linearization boxes or that the local model bounds are not valid over the full certified ellipsoid. See box violation columns.

## Scaling sweep
- Largest median volume: `xy_10_angles_1`.
- Largest own-E0 stress success: `trajectory_std`.
- These are not the same configuration.

## Box sweep
- Formal candidates in box sweep: 5.

## Formal validity
No scaling sweep configuration satisfies all formal post-hoc checks.

## Recommendation
Use the scaling results as experimental evidence only if identity scaling passes and the selected configuration has zero box/input violations. Otherwise, large physical volumes are suspicious and should not be presented as formal certificate improvements.
