# Nonlinear BM vs MPPI benchmark

- Dynamics, domain, obstacle, initial ellipsoid, target ellipsoid and input bounds are reused from `problems/non_linear.jl`.
- The finite horizon is overridden locally inside the benchmark code; the original problem file is left untouched.
- Both planners start from the center of the common initial ellipsoid for nominal trajectory generation.
- MPPI uses the true nonlinear discrete rollout and a neutral deterministic warm start, without abstraction or reference tracking.
- BM cannot directly optimize the repository nonlinear symbolic system, so it is run on a benchmark-local PWA/hybrid surrogate only.
- The BM surrogate uses the same input bounds and horizon, local linear maps on a coarse corridor decomposition, a conservative obstacle outer box, an inner target box contained in the true target ellipsoid, and mode-dependent constant stage costs so it remains solvable with the MILP solvers available in this environment.
- Every BM control sequence is re-simulated on the true nonlinear dynamics before common evaluation and certification.
- Both methods are evaluated with the same true-rollout cost and the same ellipsoidal backward certifier configuration.
- Ellipsoid areas are computed with `UT.get_volume`, which matches the repository ellipsoid convention `(x-c)'P(x-c) <= 1`.
- Because the certifier works backward from the terminal state, common-prefix ellipsoid metrics compare the common terminal-nearest certified ellipsoids available for both methods.

## Default reproducibility settings

- horizons = (6, 8, 10, 12, 14, 16, 18, 20)
- rng_seed = 1234
- mppi_nsamples = 1000
- mppi_niter = 12
- mppi_lambda = 5.0
