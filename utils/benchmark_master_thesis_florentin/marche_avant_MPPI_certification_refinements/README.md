# Marche avant MPPI - state-scaling certification

This benchmark compares the historical ellipsoidal backward certifier with the same certifier solved in normalized state coordinates.

Only two methods are compared:

- `baseline`: no state scaling;
- `scaled`: LMI solved with `state_scaling`, then ellipsoids and feedbacks converted back to physical coordinates.

The Monte Carlo stress test is a sanity check, not a replacement for the LMI certificate. The formal certificate is still the feasibility of the backward LMI chain.

## Step 1: generate MPPI trajectory

```bash
julia --project=. utils/benchmark_master_thesis_florentin/marche_avant_MPPI_certification_refinements/01_generate_and_save_mppi_trajectory.jl
```

This runs MPPI once and saves:

- `results/saved_mppi_trajectory.jld2`

## Step 2: run baseline vs scaled

```bash
julia --project=. utils/benchmark_master_thesis_florentin/marche_avant_MPPI_certification_refinements/02_run_baseline_vs_scaled.jl
```

The script reloads the saved trajectory, certifies it twice, and runs a paired Monte Carlo stress test on samples drawn from the baseline initial ellipsoid.

Outputs:

- `results/baseline_vs_scaled_summary.csv`
- `results/baseline_vs_scaled_detailed.csv`
- `results/baseline_result.jld2`
- `results/scaled_result.jld2`
- `results/statistical_stress_test_summary.csv`
- `results/statistical_stress_test_samples.csv`
- `results/statistical_comparison.csv`

The volumes are physical volumes after reconverting scaled certificates to physical coordinates.

## Step 3: plot

```bash
julia --project=. utils/benchmark_master_thesis_florentin/marche_avant_MPPI_certification_refinements/03_plot_baseline_vs_scaled.jl
```

Plots:

- `results/plots/volume_vs_transition_baseline_scaled.png`
- `results/plots/initial_ellipsoid_comparison.png`
- `results/plots/stress_test_success_rates.png`
- `results/plots/left_chain_rates.png`
- `results/plots/sampled_trajectories_xy.png`

The stress test samples points, replays the certified feedback chain on the concrete dynamics, and reports whether rollouts stay inside the ellipsoid chain and reach the target. It is intended to catch implementation or numerical inconsistencies, not to prove robustness.
