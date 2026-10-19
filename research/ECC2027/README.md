# Symbolic Control of Velocity-Controlled Mechanical Systems

Reproducibility capsule of the ECC 2027 paper by R. Bertollo, J. Calbert and R. M. Jungers.

It regenerates every number, table and plot of the paper (Tables I–II, Figs. 4–5), and gives the
parameters that the paper leaves out for space: masses, grids, command sets, torque and force
steps, and the gains of the PI loop of Section VI-B.

## Run

On CodeOcean, press *Reproducible Run*: it runs `./run`, and the results appear in `/results`.
Locally, from this folder inside the Dionysos.jl repository:

```bash
./run                                          # instantiate, then run everything
julia --project=. --threads=auto run_all.jl    # the same, once the environment is instantiated
```

An experiment can also run alone, e.g. `julia --project=. --threads=auto experiments/E1_pendulum.jl`.
E2 builds its abstractions in parallel, hence `--threads=auto`.

| Variable | Default | Effect |
|---|---|---|
| `SAMPLES=n` | 3 | report the median wall time of `n` timed rounds (the force model (T) of E2 is timed once) |
| `FIGURES=0` | on | skip the plots and animations |
| `RESULTS_DIR=path` | `results/` | where the results go (`/results` on CodeOcean) |

A full run takes about 20 minutes, most of it for the abstraction of the cart-pole force model
(T). The timings of the paper come from a 14-thread laptop (`results/environment.txt`); fewer
cores take longer.

## From the paper to the code

| Paper | Script | Outputs in `results/` |
|---|---|---|
| Section VI-A, pendulum swing-up; Table I, Fig. 4 | `experiments/E1_pendulum.jl`, models in Dionysos's `problems/Pendulum/` | `params_e1.md`, `table_e1.tex`, `fig_e1.pdf`, `e1_{T,F}.gif`, numbers `e1.*` |
| Section VI-B, cart-pole swing-up; Table II, Fig. 5 | `experiments/E2_cartpole.jl`, models in Dionysos's `problems/CartPendulum/point_mass_cart_pendulum.jl` | `params_e2.md`, `table_e2.tex`, `fig_e2.pdf`, `e2_F.gif` ((T) is unsolved), numbers `e2.*` |

Every abstraction and synthesized controller comes from Dionysos's uniform grid abstraction
solver; the scripts set the problems, run the solver and replay its controllers. The numbers,
tables and plots of the paper are copies of `numbers.tex`, `table_e*.tex` and `fig_e*.pdf` from
the final run, which a new run rewrites for comparison. The `fig_e*.png` previews and the `.gif`
animations are not in the paper.

## The control laws and their replay

Each synthesis returns a feedback on the cells of its grid: at each sampling instant `kτ`, the
controller reads the cell `p` of the state and returns the input held over `[kτ, (k+1)τ)`.

- (T): a torque or force, `u_k = C(p)`.
- (F): a velocity command, `u_k = C(p, u_{k-1})` from `u_{-1} = 0`. The last command is the
  controller's memory; the impulse budget decides which commands may follow it.

In the code, the controller is `MOI.get(optimizer, MOI.RawOptimizerAttribute("concrete_controller"))`,
driven by `initial_state`, `output_control` and `update_state` of `Dionysos.System`. A run ends in
the target; in E1 the final command 0 must also be admissible, so that the pendulum can stop.

The scripts replay each controller on the continuous-time model, not on the abstraction:

- (T), `closed_loop_torque` in E1 and `replay_T` in E2: the input is held over the period and the
  dynamics are integrated by Runge-Kutta (RK4).
- (F) with the ideal servo of Assumption "Velocity servos", `replay` in E1 and `replay_F` in E2:
  the actuated velocity jumps from `u_{k-1}` to `u_k` at `kτ` and stays at `u_k`; the rest of the
  state is integrated by RK4 (E2, 20 substeps per period) or exactly (E1, θ moves by `τ u_k`). The
  replay computes what the servo must deliver: at each instant the impulse `α(q)(u_k − u_{k-1})`,
  which the budget bounds (column "impulse" of Tables I–II), and between instants the hold torque
  or force, which nothing bounds (column "torque" or "force").
- (F) with a force-limited servo, `replay_servo` in E2 (Section VI-B, "Two time levels"): the ideal
  servo ignores the duration of a real servo's ramp and its tracking error. Here, on the force
  model, a PI loop on the cart velocity (`velocity_loop` in the problem file) sets the force every
  1 ms, clamped to the 20 N limit of (T), with anti-windup, while the controller of (F) sets its
  setpoint every τ = 0.05 s. Nothing certifies this loop; the replay measures its gap to the
  idealization, from the start of the paper (green in Fig. 5, whose middle panel shows the two time
  levels) and from every cell where the controller starts at rest (numbers `e2.F.servo.*`).

## Parameters

The scripts write the `params_e*.md` files from the constants they compute with, so the files
cannot drift from the results:

- `params_e1.md`: pendulum physics, torque limit, impulse budget, sampling period, speed limit,
  target; the grid and commands of (F); the grid, torque step and arrival speed of (T).
- `params_e2.md`: cart and pole masses and length, force limit, impulse budget, speed limit,
  sampling period, target; the grid, commands and budget steps of (F); the grid and forces of (T);
  the PI loop of the servo replay.

## Checks

A run stops with an error, or warns, when a claim of the paper fails on its results:

- every replay of (F) must reach the target within the impulse budget (E1, E2), that of (T) must
  reach the target (E1), and a solved (T) must stay within its force limit (E2);
- the steps of (F) from the bottom must be the fewest possible, against an independent
  breadth-first search on the augmented abstraction (E1, optimality in Section VI-A);
- (F) must reach the target from every cell where its controller starts at rest, with the ideal
  servo and with the PI loop (E2, a warning).

## Environment

- Julia 1.12 (the results come from 1.12.5); `Project.toml` and `Manifest.toml` pin the packages.
- Dionysos is loaded from source two folders up (`path = "../.."` in the manifest), with the
  problem definitions of its `problems/` folder: locally, the enclosing Dionysos.jl repository; on
  CodeOcean, the clone that `environment/postInstall` makes at the commit of the paper's results,
  installing and precompiling the packages when the environment is built, so that a Reproducible
  Run needs no network.
- Each run records the date, the Dionysos commit, the Julia version and the machine in
  `environment.txt`, and the Julia version and CPU in `numbers.tex`.

## Layout

```
run, run_all.jl           entry points
common.jl                 shared helpers: results folder, numbers, tables, environment
Project.toml, Manifest.toml
experiments/              one script per experiment (see the table above)
environment/postInstall   CodeOcean build step: Julia, Dionysos.jl at a fixed commit, packages
results/                  everything the paper quotes, plus previews, regenerated by each run
```
