# Symbolic Control of Velocity-Controlled Mechanical Systems

Reproducibility capsule of the ECC 2027 paper by R. Bertollo, J. Calbert and R. M. Jungers.

It regenerates the numbers, tables and plots of the paper (Table I, Figs. 4–5), which compare two
approaches (Fig. 1 of the paper): the baseline (B) abstracts the torque or force model, and the
proposed approach (P) abstracts the reduced model under the impulse budget. It also gives the
parameters left out of the paper: masses, grids, command sets, torque and force steps, and the
gains of the PI loop of Section VI-B.

## Run

On CodeOcean, press *Reproducible Run*. It runs `./run` and writes the results to `/results`.
Locally, from this folder inside the Dionysos.jl repository:

```bash
./run                                          # instantiate, then run everything
julia --project=. --threads=auto run_all.jl    # the same, once instantiated
julia --project=. --threads=auto experiments/E1_pendulum.jl   # one experiment
```

| Variable | Default | Effect |
|---|---|---|
| `SAMPLES=n` | 3 | median wall time of `n` rounds (the baseline (B) of E2 is timed once) |
| `FIGURES=0` | on | skip the plots and animations |
| `RESULTS_DIR=path` | `results/` | output folder (`/results` on CodeOcean) |

A full run takes about 15 minutes, mostly for the abstraction of the cart-pole force model by
(B), which is built in parallel. The timings in the paper come from a 14-thread laptop
(`results/environment.txt`).

## From the paper to the code

| Paper | Script | Outputs in `results/` |
|---|---|---|
| Section VI-A, pendulum; Table I (top), Fig. 4 | `experiments/E1_pendulum.jl`, models in Dionysos's `problems/Pendulum/` | `params_e1.md`, `table_e1.tex`, `fig_e1.pdf`, `e1_{B,P}.gif`, numbers `e1.*` |
| Section VI-B, cart-pole; Table I (bottom), Fig. 5 | `experiments/E2_cartpole.jl`, models in Dionysos's `problems/CartPendulum/point_mass_cart_pendulum.jl` | `params_e2.md`, `table_e2.tex`, `fig_e2.pdf`, `e2_P.gif` (ideal impulses; (B) is unsolved), `e2_P_PI.gif` (inner PI loop), numbers `e2.*` |

The abstractions and controllers are computed by Dionysos's uniform grid abstraction solver. The
table of the paper stacks `table_e1.tex` and `table_e2.tex`, its plots are copies of
`fig_e*.pdf`, and its numbers come from `numbers.tex`, where `B` and `P` name the two
approaches. The `.png` previews and `.gif` animations are not used in the paper.

## Controllers and replays

At each sampling instant `kτ`, a controller reads the cell `p` of the state and returns the input
held until the next instant:

- (B): a torque or force, `u_k = C(p)`.
- (P): a velocity command, `u_k = C(p, v)`, where the memory `v = u_{k-1}` is the previous
  command and `v = 0` at the start. The impulse budget restricts which commands may follow `v`.

The controller is `MOI.get(optimizer, MOI.RawOptimizerAttribute("concrete_controller"))`, used
with `initial_state`, `output_control` and `update_state` from `Dionysos.System`. A run ends in the
target; in E1, the last command must also allow the pendulum to stop.

Each controller is replayed on the continuous-time model:

- (B) (`replay_torque` in E1, `replay_force` in E2): the input is held over each period, and the
  dynamics are integrated with RK4.
- (P) with ideal impulses (`replay_reduced` in E1 and E2): the velocity jumps to each new command,
  as in Assumption 1 (impulsive actuation), and the rest of the state is integrated with RK4 (in
  E1, θ simply moves by `τ u_k`). The replay reports the impulse `γ_k = α(q)(u_k − u_{k-1})` of
  each step, bounded by the budget, and the hold torque or force `Γ_hold` between steps, which is
  not bounded (columns `max|γ_k|` and `max|Γ|` of Table I).
- (P) with an inner PI loop (`replay_inner_loop` in E2, Section VI-B, Inner loop): a PI loop on
  the cart velocity (`velocity_loop` in the problem file) sets the force every 1 ms, clamped to
  the force limit `Γ_max` = 20 N, while the controller of (P) sets its setpoint, the command,
  every 0.05 s. This replay is not certified. It shows how close a force-limited loop comes to
  ideal impulses, from the start of the paper (the green run of Fig. 5) and from every cell where
  the controller starts at rest (numbers `e2.P.loop.*` and `e2.P.starts`).

## Parameters

The scripts write `params_e1.md` and `params_e2.md` from the constants they use:

- `params_e1.md`: pendulum physics, torque limit, impulse budget, sampling period, speed limit,
  target; grid and commands of (P); grid, torque step and arrival speed of (B).
- `params_e2.md`: masses and length, force limit, impulse budget, speed limit, sampling period,
  target; grid, commands and budget steps of (P); grid and forces of (B); gains of the PI loop.

## Checks

The run stops with an error, or warns, if a result contradicts the paper:

- the replays of (P) reach the target within the impulse budget (E1, E2), the replay of (B)
  reaches the target (E1), and a solved (B) stays within its force limit (E2);
- the steps of (P) from the bottom are the fewest possible, compared with a breadth-first search
  on the augmented abstraction (E1);
- (P) reaches the target from every cell where its controller starts at rest, with ideal
  impulses and with the inner PI loop (E2, warning).

## Environment

- Julia 1.12 (results from 1.12.5); packages pinned by `Project.toml` and `Manifest.toml`.
- Dionysos is loaded from source two folders up (`path = "../.."` in the manifest), with the
  problems of its `problems/` folder. Locally this is the enclosing repository. On CodeOcean,
  `environment/postInstall` clones Dionysos.jl at a fixed commit and installs the packages when the
  environment is built, so a run needs no network.
- Each run records the date, Dionysos commit, Julia version and machine in `environment.txt`.

## Layout

```
run, run_all.jl           entry points
common.jl                 shared helpers
Project.toml, Manifest.toml
experiments/              one script per experiment
environment/postInstall   CodeOcean build step
results/                  outputs, regenerated by each run
```
