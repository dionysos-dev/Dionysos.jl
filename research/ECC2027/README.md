# Symbolic Control of Velocity-Controlled Mechanical Systems

Reproducibility capsule of the ECC 2027 paper by R. Bertollo, J. Calbert and R. M. Jungers.

It regenerates the numbers, tables and plots of the paper (Table I, Figs. 4–5). It also gives
the parameters left out of the paper: masses, grids, command sets, torque and force steps, and the
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
| `SAMPLES=n` | 3 | median wall time of `n` rounds (the force model (T) of E2 is timed once) |
| `FIGURES=0` | on | skip the plots and animations |
| `RESULTS_DIR=path` | `results/` | output folder (`/results` on CodeOcean) |

A full run takes about 15 minutes, mostly for the abstraction of the cart-pole force model (T),
which is built in parallel. The timings in the paper come from a 14-thread laptop
(`results/environment.txt`).

## From the paper to the code

| Paper | Script | Outputs in `results/` |
|---|---|---|
| Section VI-A, pendulum; Table I (top), Fig. 4 | `experiments/E1_pendulum.jl`, models in Dionysos's `problems/Pendulum/` | `params_e1.md`, `table_e1.tex`, `fig_e1.pdf`, `e1_{T,R}.gif`, numbers `e1.*` |
| Section VI-B, cart-pole; Table I (bottom), Fig. 5 | `experiments/E2_cartpole.jl`, models in Dionysos's `problems/CartPendulum/point_mass_cart_pendulum.jl` | `params_e2.md`, `table_e2.tex`, `fig_e2.pdf`, `e2_R.gif` (ideal impulses; (T) is unsolved), `e2_R_PI.gif` (inner PI loop), numbers `e2.*` |

The abstractions and controllers are computed by Dionysos's uniform grid abstraction solver. The
table of the paper stacks `table_e1.tex` and `table_e2.tex`, its plots are copies of
`fig_e*.pdf`, and its numbers come from `numbers.tex`, where `T` and `R` name the two
approaches. The `.png` previews and `.gif` animations are not used in the paper.

## Controllers and replays

At each sampling instant `kτ`, a controller reads the cell `p` of the state and returns the input
held until the next instant:

- (T): a torque or force, `u_k = C(p)`.
- (R): a velocity command, `u_k = C(p, u_{k-1})`, with `u_{-1} = 0`. The impulse budget restricts
  which commands may follow `u_{k-1}`.

The controller is `MOI.get(optimizer, MOI.RawOptimizerAttribute("concrete_controller"))`, used
with `initial_state`, `output_control` and `update_state` from `Dionysos.System`. A run ends in the
target; in E1, the last command must also allow the pendulum to stop.

Each controller is replayed on the continuous-time model:

- (T) (`closed_loop_torque` in E1, `replay_T` in E2): the input is held over each period, and the
  dynamics are integrated with RK4.
- (R) with ideal impulses (`replay` in E1, `replay_R` in E2): the velocity jumps to each new
  command, as in Assumption 1 (impulsive actuation), and the rest of the state is integrated with RK4
  (in E1, θ simply moves by `τ u_k`). The replay reports the impulse `α(q)(u_k − u_{k-1})` of
  each step, bounded by the budget, and the hold torque or force between steps, which is not
  bounded (columns "impulse" and "force" of Table I).
- (R) with an inner PI loop (`replay_inner_loop` in E2, Section VI-B, Inner loop): a PI loop on
  the cart velocity (`velocity_loop` in the problem file) sets the force every 1 ms, clamped to
  the 20 N of (T), while the controller of (R) sets its setpoint, the command, every 0.05 s. This
  replay is not certified; it shows how close a force-limited loop comes to ideal impulses. It
  runs from the start used in the paper (green in Fig. 5) and from every cell where the
  controller starts at rest (numbers `e2.R.loop.*`).

## Parameters

The scripts write `params_e1.md` and `params_e2.md` from the constants they use:

- `params_e1.md`: pendulum physics, torque limit, impulse budget, sampling period, speed limit,
  target; grid and commands of (R); grid, torque step and arrival speed of (T).
- `params_e2.md`: masses and length, force limit, impulse budget, speed limit, sampling period,
  target; grid, commands and budget steps of (R); grid and forces of (T); gains of the PI loop.

## Checks

The run stops with an error, or warns, if a result contradicts the paper:

- the replays of (R) reach the target within the impulse budget (E1, E2), the replay of (T)
  reaches the target (E1), and a solved (T) stays within its force limit (E2);
- the steps of (R) from the bottom are the fewest possible, compared with a breadth-first search
  on the augmented abstraction (E1);
- (R) reaches the target from every cell where its controller starts at rest, with ideal
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
