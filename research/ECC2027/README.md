# Symbolic Control of Velocity-Controlled Mechanical Systems

Reproducibility capsule of the ECC 2027 paper by R. Bertollo, J. Calbert and R. M. Jungers.

This folder regenerates every number, table and plot of the paper (Tables I–II, Figs. 4–5). It
also holds the parameters that the paper leaves out for space: masses, grids, command sets, and
torque and force steps.

## Run

On CodeOcean, press *Reproducible Run*: it runs `./run`, and the results appear in `/results`.
Locally, from this folder inside the Dionysos.jl repository:

```bash
./run                                          # instantiate, then run everything
julia --project=. --threads=auto run_all.jl    # the same, once the environment is instantiated
```

Each experiment can also run alone, for example
`julia --project=. --threads=auto experiments/E1_pendulum.jl`. The cart-pole abstractions are
built in parallel, so `--threads=auto` matters for E2.

Switches, read from the environment:

| Variable | Default | Effect |
|---|---|---|
| `SAMPLES=n` | 3 | report the median wall time of `n` timed rounds, except for the force model (T) of E2, timed once |
| `FIGURES=0` | on | skip the plots and animations |
| `RESULTS_DIR=path` | `results/` | where the results go (`/results` on CodeOcean) |

A full run takes about 20 minutes once the packages are installed, most of it for the abstraction
of the cart-pole force model (T), the largest of Table II. Timings in the paper come from the laptop
named in `results/numbers.tex` (`env.cpu`), with 14 threads (`results/environment.txt`); on a
machine with fewer cores the run takes longer.

## From the paper to the code

| Paper | Script | Outputs in `results/` |
|---|---|---|
| Section VI-A, pendulum swing-up; Table I, Fig. 4 | `experiments/E1_pendulum.jl`, models in Dionysos's `problems/Pendulum/` | `params_e1.md`, `table_e1.tex`, `fig_e1.pdf`, `e1_{T,F}.gif`, numbers `e1.*` |
| Section VI-B, cart-pole swing-up; Table II, Fig. 5 | `experiments/E2_cartpole.jl`, models in Dionysos's `problems/CartPendulum/point_mass_cart_pendulum.jl` | `params_e2.md`, `table_e2.tex`, `fig_e2.pdf`, `e2_F.gif` (no animation for (T), which fails), numbers `e2.*` |

Every abstraction and controller is computed by Dionysos's uniform grid abstraction solver;
the scripts set the problem, run the solver, and replay its controllers.

The numbers, the two tables and the two plots of the paper are copies of `results/numbers.tex`,
`table_e1.tex`, `table_e2.tex`, `fig_e1.pdf` and `fig_e2.pdf` from the final run. A new run
rewrites these files, so they can be compared with the paper. `fig_e*.png` are previews of the
plots, and the `.gif` files animate the certified swing-ups; the paper does not use them.

## Parameters

The `params_*.md` files are the full specification behind the shortened text of Section VI.
They are written by the scripts from the constants used in the computation, so they cannot
drift from the results:

- `params_e1.md`: pendulum physics, torque limit, impulse budget, sampling period, speed limit,
  target; the grid and commands of (F); the grid, torque step and arrival speed of (T).
- `params_e2.md`: cart and pole masses and length, force limit, impulse budget, speed limit,
  sampling period, target; the grid, commands and budget steps of (F); the grid and forces of (T).

## Checks performed by the run

A run stops with an error, or warns, if a claim of the paper fails on its results:

- every controller is replayed on the continuous-time dynamics; the controller of (F) must reach
  the target within the impulse budget (E1, E2), that of (T) must reach the target (E1), and any
  solved (T) must stay within its force limit (E2);
- the steps of (F) from the bottom must be the fewest possible, checked against an independent
  breadth-first search on the augmented abstraction (E1, Completeness paragraph);
- every cell must be controllable from rest (E1).

## Environment

- Julia 1.12 (the results come from 1.12.5); the packages are pinned by `Project.toml` and
  `Manifest.toml`.
- Dionysos is loaded from source, two folders up (`path = "../.."` in the manifest), and the
  problem definitions from its `problems/` folder. Locally, this is the enclosing Dionysos.jl
  repository. On CodeOcean, `environment/postInstall` clones Dionysos.jl at the commit of the
  paper's results into `/opt/Dionysos.jl`, where this capsule's manifest finds the clone itself,
  then installs and precompiles the packages once, when the environment is built; a Reproducible
  Run therefore needs no network.
- Each run records the date, the Dionysos commit, the Julia version and the machine in
  `environment.txt`, and the Julia version and CPU in `numbers.tex`, next to the results.

## Layout

```
run, run_all.jl           entry points
common.jl                 shared helpers: results folder, numbers, tables, environment
Project.toml, Manifest.toml
experiments/              one script per experiment (see the table above)
environment/postInstall   CodeOcean build step: Julia, Dionysos.jl at a fixed commit, packages
results/                  everything the paper quotes, plus previews, regenerated by each run
```
