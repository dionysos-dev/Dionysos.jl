# Data-driven abstraction: successors from simulated sample points

This document describes `SampledAbstraction.jl` in `problems/BipedRobot/4D_model_underactuated/`
(module `SampledAbstraction`): the "data-driven" way of computing the successors of a cell, to compare
with center simulation, the global growth bound, the local growth bound and the certified RK4
successors. Nothing in it is specific to the biped: it takes any black-box system, in continuous or
discrete time, in any dimension, and was tested on the simple and double pendulum.

| File | Module | One-line description |
| :--- | :--- | :--- |
| `SampledAbstraction.jl` | `SampledAbstraction` | The successors of a cell are the cells of the images of `N` simulated points of it. The points are random (`:random`) or a regular design (`:uniform`). Cache files, Slurm chunks and provenance as for the other abstractions. |
| `tests/test_sampled_abstraction.jl` | | 480 checks on the simple and double pendulum. |

The other modules are described in
[`symbolic_multibody_and_local_abstraction.md`](symbolic_multibody_and_local_abstraction.md); the
experiment that uses the biped in
[`biped_4d_normal_form_experiment.md`](biped_4d_normal_form_experiment.md).

**Contents**

1. [The idea, and what it is not](#1-the-idea-and-what-it-is-not)
2. [Using it](#2-using-it)
3. [Choosing the points: `:random` and `:uniform`](#3-choosing-the-points-random-and-uniform)
4. [Reference](#4-reference)
5. [How it plugs into Dionysos](#5-how-it-plugs-into-dionysos)
6. [Tests and known limitations](#6-tests-and-known-limitations)

---

## 1. The idea, and what it is not

### 1.1 The idea

For a cell `R` of the state grid and an input `u`, the abstraction needs the set of cells that the
system can reach from `R` in one sampling period. The other methods **bound** that set: a growth bound
inflates the image of the center by a radius, the certified RK4 encloses the image with interval
arithmetic. Here nothing is bounded:

1. take `N` points `x₁ … x_N` of the cell `R` (§3);
2. simulate each of them with the RK4 map of the system, `yᵢ = F(xᵢ, u)`;
3. find the cell that each `yᵢ` belongs to (the quantizer of the state grid, the one that Dionysos uses
   to say which abstract state a concrete state is in);
4. the successors of `(R, u)` are the **union** of those cells.

If one of the `yᵢ` is outside the mapped domain (or is not a finite number), the input is disabled in
that cell: no transition at all, as for every over-approximation in Dionysos.

Only the ability to **evaluate** the dynamics is needed: no symbolic expression, no Jacobian, no
interval arithmetic. That is the point of the method, and its price: the result is not guaranteed.

### 1.2 What it is not

* **It is not certified.** The union of the cells of `N` images need not contain the image of every
  point of `R`. It is an over-approximation in the empirical sense only; the more points, the closer
  to the true successors (§3.3 measures how close). With random points, the probability that a
  further random point of the cell leaves the set is what a probabilistic certificate would bound;
  computing such a bound is **not done** here.
* **It is not the set of cells that the image touches.** The union of the cells *of the images of
  the points* is not the union of the cells met by the convex hull of those images: a cell between
  two images is not listed. The other methods list it.

### 1.3 The five methods side by side

All five are built by the same optimizer of Dionysos (`UniformGridAbstraction`) on the same grids, and
return the same kind of optimizer with the abstraction cached, so the control synthesis that follows
does not change.

| Method | How it is selected | What it needs | Guarantee |
| :--- | :--- | :--- | :--- |
| Center simulation | `approx_mode = CENTER_SIMULATION` | the dynamics as a function | none (one point) |
| Global growth bound | `approx_mode = GROWTH`, `jacobian_bound = u -> L` | a bound `L` of the Jacobian over the state space | sound for the ODE |
| Local growth bound | `GROWTH`, `jacobian_bound = cellwise_bound(SymbolicDynamics(system))` | symbolic dynamics | sound for the ODE |
| Certified RK4 | `build_abstraction(CertifiedRK4GrowthBound(...))` | symbolic, branch-free dynamics | sound for the Float64 RK4 map |
| **Sampled** | `build_sampled_abstraction(...; sampling = :random / :uniform)` | the dynamics as a function | none; empirical (§1.2) |

Per pair `(cell, input)`, measured on the grids below with one thread on a laptop (indicative: they
depend on the machine and on the grid):

| | Simple pendulum | Double pendulum |
| :--- | ---: | ---: |
| cells × inputs | 2 989 × 9 | 9 801 × 3 |
| center simulation | 0.1 µs | 0.7 µs |
| global growth bound | 0.3 µs | 445 µs (no transition survives) |
| local growth bound | 4.2 µs | 14 µs |
| certified RK4 | 73 µs | 467 µs |
| sampled, random, N = 16 | 2.1 µs | 11 µs |
| sampled, random, N = 64 | 7.3 µs | 54 µs |

The simple pendulum has `ẍ = −g sin θ + u`, a grid of `(0.1, 0.2)` and inputs every 1.5, the double
pendulum is `problems/Pendulum/double_pendulum.jl` on a grid of `(0.5, 0.5, 1, 1)` with inputs
`{−11, 0, 11}`, both with a sampling period of 0.1 and 0.05. The cost of the sampled method is `N`
simulations and `N` cell lookups, so it grows with `N` and with the cost of the dynamics. On the
double pendulum the global growth bound takes long because its radius is larger than the domain, so
that every image box is enumerated cell by cell before the input is found to leave it.

### 1.4 What Dionysos already has

Dionysos has `approx_mode = RANDOM_SIMULATION` with `n_samples = N`, which is the `:random` version of
this method. This module exists because that mode cannot serve a comparison of methods or a
computation on a cluster:

| | `RANDOM_SIMULATION` | `SampledAbstraction` |
| :--- | :--- | :--- |
| Points | drawn from the global random generator: no two runs agree, and neither do two numbers of threads or two Slurm chunks | a function of the seed and of the cell only: the same in every run, order, thread, chunk and Julia version |
| Same points for every input | no (redrawn for each pair) | yes |
| Cost of drawing 16 points of a 2-D cell | about 8.5 µs (`LazySets.sample`) | about 20 ns |
| Cache files and chunks | no record of how they were computed | provenance: `N`, seed, sampler, period, dynamics (§4.5) |
| Regular design of the points | no | `sampling = :uniform` (§3.2) |

Both give the same abstraction up to the randomness: on the simple pendulum with 16 points their
numbers of transitions differ by less than 1%.

---

## 2. Using it

### 2.1 Loading

```julia
using Symbolics, IntervalArithmetic      # loaded by BipedAbstraction (certified RK4), and the interval extension of Dionysos
import Dionysos, MathOptInterface as MOI
using StaticArrays

dir = joinpath(dirname(dirname(pathof(Dionysos))), "problems", "BipedRobot", "4D_model_underactuated")
include(joinpath(dir, "SampledAbstraction.jl"))      # also includes BipedAbstraction.jl, unless it is loaded
using .SampledAbstraction
import .BipedAbstraction as BA                       # the other builder, for the comparison
const MP, SY, ST = Dionysos.Mapping, Dionysos.Symbolic, Dionysos.System
const UGA = Dionysos.Optim.Abstraction.UniformGridAbstraction
```

`SampledAbstraction.jl` builds on `BipedAbstraction.build_abstraction`: if `BipedAbstraction` is
already defined where the file is included, that copy is used (it must be one copy, because the
provenance of a sampled abstraction is registered in it); otherwise the file includes
`BipedAbstraction.jl` itself. The packages are those of `BipedAbstraction` (`Dionysos`, `Symbolics`,
`IntervalArithmetic`, `CoreMath`, `MathematicalSystems`, `LazySets`, `StaticArrays`,
`MathOptInterface`, `JLD2`); the sampling itself needs none of the symbolic ones.

### 2.2 One abstraction

```julia
include(joinpath(dirname(dirname(pathof(Dionysos))), "problems", "Pendulum", "simple_pendulum.jl"))
system = SimplePendulum.system()                       # ẋ = (ω, −g sin θ + u), continuous time

state_grid = MP.GridFree(SVector(0.0, 0.0), SVector(0.1, 0.2))
input_grid = MP.GridFree(SVector(0.0), SVector(1.5))

opt = build_sampled_abstraction(
    system; state_grid, input_grid,
    time_step = 0.1,                    # continuous time only: the sampling period
    nsamples = 16, seed = 1,            # 16 random points per cell
)
abstract_system = MOI.get(opt, MOI.RawOptimizerAttribute("abstract_system"))
SY.get_n_transitions(abstract_system)   # 77 209
```

The call is `BipedAbstraction.build_abstraction` for this approximation: the same grids, the same
optimizer attributes (`print_level`, `execution_backend`, `state_input_filter`,
`use_implicit_mapping`, ...), the same optimizer returned with the abstraction cached, so that a
control problem is set on it afterwards exactly as for the other methods. A **discrete-time** system
is given the same way, without `time_step`. The concrete system for simulation and synthesis is the RK4
map that was simulated, `ST.get_system(approx)` (the optimizer's `discrete_time_system`): the same map
that `CENTER_SIMULATION` uses, so one point at the center gives that abstraction exactly.

The same with the regular design, 17 points being the center and the 16 corners of a four-dimensional
cell (for the pendulum, 5 points are the center and the 4 corners):

```julia
opt = build_sampled_abstraction(
    system; state_grid, input_grid, time_step = 0.1, nsamples = 5, sampling = :uniform,
)
```

### 2.3 The five methods on the same grids

```julia
grids = (; state_grid, input_grid)
tstep = 0.1
include(joinpath(dir, "LocalGrowthBound.jl")); using .LocalGrowthBound: SymbolicDynamics, cellwise_bound

center    = BA.build_abstraction(system; grids..., time_step = tstep, approx_mode = UGA.CENTER_SIMULATION)
global_gb = BA.build_abstraction(system; grids..., time_step = tstep, approx_mode = UGA.GROWTH,
                                 jacobian_bound = SimplePendulum.jacobian_bound())
local_gb  = BA.build_abstraction(system; grids..., time_step = tstep, approx_mode = UGA.GROWTH,
                                 jacobian_bound = cellwise_bound(SymbolicDynamics(system)))
certified = BA.build_abstraction(BA.CertifiedRK4GrowthBound(system, tstep); grids...)
sampled   = build_sampled_abstraction(system; grids..., time_step = tstep, nsamples = 64, seed = 1)
```

On this grid the sampled abstraction lies inside the three over-approximations: wherever one of them
enables an input, the sampled one does too, and its successors are among theirs (the tests check it
with 64 random points and with the 25 points of the uniform lattice). The numbers of transitions go
from the center simulation, through the sampled ones, to the certified RK4, the local and the global
growth bound.

### 2.4 Cache file and Slurm array

They are those of `BipedAbstraction` (see its documentation), with the provenance of the sampled
approximation stored and checked:

```julia
info = (; tstep, h = [0.1, 0.2])
# a cache file: computed once, loaded afterwards
build_sampled_abstraction(system; grids..., time_step = tstep, nsamples = 16, seed = 1,
                          abstraction_file = "pendulum_sampled.jld2", abstraction_info = info)

# a Slurm array: one task per chunk, then the merge
backend = SY.SlurmArrayBackend(nchunks, chunk_id, outdir, :contiguous, true)
build_sampled_abstraction(system; grids..., time_step = tstep, nsamples = 16, seed = 1,
                          execution_backend = backend)                       # task `chunk_id`
build_sampled_abstraction(system; grids..., time_step = tstep, nsamples = 16, seed = 1,
                          merge_chunks = (; outdir, nchunks),
                          abstraction_file = "pendulum_sampled.jld2", abstraction_info = info)
```

A file or chunk computed with another `nsamples`, `sampling`, `seed` (when the sampling is random),
sampling period, number of RK4 substeps or dynamics is **refused**, whatever `abstraction_info` says
(§4.5); so is a file computed by `RANDOM_SIMULATION`, which has no record, and a growth bound's.
Because the points of a cell depend only on the seed and the cell, the chunks of an array agree with
the whole computation.

---

## 3. Choosing the points: `:random` and `:uniform`

The flag `sampling` of `SampledSimulation` and `build_sampled_abstraction` chooses the points of a
cell: `:random` (the default), `:uniform`, or a sampler of your own (§4.3). Both are given as a name
(a `Symbol` or a `String`). `nsamples` is the number of points; with `:uniform` it includes the
center.

### 3.1 `:random`: `RandomSampler`

`N` points independent and uniformly distributed in the box of the cell, a Monte Carlo sample. It is
the choice for which a probabilistic statement can be made: the points are independent and uniform,
which is the hypothesis such a statement needs.

The points are **not** drawn from a stateful random number generator. Point `k` of a cell is a pure
function of the seed, the coordinates of the center of the cell and `k`, through a counter-based
generator (SplitMix64, started from a FNV-1a hash of the seed and of the center). Hence:

* the same points in every run, in any order of the cells, with any number of threads, in any chunk of
  a Slurm array, in any Julia version and on any machine;
* the same points for every input of a cell;
* the points for `N` are the first `N` of those for more (so the abstraction can only grow with `N`,
  except that an input is disabled as soon as one point escapes);
* about 20 ns for the 16 points of a two-dimensional cell.

### 3.2 `:uniform`: `UniformSampler`

A regular, deterministic design, the same in every cell (relative to the cell): **no randomness and no
seed**, so there is no probability to speak of. The points are, in this order:

1. **the center** of the cell;
2. **the `2ⁿ` corners** (`n` is the number of states);
3. the points of **finer and finer lattices**: those of the lattice with 3 points per side
   (`{−1, 0, 1}ⁿ` in cell units) that are not yet there, then those of the lattice with 5, then 9, 17...
   points per side. Each level is completed before the next one starts.

Hence the numbers of points for which the design is a regular lattice (the center being a lattice point
from the second level on):

| `n` | center | + corners | 3 per side | 5 per side | 9 per side |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 2 | 1 | **5** | 9 | 25 | 81 |
| 3 | 1 | **9** | 27 | 125 | 729 |
| 4 | 1 | **17** | 81 | 625 | 6 561 |

`nsamples = 1` is the center simulation. The numbers in between end a level halfway: the points of
that level are then taken **one at a time, each the farthest from all those already chosen** (the
Euclidean distance, in units of the half-width of the cell in each coordinate; ties are broken by the distance to the points chosen in the same
level, then by the order of the lattice), so a partial level is spread over the cell instead of
being a corner of it. In two dimensions:

```
N = 1   (0, 0)                              the center
N = 2…5 (−1,−1) (1, 1) (1,−1) (−1, 1)       the corners: a diagonal pair first
N = 6…9 (0,−1) (0, 1) (−1, 0) (1, 0)        the midpoints of the sides: opposite pairs first
N = 10… the points of the 5×5 lattice that are not there yet
```

The design for `N` points is the first `N` of that for more, like the random one. The design is
computed once, when the `SampledSimulation` is made, and not per cell.

**Limits.** The lattice of a level is listed, so a design is refused (`ArgumentError`, which says to use
`sampling = :random`) if it needs a lattice of more than 5 million points or more than about
`2·10⁹` coordinate comparisons to choose the points: for instance 100 000 points in ten dimensions. In
four dimensions `N = 6 561` takes 0.2 s.

**Corners.** The corners lie on the boundary of the cell, which belongs to its neighbors too. Their
successors are those of states of the *closed* cell, which is what an over-approximation of a cell is.
A state exactly on the boundary of two cells is assigned by the quantizer of the grid (ties to even),
so its source cell may be the neighbor, which does not matter for the abstraction. An *image* exactly
on a boundary is assigned the same way.

### 3.3 Which one, and how many points

The quantity to look at is how often a **fresh random point** of a cell goes to a cell that the
abstraction does not list (it is what an empirical coverage means; it says nothing about states
outside the cells). Measured on 3 000 pairs of the grids of §1.3, with 300 fresh points each (400 for
the random-only series below):

| | `N` | random | uniform |
| :--- | ---: | ---: | ---: |
| simple pendulum | 5 | 0.160 | 0.019 |
| | 9 | 0.080 | 0.001 |
| | 25 | 0.019 | 0.000 |
| | 81 | 0.003 | 0.000 |
| double pendulum | 17 | 0.160 | 0.050 |
| | 81 | 0.033 | 0.001 |

and, for the random points only, for `N = 1, 4, 16, 64, 256`: 0.57, 0.20, 0.035, 0.005, 0.0006
(simple pendulum) and, for `N = 1, 4, 16, 64`: 0.77, 0.45, 0.17, 0.04 (double pendulum). Four
observations:

* The regular design covers much better than random points for the same `N`, because the corners are
  where the extreme successors of a smooth map are.
* Its abstraction is **larger** (more transitions, for the same reason: 243 217 against 119 831 on
  the double pendulum for `N = 17`) and has **fewer enabled pairs** (15 583 against 18 891), because a
  corner leaves the domain more often than a random point.
* The random miss rate is below `1/(N+1)` from `N = 16` on for the simple pendulum (0.035 against
  0.059), but above it for the double pendulum (0.17 against 0.059 for `N = 16`, 0.04 against 0.015 for
  `N = 64`): in four dimensions the successors of a cell spread over many more cells, and the number
  of points that a given coverage needs grows with the number of target cells, so a certificate has to
  account for it.
* The cost per pair is the same for both, `N` simulations.

---

## 4. Reference

### 4.1 `build_sampled_abstraction(system; state_grid, input_grid, nsamples, sampling = :random, seed = 0, include_center = false, time_step = nothing, nsystem = 5, attributes...) -> optimizer`, and `build_sampled_abstraction(approx::SampledSimulation; state_grid, input_grid, attributes...)`

The transition system of `system` (a `ConstrainedBlackBoxControlContinuousSystem`, with its
`time_step`, or a `ConstrainedBlackBoxControlDiscreteSystem`, without it), by the successors of §1.1.
`attributes` are forwarded to `BipedAbstraction.build_abstraction`, hence to the optimizer, with the
`abstraction_file`, `abstraction_info` and `merge_chunks` of the cache and of the Slurm chunks.
`overapproximation_map`, `approx_mode`, `approx`, `time_step` (with a `SampledSimulation`) and
`nsystem` are set by the function and cannot be overridden (`ArgumentError`). `time_step` is an error
for a discrete-time system and required for a continuous-time one.

### 4.2 `SampledSimulation(system, tstep; nsamples, sampling = :random, seed = 0, include_center = false, nsystem = 5)`, `SampledSimulation(system; ...)`

The approximation itself, a `DiscreteTimeSystemOverApproximation` of Dionysos (**in the empirical
sense only**). For a continuous-time system, `ST.discretize_continuous_system` gives the RK4 map
(`nsystem` substeps of `tstep`, as Dionysos discretizes every continuous system); `ST.get_system` returns
that discrete-time system, `ST.get_over_approximation_map` the function that the `USER_DEFINED` mode of
the optimizer takes, and `ST.reach_set(approx, cell, u, nothing)` the `SampledImage` of a cell.

* `nsamples`: the number of points (for `:uniform`, at least 1: the first is the center).
* `include_center = true`: the center is simulated in addition to the `nsamples` points
  (`nsamples = 0, include_center = true` is the center simulation). `:uniform` refuses it, since it has
  the center already.
* `seed`: for `:random` only; ignored (and not recorded) by `:uniform`.
* `sampling`: `:random`, `:uniform`, or a sampler (§4.3). Anything else is an `ArgumentError`.
* The state dimension is a type parameter of the approximation, so that the points of a cell are
  `SVector`s of known size.

### 4.3 `RandomSampler()`, `UniformSampler()`, and writing a sampler

A sampler is any callable

```julia
sampler(cell, nsamples::Int, key::UInt64) -> AbstractVector of points
```

where `cell` is a `LazySets.AbstractHyperrectangle` and `key` identifies the seed and the cell (a
`UInt64` from the hash above). To be reproducible it must be a function of `(cell, nsamples, key)`
alone. Its `string` is stored in the provenance, so it must identify it (a callable struct with a
`Base.show` is enough; the name of a closure may change from one session to the next). Three optional
functions
describe a sampler to the approximation:

* `SampledAbstraction.prepare(sampler, n, nsamples) -> sampler`: called once with the state dimension
  and the number of points, to compute what depends on nothing else (the design of `UniformSampler`).
  The default returns the sampler.
* `SampledAbstraction._uses_seed(sampler)`: `true` by default; `false` means that the seed is not part
  of the provenance.
* `SampledAbstraction._contains_center(sampler)`: `false` by default; `true` means that
  `include_center` would simulate the center twice and is refused.

### 4.4 `SampledImage(points)`

The set that the over-approximation map returns: a `LazySets.LazySet` that holds the images of the
points (`dim`, and `low`/`high`, their bounding box). Its cells are looked up by the method
`MP.get_states_from_set_strict(mapping::MP.GridMapping, ::SampledImage, mode)`: the cell of each
point with `MP.get_state_by_coord`, each cell once. `allin = false` if a point has no cell in the
mapped domain or is not a finite number (a `NaN`, an infinity, a number beyond the integers of the
quantizer), and Dionysos then lists no transition for that cell and input.

### 4.5 Provenance: `get_provenance(approx)`, `KERNEL_VERSION`

`get_provenance(approx)` is the plain data that identifies how the successors are computed:

```julia
(; approximation = "SampledSimulation", kernel_version, nsamples, seed,   # seed = nothing for :uniform
   include_center, sampler,                                                # sampler = string, "RandomSampler()"
   tstep, nsystem,                                                         # nothing for a discrete-time system
   dynamics)                                                               # fingerprint, see below
```

`BipedAbstraction` stores it (with a fingerprint of the cells and inputs behind the ids) in cache
files and in a small file next to each chunk, and refuses anything that carries another record or none.

* **`dynamics`** fingerprints a black-box system, which has no expression to print, by its **values**:
  the discrete-time map at 16 states and inputs, always the same ones, in the box of the state set
  and of the input set (`[−1, 1]` if a set is unbounded), rounded to ten significant digits so that a
  last-bit difference between two machines does not tell two runs apart. A system written
  differently that computes the same function has the same fingerprint, one whose parameters differ
  does not; a map that is changed only far from those points is not detected.
* **`KERNEL_VERSION`** (currently `1`) must be raised whenever a change can alter the transitions that
  a given provenance returns: the generator, the way a cell is keyed, the uniform design. The test pins
  the dynamics fingerprint of the simple pendulum and the stream key of a cell, to notice an accidental
  change.

---

## 5. How it plugs into Dionysos

Nothing of Dionysos was changed. It uses the extension points that exist already, as `CertifiedRK4`
does:

```
SampledSimulation  ──ST.get_over_approximation_map──►  (cell, u) -> SampledImage
        │                                                      │
        │ USER_DEFINED mode of the UniformGridAbstraction      │ MP.get_states_from_set_strict(GridMapping, SampledImage, ·)
        ▼                                                      ▼
BipedAbstraction.build_abstraction ───► optimizer ──► collect_abstract_transitions!  ──►  cells of the points, if all are
(cache file, chunks, provenance)         (SY kernel for over-approximations)               in the retained states
```

* The kernel for over-approximations of `SY` calls `ST.reach_set`, whose default calls the map; the map
  returns a `SampledImage`; `compute_abstract_transitions_from_set!` asks the mapping for its cells and
  refuses the pair if any is outside the retained states.
* The cell of a point is the one the quantizer of the grid gives, the same lookup as in
  `CENTER_SIMULATION`; the two give the same transitions for one point at the center (tested).
* A point whose image leaves the domain disables the input for the cell: the `RANDOM_SIMULATION` rule
  as well.
* Execution backends (sequential, threaded, distributed, Slurm array) work unchanged; the result does
  not depend on them (tested for the threaded one and for the chunks).

---

## 6. Tests and known limitations

```
julia --project=<env> problems/BipedRobot/4D_model_underactuated/tests/test_sampled_abstraction.jl      # 480 checks, about 30 s
```

`<env>` needs the packages of §2.1 and `Plots` (the pendulum problems of `problems/Pendulum` load it).
The test checks, on the simple pendulum (a grid of `(0.3, 0.5)`, 361 cells and 5 inputs) and the
double pendulum (the grid of §1.3, 9 801 cells):

* the generator against the reference outputs of SplitMix64; the uniformity, independence and
  reproducibility of the random points (mean, variance, halves and correlations within five
  standard deviations; prefixes; the points inside the cell);
* the uniform design: `N = 17` in four dimensions is the center and the 16 corners, the regular
  lattices for `3ⁿ`, `5ⁿ`, `9ⁿ`, prefixes, distinct points, opposite pairs, the farthest-first property
  by brute force, the refusal of a design that is too large;
* the lookup of an image (explicit and implicit mappings; a point outside the domain, a `NaN`, an
  infinity, `10³⁰⁰`);
* the abstractions: reproducible, equal under the threaded backend, monotone in `N`, equal to
  `CENTER_SIMULATION` for one point at the center (uniform `N = 1`, or random with
  `nsamples = 0, include_center = true`), equal to a by-hand simulation of the corners, close to
  `RANDOM_SIMULATION` in size, inside the global and local growth bounds and the certified RK4,
  the same for a discrete-time system;
* the provenance, the cache file and the chunks: reuse without recomputation, refusal for every other
  option, of a file written by `RANDOM_SIMULATION`, and of a chunk of another seed; a uniform file
  does not depend on the seed;
* the arguments (a missing or surplus `time_step`, an unknown `sampling`, `include_center` with
  `:uniform`, ...).

Deliberately breaking the module in a dozen ways (points not simulated, a wider cell, a seed left out
of the provenance, an escaping image ignored, the center not simulated, corners moved inward, the
tie-break removed, ...) made the tests fail each time.

**Known limitations**

* **Not certified** (§1.2): the abstraction can miss a successor; no probability is computed.
* **Cells between images are missed** (§1.2). A cell reached only by points of the cell that no
  sample is near is not listed, however regular the design.
* **The cost is `N` simulations per pair**, so a four-dimensional model with a costly map and a large
  `N` is a long computation; the Slurm array applies as for the other methods. For the biped, which
  has not been run with this method yet, expect the cost per pair to be that of its Float64 RK4 map
  times `N` (the certified interval computation costs about 0.5 ms per pair).
* **Corners escape more** (§3.3): the uniform design disables more pairs than random points do.
* **Reproducibility is that of the points**, not of the floating-point arithmetic of the dynamics:
  another machine may round an image differently, which changes the cell only for an image within
  rounding error of a cell boundary.
* **The dynamics fingerprint is behavioral** (§4.5): a change that does not show at its 16 points is
  not detected; the `abstraction_info` you pass should identify what the provenance cannot see (the
  specification, `state_input_filter`, ...).
* **Periodic mappings** have not been tried: the lookup is the `MP.get_state_by_coord` of Dionysos, but
  the method for a `SampledImage` is written for `MP.GridMapping` and was only tested with explicit and
  implicit grid mappings.
* **`seed` is an `Int`**; the points of a cell are hashed from it and the center of the cell, so
  changing the grid origin changes the points.
