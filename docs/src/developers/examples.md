# Adding an example

Examples are the part of the documentation people actually read, and the part most likely to drift.
This page is the contract they follow. It exists so that adding one is a copy-paste rather than a
design exercise, and so that ten pages written by different people still read as one document.

## Where it goes

```
docs/src/examples/
  getting_started.jl        the guided first contact (pendulum) — not a normal example
  jump/                     "Examples" — problems stated through the JuMP front-end
  solvers/                  "Solver families" — one page per algorithm, driven through MOI
```

The two tiers answer different questions. `jump/` shows **how you state a problem**; `solvers/`
shows **how many algorithms will accept it**. An example belongs in `jump/` if the front-end can
express it, and in `solvers/` if it exists to demonstrate a particular optimizer — including when
that optimizer needs inputs the front-end cannot build (PWA systems, ellipsoid templates: see
[`src/wrapper/README.md`](https://github.com/dionysos-dev/Dionysos.jl/blob/master/src/wrapper/README.md) §11).

Drop a `.jl` file in either folder and it appears in the navigation — no edit to `docs/make.jl` is
required. To place it in a particular position, add its basename to `ORDER` in `docs/make.jl`;
files not listed there are appended alphabetically, and the build prints a notice saying so.

**Filenames are lowercase, `snake_case`, no spaces** — they become the page URL.

## The title line

The first line is the title, and it carries two things at once:

```julia
# # DC-DC converter: keeping a switching converter inside a safe band
```

The sidebar shows everything **before the colon** ("DC-DC converter"); the page heading shows the
whole line. Keep the left half short enough to fit a navigation entry.

## Skeleton

```julia
# # <Title>: <what is achieved>
#
# | | |
# |:--|:--|
# | **System**        | 3-D continuous, nonlinear |
# | **Specification** | reach-avoid |
# | **Solver**        | uniform grid abstraction |
#
# <2–5 sentences: the physical system, its dynamics in LaTeX, the control objective.>

using StaticArrays, JuMP, Plots
using Symbolics, MathOptSymbolicAD          # only when the dynamics are JuMP expressions

using Dionysos
const DI = Dionysos
const UT = DI.Utils
const MP = DI.Mapping
const AB = DI.Optim.Abstraction

# ## The model
model = Model(Dionysos.Optimizer)
# variables → dynamics → specification, in that order

# ## The abstraction
# jacobian_bound / growthbound_map, time_step, approx_mode, state_grid, input_grid

# ## Solving
optimize!(model)

# ## Results
termination_status(model)

# ## Closed loop
trajectory = Dionysos.simulate(model, x0; nsteps = 100)

# ## Visualisation
# one static figure, then one animated dashboard
```

The three-row table at the top is what makes the examples index navigable: it feeds the capability
matrix, and its **System** row uses the vocabulary of the [manual](@ref Overview) — *continuous*,
*discrete-time*, *switched*, *hybrid*.

## Rules

1. **Write the model inline.** Never `include` it from `problems/`. An example has to read
   top-to-bottom as one story — that is the entire point of the front-end. (`problems/` keeps serving
   the driver scripts in `examples/` and `research/`; that is a different audience.)

2. **Reuse the drawing, not the model.** The `system_plot!` closure in `problems/<Name>/` is genuinely
   shared code and belongs there. Duplicating a pendulum renderer into the documentation would rot.

3. **Simulate with `Dionysos.simulate`.** Not `get_attribute(model, "discrete_time_system")` plus
   `ST.get_closed_loop_trajectory`: `simulate` already derives the stopping criterion from the
   specification — reaching the target, or leaving the safe set — including the reach-avoid case.

4. **Spend the comment budget on choices, not on mechanics.** Prose at section boundaries, plus one
   sentence wherever a *decision* is non-obvious: why this grid step, why this growth bound, why
   `INNER`. A comment that restates the next line (`# Define the control variables`) is noise; delete
   it.

5. **Declare the aliases you use, and only those.** `const SY = DI.Symbolic` in a file that never
   says `SY.` is a small lie about what the example needs.

6. **Assert the result.** Documenter fails the build when an example *throws*, so the documentation
   build already checks that every page still runs. What it cannot check is whether the example still
   does what its prose claims — a page keeps building perfectly well while synthesizing a controller
   that no longer reaches its target.

   Close that gap with `#src` lines. Literate strips them from both the published page and the
   generated script, so they do not run during the documentation build; they run when the example is
   executed **as a file**, which is how you check an example you have just edited:

   ```julia
   using Test                                                        #src
   @test is_solved_and_feasible(model)                               #src
   @test last(ST.states(trajectory)) ∈ concrete_problem.target_set   #src
   ```

   ```
   julia --project=docs docs/src/examples/jump/path_planning.jl
   ```

   Run that before opening a pull request that touches an example.

   **Each `#src` line must be a complete statement on its own line.** Literate strips lines
   *individually*, so an assertion split across several lines leaves the un-marked lines behind and
   the generated page fails to parse — with a syntax error pointing at the markdown, not at your
   example. Note that the formatter will happily split a long call for you, so keep these lines
   short enough to survive it:

   ```julia
   T_end, mode_end = last(ST.states(trajectory)), last(ST.modes(trajectory))   #src
   @test PR.satisfies(concrete_problem.target_set, T_end, mode_end)           #src
   ```

7. **End with an animation.** See below.

8. **Suppress results you are not showing on purpose.** Documenter prints the value of each block's
   last expression into the page. A controller or a value function over a fine grid is *hundreds of
   kilobytes* of text, and Documenter hard-errors above 200 KiB per page — so a missing semicolon
   fails the build with a size error that says nothing about which line caused it:

   ```julia
   concrete_controller = get_attribute(model, "concrete_controller");   # note the ;
   ```

   Display an object only when a reader learns something from it — `termination_status(model)`, or
   `concrete_problem.system` to show what the model lowered to.

9. **No timing `println` blocks.** They date instantly and say nothing a reader can act on.

## Visualisation

Every example closes with a static figure and an animated dashboard. The dashboard comes from
[`animate_trajectory_dashboard`](@ref Dionysos.animate_trajectory_dashboard), which pairs a drawing of
the physical system with the state, input and (for hybrid models) mode channels of the trajectory:

```julia
anim = Dionysos.animate_trajectory_dashboard(
    system_plot!, trajectory;
    xdims = (1, 2), udims = (1,), Δt = 0.1, frame_step = 2,
)
gif(anim; fps = 10)
```

**The `gif(anim)` wrapper is required.** With no `filename`, the dashboard returns a
`Plots.Animation`, which defines no `show` method for any MIME and therefore displays nothing at all.
Only the `AnimatedGif` that `gif()` returns carries the methods Documenter needs. Passing `filename`
to the dashboard instead is also wrong here: that returns a `String`.

GIF rather than MP4 is deliberate. An `AnimatedGif` holding a `.gif` is showable as `image/gif`, so
Documenter writes it to a file next to the page and links it, keeping the page HTML small. An `.mp4`
only offers `text/html`, so the whole base64 video is inlined into the page and counts against
Documenter's 200 KiB per-page limit — which is a hard build failure, not a warning.

Keep an animation **under about 2 MB**. The knobs, cheapest first: `frame_step` (render every n-th
frame), the trajectory length, `fps`, and the panel `size`.

## Build cost

Every example executes on every documentation build. There is no per-page time limit — one chosen
without measurement would quietly force the grids, and therefore every published figure, to change —
but the cost is real and worth minding. Prefer the levers that change nothing a reader sees:
`frame_step`, not solving the same problem twice on one page, a shorter closed loop. Coarsening a grid
changes the controller and every figure with it, so if an example needs that, say so in its prose.

## Checklist

- [ ] `snake_case.jl` in `jump/` or `solvers/`
- [ ] title line with a short left half, and the three-row summary table
- [ ] sections in order: The model / The abstraction / Solving / Results / Closed loop / Visualisation
- [ ] model written inline; only `system_plot!` reused from `problems/`
- [ ] `Dionysos.simulate` for the closed loop
- [ ] `#src` assertions on the result
- [ ] static figure **and** `gif(anim)` dashboard, under ~2 MB
- [ ] `julia --project=docs docs/src/examples/<tier>/<name>.jl` runs standalone
- [ ] a row added to the capability matrix on the examples index
