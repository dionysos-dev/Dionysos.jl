# Plan — JuMP wrapper refactor & extension

Tracking issue: **#589** (umbrella for #510 custom dynamics, #511 fit the JC refactoring, #512 full
hybrid support). Prior art: **PR #395** (`ja/ft/hybridsystem`, Nov 2024, stale/dirty) — its *idea* is
kept, its *implementation* is not (§3).

User guide for the target DSL: [`src/wrapper/README.md`](src/wrapper/README.md).

Status: **all phases executed** (branch `jc_issue_4`), including the optional LTL layer.
Revision 2 — every JuMP mechanism claimed here was executed against JuMP 1.30.1 (§9), and
revision 1 was adversarially reviewed; the defects it found are fixed in place and recorded in §10.
Execution notes, including five further defects the work surfaced and one design decision taken
*against* this plan, are in §11.

---

## 1. What the wrapper is, and what it is for

Dionysos' pitch is *describe the system → pick the specification → pick a solver → get a certified
controller*. The **JuMP/MOI front-end is the surface that pitch is delivered on**: it is the only
place where a user who is not a Dionysos developer writes a control problem.

Today the front-end reaches a small fraction of what the library underneath can do. This document is
about closing that gap without creating a second, divergent modelling language.

Two non-negotiable design rules:

* **R1 — the front-end owns no semantics.** It is a *compiler*: JuMP model → `(system, problem)` →
  existing solver. If a concept cannot be expressed as a `MathematicalSystems`/`HybridSystems`
  object plus a `PR.ProblemType`, it does not belong in the wrapper; it belongs in `Problem`/`System`
  first. Where the library is genuinely missing something, §10 says so explicitly rather than
  papering over it in the DSL.
* **R2 — infer, don't declare.** The current wrapper's one good idea is "a variable is a state iff
  it has dynamics". Generalise that principle rather than replacing it with ceremony.

---

## 2. Audit of the current wrapper

Source: `src/MOI_wrapper.jl` (503 lines), `include`d **only** from
`ext/DionysosMathOptSymbolicADExt.jl`, which defines `Optimizer() = SymbolicsOptimizer()`.
`Dionysos.Optimizer` in `src/Dionysos.jl` is an empty stub function.

| # | Limitation | Evidence |
| :- | :--- | :--- |
| L1 | **Lives in an extension.** Nothing in the front-end can be loaded, tested, or documented without `using Symbolics, MathOptSymbolicAD`. | `ext/DionysosMathOptSymbolicADExt.jl:15` |
| L2 | **Time type is a 3-valued enum** `UNKNOWN/CONTINUOUS/DISCRETE`. No hybrid, no clock. | `MOI_wrapper.jl:8` |
| L3 | **Two variable roles only.** `variable_type` = `INPUT` unless the variable has dynamics. `MODE` is declared in the enum and never produced. | `MOI_wrapper.jl:412-418` |
| L4 | **An unused variable silently becomes an input.** A typo'd or leftover `@variable` inflates the input dimension; if unbounded, `_U_` becomes an infinite box and the input grid explodes. | `MOI_wrapper.jl:426-430`, `:361` |
| L5 | **One specification: reach-avoid.** `problem(model)` unconditionally builds `OptimalControlProblem(sys, _I_, _T_, nothing, nothing, Infinity())`. | `MOI_wrapper.jl:388-410` |
| L6 | **Sets are boxes only**, assembled coordinate-by-coordinate from `MOI.Interval`s — even though `MP.get_states_from_set` accepts any bounded `LazySet`. | `MOI_wrapper.jl:391-398` |
| L7 | **Unconstrained `final`/`start` coordinates default to `±Inf`**, silently producing an unbounded target. Path planning works around it with `final(x[3]) in MOI.Interval(-100.0, 100.0)`. | `MOI_wrapper.jl:77-78`, `docs/…/Path planning.jl:69` |
| L8 | **No horizon.** `time` is hardwired to `Infinity()`. | `MOI_wrapper.jl:407` |
| L9 | **`@objective` is accepted and silently discarded**; `state_cost`/`transition_cost` are always `nothing`. | `MOI_wrapper.jl:249-265` vs `:401-408` |
| L10 | **The solver is hardcoded** to `UniformGridAbstraction.Optimizer` in the constructor. | `MOI_wrapper.jl:30`, `:473-481` |
| L11 | **Parsing and lowering are fused** — an `if/elseif` chain over `ScalarNonlinearFunction` heads ending in `dump(func); error("Unsupported")`. No inspectable intermediate model. | `MOI_wrapper.jl:141-218` |
| L12 | **No reset-map type in the library.** **11 near-identical copies** of `struct …ResetMap <: MS.AbstractMap` across `problems/{Thermostat,FlowShopScheduling}/*`, `research/Settings/`, `test/optim/HybridSystemAbstraction/*`. | grep `<: MS.AbstractMap` |
| L13 | **No solution status.** `termination_status(model)` / `is_solved_and_feasible(model)` are unimplemented, although every control solver carries a `success::Bool`. The front-end feels foreign to a JuMP user. | no `MOI.TerminationStatus` in `src/optim/continuous_systems/` or the wrapper |
| L14 | **No simulation path.** Every example drops into `ST.get_closed_loop_trajectory(get_attribute(model, "discrete_time_system"), …)`; the hybrid equivalent has a *different signature* returning a raw tuple. | all three docs examples |
| L15 | **Only 4 consumers exist in the whole repo** — `test/optim/jump_frontend.jl` plus three docs examples. Nothing in `problems/` uses the wrapper. This is both why a refactor is safe and why the safety net must be built first. | grep `Dionysos.Optimizer` |

Preserve: `∂`/`Δ`/`final`/`start` as `JuMP.NonlinearOperator`s, the `x ∉ HyperRectangle` obstacle
syntax (documented `OuterSet` + `parse_constraint_call` piracy), the Symbolics→
`RuntimeGeneratedFunction` codegen with CSE, and transparent raw-attribute forwarding.

---

## 3. What to take from PR #395, and what to drop

**Take — blegat's transition proposal** (2024-11-13), which johnaoga accepted: a transition is a
*first-class object* you attach guards and reset maps to with `@constraint`, plus `add_transition!`
do-block sugar. This is right: it keeps guard and reset linked even when several transitions share a
`(source, target)` pair.

**Drop — the implementation.** The PR: re-introduces the very ambiguity the proposal fixes by keying
transitions `Dict{String, Transition}` on `"$src->$dst"`; pattern-matches ~10 bespoke predicates on
nonlinear expression heads; **hardcodes the thermostat inside the library** (`_hybrid_system` builds
two fixed affine systems and a literal
`guard_constraints = [HyperRectangle([19.0],[Inf]), HyperRectangle([-Inf],[21.0])]`); hardcodes a
Bemporad–Morari + OSQP/HiGHS/Ipopt/Pavito stack that Dionysos does not depend on; pirates
`JuMP._valid_model` (unnecessary — see §4.5) and `Base.rem`; leaves `println` tracing in the hot
path; and predates the entire Problem/Mapping/System refactor (issue #511).

---

## 4. Target architecture

### 4.1 The one structural change: an explicit IR

```
JuMP model ──parse──▶ ModelIR ──lower──▶ (system, problem) ──select──▶ solver ──▶ controller
             (pure)   (data)   (pure)                       (dispatch)
```

`ModelIR` is a plain, dependency-free description of what the user wrote: variables and roles,
modes, transitions, dynamics *expressions*, specification entries, obstacles, attributes. It holds
no `Symbolics` objects and can be printed.

This split buys, in order of importance:

1. **The wrapper moves to `src/`.** Only *compiling a dynamics expression into a callable* needs
   Symbolics; parsing does not. Fixes L1 and makes the front-end unit-testable in the base test env.
2. **Validation happens once, on a complete model.** Errors can name the offending variable and say
   what was expected, replacing `dump(func); error("Unsupported")`.
3. **New features stop touching `add_constraint`.** A new specification kind is one IR field plus
   one lowering method.

### 4.2 File layout

```
src/wrapper/
  README.md               user guide (already written)
  wrapper.jl              module + includes; loaded from src/Dionysos.jl
  variables.jl            VariableInfo, roles, inference
  model_ir.jl             ModelIR + Base.show
  modes.jl                Mode       <: JuMP.AbstractModel
  transitions.jl          Transition <: JuMP.AbstractModel, add_transition!
  operators.jl            ∂ Δ final start; Start/Final/Always/EventuallyAlways/Label sets; ∉ piracy
  specification.jl        spec IR + problem-type inference (+ optional formula normaliser, P9)
  parse.jl                MOI.add_constraint / MOI.set → ModelIR     ← the only pattern matching
  lower_system.jl         ModelIR → MathematicalSystems / HybridSystems
  lower_problem.jl        ModelIR → Problem.ProblemType
  solver_selection.jl     select_solver, supports_problem
  dynamics_backend.jl     AbstractDynamicsBackend + evaluator / user-function backends
  results.jl              TerminationStatus / PrimalStatus / simulate
  optimizer.jl            Dionysos.Optimizer <: MOI.AbstractOptimizer

ext/DionysosMathOptSymbolicADExt.jl   ← shrinks to: SymbolicADBackend only
```

`src/MOI_wrapper.jl` is deleted at the end of Phase 1.

### 4.3 Variable roles (L3, L4)

Roles: `STATE`, `CLOCK`, `INPUT`, and — later, optional — `PARAMETER`, `DISTURBANCE`.

> **There is no `MODE` role.** Revision 1 carried one, inherited from the current wrapper's
> never-produced enum value and PR #395's `mode_variable`. The §4.5 design has no mode *variable*
> at all: `@mode` returns a scope object. Dropped.

Inference runs **once at `optimize!`**, over the complete model:

| Rule | Condition | Role | Phase |
| :- | :--- | :--- | :-: |
| I1 | appears on the LHS of `∂(·) == …` / `Δ(·) == …` in **any** mode | `STATE` | 1 |
| I2 | its dynamics is exactly `∂(t) == 1` or `∂(t) == 0`, and it appears in no other dynamics RHS | `CLOCK` (a `STATE` specialisation) | 6 |
| I3 | appears on some dynamics RHS, none of the above | `INPUT` | 1 |
| I4 | appears nowhere | **error** naming the variable | 1 |
| I5 | JuMP-fixed (`@variable(model, p == 3.0)`), appears only on RHSs | `PARAMETER` (inlined, not gridded) | later |
| I6 | — | `DISTURBANCE`, via explicit `set_role!` only | later |

**Rule I4 is the valuable one** and the fix for L4: today a leftover variable silently enlarges the
input grid. It is ~10 lines and ships in Phase 1.

I5/I6 are deliberately deferred: no issue or user report asks for them, and the escape hatch
(`set_role!`) covers the cases inference cannot see (e.g. dynamics supplied as a raw Julia function).

**Ordering.** State/input vector order remains declaration order, but becomes *inspectable*:
`get_attribute(model, "state_variables")` returns the ordered `VariableRef`s.

### 4.4 Specifications (L5–L8)

Three JuMP mechanisms, three meanings, never overloaded:

| Mechanism | Carries | Cardinality |
| :--- | :--- | :--- |
| `@constraint` | the specification **sets** | many — they accumulate |
| `@objective` | the **cost** (`state_cost`/`transition_cost`) | exactly one |
| `@specification` *(P9, optional)* | a temporal **formula** over named sets | exactly one |

#### The sugar layer (Phase 3 — the main deliverable)

**Two syntactic forms, and the choice is forced by JuMP, not by taste.** `final(x[1]) in
MOI.Interval(a, b)` works because `final` is a `JuMP.NonlinearOperator` and `x[1]` is an
`AbstractJuMPScalar`. The same operator on a *whole vector* does **not** work: with no JuMP scalar
among its arguments `NonlinearOperator` evaluates the underlying function, and
`final(x::Vector{VariableRef})` raises `MethodError: no method matching _final(::Vector{VariableRef})`
(verified). So the vector form wraps the **set**, exactly as the obstacle syntax already does with
`OuterSet`:

| Syntax | Meaning | Lowers to |
| :--- | :--- | :--- |
| `@constraint(model, x in Start(S))` | initial set | `initial_set` |
| `@variable(model, x, start = v)` | singleton initial set (unchanged) | `initial_set` |
| `@constraint(model, x in Final(S))` | ◇ reach | `target_set` |
| `@constraint(model, x in Always(S))` | □ invariant | `safe_set` — but see below |
| `@constraint(model, x in EventuallyAlways(S))` | ◇□ reach-and-stay | `target_set` + `safe_set` |
| `@constraint(model, x ∉ O)` | obstacle (unchanged) | carved out of `X` |
| `@constraint(model, final(x[i]) in MOI.Interval(a, b))` | per-coordinate reach (unchanged) | `target_set` |
| `@constraint(model, start(x[i]) in MOI.Interval(a, b))` | per-coordinate initial (unchanged) | `initial_set` |

`Start`/`Final`/`Always`/`EventuallyAlways` are thin `MOI.AbstractVectorSet` wrappers around any
bounded `LazySet` — box, ellipsoid, ball, polytope, zonotope. Verified: `x in Final(S)` builds
`VectorOfVariables`-in-`Final{S}`, slices included, needing **no new piracy** (unlike `∉`, `in`
already parses).

**Problem type is inferred from which markers appear** — the same spirit as rule I1:

| `Final` | `Always` | `EventuallyAlways` | ⇒ problem |
| :-: | :-: | :-: | :--- |
| ✓ | | | `OptimalControlProblem` (reach) |
| ✓ | ✓ | | `OptimalControlProblem`, avoid-set folded into `X` — see the asymmetry below |
| | ✓ | | `SafetyProblem` |
| | opt | ✓ | `ReachAndStayProblem` |
| | | | `AlternatingSimulationProblem` — *build the abstraction only* |

The empty row matters: "just build me the abstraction and give it back" is impossible from JuMP
today and is the natural entry point for exploring a new system. All of it is overridable via
`set_attribute(model, "problem_type", PR.SafetyProblem)`.

#### The `Always` / `∉` asymmetry — a library fact, not a DSL choice

`x ∉ O` and `x in Always(complement(O))` are logically equivalent but lower differently, and this
must be documented rather than hidden:

* `∉` restricts the **system**: the region is carved out of `X`, so those cells are never created.
  Fewer states, faster abstraction. This is the established SCOTS-style behaviour and today's.
* `Always(S)` constrains the **specification**: states outside `S` stay representable so synthesis
  can reason about leaving them.

The split maps exactly onto the top-level *system vs specification* architecture, so the README
teaches `∉` under "system" and `Always` under "specification".

The reach-avoid row is forced by the library: **`OptimalControlProblem` has no `safe_set` field**
(`system, initial_set, target_set, state_cost, transition_cost, time`). With both `Final` and
`Always` present there is nowhere else to put the safe set, so it folds into `X`. Recorded as a
library gap in §10.

#### Horizon (L8) — the wrapper converts, because nothing else does

Revision 1 claimed `PR.discretize_time` handles this. **It does not**: `discretize_problem` is never
called on the abstraction path — its only callers are
`src/optim/trajectory_generators/composite_trajectory_generator.jl:66` and tests. Every abstraction
solver passes `concrete_problem.time` **verbatim** into the abstract problem
(`UniformGridAbstraction/{optimal_control,safety,reach_and_stay}_problem.jl`), so a horizon in
seconds would be silently consumed as a step count.

Design, therefore:

* **Continuous-time model (`∂`)** — `set_attribute(model, "horizon", 10.0)` is in **seconds**, and
  the *wrapper* converts it to steps with `PR.discretize_time(h, time_step; round_up = PR.horizon_round_up(problem))`
  before constructing the problem. The wrapper knows both `time_step` and the problem type, so it
  can pick the rounding direction correctly.
* **Discrete-time model (`Δ`)** — `"horizon"` is a **step count**, passed through.
* Default `PR.Infinity()` in both cases.
* **Documented caveat**: `SafetyProblem.time` is currently *ignored* by the solver — `MOI.optimize!`
  in `src/optim/discrete_systems/safety_problem.jl` computes a maximal invariant set (an
  infinite-horizon notion) and never reads `time`. Setting `"horizon"` on a safety problem is inert
  today. §10.

#### Costs (L9) — `@objective` is the cost, and only the cost

`@objective` *cannot* carry a specification even if we wanted it to: JuMP requires an
`MOI.AbstractScalarFunction`, and `@objective(model, Min, "F(goal)")` fails with *"The objective
function `F(goal)` is not supported by JuMP"* (verified). Decision settled with the user: specs stay
constraints, `@objective` models only an additional cost function.

```julia
@objective(model, Min, sum(u[i]^2 for i in 1:2))       # → transition_cost
```

*Phase 3*: raise a clear error instead of today's silent discard.
*Later*: map the separable `Min Σ c(x, u)` form onto `transition_cost`.

#### Arbitrary sets (L6) and the `±Inf` fix (L7)

The per-coordinate `MOI.Interval` form stays supported and is assembled into a box **using each
variable's own declared bounds for the coordinates the user did not constrain** rather than `±Inf`.
*This is the one deliberate behaviour change in the plan*; the Path-planning example's
`final(x[3]) in MOI.Interval(-100.0, 100.0)` line becomes deletable.

### 4.5 Modes and transitions (L2, #512)

```julia
off = @mode(model, off)
on  = @mode(model, on)

@constraint(off, ∂(T) == -α * (T - Ta))          # per-mode dynamics
@constraint(on,  ∂(T) == -α * (T - Ta) + β * u)
@constraint(off, u == 0)                         # per-mode input set
@constraint(on, 0.2 <= u <= 1)

add_transition!(model, off => on) do t
    @constraint(t, T <= 19)                      # guard
end
add_transition!(model, on => off) do t
    @constraint(t, T >= 21)
    @constraint(t, Δ(T) == T)                    # reset map (identity ⇒ omittable)
end
```

* **`Mode` and `Transition` subtype `JuMP.AbstractModel`**, so `@constraint(mode, …)` works with
  **zero piracy** — JuMP already defines `_valid_model(::AbstractModel, ::Any) = nothing`
  (`JuMP/src/macros.jl:362`), `model_convert(::AbstractModel, …)` and
  `value_type(::Type{<:AbstractModel}) = Float64`. Only `JuMP.add_constraint` must be implemented.
  Strictly better than PR #395's `JuMP._valid_model` piracy.

  **Verified** with a throwaway `struct Scope <: JuMP.AbstractModel`; `JuMP.moi_function` resolves
  the *parent* model's variable refs without complaint:

  | written on the scope | `moi_function` | `moi_set` |
  | :--- | :--- | :--- |
  | `@constraint(off, ∂(T) == -α*(T-Ta) + β*u)` | `ScalarNonlinearFunction`, `head = :-`, `args[1].head = :∂` | `EqualTo(0.0)` |
  | `@constraint(off, T <= 19)` | `ScalarAffineFunction` | `LessThan(19.0)` |
  | `@constraint(off, 0.2 <= u <= 1.0)` | `ScalarAffineFunction` | `Interval(0.2, 1.0)` |

  One extra method is required for the vector spec form: `Base.broadcastable(m::Mode) = Ref(m)`,
  because `model_convert(::AbstractModel, ::VectorConstraint)` broadcasts over the model
  (`JuMP/src/macros.jl:342`) and JuMP ships `broadcastable` for `GenericModel` only
  (`JuMP/src/JuMP.jl:415`). A method on our own type — not piracy.

* **Guard vs reset is decided by MOI function *type***, not by pattern-matching: a
  `ScalarAffineFunction` (no `∂`/`Δ` head) is a *guard*; a `ScalarNonlinearFunction` whose LHS head
  is `∂`/`Δ` is a *reset map*. The same rule separates per-mode bounds from per-mode dynamics. None
  of PR #395's ten head-pattern predicates are needed.
* **Transitions are objects in a `Vector`, not string keys**, so several may share `(source, target)`.
* **Per-mode bounds are just scoped bound constraints** — this is what lets the wrapper express the
  thermostat's per-mode input sets, which it cannot today.
* Switching type is a kwarg, default `AutonomousSwitching()` (as blegat asked).
* `time_type` becomes **two orthogonal axes** — `CONTINUOUS`/`DISCRETE` from `∂` vs `Δ`, and
  `MONOLITHIC`/`HYBRID` from the presence of modes. PR #395 conflated them into a single
  `@enum TimeType … HYBRID`, which structurally cannot express a *discrete-time* hybrid system.

#### Library prerequisite (L12): `GuardedResetMap`, and it must be arity-aware

Guards and resets in `HybridSystemAbstraction` operate on the **augmented** state, and the arity
depends on whether the mode is clock-lifted
(`src/optim/hybrid_systems/HybridSystemAbstraction/hybrid_system_abstraction.jl:156-183`):

| mode kind | reset called as | guard set lives over |
| :--- | :--- | :--- |
| clock-lifted | `MS.apply(reset, vcat(x, t))`, result split back into `(x, t)` | `(x, t)` — e.g. `UT.box([0.2, 0.0], [1.0, 2.0])` |
| time-free | `MS.apply(reset, x)` | `x` |

So a user's `@constraint(tr, T <= 19)` is a guard over `x` alone, and the lowering must **extend it
across the clock domain** to build the `(x, t)` guard box. `GuardedResetMap` therefore carries the
mode's arity, and mixed timed/plain hybrids are rejected with a clear message — this is already a
documented open limitation of the hybrid solver, not something the wrapper can paper over.

```julia
struct GuardedResetMap{G, F} <: MS.AbstractMap   # MS.stateset → guard, MS.apply → reset
```

with an identity-reset default. This also deletes the eleven copy-pasted `…ResetMap` structs, and
is worth doing regardless of the wrapper.

**Guard representation.** Affine guards become half-spaces → a polyhedron/box. Nonlinear guards
error, because the hybrid abstraction cannot consume them either.

### 4.6 Clocks and timed hybrids — free, via rule I2

`HybridSystemAbstraction` distinguishes a *plain* mode (`(x, k)`) from a *clock-lifted* one
(`ST.VectorContinuousSystem([phys, clock])`, `(x, t, k)`), and `SY.ClockAbstraction` accepts exactly
`A = I` (running) or `A = 0` (frozen). That maps onto existing syntax with **no new operators**:

```julia
@variable(model, 0 <= t <= 50)
@constraint(off, ∂(t) == 1)      # clock runs   → A = I
@constraint(on,  ∂(t) == 0)      # clock frozen → A = 0
```

A time window on a target is a reach constraint on the clock (scalar, so the per-coordinate form
applies):

```julia
@constraint(on, final(T) in MOI.Interval(21.0, 23.0))
@constraint(on, final(t) in MOI.Interval(15.0, 50.0))   # ⇒ PR.TimedSpec(base, 15.0, 50.0)
```

Lowering: a reach constraint on a `CLOCK` variable inside a mode becomes the `TimedSpec` window
rather than a state coordinate; per-mode specs collect into `PR.HybridSpec`. This reproduces
`PR.hybrid_reach_spec(Xs, Ts, Ns)` from ordinary constraints.

### 4.7 Dynamics backends (#510)

```julia
abstract type AbstractDynamicsBackend end
compile_dynamics(backend, ir, mode) -> (f, statedim, inputdim)
```

| Backend | Needs | Notes |
| :--- | :--- | :--- |
| `SymbolicADBackend()` | Symbolics + MathOptSymbolicAD | today's path (fused `RuntimeGeneratedFunction`, CSE). **The production backend**; default when loaded. |
| `NonlinearEvaluatorBackend()` | nothing | `MOI.eval_constraint` on an `MOI.Nonlinear.Evaluator`. **Test/fallback only** — interpreted and allocating, called millions of times during abstraction, so expect it to be far slower. Its purpose is making the front-end loadable and testable with no weak deps, not production use. |
| `UserDynamicsBackend(f)` | nothing | **issue #510**: `set_attribute(model, "dynamics", f)`, also per-mode. Bounds and roles still come from the model. |

**Free win:** `ST.compute_jacobian_bound` already exists in `DionysosSymbolicsExt` (traces the
dynamics, bounds the Jacobian over `X` with interval arithmetic) and
`ContinuousTimeGrowthBound(system; jacobian_bound = nothing)` already falls back to it. With the
symbolic backend active the wrapper should **stop requiring a hand-written `jacobian_bound`** —
removing a correctness footgun from every `GROWTH` example.

### 4.8 Solver selection (L10)

```julia
model = Model(Dionysos.Optimizer)                                      # auto-select
set_attribute(model, "solver", AB.UniformGridAbstraction.Optimizer)    # explicit
```

* `select_solver(system, problem)` — one method per case, extensible like the existing
  `control_solver_for`. `HybridSystem` ⇒ `HybridSystemAbstraction.Optimizer`; otherwise
  `UniformGridAbstraction.Optimizer`.
* `supports_problem(::Type{<:MOI.AbstractOptimizer}, ::Type{<:PR.ProblemType})` — a trait each
  family declares, so an unsupported combination fails **at the front-end**, naming both.
* The solver is chosen *after* the model is known, so raw attributes are **buffered** in insertion
  order and replayed once it is instantiated. This also fixes a latent problem today: attributes
  belonging to a control sub-solver that does not exist yet.
* **Per-mode attributes**: `set_attribute(off, "state_grid", g1)` fills that mode's entry of
  `optimizer_kwargs_dict`, removing the parallel `optimizer_list`/`optimizer_kwargs_dict` vectors
  users hand-build today (compare `hybrid_system_abstraction.jl:77-102`). Model-level attributes
  propagate to every mode unless the mode overrides them.

**Honest scope**: this wires **2 of ~10** solver families. The rest (BemporadMorari, BranchAndBound,
LazyEllipsoids, UniformEllipsoid, PCLF, trajectory generators/certifiers) need inputs the DSL cannot
express — PWA systems, ellipsoid templates, observation regions, RRT parameters. The dispatcher is
extensible; the coverage is not comprehensive, and §8 says so.

### 4.9 Results surface (L13, L14) — making it feel like JuMP

Every control solver already carries a `success::Bool`
(`OptimizerOptimalControlProblem`, `OptimizerSafetyProblem`, the generators/certifiers via
`get_success`). Expose it through the standard MOI attributes:

| `success` | `MOI.TerminationStatus` | `MOI.PrimalStatus` | `MOI.ResultCount` |
| :--- | :--- | :--- | :-: |
| not yet run | `OPTIMIZE_NOT_CALLED` | `NO_SOLUTION` | 0 |
| `true` | `OPTIMAL` | `FEASIBLE_POINT` | 1 |
| `false` | `LOCALLY_INFEASIBLE` | `NO_SOLUTION` | 0 |
| abstraction-only problem | `OPTIMAL` | `NO_SOLUTION` | 0 |

`LOCALLY_INFEASIBLE` rather than `INFEASIBLE` is deliberate and semantically load-bearing: the
abstraction is **sound but not complete**, so failure means "no controller exists *on this
abstraction*" — a finer grid may succeed. Claiming `INFEASIBLE` would assert something the method
cannot prove. `MOI.RawStatusString` carries the human-readable reason.

This makes `termination_status(model)` and `is_solved_and_feasible(model)` work as a JuMP user
expects, for ~20 lines.

**Simulation (L14).** One uniform entry point, hiding the continuous/hybrid signature split:

```julia
traj = Dionysos.simulate(model, x0; nsteps = 100)      # → ST.Trajectory
```

It reads the discrete-time system (or hybrid system) and the concrete controller off the solved
model, dispatches to the right closed-loop engine, and returns the unified channelled
`ST.Trajectory`. Default stopping criterion derived from the problem
(`PR.trajectory_success`-style), overridable via `stopping =`.

---

## 5. Worked example — the thermostat, end to end

The deliverable PR #395 never reached. Compare with
`problems/Thermostat/thermostat_hybrid_time_system.jl` (~90 lines of hand-built `HybridSystem`) plus
the 25-line optimizer-configuration block in the hybrid tests.

```julia
using Dionysos, JuMP, StaticArrays

Ta, α, β = 18.0, 0.1, 2.0

model = Model(Dionysos.Optimizer)
@variable(model, 17 <= T <= 25, start = 21.5)
@variable(model, 0 <= u <= 1)
@variable(model, 0 <= t <= 50)                       # clock (rule I2)

off = @mode(model, off)
on  = @mode(model, on)

@constraint(off, ∂(T) == -α * (T - Ta))
@constraint(on,  ∂(T) == -α * (T - Ta) + β * u)
@constraint(off, ∂(t) == 1)
@constraint(on,  ∂(t) == 1)

@constraint(off, u == 0)                             # per-mode input set
@constraint(on, 0.2 <= u <= 1)

add_transition!(model, off => on) do tr; @constraint(tr, T <= 19) end
add_transition!(model, on => off) do tr; @constraint(tr, T >= 21) end

@constraint(off, final(T) in MOI.Interval(21.0, 23.0))    # ⇒ HybridSpec + TimedSpec
@constraint(on,  final(T) in MOI.Interval(21.0, 23.0))
@constraint(off, final(t) in MOI.Interval(15.0, 50.0))
@constraint(on,  final(t) in MOI.Interval(15.0, 50.0))

set_attribute(model, "state_grid", MP.GridFree(SVector(0.0), SVector(0.1)))   # → all modes
set_attribute(model, "input_grid", MP.GridFree(SVector(0.0), SVector(0.1)))
set_attribute(model, "time_step", 0.1)

optimize!(model)                                     # → HybridSystemAbstraction, auto-selected

@assert is_solved_and_feasible(model)
controller = get_attribute(model, "concrete_controller")
traj       = Dionysos.simulate(model, (SVector(21.5), 0.0, 1); nsteps = 200)
```

Every construct maps to something that already exists downstream: `HybridSystem`,
`VectorContinuousSystem([phys, clock])`, `GuardedResetMap`, `PR.HybridSpec`/`TimedSpec`,
`PR.OptimalControlProblem`, `HybridSystemAbstraction.Optimizer`. **R1 holds: the wrapper adds no
semantics.**

---

## 6. Phasing

Each phase is independently gated and committable (`[REF|ADD] wrapper: …`).

### The safety net must be built first

Revision 1 proposed gating on "`problems/*` and the docs examples". **`problems/` does not use the
wrapper at all** (L15). The real net is *one test file with three testsets*, plus three docs
examples that only execute during a docs build. That is not enough to gate a rewrite of the only
user-facing surface — hence Phase 0.

| # | Phase | Content | Gate |
| :- | :--- | :--- | :--- |
| **0** | **Characterisation** | `test/wrapper/current_behaviour.jl`: pin today's wrapper end-to-end — `∂`/`Δ`/`final`/`start`/`∉`, attribute forwarding, every error path, and the *bugs* (L4, L7, L9) explicitly marked as "current behaviour, changes in Phase N". | new tests pass against the **unmodified** wrapper |
| **1** | **IR + move to `src/`** | `ModelIR`; parse/lower split; `Dionysos.Optimizer` becomes a type in `src/wrapper/`; extension keeps only the symbolic backend; roles I1/I3 + the I4 unused-variable error. | Phase-0 suite green except the I4 case; Path planning / Unicycle / Simple pendulum unchanged |
| **2** | **Results surface** | `TerminationStatus`/`PrimalStatus`/`ResultCount`/`RawStatusString` (§4.9), `Dionysos.simulate` | `test/wrapper/results.jl`; the docs examples' trajectory blocks shrink to one call |
| **3** | **Specifications** | marker sets, problem-type inference, horizon conversion, L7 fix, `@objective` error | `test/wrapper/specifications.jl`; one safety + one reach-and-stay end-to-end |
| **4** | **Solver selection** | `"solver"` attribute, buffered attributes, `select_solver`, `supports_problem` | `test/wrapper/solver_selection.jl` |
| **5** | **Modes & transitions** | `Mode`/`Transition`, `@mode`, `add_transition!`, per-mode bounds/dynamics/attributes, arity-aware `ST.GuardedResetMap`, lowering to `HybridSystem` | `test/wrapper/hybrid.jl` reproduces the hand-built system in `hybrid_system_abstraction.jl` |
| **6** | **Clocks** | rule I2, `TimedSpec` from `final(clock)`, `HybridSpec` assembly | thermostat JuMP model matches `problems/Thermostat/*` |
| **7** | **Dynamics backends** | `NonlinearEvaluatorBackend`, `UserDynamicsBackend` (#510), automatic `jacobian_bound` | a **separate `julia -e` process** (see below) |
| **8** | **Docs & examples** | Literate `Thermostat.jl` (the #395 deliverable), a safety example; `docs/src/reference/Wrapper.md`; delete the `Interval(-100, 100)` workaround | docs build clean (it errors on undocumented exported symbols) |
| **9** | **LTL layer** *(optional)* | `Label` sets named via `MOI.ConstraintName`, `@specification`, formula normaliser → `CoSafeLTLProblem` | the ToyProblem co-safe suite driven from JuMP |

**Phase 7's gate needs its own process.** `test/runtests.jl:85-95` `include`s every file into a
**single** Julia session, and `jump_frontend.jl` does `using Symbolics, MathOptSymbolicAD` — once any
file loads the extension it is live for the rest of the run. "Green with no Symbolics" must be a
separate `julia --project=test -e '...'` invocation, wired into `runtests.jl` as a subprocess check.

**Dependencies:** 0 → 1 → {2, 3, 4, 7}; 4 → 5 → 6; 3 → 9; everything → 8.
The spine is **1, 4, 5**; 3 is what was explicitly asked for; 7 covers #510. Phases 2, 6, 8, 9 can
slip without blocking the rest.

**Design note for Phase 3.** Specifications written on a *mode* become `HybridSpec` entries (§4.6),
which only exist from Phase 5. Give the spec IR an **optional mode tag from the start**, so Phase 5
populates it rather than rewriting the lowering.

---

## 7. Decisions

| # | Question | Resolution |
| :- | :--- | :--- |
| D1 | Marker names: `Final`/`Always`/`EventuallyAlways` or `Reach`/`Safe`/`ReachAndStay`? | The former — `Final` already exists as the `final` operator's name and the others read as LTL. The `Reach…` family names *problems*, i.e. what they lower to, not what they say. |
| D2 | `@mode(model, off)` macro vs `add_mode!(model, :off)`? | Both; the macro registers the name in the object dictionary (JuMP-idiomatic). |
| D3 | Specification in `@objective` rather than `@constraint`? | **Settled with the user — no.** JuMP rejects a non-`AbstractScalarFunction` objective outright (verified). Specs stay constraints; `@objective` models only an additional cost. |
| D4 | Break the `final`-defaults-to-`±Inf` behaviour (L7)? | Yes — it is a bug and the workaround is visible in the flagship example. Phase 0 pins the old behaviour first so the change is deliberate. |
| D5 | Keep PR #395's `MOI.Integer`/`ZeroOne` support? | No — Bemporad–Morari-specific, out of scope; error clearly. |
| D6 | Keep the `∉` type piracy? | Yes, unchanged and already documented as deliberate. Add **no** new piracy. |
| D7 | Horizon in seconds or steps? | **Seconds for `∂` models** (wrapper converts via `PR.discretize_time`), **steps for `Δ` models**. Forced by the fact that no solver calls `discretize_problem`. |
| D8 | Build the LTL formula layer? | Yes, but **last and optional** (Phase 9). It is what makes the front-end general rather than a fixed menu, and the AP names come free from `MOI.ConstraintName` — but nothing else depends on it. An anonymous `Label` constraint must be a hard error. |
| D9 | Implement `PARAMETER`/`DISTURBANCE` roles? | **Deferred.** No issue or user report asks for them; `set_role!` covers the gap. Only rule I4 (unused-variable error) ships early. |

---

## 8. Explicit non-goals

* Not a general modelling language — no algebraic constraints, no DAEs, no free-final-time.
* Not a replacement for the direct-MOI entry style; both stay supported.
* **Not comprehensive solver coverage.** The dispatcher is extensible, but this plan wires
  **2 of ~10** families (§4.8); the rest need inputs the DSL cannot express.
* No new solver algorithms; the wrapper only *reaches* existing ones.
* No dependency creep — Spot/CDDLib/HiGHS and friends stay weak deps, and the wrapper imports none
  of them.

---

## 9. Verification log

Every JuMP mechanism claimed above was executed against the pinned toolchain
(**JuMP 1.30.1, MOI 1.51.1**, `--project=test`), not assumed.

| Claim | Result |
| :--- | :--- |
| A bare `struct <: JuMP.AbstractModel` supports `@constraint` with only `JuMP.add_constraint` defined — no `_valid_model`/`model_convert`/`value_type` piracy | **holds** |
| Per-mode dynamics, affine guards, and two-sided per-mode bounds all parse on such a scope; `JuMP.moi_function` resolves the *parent* model's variable refs | **holds** (table in §4.5) |
| Guard vs reset separable without head-pattern predicates | **holds** — they differ by MOI function *type* |
| `final(x)` on a whole `Vector{VariableRef}` builds a spec expression | **FAILS** — `MethodError: no method matching _final(::Vector{VariableRef})`. Drove the set-wrapper design of §4.4. |
| `@constraint(scope, x in Final(S))` with `Final <: MOI.AbstractVectorSet` builds `VectorOfVariables`-in-`Final{S}`, slices included | **holds** |
| A custom `AbstractModel` needs `Base.broadcastable(m) = Ref(m)` for vector constraints | **required** — `model_convert(::AbstractModel, ::VectorConstraint)` broadcasts over the model |
| `@objective` can carry a temporal formula | **FAILS** — *"The objective function `F(goal)` is not supported by JuMP"*. Settled D3. |
| `@constraint(model, goal, x in Label(S))` forwards the JuMP constraint **name** to the optimizer | **holds** — via `MOI.ConstraintName`; AP symbols come free from JuMP's naming |
| A model-level formula rides a raw attribute with the right cardinality | **holds** — single-valued, non-accumulating |
| `PR.discretize_time` is applied on the abstraction path | **FALSE** — `discretize_problem` is never called there; solvers pass `time` verbatim. Drove D7. |
| `SafetyProblem.time` is honoured by the safety solver | **FALSE** — never read; the solver computes a maximal invariant set. §10. |
| `test/runtests.jl` isolates test files enough to check "no Symbolics loaded" | **FALSE** — single process, `include` loop. Drove the Phase-7 gate change. |

The probes were throwaway scripts; the tables are the record. They are cheap to reconstruct — a
`Scope <: JuMP.AbstractModel` holding a `Vector{Any}`, an `add_constraint` that pushes
`(JuMP.moi_function(JuMP.jump_function(con)), JuMP.moi_set(con))`, and one `@constraint` per row.
Worth redoing if the `JuMP = "1"` compat bound is ever widened, since Phase 5 rests entirely on
these behaviours.

---

## 10. Library gaps this surfaces

The wrapper is a compiler (R1), so where the DSL feels awkward the cause is usually downstream.
These are recorded rather than worked around, and each is a candidate standalone fix:

| Gap | Consequence for the wrapper | Suggested fix |
| :--- | :--- | :--- |
| ~~`OptimalControlProblem` has no `safe_set` field~~ **fixed** | reach-avoid had to fold the `Always` set into `X`, so `Always` meant two different things depending on context (§4.4) | done: `safe_set` added as an optional last field, honoured by the reachability fixed points; see below |
| `discretize_problem` is never called by the abstraction solvers | the wrapper must convert seconds→steps itself (D7); the same trap awaits any direct-MOI user | call it in the abstraction pipeline, or rename `time` to make the unit explicit |
| `SafetyProblem.time` is ignored by the solver | `"horizon"` is inert for safety problems | implement finite-horizon safety, or type the field as `Infinity` |
| No reset-map type; 11 copy-pasted `…ResetMap` structs | Phase 5 must add `ST.GuardedResetMap` first | add it and sweep `problems/`, `research/`, `test/` |
| Reset/guard arity differs between clock-lifted and time-free modes; mixed hybrids unsupported | the wrapper must extend `x`-guards across the clock domain and reject mixed models | unify the reset-map calling convention in `HybridSystemAbstraction` |
| No `MOI.TerminationStatus` on any abstraction solver | Phase 2 synthesises it from `success::Bool` at the wrapper level | push the mapping down into the solvers so direct-MOI users get it too |

---

## 11. Execution log

Phases 0–8 are done on branch `jc_issue_4`, each formatted, gated on the `--fast` suite (59 files)
plus the golden regression and the slow hybrid suites, and committed separately.

| Commit | Phase | What landed |
| :--- | :--- | :--- |
| `b9331af1` | — | `FIX optim: no success when the initial set is empty` (see below) |
| `f961e894` | 0 + 1 + 2 | `ModelIR`, parse/lower split, wrapper moved to `src/wrapper/`, extension reduced to the symbolic backend, rule I4, solution status, `simulate`, `test/wrapper/lowering.jl` |
| `816610e1` | — | this plan and `src/wrapper/README.md` |
| `e4271f81` | 3 | spec markers, problem-type inference, horizon, L7 fix, `@objective` error, `test/wrapper/specifications.jl` |
| `c320ddc2` | 4 | `"solver"` attribute, attribute replay, `select_solver`/`supports_problem`, `test/wrapper/solver_selection.jl` |
| `18475055` | 5 | `ST.GuardedResetMap` + `test/system/reset_map.jl` |
| `b7b4862a` | 5 | `Mode`/`Transition` as scopes, `@mode`, `add_transition!`, per-mode bounds and options, hybrid lowering, `test/wrapper/hybrid.jl` |
| `dc8b9a38` | 6 | clocks (rule I2), `(x, t, mode)` states, time-windowed specs, `test/wrapper/clock.jl` |
| `a2555969` | 7 | supplied dynamics (#510), `set_role!`, `NonlinearEvaluatorBackend`, `test/wrapper/{dynamics,no_symbolics}.jl` |
| `3a5e38c8` | 8 | options moved off the IR (see below), `lower` entry point |
| `2135bc76` | 8 | `docs/src/examples/solvers/Thermostat.jl`, `docs/src/reference/Wrapper.md`, Path-planning workaround deleted |

Phase 2 is folded into `f961e894`: it touches the same files and was developed and gated with them.

### Five defects the execution surfaced

None of these were visible from reading the code; each was found by a test.

1. **L7 is a crash, not a silent bug.** An unconstrained `final` coordinate produced the interval
   `(-Inf, Inf)`, from which `UT.box` built a NaN centre and radius and threw
   `AssertionError: radius must be nonnegative` from inside LazySets — with no mention of the
   model. So the Path-planning example's `final(x[3]) in MOI.Interval(-100.0, 100.0)` was load
   bearing, and any multi-state model had to constrain *every* coordinate. Fixed in Phase 3.

2. **Success was reported vacuously.** `optimizer.success = all(q -> q in controllable_set,
   problem.initial_set)` is `true` for an **empty** initial set, so a model whose initial region is
   not represented in the abstraction reported success without verifying anything — reaching the
   user as a spurious `MOI.OPTIMAL` once Phase 2 exposed the status. Fixed with
   `DiscreteSystems.covers_initial_set`, applied to the optimal-control, safety and reach-and-stay
   solvers.

3. **Attaching an abstraction-only problem discarded the configuration.**
   `set_concrete_problem!(::UniformGridAbstraction.Optimizer, ::AlternatingSimulationProblem)`
   replaced the abstraction sub-solver with a fresh one, carrying over only four fields — so
   `state_grid`, `input_grid`, `time_step` and `approx_mode` set *before* the problem were
   silently lost. It now keeps the configured solver and clears just its cached results. This is a
   trap for direct-MOI users too, not only the wrapper.

4. **Every attribute-carried option was dropped in automatic mode.** `MOI.empty!` clears the
   model, and JuMP calls it on the optimizer before copying the model in — so the horizon,
   supplied dynamics, declared roles and per-mode grids, all of which the wrapper stored in the
   `ModelIR`, were silently wiped on every `optimize!`. The identical model built with
   `direct_model` worked, which is why every wrapper test missed it: they all used
   `direct_model` for speed, and the bug lives only in `Model(Dionysos.Optimizer)` — the
   documented form. **It was the thermostat example that found it**, on its first run. Options
   now live on the `Optimizer` and are folded into the IR by `_apply_options!`. The fix also
   produced `lower`, the first half of `optimize!` exposed on its own, which the test
   helpers now share so they cannot drift from what really runs.

5. **`Δ` on a transition is not a time-domain declaration.** A reset map is an instantaneous
   jump, so `Δ(x) == …` written on a transition says nothing about whether the *modes* evolve in
   continuous or discrete time. Setting the domain from it made `∂` in a mode collide with `Δ`
   in a transition — in code whose own docstring already said the right thing.

### Testing lessons

- `direct_model` is faster and more precise for inspecting the lowering, but it bypasses the
  caching layer where defect 4 lived. `test/wrapper/lowering.jl` now carries one automatic-mode
  test covering horizon, declared roles and supplied dynamics together.
- Expression shapes must be checked, not assumed: `∂(t) == 1` arrives as the expression `+(1.0)`,
  not the number `1`, because JuMP wraps a literal right-hand side.

### One decision taken against this plan: no formula normaliser

§4.4 specified a **normaliser** that would match a temporal formula against the patterns with
specialised solvers — `F(s)` → `OptimalControlProblem`, `G(s)` → `SafetyProblem`, `F(G(s))` →
`ReachAndStayProblem` — and reroute it, falling back to `CoSafeLTLProblem`. **It was not built.**

A formula attached with `@specification` always lowers to `CoSafeLTLProblem`; the markers of
§4.4's sugar layer remain the only way to reach the specialised fixed points. The reason is that
recognising formula shapes reliably is a rewriting problem — `F(a & true)`, `!G(!a)` and `F(a)`
are the same specification — and a normaliser that gets it *nearly* right silently changes which
algorithm runs, which is exactly the class of bug this plan spent its effort eliminating. The two
layers are therefore presented as a genuine choice in `src/wrapper/README.md` §5.3 rather than as
a fast path and a slow path over one representation.

Rerouting can be added later without touching the DSL: both layers already land in
`build_problem`, so a normaliser would be one function between them.

### Follow-up: `safe_set` on `OptimalControlProblem` (§10, gap 1 — done)

`OptimalControlProblem` gained an optional `safe_set`, closing the asymmetry where `Always`
became a genuine safe set in a pure safety model but was folded into `X` — behaving like `∉` —
as soon as a `Final` set was written next to it.

Shape of the change:

* **Field last, defaulted.** `safe_set` is the seventh field, `nothing` by default, so the 70
  existing positional call sites are untouched and `remake` reaches it for free. Inserting it
  next to `target_set` — where it reads better — would have rewritten all of them.
* **Semantics: `safe U target`.** The target is intersected with the safe set, so a target state
  outside it does not count as reached. This is the sound reading of `◻ safe ∧ ◇ target`, and it
  is what folding into `X` used to do implicitly, since unsafe cells simply did not exist.
* **Enforced in the fixed points, not the lifting.** Both reachability fixed points
  (`compute_worst_case_uniform_cost_controller` and the Dijkstra `compute_optimal_controller`)
  refuse to admit an unsafe state into the controllable set. Worst-case avoidance then follows
  for free from the existing counter: an input that *may* land in an unsafe state never sees its
  counter reach zero, so it is never selected.
* **Solvers that cannot honour it must say so.** `PR.check_safe_set_supported` raises rather than
  let `BemporadMorari` or `BranchAndBound` return a controller certified against a weaker
  specification than the one asked for. Silent weakening is the failure mode this whole plan has
  been chasing (§11, defects 2 and 3).

It also fixed a live bug in the hybrid path: `build_hybrid_problem` **dropped** the `Always` set
whenever a `Final` set was present, because there was nowhere to put it.

### Remaining

Nothing in this plan. Open follow-ups are the remaining library gaps of §10, the
`PARAMETER`/`DISTURBANCE` roles deferred by D9, and mapping `@objective` onto `transition_cost`
(D3) — whose error message currently names an attribute that does not exist.
