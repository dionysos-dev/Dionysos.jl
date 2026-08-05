# `Dionysos.Wrapper` — the JuMP front-end

> **Status: design document.** The code described here does not exist yet; this folder currently
> holds only this README. The implementation plan (phases, gates, decisions) is [`plan.md`](../../plan.md)
> at the repo root. Today's front-end lives in [`src/MOI_wrapper.jl`](../MOI_wrapper.jl) and is
> loaded from [`ext/DionysosMathOptSymbolicADExt.jl`](../../ext/DionysosMathOptSymbolicADExt.jl).
> Sections marked **(today)** already work; everything else is the target.

This is the surface where a user who is not a Dionysos developer writes a control problem. You
describe a **system**, state a **specification**, optionally pick a **solver**, and get back a
correct-by-construction **controller**.

---

## 1. Thirty seconds

```julia
using Dionysos, JuMP, StaticArrays
using Symbolics, MathOptSymbolicAD          # enables the fast dynamics backend

model = Model(Dionysos.Optimizer)

@variable(model, -1.0 <= x[i = 1:2] <= 1.0, start = -0.75)   # states + where we start
@variable(model, -1.0 <= u[1:2] <= 1.0)                      # inputs

@constraint(model, ∂(x[1]) == u[1])                          # ẋ = u
@constraint(model, ∂(x[2]) == u[2])

@constraint(model, x in Final(UT.box([-0.5, -0.5], [0.5, 0.5])))   # ◇ reach this box

set_attribute(model, "time_step", 1.0)
set_attribute(model, "state_grid", MP.GridFree(SVector(0.0, 0.0), SVector(0.25, 0.25)))
set_attribute(model, "input_grid", MP.GridFree(SVector(0.0, 0.0), SVector(0.5, 0.5)))

optimize!(model)

@assert is_solved_and_feasible(model)
controller = get_attribute(model, "concrete_controller")
traj       = Dionysos.simulate(model, SVector(-0.75, -0.75); nsteps = 100)
```

That is the whole contract. Everything below is detail.

---

## 2. The mental model

A Dionysos task is a triple:

```
     system  𝒮          +      specification  Σ       →     solver  𝒪      →   controller
  ∂/Δ + bounds + modes      Start/Final/Always/LTL       auto or explicit     + certificate
```

The wrapper is a **compiler**, not a solver. It lowers your JuMP model into a
`MathematicalSystems`/`HybridSystems` object plus a `Dionysos.Problem.ProblemType`, then hands both
to an existing optimizer. It never invents semantics: if something cannot be expressed as one of
those two objects, it does not belong here.

Three JuMP mechanisms, three distinct meanings — never overloaded:

| JuMP mechanism | Carries | Cardinality |
| :--- | :--- | :--- |
| `@variable` + `@constraint` | the system and the specification *sets* | many, they accumulate |
| `@specification` / `"specification"` attribute | the temporal **formula** combining them | exactly one |
| `@objective` | the **cost** to minimise (not the specification) | exactly one |

---

## 3. Variables — you declare, we infer the role

You never say "this is a state" and "this is an input". The role is **inferred from how the
variable is used**, and this is deliberately the wrapper's central idea.

| Rule | If the variable… | …its role is |
| :-: | :--- | :--- |
| I1 | appears on the left of `∂(·) == …` or `Δ(·) == …` | `STATE` |
| I2 | has dynamics exactly `∂(t) == 1` or `∂(t) == 0`, and appears in no other right-hand side | `CLOCK` |
| I3 | appears on some right-hand side, none of the above | `INPUT` |
| I4 | appears nowhere | **error** |
| I5 | is JuMP-fixed (`@variable(model, p == 3.0)`) and appears only on right-hand sides | `PARAMETER` — *planned, not in the first release* |

Rule I4 matters most: **(today)** a leftover or mistyped variable silently becomes an input,
enlarging the input grid — and if it is unbounded, the abstraction explodes. It becomes a hard error
naming the variable.

There is **no mode "role"** — a mode is a scope you write constraints on (§6), not a variable.

Override when inference cannot see the truth (e.g. dynamics supplied as a raw Julia function):

```julia
set_role!(w, Dionysos.DISTURBANCE)
```

**Ordering.** The state vector `x[i]` follows *declaration order*. Check it with
`get_attribute(model, "state_variables")` rather than guessing.

---

## 4. Dynamics

Continuous time uses `∂`, discrete time uses `Δ`. Mixing them in one model is an error. **(today)**

```julia
@constraint(model, ∂(x[1]) == u[1] * cos(x[3]))     # ẋ₁ = u₁cos(x₃)
@constraint(model, Δ(x[1]) == x[1] + u[1])          # x₁⁺ = x₁ + u₁
```

Bounds come from the variable declaration and become the state set `X` and input set `U`:

```julia
@variable(model, 0.0 <= x[1:3] <= 4.0)    # X
@variable(model, -1 <= u[1:2] <= 1)       # U
```

Obstacles are carved out of `X`: **(today)**

```julia
@constraint(model, x[1:2] ∉ MOI.HyperRectangle([1.0, 0.0], [1.2, 9.0]))
```

You can also hand over a Julia function instead of writing the equations
(issue [#510](https://github.com/dionysos-dev/Dionysos.jl/issues/510)):

```julia
set_attribute(model, "dynamics", (x, u) -> [x[2], -sin(x[1]) + u[1]])
```

---

## 5. Specifications — two layers

This is where the wrapper differs most from a plain optimisation model, so it is worth being
explicit about the design.

### 5.1 The general layer: a temporal formula over named sets

> **Planned last** (Phase 9). If you only need reach, safety, reach-and-stay, or reach-avoid, §5.2
> covers you and ships first.

The truly general way to state a control specification is a **temporal-logic formula over atomic
propositions**, where each proposition names a region of the state space. The proposition names come
from **JuMP's own constraint naming** — no separate registration call:

```julia
@constraint(model, goal,   x in Label(UT.box([3.0, 0.3], [3.6, 0.8])))
@constraint(model, hazard, x in Label(obstacle_region))

@specification(model, "F(goal) & G(!hazard)")     # ◇goal ∧ □¬hazard
```

The constraint name *is* the atomic proposition. `@constraint(model, goal, …)` already registers
`goal` in JuMP's object dictionary and forwards the name to the optimizer through
`MOI.ConstraintName` — verified working, so this costs no new syntax at all.

With Spot loaded you can pass a formula object instead of a string:

```julia
using Spot
@specification(model, ltl"F(goal) & G(!hazard)")
```

### 5.2 The sugar layer: named patterns for the common cases

Most control problems are one of four shapes, and those shapes have **dedicated fixed-point
algorithms that are far faster than going through an automaton product**. So the common patterns get
their own spelling, which builds the *same* internal formula:

| Sugar | Formula | Meaning |
| :--- | :--- | :--- |
| `@constraint(model, x in Start(S))` | — | initial set |
| `@variable(model, x, start = v)` **(today)** | — | singleton initial set |
| `@constraint(model, x in Final(S))` | `F(s)` | ◇ eventually reach `S` |
| `@constraint(model, x in Always(S))` | `G(s)` | □ always stay in `S` |
| `@constraint(model, x in EventuallyAlways(S))` | `F(G(s))` | ◇□ reach `S` and stay |
| `@constraint(model, x ∉ O)` **(today)** | `G(!o)` | avoid `O` (carved out of `X`) |
| `@constraint(model, final(x[i]) in MOI.Interval(a, b))` **(today)** | `F(s)` | per-coordinate reach |
| `@constraint(model, start(x[i]) in MOI.Interval(a, b))` **(today)** | — | per-coordinate initial |

`S` is **any bounded `LazySet`** — a box, an ellipsoid, a ball, a polytope, a zonotope. The
discretisation layer (`MP.get_states_from_set`) has always accepted these; only the front-end
restricted you to boxes.

### 5.3 How the two layers meet

Both layers produce one formula. A **normaliser** then matches it against the patterns the library
has specialised solvers for, and routes to the fastest one:

| Normalised formula | Problem built | Why specialised |
| :--- | :--- | :--- |
| `F(s)` | `OptimalControlProblem` | backward fixed point on the pre-image |
| `F(s) & G(¬o)` | `OptimalControlProblem` (reach-avoid) | same, with `¬o` folded into `X` |
| `G(s)` | `SafetyProblem` | maximal controlled-invariant set |
| `F(G(s))` | `ReachAndStayProblem` | invariance then reachability |
| any other **co-safe** formula | `CoSafeLTLProblem` | product with a deterministic automaton |
| *(no formula at all)* | `AlternatingSimulationProblem` | just build the abstraction and return it |
| anything else | **error**, naming the unsupported operator | — |

That last row is the honest part, and it is why the general layer is *LTL* but the general *solver*
is co-safe only. **Not every specification is co-safe.** Co-safe means "satisfied by some finite
prefix" — good for `F`, `U`, `X`. But `G(s)` is never satisfied by a finite prefix (it is the dual,
a *safety* property), and `F(G(s))` is neither. The wrapper therefore recognises `G` and `F(G)`
as patterns and sends them to their dedicated solvers rather than to the automaton machinery, which
genuinely cannot take them.

This mirrors the `Problem` module exactly: `CoSafeLTLProblem` is the general type, and
`OptimalControlProblem`/`SafetyProblem`/`ReachAndStayProblem` exist because the specialised
algorithms are worth having.

### 5.4 Why not `@objective`?

It reads like the right home for a goal, but JuMP forbids it: `@objective` requires an
`MOI.AbstractScalarFunction`, and a temporal formula is not one —
`@objective(model, Min, "F(goal)")` fails with *"The objective function `F(goal)` is not supported
by JuMP"* (verified). The instinct behind the question is still correct, though: a specification is
**one, model-level, non-accumulating** declaration, exactly like an objective and unlike a
constraint. That cardinality is why the *formula* is a single `@specification`, while the *sets* it
names — which really do accumulate — are constraints.

`@objective` is reserved for what it actually means: the **cost**.

```julia
@objective(model, Min, sum(u[i]^2 for i in 1:2))    # → transition_cost
```

**(today)** `@objective` is accepted and *silently discarded*. It will raise a clear error until the
cost mapping lands.

### 5.5 Avoiding a region: `∉` and `Always` are *not* the same

`x ∉ O` and `x in Always(complement(O))` say the same thing logically but mean different things
here, and the difference is worth understanding because it affects both speed and what the
controller can promise:

* **`x ∉ O` restricts the system.** The region is carved out of the state space `X`, so those cells
  are never created. Fewer states, faster abstraction. Use it for walls, obstacles, and anywhere the
  model simply is not valid.
* **`x in Always(S)` constrains the specification.** States outside `S` stay representable, so
  synthesis can reason about leaving them and actively prevent it. Use it when "stay inside `S`" is
  the *goal*.

Rule of thumb: `∉` describes the **world**, `Always` describes the **requirement**.

One asymmetry to know about: `Dionysos.Problem.OptimalControlProblem` has no safe-set field, so in a
*reach-avoid* model (`Final` together with `Always`) the `Always` set is folded into `X` — it behaves
like `∉`. In a pure safety model (`Always` alone) it becomes a genuine `safe_set`.

### 5.6 Horizon

```julia
set_attribute(model, "horizon", 10.0)     # default: PR.Infinity()
```

* **Continuous-time model (`∂`)** — the horizon is in **seconds**, and is converted to a step count
  using your `time_step`, rounded in the direction the specification requires ("within at most `T`"
  rounds down, "for at least `T`" rounds up).
* **Discrete-time model (`Δ`)** — the horizon is a **step count**, used as given.

> **Caveat:** a horizon on a *safety* problem currently has no effect. The safety solver computes a
> maximal controlled-invariant set, which is an infinite-horizon notion, and never reads the field.

---

## 6. Hybrid systems — modes and transitions

A **mode** is a scope you attach dynamics, bounds, and specifications to. A **transition** is a
scope you attach a guard and a reset map to. Both are ordinary JuMP models, so `@constraint` works
on them unchanged.

```julia
model = Model(Dionysos.Optimizer)
@variable(model, 17 <= T <= 25)
@variable(model, 0 <= u <= 1)

off = @mode(model, off)
on  = @mode(model, on)

@constraint(off, ∂(T) == -α * (T - Ta))              # per-mode dynamics
@constraint(on,  ∂(T) == -α * (T - Ta) + β * u)

@constraint(off, u == 0)                             # per-mode input set
@constraint(on, 0.2 <= u <= 1)

add_transition!(model, off => on) do t
    @constraint(t, T <= 19)                          # guard
end

add_transition!(model, on => off) do t
    @constraint(t, T >= 21)                          # guard
    @constraint(t, Δ(T) == T)                        # reset map (identity ⇒ omittable)
end
```

**Guard or reset?** Decided by the *shape* of what you wrote, with no ambiguity: a constraint
containing `∂`/`Δ` is a **reset map**; anything else is a **guard**. The same rule separates
per-mode dynamics from per-mode bounds.

**Several transitions between the same two modes** are fine — a transition is an object, not a
`"src->dst"` key, so their guards and resets never get mixed up.

Switching type, when you need it:

```julia
add_transition!(model, off => on; switching = HybridSystems.ControlledSwitching()) do t … end
```

Solver options are set **per mode**, which is what the hybrid abstraction wants:

```julia
for m in (off, on)
    set_attribute(m, "state_grid", MP.GridFree(SVector(0.0), SVector(0.1)))
    set_attribute(m, "input_grid", MP.GridFree(SVector(0.0), SVector(0.1)))
    set_attribute(m, "time_step", 0.1)
end
```

### Clocks and timed specifications

A clock is just a variable whose dynamics is `∂(t) == 1` (running) or `∂(t) == 0` (frozen) —
rule I2. No new syntax:

```julia
@variable(model, 0 <= t <= 50)
@constraint(off, ∂(t) == 1)          # clock runs while off
@constraint(on,  ∂(t) == 0)          # clock frozen while on
```

A time window on a target is then simply a reach constraint on the clock:

```julia
@constraint(on, final(T) in MOI.Interval(21.0, 23.0))
@constraint(on, final(t) in MOI.Interval(15.0, 50.0))   # … within t ∈ [15, 50]
```

Specifications written on a mode become per-mode specifications, and with a clock they become
time-windowed ones — reproducing `PR.hybrid_reach_spec(Xs, Ts, Ns)` from ordinary constraints.

---

## 7. Choosing a solver

By default the wrapper picks one from the shape of your problem: a model with modes goes to the
hybrid abstraction, everything else to the uniform grid abstraction. Those two are the families
reachable from JuMP — see the last entry in §11 for why the others are not. Override explicitly:

```julia
set_attribute(model, "solver", AB.UniformGridAbstraction.Optimizer)
```

Every other `set_attribute` is forwarded to the chosen solver and its sub-solvers, so all the
options documented on those optimizers (`approx_mode`, `jacobian_bound`, `print_level`,
`execution_backend`, …) are reachable from JuMP. An attribute that no solver recognises raises
`MOI.UnsupportedAttribute` naming it, rather than being silently dropped.

An unsupported *combination* (a solver family that cannot handle your problem type) is reported
here, naming both — not deep inside a sub-solver.

### Growth bounds come for free

**(today)** the `GROWTH` approximation mode requires you to hand-write a `jacobian_bound`. With the
symbolic backend loaded this is derived automatically by tracing the dynamics and bounding the
Jacobian over `X` with interval arithmetic. Supply one only to override.

---

## 8. Reading the results

### Did it work?

```julia
optimize!(model)

termination_status(model)          # OPTIMAL | LOCALLY_INFEASIBLE | OPTIMIZE_NOT_CALLED
is_solved_and_feasible(model)      # true ⇒ you have a certified controller
raw_status(model)                  # human-readable reason
```

`LOCALLY_INFEASIBLE` — not `INFEASIBLE` — is what you get on failure, and the distinction is the
whole point of abstraction-based control: the abstraction is **sound but not complete**, so "no
controller was found" means *no controller exists on this abstraction*. A finer grid, a smaller time
step, or a different `approx_mode` may well succeed. Nothing here proves your problem is unsolvable.

### The controller and the certificate

```julia
controller       = get_attribute(model, "concrete_controller")
abstract_system  = get_attribute(model, "abstract_system")
concrete_problem = get_attribute(model, "concrete_problem")
value_function   = get_attribute(model, "abstract_value_function")   # the certificate

MOI.get(model, MOI.RawOptimizerAttribute("abstraction_construction_time_sec"))
MOI.get(model, MOI.SolveTimeSec())
```

### Running it

```julia
traj = Dionysos.simulate(model, x0; nsteps = 100)
```

One call for every kind of model — it finds the discrete-time (or hybrid) system and the controller
on the solved model and picks the right closed-loop engine. For a hybrid model `x0` is the augmented
state, `(x, mode)` or `(x, t, mode)`. The result is a channelled `Dionysos.System.Trajectory`; read
it with `ST.states`, `ST.inputs`, `ST.times`, `ST.modes`, plot it with `plot(traj)`, or animate it
with `Dionysos.animate_trajectory_dashboard`.

Stopping is derived from your specification by default — pass `stopping = x -> …` to override.

---

## 9. Folder structure

```
src/wrapper/
  wrapper.jl              module + includes; loaded from src/Dionysos.jl
  variables.jl            VariableInfo, roles, the I1–I6 inference
  model_ir.jl             ModelIR — the parsed model, printable, dependency-free
  modes.jl                Mode       <: JuMP.AbstractModel
  transitions.jl          Transition <: JuMP.AbstractModel, add_transition!
  operators.jl            ∂ Δ final start; Start/Final/Always/EventuallyAlways/Label sets; ∉
  specification.jl        spec IR + problem-type inference (+ the §5.3 normaliser)
  parse.jl                MOI.add_constraint / MOI.set → ModelIR      ← the only pattern matching
  lower_system.jl         ModelIR → MathematicalSystems / HybridSystems
  lower_problem.jl        ModelIR → Problem.ProblemType
  solver_selection.jl     select_solver, supports_problem
  dynamics_backend.jl     AbstractDynamicsBackend + evaluator / user-function backends
  results.jl              TerminationStatus / PrimalStatus / simulate
  optimizer.jl            Dionysos.Optimizer <: MOI.AbstractOptimizer

ext/DionysosMathOptSymbolicADExt.jl    SymbolicADBackend only (codegen + automatic Jacobian bound)
```

The pipeline is strictly one-directional:

```
JuMP model ──parse──▶ ModelIR ──lower──▶ (system, problem) ──select──▶ solver ──▶ controller
             (pure)   (data)   (pure)                       (dispatch)
```

`ModelIR` is the reason the wrapper lives in `src/` and not in an extension: parsing needs no
Symbolics, only *compiling a dynamics expression into a callable* does. That is what makes the
front-end loadable, testable, and documentable without any optional dependency.

---

## 10. Extending it

| To add… | Touch |
| :--- | :--- |
| a specification kind | one marker set in `operators.jl` + one normaliser row in `specification.jl` |
| a problem type | one `lower_problem.jl` method + one `supports_problem` declaration |
| a solver family | one `select_solver` method + one `supports_problem` declaration |
| a dynamics source | one `AbstractDynamicsBackend` subtype + `compile_dynamics` |
| a system class | one `lower_system.jl` method |

Nothing in that table requires touching `parse.jl`. That is the point of the IR: **(today)** every
new feature grows an `if/elseif` chain inside `MOI.add_constraint` that ends in
`dump(func); error("Unsupported")`.

Two rules for contributors:

* **The wrapper owns no semantics.** If a concept has no `MathematicalSystems`/`HybridSystems` +
  `ProblemType` representation, add it to `System`/`Problem` first, then expose it here.
* **No new type piracy.** The `∉` parsing rule is a documented, deliberate exception because JuMP
  offers no extension point for it. `Mode`/`Transition` subtype `JuMP.AbstractModel`, which gives
  `@constraint` support with no piracy at all.

---

## 11. Gotchas

* **Mixing `∂` and `Δ`** in one model is an error — pick continuous or discrete time.
* **Variable order is declaration order.** Interleaving state and input declarations reorders the
  vectors. Check with `get_attribute(model, "state_variables")`.
* **`final(x)` on a whole vector does not work.** `final` is a `JuMP.NonlinearOperator`, which only
  builds an expression when one of its arguments is a JuMP scalar; on a `Vector{VariableRef}` it
  raises `MethodError: no method matching _final(::Vector{VariableRef})`. Use the vector form
  `x in Final(S)`, or the scalar form `final(x[i]) in MOI.Interval(a, b)`.
* **`start = v` gives a singleton initial set**, not a region. Use `x in Start(S)` for a set.
* **(today) unconstrained `final` coordinates default to `±Inf`**, silently producing an unbounded
  target — which is why the Path planning example carries a
  `final(x[3]) in MOI.Interval(-100.0, 100.0)` line. They will default to the variable's own
  declared bounds instead. This is the one intentional behaviour change in the plan.
* **Nonlinear guards are rejected.** Guards must be affine so they can become a polyhedron; the
  hybrid abstraction cannot consume anything else either.
* **A horizon on a safety problem does nothing** — the safety solver is infinite-horizon by
  construction (§5.6).
* **Mixing clock-lifted and time-free modes in one hybrid model is rejected.** Guards and reset maps
  are applied to the augmented state, whose arity differs between the two (`[x; t]` vs `x`), so a
  mixed model has no consistent reset convention. Give every mode a clock, or none.
* **`Always` and `∉` are not interchangeable** — see §5.5.
* **Integer / binary variables are not supported.** They belong to the MIQP solvers, which the
  wrapper does not target.
* **Only two solver families are reachable from JuMP** — the uniform grid abstraction and the hybrid
  abstraction. The others (Bemporad–Morari, branch-and-bound, the ellipsoidal families, the PCLF
  bisimulation quotient, trajectory generators and certifiers) need inputs this DSL cannot express,
  such as PWA systems, ellipsoid templates, or observation regions. Use the direct-MOI entry style
  for those.

---

## 12. See also

* [`plan.md`](../../plan.md) — the implementation plan: audit, phases, gates, open decisions.
* [`docs/src/examples/solvers/`](../../docs/src/examples/solvers/) — runnable examples.
* [`docs/src/manual/abstraction-based-control.md`](../../docs/src/manual/abstraction-based-control.md) — what the solvers actually do.
* [`docs/src/developers/conventions.md`](../../docs/src/developers/conventions.md) — house style.
