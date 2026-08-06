# `Dionysos.Wrapper` — the JuMP front-end

This is where you write a control problem. You describe a **system**, state a **specification**,
optionally pick a **solver**, and get back a correct-by-construction **controller**.

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

Your JuMP model is compiled into a `MathematicalSystems`/`HybridSystems` object plus a
`Dionysos.Problem.ProblemType`, which an existing solver then runs.

Three JuMP mechanisms, three distinct meanings:

| JuMP mechanism | Carries | How many |
| :--- | :--- | :--- |
| `@variable` + `@constraint` | the system and the specification *sets* | many, they accumulate |
| `@specification` | the temporal **formula** combining them | exactly one |
| `@objective` | the **cost** to minimise, not the specification | exactly one |

---

## 3. Variables — you declare, the role is inferred

You never say "this is a state" and "this is an input":

| Rule | If the variable… | …its role is |
| :-: | :--- | :--- |
| I1 | appears on the left of `∂(·) == …` or `Δ(·) == …` | `STATE` |
| I2 | has dynamics exactly `∂(t) == 1` or `∂(t) == 0`, and appears in no other right-hand side | `CLOCK` |
| I3 | appears on some right-hand side, none of the above | `INPUT` |
| I4 | appears nowhere | **error**, naming the variable |

Rule I4 catches a leftover or mistyped variable, which would otherwise be discretized as an input
and enlarge the abstraction for nothing.

There is **no mode "role"** — a mode is a scope you write constraints on (§6), not a variable.
`PARAMETER` and `DISTURBANCE` roles are not supported yet.

Declare a role yourself when there is nothing to infer from — which is the case when the dynamics
are a Julia function rather than equations:

```julia
set_role!(x, Dionysos.STATE)      # a single variable or a whole array
```

Takes `Dionysos.STATE`, `Dionysos.INPUT` or `Dionysos.CLOCK`. Anything left undeclared is an
input.

**Ordering.** The state vector `x[i]` follows *declaration order*. Check it with
`get_attribute(model, "state_variables")` rather than guessing.

---

## 4. Dynamics

Continuous time uses `∂`, discrete time uses `Δ`. Mixing them in one model is an error.

```julia
@constraint(model, ∂(x[1]) == u[1] * cos(x[3]))     # ẋ₁ = u₁cos(x₃)
@constraint(model, Δ(x[1]) == x[1] + u[1])          # x₁⁺ = x₁ + u₁
```

Bounds come from the variable declaration and become the state set `X` and input set `U`:

```julia
@variable(model, 0.0 <= x[1:3] <= 4.0)    # X
@variable(model, -1 <= u[1:2] <= 1)       # U
```

Obstacles are carved out of `X`:

```julia
@constraint(model, x[1:2] ∉ MOI.HyperRectangle([1.0, 0.0], [1.2, 9.0]))
@constraint(model, x[1:2] ∉ LazySets.Ball2([0.0, 0.0], 0.5))
```

An `MOI.HyperRectangle` may be written over a subset of the coordinates and spans the variable
bounds on the rest — `x[1:1] ∉ MOI.HyperRectangle([1.0], [1.2])` is a full-height wall. Any other
bounded `LazySet` is taken as written and must span the whole state vector.

The left-hand side must be a **vector**, even for a single coordinate: `x[1] ∉ …` is a scalar
function against a vector set, which JuMP rejects with *"we don't recognize x[1] as a valid JuMP
function"*. Write `x[1:1]`, or `[x1]` for a scalar variable.

### Dynamics as a Julia function

For a simulator, a lookup table, or anything else that is not a JuMP expression:

```julia
set_role!(x, Dionysos.STATE)          # nothing to infer from, so name the states
set_attribute(model, "dynamics", (x, u) -> [x[2], -sin(x[1]) + u[1]])
```

Bounds, specifications and solver options are unaffected — only where the equations come from
changes. A mode can carry its own: `set_attribute(off, "dynamics", f_off)`.

### Choosing how equations are compiled

```julia
set_attribute(model, "dynamics_backend", Dionysos.Wrapper.NonlinearEvaluatorBackend())
```

`SymbolicADBackend` (the default) traces the equations with Symbolics into one fused function, and
is what you want for real problems. `NonlinearEvaluatorBackend` evaluates them through
`MOI.Nonlinear` with no optional dependency — much slower, useful as a fallback.

---

## 5. Specifications

There are two ways to say what the controller must achieve: named patterns for the four common
shapes, and a temporal formula for anything else.

### 5.1 Named patterns

| Write | Meaning |
| :--- | :--- |
| `@constraint(model, x in Start(S))` | start in `S` |
| `@variable(model, x, start = v)` | start at the single point `v` |
| `@constraint(model, x in Final(S))` | ◇ eventually reach `S` |
| `@constraint(model, x in Always(S))` | □ never leave `S` |
| `@constraint(model, x in EventuallyAlways(S))` | ◇□ reach `S` and stay |
| `@constraint(model, x ∉ O)` | avoid `O` (carved out of `X`) — a box or any bounded `LazySet` |
| `@constraint(model, final(x[i]) in MOI.Interval(a, b))` | reach, one coordinate at a time |
| `@constraint(model, start(x[i]) in MOI.Interval(a, b))` | start, one coordinate at a time |

`S` is **any bounded `LazySet`** — a box, an ellipsoid, a ball, a polytope, a zonotope.

A coordinate with no `final` constraint falls back to that variable's own bounds, so leaving one
unconstrained means "any value counts as reaching the target".

### 5.2 A temporal formula

Name regions of the state space, then write a formula over them. The atomic proposition is the
**constraint's name**:

```julia
@constraint(model, goal,   x in Label(UT.box([3.0, 0.3], [3.6, 0.8])))
@constraint(model, hazard, x in Label(obstacle_region; semantics = MP.OUTER))

@specification(model, "F(goal) & G(!hazard)")     # ◇goal ∧ □¬hazard
```

An anonymous `Label` is an error — the formula would have nothing to refer to.

`semantics` says how a region is discretized: `INNER` (the default) keeps only cells lying fully
inside it, which is the conservative reading for something you must reach; `OUTER` keeps every cell
that touches it, the conservative reading for something you must avoid.

With Spot loaded you can pass a formula object instead of a string:

```julia
using Spot
@specification(model, ltl"F(goal) & G(!hazard)")
```

A hand-written monitor — any `Optim.DiscreteSystems.AbstractSpecStepper` — is accepted too.

### 5.3 Which one to use

The two forms reach different solvers:

| What you wrote | Problem built | How it is solved |
| :--- | :--- | :--- |
| `Final(S)` | `OptimalControlProblem` | backward fixed point on the pre-image |
| `Final(S)` + `Always(S')` | `OptimalControlProblem` (reach-avoid) | same, restricted to the `safe_set` `S'` |
| `Always(S)` | `SafetyProblem` | maximal controlled-invariant set |
| `EventuallyAlways(S)` | `ReachAndStayProblem` | invariance, then reachability |
| a formula over `Label`s | `CoSafeLTLProblem` | product with a deterministic automaton |
| nothing at all | `AlternatingSimulationProblem` | build the abstraction and return it |

Use a named pattern when your specification is one of the four shapes: they reach dedicated
fixed-point algorithms that are much faster than an automaton product. Use a formula when it is
not. Writing a formula **takes precedence** if you write both.

Note that a formula is always solved by the co-safe machinery, which handles properties satisfied
by a finite prefix (`F`, `U`, `X`). `G` and `F(G)` are not of that kind — express those with
`Always` and `EventuallyAlways`.

### 5.4 Costs

`@objective` carries the cost, never the specification:

```julia
@objective(model, Min, sum(u[i]^2 for i in 1:2))
```

It is not implemented yet, and raises an error rather than being ignored. A control cost is the
`transition_cost` of `Problem.OptimalControlProblem`; until the front-end can set it, build that
problem directly and hand it to a solver through the direct-MOI entry style.

### 5.5 Avoiding a region: `∉` and `Always` differ

`x ∉ O` and `x in Always(complement(O))` say the same thing logically, but do different things:

* **`x ∉ O` restricts the system.** The region is carved out of the state space `X`, so those cells
  are never created. Fewer states, faster abstraction. Use it for walls, obstacles, and anywhere
  the model simply is not valid.
* **`x in Always(S)` constrains the specification.** States outside `S` stay representable, so
  synthesis can reason about leaving them and actively prevent it. Use it when "stay inside `S`" is
  the *goal*.

Rule of thumb: `∉` describes the **world**, `Always` describes the **requirement**.

`Always` means the same thing everywhere it appears: it becomes the `safe_set` of the problem —
`SafetyProblem` on its own, `OptimalControlProblem` when written together with a `Final` set. The
reach-avoid reading is `safe U target`: every state up to and including the one that reaches the
target must be safe, and what happens afterwards is not constrained.

### 5.6 Horizon

```julia
set_attribute(model, "horizon", 10.0)     # default: infinite
```

* **Continuous-time model (`∂`)** — the horizon is in **seconds**, converted to a step count using
  your `time_step`, rounded in the direction the specification requires ("within at most `T`"
  rounds down, "for at least `T`" rounds up).
* **Discrete-time model (`Δ`)** — the horizon is a **step count**, used as given.

> **Caveat:** a horizon on a *safety* problem has no effect. The safety solver computes a maximal
> controlled-invariant set, which is an infinite-horizon notion, and never reads the field.

---

## 6. Hybrid systems — modes and transitions

A **mode** is a scope you attach dynamics, bounds and specifications to. A **transition** is a
scope you attach a guard and a reset map to. Both accept ordinary `@constraint`s.

`@mode` binds the name in your scope and registers it in the model, exactly like `@variable` —
write `@mode(model, off)`, not `off = @mode(model, off)`.

```julia
model = Model(Dionysos.Optimizer)
@variable(model, 17 <= T <= 25)
@variable(model, 0 <= u <= 1)

@mode(model, off)
@mode(model, on)

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

**Guard or reset?** A constraint containing `∂`/`Δ` is a **reset map**; anything else is a
**guard**. The same rule separates per-mode dynamics from per-mode bounds.

**Several transitions between the same two modes** are fine — each is its own object, so their
guards and resets never get mixed up.

**A guard need not be a box.** Everything written on the transition intersects:

```julia
add_transition!(model, a => b) do t
    @constraint(t, x >= 1)                              # per-coordinate bound
    @constraint(t, x + y <= 3)                          # half-space
    @constraint(t, [x, y] in Guard(LazySets.Ball2([0.0, 0.0], 1.0)))   # any bounded LazySet
end
```

A `Guard` set must span the state vector, in declaration order — same rule as `Final`/`Always`.
The marker is needed because JuMP cannot tell a bare set apart from a bound. Guards must stay
affine or set-valued; a *nonlinear* guard is rejected.

On a **clocked** model a guard must be a plain box: a non-box guard would have to be extruded
across the time axis, and the discretisation cannot enumerate that product.

Only `HybridSystems.AutonomousSwitching` (the default) is supported — the switch is taken
because the state entered the guard. `ControlledSwitching` is rejected rather than accepted and
ignored, since the abstraction builds switch transitions from the guard alone. Model a
controller-chosen switch by writing the choice into the guard.

Solver options are set **per mode**, since each mode is abstracted by its own sub-solver:

```julia
for m in (off, on)
    set_attribute(m, "state_grid", MP.GridFree(SVector(0.0), SVector(0.1)))
    set_attribute(m, "input_grid", MP.GridFree(SVector(0.0), SVector(0.1)))
    set_attribute(m, "time_step", 0.1)
end
```

Options set on the model itself apply to every mode unless a mode overrides them.

### Clocks and timed specifications

A clock is a variable whose dynamics is `∂(t) == 1` (running) or `∂(t) == 0` (frozen) — rule I2.
No new syntax:

```julia
@variable(model, 0 <= t <= 50)
@constraint(off, ∂(t) == 1)          # clock runs while off
@constraint(on,  ∂(t) == 0)          # clock frozen while on
```

A time window on a target is then a reach constraint on the clock:

```julia
@constraint(on, final(T) in MOI.Interval(21.0, 23.0))
@constraint(on, final(t) in MOI.Interval(15.0, 50.0))   # … within t ∈ [15, 50]
```

Specifications written on a mode apply in that mode only; with a clock they become time-windowed.

---

## 7. Choosing a solver

A model with modes goes to the hybrid abstraction, everything else to the uniform grid. Override
it explicitly:

```julia
set_attribute(model, "solver", AB.UniformGridAbstraction.Optimizer)
```

Every other `set_attribute` is forwarded to the chosen solver and its sub-solvers, so all the
options documented on those optimizers (`approx_mode`, `jacobian_bound`, `print_level`,
`execution_backend`, …) are reachable. An attribute no solver recognises raises
`MOI.UnsupportedAttribute` naming it, rather than being silently dropped. A solver that cannot
handle your problem type is reported here, naming both.

### Growth bounds come for free

The `GROWTH` approximation mode does not need a hand-written `jacobian_bound`: with Symbolics and
IntervalArithmetic loaded, one is derived by tracing the dynamics and bounding the Jacobian over
`X`. Supply your own only to override it.

---

## 8. Reading the results

### Did it work?

```julia
optimize!(model)

termination_status(model)          # OPTIMAL | LOCALLY_INFEASIBLE | OPTIMIZE_NOT_CALLED
is_solved_and_feasible(model)      # true ⇒ you have a certified controller
raw_status(model)                  # human-readable reason
```

Failure is `LOCALLY_INFEASIBLE`, not `INFEASIBLE`, and the distinction matters: the abstraction is
**sound but not complete**, so "no controller was found" means no controller exists *on this
abstraction*. A finer grid, a smaller time step, or a different `approx_mode` may well succeed.
Nothing here proves your problem is unsolvable.

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

One call for every kind of model. For a hybrid model `x0` is the augmented state, `(x, mode)` or
`(x, t, mode)`. The result is a channelled `Dionysos.System.Trajectory`; read it with `ST.states`,
`ST.inputs`, `ST.times`, `ST.modes`, plot it with `plot(traj)`, or animate it with
`Dionysos.animate_trajectory_dashboard`.

Stopping is derived from your specification by default — pass `stopping = x -> …` to override.

---

## 9. Folder structure

```
src/wrapper/
  wrapper.jl              module + includes; loaded from src/Dionysos.jl
  variables.jl            VariableInfo, roles, the role inference, set_role!
  specification.jl        specification markers, Label, @specification
  model_ir.jl             ModelIR — the parsed model, printable, dependency-free
  modes.jl                Mode and Transition scopes, @mode, add_transition!
  clock.jl                clock detection and the per-mode clock subsystem
  operators.jl            ∂ Δ final start, and the ∉ obstacle syntax
  parse.jl                MOI constraints → ModelIR
  parse_scoped.jl         the same, for constraints written on a mode or transition
  lower_system.jl         ModelIR → MathematicalSystems
  lower_problem.jl        ModelIR → Problem.ProblemType
  lower_hybrid.jl         ModelIR with modes → HybridSystems.HybridSystem
  solver_selection.jl     select_solver, supports_problem
  dynamics_backend.jl     AbstractDynamicsBackend + the evaluator backend
  results.jl              solution status and simulate
  optimizer.jl            Dionysos.Optimizer <: MOI.AbstractOptimizer

ext/DionysosMathOptSymbolicADExt.jl    SymbolicADBackend (codegen + automatic Jacobian bound)
```

The pipeline is one-directional:

```
JuMP model ──parse──▶ ModelIR ──lower──▶ (system, problem) ──select──▶ solver ──▶ controller
```

`ModelIR` is why the front-end lives in `src/` rather than in an extension: parsing needs no
Symbolics, only *compiling a dynamics expression into a callable* does. That is what makes it
loadable, testable and documentable without any optional dependency.

`Dionysos.Wrapper.lower(model)` runs the first half on its own, if you want to inspect the
`(system, problem)` pair without paying for an abstraction.

---

## 10. Extending it

| To add… | Touch |
| :--- | :--- |
| a specification kind | one marker set in `specification.jl` + one branch in `lower_problem.jl` |
| a problem type | one `lower_problem.jl` method + one `supports_problem` declaration |
| a solver family | one `select_solver` method + one `supports_problem` declaration |
| a dynamics source | one `AbstractDynamicsBackend` subtype + `compile_dynamics` |
| a system class | one `lower_system.jl` method |

Nothing in that table requires touching `parse.jl`.

Two rules for contributors:

* **The wrapper owns no semantics.** If a concept has no `MathematicalSystems`/`HybridSystems` +
  `ProblemType` representation, add it to `System`/`Problem` first, then expose it here.
* **No new type piracy.** The `∉` parsing rule is a deliberate, documented exception because JuMP
  offers no extension point for it. `Mode`/`Transition` subtype `JuMP.AbstractModel`, which gives
  `@constraint` support with no piracy at all.

---

## 11. Gotchas

* **Mixing `∂` and `Δ`** in one model is an error — pick continuous or discrete time.
* **Variable order is declaration order.** Interleaving state and input declarations reorders the
  vectors. Check with `get_attribute(model, "state_variables")`.
* **`final(x)` on a whole vector does not work.** `final` is a `JuMP.NonlinearOperator`, which only
  builds an expression when one of its arguments is a JuMP scalar. Use the vector form
  `x in Final(S)`, or the scalar form `final(x[i]) in MOI.Interval(a, b)`.
* **`start = v` gives a single starting point**, not a region. Use `x in Start(S)` for a set.
* **Nonlinear guards are rejected.** A guard is built from affine constraints and `Guard` sets;
  an arbitrary nonlinear expression has no set representation the abstraction can enumerate.
* **A mode's own bounds stay per-coordinate.** `@constraint(mode, x + y <= 1)` is rejected: a
  mode's state and input sets are boxes the abstraction discretizes axis by axis. The same
  constraint *on a transition* is fine — it becomes a half-space of the guard.
* **A horizon on a safety problem does nothing** — see §5.6.
* **Mixing clock-lifted and time-free modes in one hybrid model is rejected.** Guards and reset maps
  apply to the augmented state, whose arity differs between the two (`[x; t]` vs `x`), so a mixed
  model has no consistent reset convention. Give every mode a clock, or none.
* **`Always` and `∉` are not interchangeable** — see §5.5.
* **Integer / binary variables are not supported.** They belong to the MIQP solvers, which this
  front-end does not target.
* **Only two solver families are reachable from JuMP** — the uniform grid abstraction and the hybrid
  abstraction. The others (Bemporad–Morari, branch-and-bound, the ellipsoidal families, the PCLF
  bisimulation quotient, trajectory generators and certifiers) need inputs this DSL cannot express,
  such as PWA systems, ellipsoid templates, or observation regions. Use the direct-MOI entry style
  for those.

---

## 12. See also

* [`docs/src/examples/jump/`](../../docs/src/examples/jump/) — runnable front-end examples.
* [`docs/src/manual/abstraction-based-control.md`](../../docs/src/manual/abstraction-based-control.md) — what the solvers actually do.
* [`docs/src/developers/conventions.md`](../../docs/src/developers/conventions.md) — house style.
