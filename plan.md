# Robust control in Dionysos

A plan for making **robust control** — a controller that must work against an adversarial
environment — the general case every solver in the toolbox solves, and for obtaining the two
questions Dionysos users actually ask as its degenerate cases.

The general problem is the game

> `∃u ∀w` : does the controller have a choice `u` such that **every** environment choice `w` yields
> a trajectory satisfying `φ`?

and the two extremes are:

| | `U` (controller) | `W` (environment) | answer |
| :--- | :--- | :--- | :--- |
| **pure synthesis** | non-trivial | absent | a controller — this is Dionysos today |
| **robust synthesis** | non-trivial | non-trivial | a controller valid against every `w` |
| **pure verification** | absent | non-trivial | a set, and on failure a counterexample |

Robust control is the subject. Pure synthesis and pure verification are not two features to build:
they are what the same code does when one side of the game is empty. **Everything below is written
so that no solver ever learns which case it is in.**

Three invariants bound every design choice below, and each later section must answer to them:

- **I1 — pure synthesis does not move.** Not "stays fast": the same machine code on the same
  path. No wrapper, no extra type parameter, no branch is added to the path a no-`W` problem
  takes, and the existing goldens stay byte-identical. Enforced by a benchmark, not by intent.
- **I2 — the abstraction cache survives.** "The abstraction is cached: switching the control
  task on the same system does not recompute it" is the framework's keystone promise, and
  synthesis-then-verification on one system is exactly such a task switch. Both questions must
  be answerable from **one** abstraction.
- **I3 — one concept, one mechanism.** The environment is folded into successor sets — never
  into the alphabet, never into the solvers. Anything that needs a second mechanism for the
  same concept is redesigned until it does not.

Written 2026-09-01 on `jc_pclf_abstraction`. This supersedes nothing: the previous repo-root
`plan.md` was the v0.2 refactor roadmap, removed in `9a1a9b739` and recoverable with
`git show fb2266793:plan.md`.

---

## 1. Why this framing

### The modelling observation it rests on

For a switched system the mode signal is *one* mathematical object, but it plays two roles:

| | the mode signal is | quantifier | the answer is |
| :--- | :--- | :--- | :--- |
| synthesis | a control input the designer picks | `∃` | a controller |
| verification | a disturbance the environment picks | `∀` | a set + counterexample |

So "for all inputs" is really "for all **disturbances**": under `⊨_∀` the signal is not an input at
all, it is uncontrolled. **This is the design lever.** The quantifier is not a property of the
specification — `SafetyProblem` means the same thing in both worlds — it is a property of *who owns
the signal*, and it therefore belongs in the system description rather than in the `ProblemType`.

Building only verification would be less code now and the wrong shape later: it would add a second
`∀` path beside the existing `∃` one, and the mixed case — the one an engineer actually wants,
a controller that survives disturbance — would then have nowhere to live. Building the game gets all
three for one implementation, and the two extremes serve as its regression tests.

### The example that validates it

Gol, Ding, Lazar & Belta (arXiv:1208.5471) state *both* extremes on one system:

- their **Problem 5.1** is synthesis — the largest set of initial states from which *some* switching
  sequence satisfies `φ`;
- their **Problem 5.2** is verification — the set from which *every* switching sequence does.

`research/BisimulationQuotient/gol_lazar_belta_example.jl` reproduces 5.1 and says in prose that 5.2
is out of reach of the current optimizer. That pair is the ideal acceptance test for the general
code: same system, same specification, same abstraction, same solver call — only the ownership of
the switching signal differs. It also completes the reproduction of the predecessor we claim to
generalise.

**It is a demonstration, not the design target.** Nothing in §§2–4 is specific to switched systems,
to that example, or to the PCLF quotient; the example appears again only in the phasing (§5), as
what P2 must produce.

---

## 2. What exists today

**Everything is `∃`.** Every solver in `src/optim/discrete_systems/` computes a fixed point of the
*controlled* predecessor `CPre(Y) = {q : ∃u, Post(q,u) ⊆ Y}`.

The good news is that the `∀` half is **already implemented**. In
[`compute_largest_invariant_set`](src/optim/discrete_systems/safety_problem.jl#L45), `pairstable[q,u]`
is disabled as soon as *any* successor of the pair leaves the safe set — that is `∀` over successors,
i.e. over the nondeterminism the abstraction introduced. What is `∃` is only the reduction from
pairs to states:

```julia
# safety_problem.jl:90 — the state survives while SOME input still works
if safeset[source] && nsymbolslist[source] == 0
    nextunsafeset[source] = true
end
```

`nsymbolslist[q]` counts the inputs still enabled at `q`, and `q` dies when that count reaches
**zero**.

**That `∀` over successors is the whole of what `∀w` needs**, provided the environment's choices are
already inside the successor set. Which is the design of §4: fold `W` into the successor sets the
solver sees — lazily for a finite `W`, in the reachable set for a continuous one — and this loop,
unmodified, computes `∃u ∀w`. It is worth stating plainly because
the obvious alternative is to teach this loop a second threshold, and that turns out to be
unnecessary work on the wrong file.

**Disturbances are recognised, then discarded.** `MathematicalSystems` ships a full `Noisy*` family
(`NoisyConstrainedLinearControlDiscreteSystem` and 15 siblings) with `noiseset`, `noisedim` and the
`isnoisy` trait. Dionysos detects them and throws the information away:

```julia
# alternating_simulation_problem.jl:663
if concrete_system isa MS.AbstractSystem && MS.isnoisy(concrete_system) === true
    @warn("Noise is not yet accounted for in system abstraction.")
end
```

The vocabulary exists at the boundary and stops there.

---

## 3. Soundness — the part that is easy to get wrong

An abstraction sound for one quantifier need not be sound for the other. This section is why the
plan is not simply "flip the counter".

**Both quantifiers want an over-approximation, but they fail differently.** If the abstraction has
*more* behaviours than the concrete system:

- `∃` — a spurious transition can only *lose* a controller. Conservative, sound.
- `∀` — a spurious transition only *strengthens* the adversary. Conservative, sound.

Over-approximation is therefore the right direction for both. The asymmetry appears when the
abstraction has **fewer** behaviours than the concrete system: harmless-to-conservative for `∃`,
**unsound** for `∀`, because a move the real adversary has was silently removed.

### Which relation each quantifier actually needs

These are different notions and it is worth separating them, because the stronger-sounding one is
needed for the *easier*-sounding question:

| | relation needed | what is being transferred |
| :--- | :--- | :--- |
| synthesis `∃u ∀` | **alternating** simulation | a *strategy* — the controller must beat the abstraction's nondeterminism, which plain simulation says nothing about |
| verification `∀` | plain simulation (behaviour over-approximation) | a *property* — if no abstract run violates `φ`, no concrete run does |

**So verification does not need more than the existing grid abstraction gives.** The construction at
[abstraction.jl:76](src/symbolic/grid_based_symbolic_model/abstraction.jl#L76) covers the reachable
set of the *whole cell* with `MP.OUTER`, so every concrete successor's cell is an abstract successor:
that is exactly the simulation `∀` wants. Over-approximating a set of behaviours to prove a universal
property is the classical use of abstraction in model checking, and it applies here unchanged.

What `∀` is intolerant of is not the *kind* of relation but any **missing** behaviour. Four such
places exist today, all of them harmless under `∃`:

1. **Sampling approximation modes are under-approximations.** `GROWTH`, `LINEARIZED` and a
   sound `USER_DEFINED` over-approximate and are fine. `CENTER_SIMULATION` and `RANDOM_SIMULATION`
   sample the cell, so they are *under*-approximations of the reachable set: valid heuristics for
   `∃`, unsound for `∀`. Any run with a non-empty `W` must reject them by name.

2. **A dropped transition becomes a vacuous proof.** `compute_abstract_transitions_from_set!` returns
   `false` — adding *no* transitions for the pair — when the reachable set is not entirely inside the
   domain (`allin || return false`) or lands on a disallowed state. Under `∃` that merely costs the
   controller an input. Under `∀`, a state with no successors satisfies `∀u: Post(q,u) ⊆ Y`
   *vacuously* and is reported as verified. The same hazard appears as `atol` coverage loss (about
   3.7% of the domain dropped at `1e-3`, about 0.4% at `1e-4`) and as deadend states —
   `lifted_showcase.jl` reports 105 of 1143 with outgoing degree 0. Both are recorded as open in
   `research/BisimulationQuotient/README.md`.

3. **A gridded adversary is not the real adversary.** The input alphabet is a finite grid. For `∃`
   that is conservative: a controller found using grid inputs is still a controller. For `∀` it is
   unsound, because the environment is not restricted to grid points. This is the one case where the
   abstraction genuinely has to change: the disturbance must enter the reachable-set computation *as
   a set*, the way the state cell already does, rather than being enumerated. It is also exactly what
   `build_noise` currently refuses to do.
   **A switched system escapes this entirely** — its mode set is genuinely finite, so the gridded
   adversary *is* the real one. That, and not the choice of simulation relation, is why the switched
   and PCLF cases come first in this plan.

4. **Graph completeness, for the lifted PCLF abstraction.** Path-completeness licenses synthesis but
   not verification: a node with no outgoing edge for some mode removes a move the adversary has —
   an instance of hazard 2 at the graph level. `PCLF.is_complete`
   ([pclf.jl:43](src/utils/pclf.jl#L43)) already tests exactly this and its docstring already records
   why. Verification must *check* it, never assume it.

**Consequence for the API.** These are gates on any run with a non-empty `W`, not on a separate
"verification mode": a solver must refuse, loudly, on a sampling approximation mode, on an
abstraction that dropped transitions, and on a non-finite environment alphabet that was gridded.
Silent conservatism is acceptable; silent optimism is not.

---

## 4. Design

### 4.0 How the user says it

**Through the system type. Nothing new to learn.** `MathematicalSystems` already carries the
distinction, via `MS.iscontrolled` / `MS.isnoisy` and the accessors `MS.inputset` / `MS.noiseset`,
all exported and already reachable at the Dionysos boundary.

But its catalogue is **not** uniform, and the asymmetry decides the design:

- **noise-only (uncontrolled) types exist for `Linear` alone** — `NoisyLinearDiscreteSystem`,
  `NoisyConstrainedLinearDiscreteSystem` and their continuous twins. There is no Affine and no
  black-box equivalent.
- **the black-box family is noisy only in *controlled* form** —
  `NoisyConstrainedBlackBoxControlDiscreteSystem`, fields `(:f, statedim, inputdim, noisedim, X, U, W)`.

Since the JuMP front-end emits black-box types
([lower_system.jl:139](src/wrapper/lower_system.jl#L139)), a nonlinear system has *no* "just noisy"
type to build. **That is fine**: pure verification is a **singleton input set**, which by §4.1
yields an abstraction with a single input symbol — and that is exactly the verification case.

So the toolbox needs one noisy family, not three, and the three questions are three shapes of
`(U, W)`:

| `U` | `W` | abstract alphabet | the question is |
| :--- | :--- | :--- | :--- |
| non-trivial | absent | one symbol per `u` | **synthesis** — today |
| **singleton** | non-trivial | a single symbol | **verification** |
| non-trivial | non-trivial | one symbol per `u` | **robust synthesis** |

`NoisyConstrainedBlackBoxControlDiscreteSystem` covers all three rows on its own; the linear
noise-only types are a convenience for users whose system happens to be linear, and lower onto the
singleton-`U` case internally. Only the `Constrained…` variants are usable at all, since Dionysos
needs `X` — and now `W` — as sets.

**Nothing beyond the system object is ever declared.** `MS.inputset` and `MS.noiseset` are the whole
interface the abstraction reads (§4.1), and a singleton or absent side degenerates the alphabet in
the right direction on its own. The user's model *is* the declaration.

Three places already stub this out and consume nothing:

- `MS.isnoisy` is detected and the information discarded, with a warning
  ([alternating_simulation_problem.jl:663](src/optim/continuous_systems/UniformGridAbstraction/alternating_simulation_problem.jl#L663));
- the JuMP front-end has a `DISTURBANCE` variable role
  ([variables.jl:19](src/wrapper/variables.jl#L19)), exported from
  [Dionysos.jl:38](src/Dionysos.jl#L38) and referenced nowhere else;
- `build_noise` returns a zero vector of the right length.

The work is to connect them, not to invent vocabulary.

#### The JuMP front-end

Roles are *inferred* from usage ([README](src/wrapper/README.md) §3), and a disturbance is
indistinguishable from an input that way — both merely appear on a right-hand side. So it has to be
said explicitly, which is what `declared_role` and `set_role!` already exist for ("set when
inference cannot see the truth"). The model then reads:

```julia
@variable(model, -1 <= u[1:2] <= 1)          # inferred INPUT      -- the controller chooses
@variable(model, -0.05 <= w[1:2] <= 0.05)    # declared DISTURBANCE -- the environment chooses
@constraint(model, ∂(x[1]) == u[1] * cos(x[3]) + w[1])
```

and `lower_system.jl` emits `NoisyConstrainedBlackBoxControl…` instead of
`ConstrainedBlackBoxControl…` whenever any variable carries the disturbance role. A model with
disturbances and *no* inputs still lowers to that same type, with a singleton `U` synthesized — the
verification row of the table above. **Open**: the exact surface for saying it — `set_role!`, a
`@disturbance` macro, or an `in Disturbance()` set — is undecided; see §7.

#### Switched and hybrid systems

**`HybridSystems` already says this, and Dionysos should simply read it.** The package carries the
distinction as a type, citing Liberzon §1.1.3:

```julia
abstract type AbstractSwitching end
struct AutonomousSwitching <: AbstractSwitching end   # the switching signal is not the designer's
struct ControlledSwitching <: AbstractSwitching end   # the switching signal is a control input
```

and `HybridSystem` stores it in a `switchings` field, indexed by mode label and exposed as the
system's fourth type parameter. The mapping onto §4.1 is immediate and needs no Dionysos vocabulary
at all:

| `switchings` says | the mode is | the abstraction | the question |
| :--- | :--- | :--- | :--- |
| `ControlledSwitching` | `U` | one symbol per mode | synthesis — GLB Problem 5.1 |
| `AutonomousSwitching` | `W` | modes folded into one symbol | verification — GLB Problem 5.2 |

A Dionysos wrapper type and a raw `"switching"` optimizer attribute were both considered and
rejected. The attribute would put half the model in the model and half in the solver configuration,
contradicting this section's whole argument, and would default silently — the silent optimism §3
forbids. The wrapper would duplicate a distinction the ecosystem already draws.

**One migration hazard, and it is not small.**
`HybridSystems.discreteswitchedsystem` hardcodes `Fill(AutonomousSwitching(), 1)`, while every
switched-system solver in Dionysos today treats the mode as a *control* input. So the label already
present on every switched system in this repository says "not the designer's", and the code assumes
the opposite. Honouring the field naively would silently turn every existing example into a
verification query.

Two consequences for the work:

- **P2 must migrate, not merely read.** Either the examples construct their systems with
  `ControlledSwitching`, or Dionysos offers a constructor that does — `discreteswitchedsystem`
  itself takes no option for it, so a thin Dionysos helper is justified *here*, as a convenience for
  building the object, not as a competing way to express the concept.
- **The reproduction script becomes a one-word diff**, which is what P2 owes: the same
  `gol_lazar_belta_problem()` built once with `ControlledSwitching` for 5.1 and once with
  `AutonomousSwitching` for 5.2.

Downstream nothing else is special: the single accessor of §4.8 is the only interface the
abstraction reads, and the switched system reaches it by the same path as a continuous one with a
disturbance.

### 4.1 Mirror the representation, not the alphabet

The instinct is to treat `U` and `W` symmetrically all the way down: a `WMapping` and a `Wset`
beside `UMapping` and `Uset`, a `w` component in the input symbol, and a quantifier consulted at
solve time. Half of that is right and the other half is both unsound and unnecessary.

**Mirror the representation.** `W` is a first-class citizen of the abstraction's inputs: a set the
user gives, carried and consumed like `X` and `U` are. But note what "like `U`" does *not* include:
`U` gets a `UMapping` because the `Mapping` layer exists to turn sets into integer alphabets, and
`W` must never have an alphabet — so **there is no `WMapping`**, by the same principle (I3) that
keeps `w` out of the symbols. `W` reaches the kernels as a set (`w_c`, `w̄`, or the box
approximation of a fancier set) and reaches the switched case as the mode labels the automaton
already has. A mapping for `W` would not be unused generality; it would be an invitation to the
exact unsoundness §3 hazard 3 describes.

**Do not mirror the alphabet.** The abstract input symbol enumerates `U` only. `W` never becomes a
symbol; it is *folded into the successor set*:

```
Post(q, u)  :=  ⋃  Post(q, u, w)
               w∈W
```

so that one abstract input symbol `u` carries every environment reaction to it.

The asymmetry is forced, and it is the difference between the two roles:

- the **controller** must be *told which* `u` to apply, so `u` needs a name in the alphabet;
- the **environment** is never instructed, so `w` needs no name — only its effect. Worse, giving it
  one is unsound: enumerating a continuum at grid points removes moves the real adversary has
  (§3, hazard 3), whereas folding `W` in as a *set* over-approximates it, which is exactly the
  direction `∀` requires.

### 4.2 Fold lazily where possible, at build time only where necessary


*Where* the fold happens is a separate decision from *whether*, and it is governed by two
constraints stated in §0: the abstraction cache must survive, and the pure-synthesis path must not
move. That yields one principle with two mechanisms:

**Finite `W` (the switched case): fold at solve time, as a view.** Per-mode transitions unioned at
the automaton level are *identical* to transitions built from a folded successor set — for a finite
`W` the fold commutes with abstraction, exactly. So the abstraction is built as today, per mode,
once; and a `FoldedAutomaton <: SY.AbstractAutomatonList` view presents it with the mode dimension
collapsed:

```julia
# One screen of code against the existing seven-method interface
# (src/symbolic/automata/automaton.jl): get_n_state unchanged, get_n_input
# reports |U|, pre/post/enum_transitions map the underlying symbol through
# symbol -> (u, w) and drop the w.
struct FoldedAutomaton{A <: AbstractAutomatonList} <: AbstractAutomatonList
    inner::A
    nu::Int      # folded alphabet size; symbol s of inner maps to u = fold(s)
end
```

Nothing is copied, nothing is rebuilt, and — decisive for the cache — **one abstraction answers
both questions.** Ask synthesis, get the automaton; ask verification, get the view over the *same*
automaton. This is what makes P2's acceptance test (§5) achievable at all: under build-time folding
the PCLF refinement, which runs per-mode pre-images
([pclf_bisimulation_quotient.jl:261](src/optim/hybrid_systems/PCLFBisimulationQuotient/pclf_bisimulation_quotient.jl#L261)),
would produce a *different, coarser quotient* for the folded system, and "same abstraction, one-word
diff" would be false by construction.

**Set-valued `W` (continuous disturbance): fold at build time, into the reachable set.** There is no
finite alphabet to view; the disturbance enters the reachable-set computation as a set, the way the
state cell already does. Concretely this is the Reissig–Weber–Rungger growth bound for perturbed
systems — `ṙ = L(u)·r + z` with the nominal at the centre of `W` — so it costs no extra reach-set
computations and never loops over `w`; P5 gives the kernel-by-kernel account. This is the only
place the abstraction builder changes.

The counter machinery is indifferent to which mechanism produced the successor set: a `(q, u)` pair
of the view is disabled by the first bad transition of *any* of its modes, which is precisely `∀w`.
Memory even improves — the counters are `nstates × |U|` instead of `nstates × |U|·|W|`.

### 4.3 Why folding loses nothing

Every fixed point in the toolbox tests one predicate — *is the whole successor set inside `Y`* — and
that predicate does not distinguish a union from a quantifier:

```
∀w : Post(q,u,w) ⊆ Y      ⟺      ⋃ Post(q,u,w) ⊆ Y
                                 w
```

The worst-case cost composes the same way — the optimal-control solver already takes the **maximum
cost over successors** ([optimal_control_problem.jl:123](src/optim/discrete_systems/optimal_control_problem.jl#L123)),
and the max over a union is the max over its parts. One genuine caveat: a stage cost that *depends
on `w`* must be maximised over the modes folded together before it is attached to the folded pair.
Whether `transition_cost` signatures permit a `w`-dependent cost today must be checked in P0, not
assumed either way.

Folding is exact for the information pattern being solved: in `∃u ∀w` the controller commits to `u`
before the environment moves, which is what a union expresses. A controller allowed to *observe* `w`
first would be playing `∀w ∃u` — a different game needing `w` in the alphabet; §7 records it as
out of scope.

### 4.4 The solvers do not change at all

This is the payoff, and it is now structural rather than aspirational: the solvers operate against
`SY.AbstractAutomatonList` already, and the view *is* an `AbstractAutomatonList`.

Look again at [`compute_largest_invariant_set`](src/optim/discrete_systems/safety_problem.jl#L45):
`pairstable[q,u]` is disabled as soon as **any** successor of the pair leaves the safe set. Over the
folded view that is already the `∀w` we need — and `nsymbolslist[q] > 0`, the surviving `∃u`, is
already the `∃u`.

> **The existing controlled predecessor over a folded automaton is the robust predecessor.**
> `{q : ∃u, Post(q,u) ⊆ Y}` over the view *is* `{q : ∃u ∀w, Post(q,u,w) ⊆ Y}`.

So no counter is regrouped, no partition object reaches a solver, no `if verification` is written
anywhere, and `src/optim/discrete_systems/` is not touched. The three cases are three *automata*
handed to the same solvers:

| the solver receives | alphabet | it computes |
| :--- | :--- | :--- |
| the automaton, as today | one symbol per `u` (or per mode) | **pure synthesis** — today, bit for bit |
| the view, `U` non-trivial | one symbol per `u` | **robust synthesis** |
| the view, `U` a singleton | a single symbol | **pure verification** |

Pure verification is the view whose alphabet has one symbol: `∃u` over a singleton is vacuous, so
what remains is exactly `∀w`. It needs no solver support because there is nothing left to choose.

### 4.5 What the controller becomes

Nothing changes. The controller stores the winning input symbol, which is a `u` in all three rows
above, and the closed-loop runtime steps with it exactly as today — the environment supplies `w` to
the *concrete* system, which was always true and was simply not modelled.

Where the single-symbol view is used, every state that survives is assigned the one symbol, the
controller carries no information, and the useful part of the result is the winning set. That is why
verification "returns a set and not a controller": not a separate code path, just a degenerate
alphabet.

The one thing worth keeping is *which* `w` produced a given transition, for counterexamples (§4.6).
The view knows it for free — its `symbol -> (u, w)` map is the record — and for build-time folding
the witness goes in the per-transition metadata that already exists
([abstraction.jl:47](src/symbolic/grid_based_symbolic_model/abstraction.jl#L47),
`Pair{TransitionKey, Any}`). An annotation, never an alphabet symbol.

### 4.6 What comes back

One result shape for the general problem, whose parts degenerate on their own:

- `winning_set` — the states from which the controller can force `φ` against every `w`. Under an
  empty `W` this is today's winning set; under an empty `U` it is the verified set. Same field.
- `controller` — the winning group per state (§4.4). Vacuous when there is nothing to choose, which
  is what makes "verification returns no controller" a consequence rather than a special case.
- `counterexample` — on failure, an environment word (finite for a co-safe violation, a lasso for a
  safety violation) together with the concrete trajectory realising it.

The counterexample is not garnish, and it is not verification-specific: a failed *robust* synthesis
needs to say which disturbance defeats the controller just as much. Replaying an abstract
counterexample concretely is also how an unsound abstraction gets caught — if it does not reproduce,
the abstraction is at fault, not the system.

### 4.7 Co-safe LTL under `∀`

*(Revised after P2: an earlier version of this section predicted a complementation argument would
be needed — it is not, and P2 proved it by running.)*

For co-safe `φ`, satisfaction is having a good prefix, and the monitor accepts exactly the good
prefixes. So "every run satisfies `φ`" ⟺ "every product run reaches the accepting set", which is
∀-reachability over the product — the machinery of §4.4, unchanged, over the folded product. The
direct route works; no complementation, no new fixed point.

What the quantifier *does* change is the standing of two assumptions the product machinery makes,
both conservative under `∃` and **silently optimistic** under `∀`:

1. **The monitor must be total.** A dropped letter drops a run; under `∃` that loses a controller,
   under `∀` it verifies a state whose violating run disappeared. The spot stepper is total today —
   missing DRA edges clamp to an explicit absorbing dead state — but the monitor is an interface,
   and the contract must be stated and gated, not inherited by luck.
2. **The accepting set must be exactly the good-prefix states.** Today it is a self-described
   heuristic (`_cosafe_done_states_dra`): a state is "done" iff every *defined* edge is a
   self-loop, which counts a missing edge — a step into the rejecting dead state — as a self-loop,
   and would classify an explicit rejecting sink as accepting. An over-large accepting set under
   `∃` yields a controller whose trajectory visibly fails `φ`; under `∀` it verifies violating
   states and nothing downstream notices.

### 4.8 Where each piece lives

The design decomposes onto the module chain `Utils → System → … → Symbolic → Optim` with **one
seam per module**, which is what makes it maintainable: each concern has exactly one home, and the
homes respect the existing dependency order.

**`System` — declaration.** One generic accessor turns the ecosystem's scattered vocabulary into a
single internal question:

```julia
"""
The environment's choice set of `system` — the modes of an autonomously switched system, the
noise set of a noisy one — or `nothing` when the environment has no move. This is the only
question the rest of the toolbox ever asks; how a system declares ownership is dispatch, here.
"""
environment_input(::MS.AbstractSystem) = nothing            # default: the environment has no move
environment_input(s::MS.NoisyConstrainedBlackBoxControlDiscreteSystem) = MS.noiseset(s)
# … one method per Noisy… family (or one via the isnoisy trait), and:
environment_input(s::HybridSystems.HybridSystem) = _modes_if_autonomous(s)   # switchings-dispatch
```

This is the standard Julia interface pattern the house style already mandates (abstract question,
per-type methods): `MS.isnoisy`/`MS.noiseset` and `HybridSystems.switchings` are *implementations*
behind it, and nothing downstream ever mentions them. A future ownership vocabulary — a stochastic
noise model, a new system package — is one new method, zero edits elsewhere. Without this seam, the
ownership read is an if-else smeared across every composing optimizer; with it, extensibility is
the dispatch table.

**`Symbolic` — mechanism.** `FoldedAutomaton` in `src/symbolic/automata/folded_automaton.jl`,
beside the three existing backends, generic over `A <: AbstractAutomatonList` so every backend and
any future one is covered. Its contract, per the conventions:

- read-only: the enumeration half of the interface is implemented; `add_transition!` and
  `Base.empty!` throw with a message saying views are immutable — never silently no-op;
- concrete parametric fields, no closures in fields: the `symbol → (u, w)` fold and its `u →
  symbols` inverse (needed by `post`) are plain stored data, type-stable;
- knows nothing about systems, ownership, or solvers — it is an automaton combinator, testable on
  hand-built automata alone.

**`Optim` — policy.** The composing optimizers ask `environment_input`, and hand their control
solver either the automaton or the view — plus the §3 gates at that hand-over. Policy is the only
part that touches both worlds, and it is a few lines per optimizer because declaration and
mechanism each already did their half.

**`Mapping` is untouched** (there is no `WMapping`, §4.1), and `Problem` is untouched (§4.0's first
principle). The one addition outside the chain is the kernel arithmetic of P5, which lives where
the kernels live, `src/system/approximation/`.

---

## 5. Phases

Each phase is independently reviewable and leaves the tree green.

The shape §4 forces: `src/optim/discrete_systems/` is not touched by any phase — its controlled
predecessor over a folded automaton is already the robust predecessor (§4.4) — and the abstraction
*builder* is not touched until P5. Everything before that is a view, a reader of existing
declarations, and gates. A phase that finds itself editing a solver has misunderstood the design.

House discipline applies to every phase and is not restated per phase: each new file mirrors the
`src/` layout, ships a standalone-runnable test wired into `TEST_FILES`, and every exported symbol
carries a docstring — the docs build errors on missing ones, so this is a gate, not a preference.

### P0 — `FoldedAutomaton` *(the whole mechanism, one file)*
- The view of §4.2: one type against the seven-method `AbstractAutomatonList` interface, plus its
  `symbol -> (u, w)` map. No abstraction, no builder change, no solver change.
- Check whether `transition_cost` can depend on the folded dimension (§4.3's caveat) and record the
  answer where the view is documented.
- **Accept**: on a hand-built two-mode automaton, the *same solver call, unmodified*, returns
  today's winning set on the bare automaton and a strictly smaller one on the view; on an automaton
  with a single mode, bare and view answers are identical. I1 holds trivially — no existing file
  changed except to include the new one.

### P1 — ownership is read once, and gated
- The `environment_input` accessor of §4.8 in `System`, with its methods for the `Noisy…` family
  and `HybridSystem` (dispatching on `switchings`). Composing optimizers call the accessor and
  nothing else; no optimizer ever mentions `isnoisy` or `AutonomousSwitching` by name.
- **Migrate first.** `discreteswitchedsystem` hardcodes `AutonomousSwitching`, so every existing
  switched example is mislabelled relative to how it is solved — 9 files today. Reading the field
  without migrating turns them all into verification queries silently. Prefer fixing this
  *upstream*: a `switching` keyword on `HybridSystems.discreteswitchedsystem` (its maintainer is a
  Dionysos author), with a thin local constructor only as the interim.
- The §3 gates land at the one hand-over point the accessor creates: completeness of the underlying
  alphabet, no deadends, no sampling approximation mode underneath.
- **Accept**: the migrated examples produce today's numbers (I1); an unmigrated construction is
  refused with a message naming the field, never solved under the silently-flipped meaning; and
  `grep -r "AutonomousSwitching\|isnoisy" src/optim/` matches only the accessor's own methods file.

### P2 — the switched case end to end *(the first result on a real abstraction)*
- Nothing new: P0's view over the PCLF quotient P1 already selects. The quotient is built once,
  per mode, exactly as today — the view collapses it at solve time, so the cache (I2) is what makes
  this phase almost empty.
- **Accept**: **Gol–Lazar–Belta Problems 5.1 and 5.2 both reproduced** from one quotient, one
  certificate, one solver call in `gol_lazar_belta_example.jl`, differing only in the `switchings`
  of the system — with the quotient *constructed once* and both questions answered from it. The
  verified set must be a subset of the synthesised one, and the two figures go side by side in the
  paper. If the diff between the two runs is more than that one word, or the quotient is built
  twice, the design has failed its own test.

### P3 — close the soundness gaps P2 exposes
- Deadend states and `atol` coverage loss, both open in the folder README. Either fix them, or have
  any run on a view compute the covered volume and report the uncovered fraction as an explicit
  caveat on the result.
- **Accept**: the winning set's volume plus the reported uncovered volume accounts for the working
  set to within a stated tolerance.

### P4 — counterexamples
- Extract the environment word from the fixed-point iteration through the view's `symbol -> (u, w)`
  map, concretise it, replay it on the concrete system.
- **Accept**: every abstract counterexample either reproduces concretely or is reported as a
  spuriousness witness. Never silently dropped.

### P5 — set-valued `W`: the continuous case *(the only builder change)*

The reference is Reissig, Weber & Rungger, *Feedback Refinement Relations for the Synthesis of
Symbolic Controllers* (IEEE TAC 2017): their growth bound is defined for the **perturbed** system
`ẋ ∈ f(x,u) + [−w̄, w̄]` from the start, and the disturbance enters it as a constant. Our
`ContinuousTimeGrowthBound`
([growth_bound.jl:174](src/system/approximation/growth_bound.jl#L174)) is exactly their
construction with that term dropped — the radius ODE is `ṙ = L(u)·r` where theirs is

```
ṙ = L(u)·r + z,        z = w̄ (additive W); in general zᵢ = sup_{x,w∈W} |fᵢ(x,u,w) − fᵢ(x,u,w_c)|
```

with the nominal trajectory simulated at the centre `w_c` of `W`.

**The consequence for cost is the point: `F(cell, u, W)` is never computed by looping over `w`.**
One nominal integration and one radius integration per `(cell, u)`, exactly as today — `w̄` is an
additive constant in an ODE that is already being integrated, and the per-input hoist
(`input_cache`, [growth_bound.jl:80](src/system/approximation/growth_bound.jl#L80)) is untouched
since the radius integration was already per `(r, u)`. For the growth-bound kernel, **robust
abstraction has the same build cost as nominal abstraction.**

Kernel by kernel:

| kernel | disturbance treatment | effort |
| :--- | :--- | :--- |
| `GROWTH` (continuous) | `+z` in the radius ODE; nominal at `w_c` | one line each |
| `GROWTH` (discrete, `x⁺ = f(x,u) + w`) | `Fr .+ w̄` on the inflated radius | one line |
| `LINEARIZED` | `w̄` absorbed into the second-order `error_map` term | care needed; derive, don't guess |
| `CENTER_/RANDOM_SIMULATION` | already banned for non-empty `W` (§3, hazard 1) | gate only |
| `USER_DEFINED` | contract documented: the map must over-approximate over all `w ∈ W` | docs + gate |

- `W` reaches the kernels as a **set, not a mapping** (§4.1): `w_c` and `w̄` from the set's centre
  and radius, a zonotope or polytope `W` through its box approximation, conservatively. Fed from
  `environment_input` (§4.8), so the kernels never ask who owns the signal — only what its set is.
- Extend `compute_jacobian_bound` in `DionysosSymbolicsExt` to trace `w` as well: interval-bound
  `L` over `x` *and* `w ∈ W`, and compute `z` by interval arithmetic — then non-additive
  disturbances `f(x,u,w)` come for free, with the same proof-not-guess character the Jacobian bound
  already has. This replaces `build_noise`'s warning
  ([alternating_simulation_problem.jl:663](src/optim/continuous_systems/UniformGridAbstraction/alternating_simulation_problem.jl#L663)).
- I1 discipline: the nominal `ContinuousTimeGrowthBound(system; jacobian_bound)` constructor is
  untouched; a noisy system selects a *different closure* at construction time, so no branch and no
  `w̄ = 0` arithmetic appears on the noise-free hot path.
- JuMP front-end vocabulary for the disturbance role (§4.0, open surface question in §7).
- **Accept**: a continuous system with a declared noise set yields a controller valid for every
  disturbance in that set, demonstrated on a benchmark where the noise-free controller fails; the
  abstraction build time with `W` is within noise of the build without it; and I1's benchmark shows
  the noise-free path unchanged.

### P6 — LTL-`∀` hardening *(no new algorithm — P2 proved the direct route)*
- Replace the accepting-set heuristic with an exact good-prefix computation: a monitor state is
  accepting iff *every* extension is accepted, decidable on the DRA with missing edges read as
  dead. Refuse, do not guess, when the check cannot be established.
- State and gate the monitor totality contract at the interface, so a user-supplied monitor
  cannot silently drop runs on a `∀` path.
- Trap tests: a DRA with an absorbing non-accepting state (the heuristic's blind spot); a
  violation that lives entirely in the dead component; a next-step formula (GLB's own
  `R₃ ⇒ X ¬R₁` shape) under both quantifiers.
- **Accept**: on a formula and system where `∃` succeeds and `∀` fails, both answers are produced,
  with a counterexample for the second; and the heuristic's blind-spot DRA is refused or answered
  correctly, never verified wrongly.

---

## 6. Decisions taken

- **The quantifier goes in the system, not the specification.** A `ProblemType` says *what* must
  hold; who owns which signal says *who* chooses. Duplicating the `ProblemType` hierarchy into
  existential and universal variants was considered and rejected: it doubles the catalogue and still
  cannot express the mixed game, which is the case that matters most.
- **The ecosystem's vocabulary wins over ours.** `MathematicalSystems` says who owns a continuous
  signal (`isnoisy` / `noiseset`) and `HybridSystems` says who owns a switching signal
  (`AutonomousSwitching` / `ControlledSwitching`). Dionysos reads both and invents neither. A
  wrapper type and a raw optimizer attribute were each considered and rejected on that ground.
- **`U` and `W` mirror each other in representation and not in the alphabet** (§4.1). Both get a
  mapping and a set; only `U` gets symbols. Full symmetry was considered and rejected on two
  grounds: enumerating a continuous `W` at grid points is unsound, and a `w` component in the symbol
  buys nothing a union does not already give (§4.2).
- **The game is the primitive; the two pure cases are degenerate abstractions, not degenerate
  algorithms.** Implementing verification on its own would be less code now and the wrong shape
  later. Adding an `∃`/`∀` flag would put a branch in every solver for something that is not a
  choice of algorithm at all.
- **No solver changes.** `src/optim/discrete_systems/` is out of scope for every phase: its
  controlled predecessor over a folded automaton already *is* the robust predecessor (§4.4).
  An earlier draft proposed regrouping the solvers' counters; unnecessary, and the measure of
  success is that those files appear in no diff at all.
- **One seam per module** (§4.8): `System` declares (the `environment_input` accessor — ownership
  vocabularies map onto one internal question by dispatch), `Symbolic` provides the mechanism (the
  read-only `FoldedAutomaton` combinator, ignorant of systems and solvers), `Optim` decides (asks
  the accessor, hands automaton or view, gates). Extending to a new ownership vocabulary is one
  accessor method; to a new automaton backend, nothing — the view is generic over the interface.
- **`W` is a set everywhere, never a mapped alphabet.** No `WMapping`: the `Mapping` layer exists
  to make integer alphabets, and `W` must not have one (§4.1) — the layering now enforces what §3
  hazard 3 forbids, instead of merely documenting it.
- **Fold as a view for finite `W`; in the reachable set for set-valued `W`** (§4.2). An earlier
  draft folded at abstraction-build time in both cases; the challenge against I1/I2 killed it —
  build-time folding recomputes the abstraction per question, and for the PCLF quotient it changes
  the quotient itself, making "one abstraction, both answers" false by construction. The view costs
  one small type against an existing interface and preserves both invariants.
- **Switched systems first, continuous disturbances last.** Not because the grid abstraction's
  relation is inadequate — it is not (§3) — but because a finite `W` folds *exactly* and lazily,
  while a continuous one needs the set-valued reachable set that does not exist yet. The ordering is
  by what can be made sound with the least new theory.
- **No result with a non-empty `W` ships without gates on the abstraction.** A missing transition
  manufactures a positive answer, so a state that lost one must be reported unverifiable rather than
  verified.

## 7. Open questions

- Is the bisimulation claim of `BisimulationQuotientProblem` *checked* anywhere, or asserted? Under
  `∀` it carries the whole soundness argument, so it needs at least a numerical check in the spirit
  of `check_pclf` ([pclf.jl:1030](src/utils/pclf.jl#L1030)).
- How should a state whose transitions were dropped be reported? "Unverifiable" is a third answer
  next to verified and refuted, and the fixed points currently have only two. Since §4.3 leaves the
  solvers untouched, this is the one thing that might yet force a change in them — and it is the
  reason to keep the gate at abstraction time if at all possible.
- The design solves `∃u ∀w`, where the controller commits before the environment moves. A controller
  permitted to *observe* `w` first plays `∀w ∃u`, a strictly easier game and a genuinely different
  information pattern; it would need `w` back in the alphabet (§4.2). Out of scope, but worth
  recording as the one thing folding forecloses.
- What is the JuMP surface for declaring a disturbance (§4.0)? `set_role!(model, w, DISTURBANCE)`
  reuses machinery that exists but is imperative and unlike the rest of the front-end, which is
  declarative; `@variable(model, w[1:2] in Disturbance())` matches JuMP idiom but needs a set type;
  a `@disturbance` macro matches `@mode` but adds vocabulary. Decide before P5, not during.
- Should a `Noisy…` system with no `W` samples — a gridded disturbance — be refused outright at the
  front end, given §3 hazard 3? Refusing is honest but makes the front-end path unusable until the
  set-valued reachable-set work lands, so the useful order may be to ship the switched case first
  and gate the continuous one.
- The 6332-versus-9677 cell discrepancy against the predecessor is unexplained
  (`research/BisimulationQuotient/README.md`). It does not block P2, whose comparison is 5.1 against
  5.2 on *our* quotient, but it would have to be settled before claiming a reproduction of their
  Problem 5.2 numbers.

---

## Execution status

Updated 2026-09-01, same session as the plan.

- **P0 — done.** `FoldedAutomaton` + `complete_with_sink` in
  `src/symbolic/automata/folded_automaton.jl`; tests in `test/symbolic/folded_automaton.jl`
  (29 assertions). The reachability testset shows the strict chain on one automaton with the
  unmodified solver: synthesis `{1,3,4}` ⊋ robust `{3,4}` ⊋ verification `{3}`. The
  `transition_cost` question of §4.3 resolved: the uniform-cost path is `w`-independent by
  construction and the general-cost path takes worst case over successors, so folding composes;
  a `w`-dependent stage cost remains future work and is not silently mis-costed (costs are per
  `(q, u)`).
- **P1 — done.** `ST.environment_input` + `ST.with_switching` in `src/system/environment.jl`;
  tests in `test/system/environment.jl` (8 assertions). The pessimistic completion doubles as the
  deadend gate: rather than refusing on empty `(state, mode)` pairs, they are routed to a losing
  sink and counted, which is sound and keeps the run alive.
- **P2 — done.** Flagship numbers: on the predecessor's own certificate and quotient (6 332
  cells), Problem 5.1 wins 4 379 cells (volume 91.88) and Problem 5.2 verifies 1 867 (volume
  51.79), verified ⊊ winning checked in the script. The two runs differ in one declaration. The seam landed exactly where §4.8 said:
  one branch in `OptimizerCoSafeLTLOnQuotient`'s `MOI.optimize!` asks `environment_input` and
  hands the sub-solver either the quotient automaton or the folded completion; translations back
  skip the sink; `solve_concrete_problem` refuses on a folded run with a message that names the
  fix. Migration done where it was owed: the three `common.jl` constructors, `lifted_showcase`'s
  local system, and the PCQ contract test all declare `ControlledSwitching`; that test now also
  runs the same problem under `AutonomousSwitching` and asserts `verified ⊊ synthesised` plus the
  controller refusal (14 assertions). `gol_lazar_belta_example.jl` asks 5.1 and 5.2 from one
  quotient, differing only in the switching declaration.
- **P3 — partially absorbed by P2**: the sink is the coverage gate for empty pairs. The `atol`
  volume accounting remains open.
- **P4 — done.** `verification_counterexample` in the PCQ optimizer: walks the folded product
  along unverified successors to a lasso (or into the pessimistic sink), names the mode behind
  every step from the completed automaton, and replays the word concretely — exact, since the
  quotient is a bisimulation. Tested in the PCQ contract test; demonstrated in the flagship on
  the paper's own initial point.
- **P6 — done** (pulled forward; no new algorithm, as §4.7 records). The absorbing-states
  heuristic in `DionysosSpotExt` is replaced by an exact good-prefix computation — two greatest
  fixed points over the emitted-label alphabet, judged against the DRA's own Rabin pairs — with
  a new two-argument `accepting_states(spec, used_labels)` seam so acceptance is relative to
  what the system can say. Conservative where it cannot certify (refuses; never guesses). Trap
  tests in `test/optim/discrete_systems/cosafe_accepting.jl`: the flagship formula still yields
  its true-sink, and `G(!a)` — which the heuristic would have declared accepting — is refused.
- **P5 — done** for the GROWTH kernel. The perturbed growth bound of RWR'17 lands as a fold at
  construction: a noisy system yields a *nominal* system plus a disturbed radius map
  (`ṙ = L(u)·r + z`, nominal at the centre of `W`), so every downstream type, `discretize` and
  the per-input hoist are untouched and the build cost is the nominal cost. `z` is read off the
  disturbance set in the additive case and must be supplied otherwise — refused, not guessed.
  The grid optimizer gains a `noise_bound` attribute and the §3 gates: CENTER/RANDOM/LINEARIZED
  refuse a declared disturbance by name, a hand-written map must own it explicitly, and
  `build_noise`'s warning is gone because every path now either folds or errors. Acceptance on
  the contraction benchmark: robust winning ⊊ nominal winning from the same solver call, and the
  fold measured exact to `w̄(1 − e^{−τ})`. Remaining within P5's scope: LINEARIZED's perturbed
  error term (derivation), tracing `w` in `compute_jacobian_bound` for non-additive systems, and
  the JuMP disturbance-role surface (§7).
