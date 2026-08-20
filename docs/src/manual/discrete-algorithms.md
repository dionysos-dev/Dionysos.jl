# Discrete synthesis algorithms

Once a system has been abstracted into a finite automaton, the control problem becomes a
graph problem. This page lists **which algorithm each discrete solver runs**, where it lives,
and what it is in the literature.

All of them are in [`src/optim/discrete_systems/`](https://github.com/dionysos-dev/Dionysos.jl/tree/master/src/optim/discrete_systems)
and consume a `Dionysos.Symbolic.AbstractAutomatonList` — so they apply equally to an
abstraction of a continuous system, to a flattened hybrid product, or to an automaton written
by hand.

Throughout, ``n`` is the number of abstract states, ``m`` the number of inputs, and ``E`` the
number of transitions. The **controllable predecessor**

```math
\mathrm{Pre}(Y) \;=\; \{\, q \;\mid\; \exists u,\ \emptyset \neq \mathrm{Post}(q,u) \subseteq Y \,\}
```

is the primitive every one of them is built from: a state is winning if *some* input keeps
*every* possible successor inside the set already won. The universal quantifier over
successors is what makes a nondeterministic abstraction sound — the controller must work
against the worst case the discretization admits.

## Summary

| Specification | Solver | Algorithm | Complexity | Reference |
| :--- | :--- | :--- | :--- | :--- |
| Safety ``\Box \mathcal{S}`` | `OptimizerSafetyProblem` | maximal controlled-invariant set: greatest fixed point of ``\mathrm{Pre}``, with a per-state counter of surviving inputs | ``O(E + nm)`` | [gradel2002automata](@cite), [tabuada2009verification](@cite) |
| Reachability ``\Diamond \mathcal{T}``, unit cost | `OptimizerOptimalControlProblem` | backward attractor by breadth-first layers; a ``(q,u)`` pair is enabled when its counter of not-yet-won successors hits zero | ``O(E + nm)`` | [gradel2002automata](@cite), [rungger2016scots](@cite) |
| Reach-avoid ``\mathcal{S} \,\mathcal{U}\, \mathcal{T}``, general cost | `OptimizerOptimalControlProblem` | Dijkstra on the AND/OR graph — the same counter turns an OR-relaxation into the worst case over successors | ``O((E + nm)\log n)`` | [dijkstra1959note](@cite), [knuth1977generalization](@cite) |
| Reach-avoid with **bounded input variation** ``d(u^-,u) \le \Delta`` | `BoundedInputVariation` | turn-restricted shortest path: Dijkstra on the **line graph**, nodes are ``(q,u)`` pairs | ``O((Em + nm^2)\log(nm))`` | [caldwell1961turn](@cite) |
| Reach-and-stay ``\Diamond\Box\,\mathcal{T}`` | `OptimizerReachAndStayProblem` | nested ``\mu Y.\,\nu Z`` fixed point, one input fixed per newly-won cell | ``O(n(E + nm))`` | [li2020robustly](@cite), [emerson1986efficient](@cite) |
| Co-safe LTL ``\varphi`` | `OptimizerCoSafeLTLProblem` | product with a deterministic monitor, then reachability on the product | ``O(\lvert Q_\varphi \rvert\,(E + nm))`` | [kupferman2001model](@cite), [baier2008principles](@cite), [duret2022spot](@cite) |

Complexities are for the dense ``n \times m`` tables the solvers allocate by default; the
reachability solvers accept `sparse_input = true`, which replaces them by dictionaries when
the input alphabet is large and mostly unused.

## Safety — the greatest fixed point

`compute_largest_invariant_set` computes the largest set ``Z \subseteq \mathcal{S}``
with ``Z \subseteq \mathrm{Pre}(Z)``, by removing states until nothing more can be removed:

```math
Z_0 = \mathcal{S}, \qquad Z_{k+1} = Z_k \cap \mathrm{Pre}(Z_k).
```

The implementation never recomputes ``\mathrm{Pre}`` from scratch. Each state carries a
counter of how many of its inputs are still admissible; when a state is lost, its predecessors
are visited once and the corresponding ``(q,u)`` pairs are disabled, and a state joins the
next removal wave exactly when its counter reaches zero. Every transition is therefore touched
once over the whole computation, which is what buys the linear bound — the standard attractor
computation for safety games [gradel2002automata](@cite).

The controller keeps **every** surviving input at each state, not one, so downstream code is
free to choose among them.

## Reachability and reach-avoid

Two implementations sit behind the same specification, chosen by whether a cost function was
given.

**Unit cost** ([`compute_worst_case_uniform_cost_controller`](@ref)) sweeps backward from the
target in breadth-first layers, and the value function is the layer index — the worst-case
number of steps to the target. The counter plays the same role as in the safety solver, but
counts *unreached* successors: input ``u`` becomes available at ``q`` only once **all** of
``\mathrm{Post}(q,u)`` is won. This is the attractor computation again, and it is why no
priority queue is needed.

**General nonnegative cost** ([`compute_optimal_controller`](@ref)) replaces the layers by a
priority queue. On a deterministic automaton this is textbook Dijkstra
[dijkstra1959note](@cite) run backwards from the target. On a nondeterministic one the counter
returns, and the value of a pair is the cost of the transition plus the **maximum** value over
its successors — Knuth's generalization of Dijkstra to superior functions on hypergraphs
[knuth1977generalization](@cite), which is exactly the min-max dynamic program worst-case
synthesis needs.

Both accept a `safe_set`, which turns ``\Diamond\mathcal{T}`` into the reach-avoid
``\mathcal{S}\,\mathcal{U}\,\mathcal{T}``: unsafe states never enter the controllable set, so
any input risking one never has its counter reach zero and is never selected.

## Bounded input variation — the line graph

Constraining *consecutive* inputs, ``d(u^-,u) \le \Delta``, is not a constraint on states, so
it cannot be expressed by removing them. The value of a state is no longer well defined
either: the cheapest input at ``q`` may be unusable given what was just played, and a locally
more expensive one may be the only way to continue.

[`compute_bounded_input_variation_controller`](@ref) therefore runs the dynamic program on the
**line graph**, whose nodes are ``(q,u)`` pairs and whose edges connect ``(q,u)`` to
``(q',u')`` when ``q \xrightarrow{u} q'`` and ``d(u,u') \le \Delta``. This is the classical
**turn-restricted shortest path**, introduced for road networks with turn penalties by
Caldwell in 1961 [caldwell1961turn](@cite); it is the same device as augmenting the state with
``\Delta u`` in model-predictive control.

The resulting controller is **dynamic**: its memory is the previously played input, and at
each step it plays the compatible input of least value. The optional `target_input` and
`initial_input` constrain the last and first input of a run — for instance forcing velocities
to ramp down before the target is entered.

!!! note "Deterministic automata only"
    The worst-case pair-graph fixed point is not implemented yet, so this solver refuses
    nondeterministic automata rather than returning an unsound controller. The exact-lattice
    abstractions it was written for are deterministic by construction.

## Reach-and-stay — a nested fixed point

``\Diamond\Box\,\mathcal{T}`` is not the composition of a reachability and a safety solve: the
robot may leave the target finitely many times, provided it always comes back. That is a
nested fixed point, the classical shape for such properties [emerson1986efficient](@cite):

```math
\mu Y.\; \mathrm{Pre}(Y) \,\cup\, \nu Z.\bigl(\mathcal{T} \cup Y\bigr) \cap \mathrm{Pre}(Z)
```

The inner greatest fixed point recomputes the invariant inside ``\mathcal{T} \cup Y`` — the
target *together with* what is already winning — which is what lets a target cell be certified
by leaving into won territory and returning.

The implementation follows Algorithm 3 of [li2020robustly](@cite), and keeps the property that
gives that paper its title: the synthesized controller is **memoryless**. That does not come
for free from the fixed point — it comes from *when* the input is chosen. Each cell receives
exactly one input, at the iteration where it first enters the winning set, pointing into the
invariant core known at that moment; later iterations enlarge the winning set but never revise
an input already fixed. A cell's decision therefore depends on the cell alone, not on how many
times the run has been through it, and the controller needs no counter to know which iterate
it is in.

The solver also offers a stricter `stay_on_first_entry` variant: the invariant core is
computed inside the target **alone** and never widened, so the run stays put from the moment
it first arrives. Its winning set is smaller, and target cells outside the core are excluded
from the reachability phase — being in the target while only able to continue by leaving it
would break the very property the variant exists to enforce.

## Co-safe LTL — product with a monitor

A co-safe formula is one every satisfying run settles in finite time
[kupferman2001model](@cite), which is what lets a finite reachability computation decide it.
`OptimizerCoSafeLTLProblem` builds the synchronous product of the abstraction with a
**deterministic monitor** over the atomic propositions, keeping only the product states
reachable from the initial set, and then runs the reachability solver on it with the monitor's
accepting states as target — the standard automata-theoretic construction
[baier2008principles](@cite).

The monitor is an interface, not a fixed type: any object answering `step`, an initial state
and a set of accepting states will do. Loading [Spot](https://spot.lre.epita.fr/) through
`DionysosSpotExt` supplies one automatically by translating an LTL formula
[duret2022spot](@cite).

Because the winning strategy depends on the monitor state and not on the abstract state alone,
the returned controller is **dynamic**. It is tabulated at construction time — the monitor is
evaluated once per (memory state, label) rather than at run time — so the controller remains
plain serializable data even when the monitor itself is a closure.

