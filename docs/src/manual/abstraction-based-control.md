# Abstraction-based control

## The problem

Designing a controller for a nonlinear system, with constraints, and *proving* that the closed loop
meets its specification is hard. Classical control gives excellent tools for stability and tracking,
but a specification like "reach this region within 20 seconds while never entering that one, and stay
there afterwards" is not naturally expressed as a gain to tune. In practice such controllers are
hand-crafted by experts and validated by simulation — which shows the absence of failure on the
trajectories tried, not on the ones that were not.

**Abstraction-based control** takes a different route: replace the system by a finite object, solve
the problem exactly on that object with graph algorithms, and carry the solution back with a
guarantee attached.

## The three steps

![Abstraction-based control.](https://github.com/dionysos-dev/Dionysos.jl/blob/master/docs/assets/abstraction.png?raw=true)

1. **Abstraction.** The state space is partitioned into **cells**, and each cell becomes one state of
   a finite automaton. For every cell and every quantized input, the automaton gets a transition to
   *every* cell the real system could land in. The specification is transposed the same way.
2. **Synthesis.** The abstract problem is a game on a finite graph, and finite graphs are something
   computers are very good at: Dijkstra, A\*, and fixed-point iterations over predecessor sets solve
   it exactly, with no local minima and no tuning.
3. **Concretization.** The abstract controller — a map from cells to quantized inputs — is turned
   back into a controller for the original system: measure the state, find its cell, apply the input
   the abstract controller assigns to it.

Throughout, everything comes in pairs, and the vocabulary is used consistently across the toolbox:

| Concrete | Abstract |
| :--- | :--- |
| the real system ``\dot x = f(x,u)`` | a finite automaton (the **symbolic model**) |
| a state ``x \in \mathbb{R}^n`` | the **cell** containing it |
| an input ``u \in \mathcal{U}`` | one of finitely many quantized inputs |
| the specification | its transposition onto cells |
| the controller you deploy | the map from cells to inputs |

## What you get: soundness

The transitions are built by **over-approximation**: a transition from cell ``c`` under input ``u``
exists whenever the real system *could* move from somewhere in ``c`` into the target cell. The
abstraction therefore admits at least every behaviour the real system has, and usually more.

That is what makes the method **correct-by-construction**. A controller that wins the abstract game
wins against a *more adversarial* opponent than reality, so it also works on the concrete system —
and this is a theorem, not a simulation campaign. The relation that formalises it is an *alternating
simulation*; when it holds in both directions the two systems are *bisimilar*.

## What it costs

Soundness is not free, and the two prices are what most of the research in Dionysos is about.

**The curse of dimensionality.** The number of cells grows exponentially with the state dimension.
Halving the grid step of a 3-D system multiplies the state count by 8 — and the transition count by
more.

**Spurious non-determinism.** Over-approximation adds behaviours the real system does not have. Push
it too far and the abstract game becomes unwinnable even though the concrete problem is perfectly
solvable. This is why a failure is reported as `LOCALLY_INFEASIBLE` rather than `INFEASIBLE`: it says
*this abstraction* admits no controller, never that no controller exists.

The two pull against each other. A finer grid means less spurious non-determinism and a better chance
of success, at exponentially more work; a coarser grid is cheap and may prove nothing. Choosing the
discretization is therefore the central modelling decision, not an implementation detail — which is
why it is written explicitly in every Dionysos model.

## The research direction: smart, lazy abstractions

The uniform grid is the obvious construction, and the wasteful one: it spends the same effort
everywhere, including in regions the controller will never visit. **Smart abstractions** design the
partition instead of fixing it a priori — using optimization-based cell shapes, or refining only
where the coarse abstraction turned out to be too coarse.

**Lazy solvers** go further and co-design the abstraction with the controller, postponing the
expensive numerical work and computing only the fragment of the abstraction that synthesis actually
needs. Both are first-class solver families in [Overview](@ref).
