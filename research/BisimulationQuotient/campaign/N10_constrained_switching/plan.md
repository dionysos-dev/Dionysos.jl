# N10 — Constrained switching: their matrices, unstable without the constraint

**Promoted to the main paper.** Planar, built on Gol–Lazar–Belta's own matrices with a single scalar
changed, and the only experiment where the predecessor's construction cannot start at all. That makes
it the clearest statement of what memory buys.

**Status:** experiments 1–3 done. Two library defects found and fixed along the way. Three *proof*
obligations remain open; they are theory, not code.

## The claim

> Take their system and scale one mode until it is no longer stable under arbitrary switching. Then
> **no common Lyapunov function exists at any complexity**, so \[gol2014finite\] has no input.
> Restrict the switching to a language on which the system *is* stable, find a path-complete
> certificate over the automaton generating it, and our construction still returns a finite
> bisimulation — and answers their specification on it.

Unlike every cost comparison in this folder, this is a separation **in kind**. Nothing about it
depends on a clock, a facet count, or a search budget.

## Why the predecessor cannot start — the crux, in one paragraph

A common Lyapunov function's condition is $V(A_m x)\le\gamma V(x)$ **for every mode at every state**,
which is literally the arbitrary-switching condition. It cannot express "decrease along admissible
words only", because it has nothing with which to remember which words are admissible. So it is
infeasible exactly when the system is unstable under arbitrary switching, **whether or not the plant
would ever execute those words**.

**The generalization, in one line: let the Lyapunov function depend on the memory state.** On the
lifted space $(x,q)$ with $q$ the automaton state, $V(x,q)=V_q(x)$ depending on $q$ is exactly a PCLF;
one that does not is exactly a common Lyapunov function. The whole gap between the two methods is that
single dependence.

**The tempting non-fix, worth dismissing explicitly.** One might keep a *common* function and enrich
the alphabet instead, treating each admissible transition $(m, q\to q')$ as a new mode. It does not
help: the common-function condition then quantifies over all the new modes at every state and is
infeasible for the same reason. **The memory must be in the Lyapunov function, not in the alphabet.**

## Two scenarios, one construction — they license different claims

The same lifted quotient serves both. What differs is **who owns the switching signal**, which is the
$\exists$/$\forall$ distinction the paper already has.

### A — the restriction is the controller's own commitment

The plant switches freely and is unstable, so no common Lyapunov function exists at any complexity.
The controller voluntarily narrows its authority — "I will never use mode 1 twice in a row" — and in
exchange gets a certificate, a finite bisimulation, and a synthesizable co-safe LTL specification.

**Sound for synthesis only.** Restricting your own choices cannot make an enforceable specification
unenforceable, so you lose completeness, never soundness. Restricting the *adversary* would be
cheating, so $\forall$ claims are not available in this reading.

**What may be claimed:** not the true controllable set, but *an* answer where there was none.

### B — the constraint belongs to the plant

The switching is genuinely restricted by something physical: a scheduling protocol, a supervisory
layer, an actuator interlock, a hardware dwell time. The restriction binds everyone, adversary
included.

**Both $\exists$ and $\forall$ are sound**, because the constrained system *is* the system rather than
a conservative under-model. The graph is problem data, so the $|S|$ copies of the domain are memory
the problem already had. `Donkers11` (networked control as switched systems with protocol-constrained
switching) is a ready-made citation.

**Lead with B** — it needs no caveat about lost completeness. A is the bonus: the same construction
lets a designer *buy* a certificate on an unconstrained plant by accepting conservatism.

## The system

$\Sigma_c = \{c A_1,\ A_2\}$ with $A_1, A_2$ exactly theirs. Measured periodic lower bound over words
of length $\le 10$:

| $c$ | $\rho(cA_1)$ | JSR lb, arbitrary | JSR lb, forbidding $11$ | |
| --: | --: | --: | --: | :-- |
| 1.00 | 0.8558 | 0.8559 | 0.8291 | their published system — stable, used unmodified in exp 1 |
| 1.10 | 0.9414 | 0.9414 | 0.8555 | still stable |
| **1.20** | **1.0270** | **1.0270** | **0.9011** | **used in exp 2** |
| 1.25 | 1.0698 | 1.0698 | 0.9120 | also works; quotient passes 25 000 cells |
| 1.50 | 1.2837 | 1.2837 | 0.9990 | in the window, no margin |
| 1.60 | 1.3693 | 1.3693 | 1.0318 | unstable even on no-$11$ |

The window is $c \in [1.17, 1.5]$. **Settled at $c = 1.20$**: the table prices instability but not the
*quotient*, and a worse constrained rate means more slices. $c=1.25$ gave a quotient past 25 000 cells
for no argumentative gain, since the margin is a **proof, not a measurement** — $\rho(1.20A_1)=1.027>1$
means $1^\omega$ diverges however thin 2.7 % looks.

**The certificate side is settled** ([`certificate_exists.jl`](certificate_exists.jl)): the
constrained language certifies at every template order tried while the unconstrained system certifies
at none, converging toward the exact 2-cycle bound $\rho(A_2\,cA_1)^{1/2}$ as the template enriches.

## Experiments

All three use the same no-$11$ automaton — $q_0 \xrightarrow{1} q_1$, $q_0 \xrightarrow{2} q_0$,
$q_1 \xrightarrow{2} q_0$ — or its dwell-time counterpart. Only the *status* of the constraint changes.

### Experiment 1 — the constraint belongs to the plant · GLB unmodified

[`glb_constraint_moves_both_ways.jl`](glb_constraint_moves_both_ways.jl) ·
[figure](fig_constraint_both_semantics.png)

Their system is stable under arbitrary switching, so their construction works here and E8's numbers
are the baseline. **This experiment validates the constrained machinery against a known answer**
rather than asserting a new one.

The prediction, and it is distinctive: constraining moves the two semantics in *opposite* directions.
$\exists$ shrinks (the controller loses options), $\forall$ grows (the adversary loses options). A
constraint moving both the same way would mean the $\exists/\forall$ machinery is wrong.

**Their method runs here but answers the wrong question.** Against a plant that cannot play `11`, a
quotient of the *unconstrained* system is unsound for synthesis (it credits the controller with words
the plant cannot execute) and needlessly conservative for verification. So scenario B is not merely a
validation: **ignoring a real constraint is unsafe, not just suboptimal.**

#### RESULT — confirmed in both directions

| language | rate | cells | $\exists$ vol | $\exists$ % | $\forall$ vol | $\forall$ % |
| :--- | --: | --: | --: | --: | --: | --: |
| unconstrained | 0.867388 | 10 611 | 186.73 | 44.2 % | 76.09 | 18.0 % |
| no $11$ | 0.830664 | 1 806 | 166.45 | **37.0 %** ↓ | 120.08 | **26.7 %** ↑ |

Volume, not cell count: the two quotients have different cells and the README's cross-certificate
warning applies. The unconstrained arm reproduces E8 on three independent figures — 10 611 cells,
8 794 certified, volume 186.73 against E8's recorded 186.728 — so the constrained arm is read against
a validated instrument.

**The first measurement was wrong, and finding out why produced the library fix below.** As first run,
$\forall$ *collapsed* to 3 cells instead of growing. That was a defect in the verification path:

| language | nodes restrict future? | completions `:language` | `:arbitrary` | $\forall$ `:language` | $\forall$ `:arbitrary` |
| :--- | :--- | --: | --: | --: | --: |
| unconstrained | no | 2 | 2 | 2 962 | 2 962 |
| no $11$ | yes | **0** | **446** | **670** | **3** |

446 is exactly the node-1 cell count. On the complete graph the two semantics agree on every column,
which is the no-op check the fix had to pass. The script now declares
`switching_semantics = :language`, because here the constraint belongs to the plant.

**The $\forall$ direction is a theorem, so this is a regression test, not a discovery.**
$W_\forall(L')\supseteq W_\forall(L)$ for $L'\subseteq L$ is set inclusion — quantifying over fewer
words is a weaker condition. A run reporting otherwise is reporting a bug, which is what happened and
how it was caught. Two other explanations (a coarser partition; differing terminal sets $D$) were
tested and rejected — $D$ differs by 8 %, 76.92 against 70.46, nowhere near enough.

### Experiment 2 — the plant is unconstrained and unstable · GLB with $1.20\,A_1$

[`quotient_on_glb_constrained.jl`](quotient_on_glb_constrained.jl) ·
[figure](fig_glb_constrained.png)

One scalar changed and $\rho(1.20A_1)=1.0270>1$: the constant word $1^\omega$ diverges, so **no common
Lyapunov function exists at any complexity** and their construction has no input at all. One
eigenvalue establishes it; no template refinement or longer search can rescue it. We then impose
no-$11$ *as a design choice*, the controller narrowing its own authority.

**Only $\exists$ is sound here** (scenario A). The certified set is a subset of the true controllable
set — *an* answer where there was none, not *the* answer.

#### RESULT

Certificate rate 0.901091; quotient of **2 138 cells** on 2 nodes, 6 slices, built in 68 s;
**1 517 cells (71.0 %) certified** for their formula. Left panel empty by construction.

**A correction worth keeping.** This script first reported **53 of 2 138**, an artefact of the co-safe
call rather than the certificate. It bypassed `synthesize_cosafe_ltl` and so carried both defects that
experiment 1 had already found: it passed `X` as the co-safe initial set (where the product's initial
*memory* is fixed by the observation there, and `X` spans all three regions), and it never set
`early_stop = false` (whose default seeds the product only from cells meeting the initial set — fatal
here, since the covered region reaches ±17 while `X` is ±10). Suppression was 29×.

**The generalizable lesson: a co-safe result that is merely small looks exactly like a correct one.**
Neither defect raises an error. Any script calling `OptimizerCoSafeLTLOnQuotient` directly needs a
point initial set, `early_stop = false`, and a validated control to read against — or it should go
through `synthesize_cosafe_ltl`, which gets all of this right.

### Experiment 3 — the opposite mechanism: minimum dwell time

[`dwell_time_example.jl`](dwell_time_example.jl) · [figure](fig_dwell_threshold.png) ·
[`plot_dwell_threshold.jl`](plot_dwell_threshold.jl)

Instability has two causes and they need **opposite** constraints:

| cause | fix | example |
| :--- | :--- | :--- |
| one bad mode: $\rho(A_1)>1$ | forbid *repeating* it (max run-length) | GLB at $c=1.20$ — exps 1 and 2 |
| bad alternation: all $\rho(A_i)<1$ but $\rho(A_1A_2)>1$ | force *staying* (min dwell $\tau$) | this experiment, $\tau^\star=3$ |

They do not transfer: a dwell time makes the first case *worse* (it forces you to sit in the diverging
mode) and forbidding repetition makes the second worse (it forces the destabilizing alternation).
"Stable modes, unstable switching" is the famous textbook phenomenon, and omitting it invites the
objection that the method handles only a contrived case.

$A_1 = \begin{psmallmatrix}0.3 & 2\\ 0 & 0.3\end{psmallmatrix}$,
$A_2 = \begin{psmallmatrix}0.3 & 0\\ -2 & 0.3\end{psmallmatrix}$: each mode contracts at 0.3 on its
own, yet $\rho(A_1A_2)^{1/2} = 1.95$, so fast alternation diverges.

#### RESULT — the threshold is discovered, not assumed

| graph | nodes | order 1 | order 2 | order 3 | |
| :--- | --: | --: | --: | --: | :-- |
| unconstrained | 1 | 2.04402 | 1.95393 | 1.95394 | no certificate |
| dwell $\tau=1$ | 2 | 2.04403 | 1.95393 | 1.95392 | constraint vacuous — reproduces the above |
| dwell $\tau=2$ | 4 | 1.14971 | 1.11405 | 1.09234 | still no certificate |
| **dwell $\tau=3$** | 6 | **0.84874** | **0.83708** | **0.82515** | **certifies** |
| dwell $\tau=4$ | 8 | 0.70558 | 0.70014 | 0.69458 | |
| dwell $\tau=5$ | 10 | 0.62271 | 0.61963 | 0.61651 | |

"The minimum dwell time is three steps" is a number the experiment *finds*, and the $\tau=1$ row is a
free sanity check: the constraint is vacuous there and the rate reproduces the unconstrained one to
six digits.

**The quotient also builds on the dwell graph**
([`quotient_on_dwell.jl`](quotient_on_dwell.jl)) — 6 nodes, 8 edges, rate 0.837, **9 452 states**,
13 slices, mean outgoing degree **1.53** with minimum 1. The dwell restriction is structural in the
abstraction, not only in the certificate: nodes $(m,k)$ with $k<\tau$ have exactly one outgoing edge,
so you *must* keep dwelling.

This answers the **code** question for obligation 1 — the construction runs on a graph
path-complete only relative to its own language, and on a bigger one than either headline experiment
uses. It does **not** discharge the theory question. Given that the certificate side also ran happily
while returning a wrong answer, that distinction is worth keeping sharp.

**No figure beyond the threshold sweep.** The shear pushes the certificate's sublevel sets out to
about ±130 against a ±4 working set, so the quotient renders illegibly; experiment 2 is the figure and
this is the table.

**"More complex constraints" needs no experiment.** The graph *is* the automaton, so the constraint may
be **any regular language**: average dwell time, round-robin scheduling, "mode 1 at most twice in any
window of five", a supervisor from a higher planning layer. Bigger automata, same machinery — one
sentence of generality, free.

## Files

| file | role |
| :--- | :--- |
| [`glb_constraint_moves_both_ways.jl`](glb_constraint_moves_both_ways.jl) | **exp 1** — both semantics, both languages; writes `fig_constraint_both_semantics` |
| [`quotient_on_glb_constrained.jl`](quotient_on_glb_constrained.jl) | **exp 2** — unstable plant, their formula; writes `fig_glb_constrained` |
| [`dwell_time_example.jl`](dwell_time_example.jl) | **exp 3** — the $\tau$ sweep and threshold table |
| [`plot_dwell_threshold.jl`](plot_dwell_threshold.jl) | exp 3's figure; writes `fig_dwell_threshold` |
| [`quotient_on_dwell.jl`](quotient_on_dwell.jl) | obligation-1 evidence: the quotient builds on a 6-node dwell graph |
| [`constrained_window.jl`](constrained_window.jl) | probe: the window of $c$ where constrained switching is the only option |
| [`certificate_exists.jl`](certificate_exists.jl) | probe: does a polyhedral PCLF exist on the constraint automaton |

Experiments 1 and 2 cache their quotients (and exp 1 its certified sets) in `*.jld2` beside the
scripts, so re-plotting costs seconds rather than minutes. The caches are gitignored. **Delete them
after changing a graph, a tolerance, `max_slices`, or the quotient struct**, or a stale quotient is
loaded against a fresh certificate.

## Library work this produced

Two genuine defects, both found because a constrained language exercises paths a De Bruijn graph never
reaches. Both are fixed and regression-tested.

1. **The PCLF bisection bracket was language-blind.** All three constructors took the lower bracket as
   $\max_m \rho(A_m)$ over *every* mode, valid under path-completeness over $\langle M\rangle^\dagger$
   where any mode may repeat for ever, and invalid for a graph forbidding repetition — so the search
   could never return anything below it and converged to the bracket floor. The fix counts only modes
   carrying a **self-loop**. Its signature was *a value that does not change when the template is
   enriched*, which this folder has now seen twice.

2. **The pessimistic completion charged the adversary for forbidden modes.** Full analysis, all six
   findings and what shipped: [`../../notes/constrained-language-review.md`](../../notes/constrained-language-review.md).
   Synthesis was never affected; only the $\forall$ path was.

## Open — the three proof obligations

Properties of the *theory*, not the experiment.

1. **Does the correctness argument need path-completeness over $\langle M\rangle^\dagger$, or only
   over $\mathcal L$?** The preliminaries already define stability over a language
   ([`../../JCVD/main.tex`](../../JCVD/main.tex) line 493) and the conclusion lists constrained
   switching as future work, so the formalism is in place. Re-read `sec:bisimulation_construction`
   wherever "for every word there is a path" is used. Empirically answered by exp 3; theoretically open.
2. **Verification needs $\mathcal L(\mathcal G)=\mathcal L$ exactly, not $\subseteq$.** If the graph
   under-generates, adversary moves are missing and $\models_\forall$ does not transfer. Synthesis is
   safe under $\subseteq$ by `thm:synthesis-invariance`. This is a *selling point* — the design rule
   the paper already owns, applied to a new setting. **Now enforced in code** as
   `switching_semantics`, which is exactly this distinction made explicit.
3. **The terminal structure is stated for complete / co-complete graphs.** `th:induced_common` and the
   $\mathcal D_1 \subseteq \mathcal D_2$ construction lean on completeness for forward invariance, and
   a constraint automaton is in general neither. `rem:general-terminal-structure` says the framework
   needs only *some* compact $\mathcal D_1 \subseteq \mathcal D_2$ with the stated property, so the
   generalization is probably available — but this is the genuine obligation.

Also open, from the review: giving `LabDigraph` an initial-node set so a certified set is projected
onto the node a run may start in rather than unioned over all nodes. Needed before any $\forall$ result
on a future-restricting graph can be quoted without a caveat.

## What would falsify it

Obligation 1 requiring path-completeness over all words, or obligation 3 not generalizing. Either
would mean the framework does not yet cover constrained switching, in which case this becomes the
extension paper after all — and the finding is worth recording either way.
