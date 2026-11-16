# N6 — The separation does not need switching: one mode, and the node is a clock

**Measured, and it changes what N5 is a statement about.** The rotation separation survives the
reduction to a *single* linear system. There is no switching, no adversary, no control input — only a
discrete-time linear map — and the single-node construction still certifies nothing where the cycle
certifies optimally.

## The result

$x_{k+1} = A x$ with $A = \rho R(2\pi/q)$. Four facets per node throughout. Single-node rate is the
best over 41 template orientations including the identity; the cycle rate is analytic and checked by
sampling.

| $q$ | $\rho$ | 1 node | $q$-cycle | verified |
| --: | --: | --: | --: | --: |
| 5 | 0.85 / 0.90 / 0.95 | 1.07089 / 1.13389 / 1.19689 | 0.85 / 0.90 / 0.95 | exact |
| 6 | 0.85 / 0.90 / 0.95 | 1.16097 / 1.22926 / 1.29755 | 0.85 / 0.90 / 0.95 | exact |
| 7 | 0.85 / 0.90 / 0.95 | 1.19438 / 1.26464 / 1.33490 | 0.85 / 0.90 / 0.95 | exact |
| 8 | 0.85 / 0.90 / 0.95 | 1.20195 / 1.27265 / 1.34336 | 0.85 / 0.90 / 0.95 | exact |

**The single node is above 1 in every case** — no certificate, no abstraction — while the cycle
attains $\rho$. The ratio single/$\rho$ reproduces the analytic gauge constants exactly: 1.2599 at
$q=5$ (the 1.260 quoted in `memory_vs_geometry.jl`), saturating at $1.4141 \approx \sqrt2$ by $q=8$,
which is a square rotated by $45°$.

**Separation threshold:** the single node certifies iff $\rho < 1/\mathrm{gauge}(q)$, that is
$\rho < 0.794$ at $q=5$ and $\rho < 0.707$ at $q=8$.

## RESULTS — a feasibility win and a cost loss, and both are informative

### The quotient builds ([`2_quotient_on_clock.jl`](2_quotient_on_clock.jl))

$q=5$, $\rho=0.9$, working set $[-6,6]^2$, two regions, formula $(\lnot R_2\,U\,R_1)$.

| | |
| :--- | :--- |
| GLB single node, 4 facets, best of 41 orientations | **1.13389 ≥ 1 — no certificate, nothing to build** |
| 5-cycle certificate | 0.9 exactly (analytic) |
| quotient | 1 694 cells, 17 slices, 28.9 s |
| cells per clock phase | 334 · 333 · 336 · 348 · 343 |
| $\exists$ and $\forall$ for $(\lnot R_2\,U\,R_1)$ | 787 cells, volume 172.499 — **identical sets** |

Three internal checks pass. $\exists \equiv \forall$ **as sets**, which only a one-mode system can
verify since neither player has a decision; cells split evenly across phases, as the rotational
symmetry demands; and completions = 10 = the deadend count, so every completion is a genuine escape
and none is language-forbidden, consistent with `restricts_future = false` for a one-letter alphabet.

### Route 1 versus route 2 ([`3_route1_vs_route2.jl`](3_route1_vs_route2.jl))

**Route 1 wins here, and the reason is structural.** Both certify at exactly $\rho$, as
`thm:graph-invariance` requires, so the guarantee is identical and only the geometry differs.

| arm | cells | slices | Σ faces | max faces | mean faces/cell | build |
| :--- | --: | --: | --: | --: | --: | --: |
| route 2 — 5 nodes × 4 facets | 1 694 | 17 | 20 817 | 80 | 12.29 | 3.00 s |
| route 1 — 1 node × 20 facets | **285** | 15 | **10 184** | 298 | 35.73 | **2.28 s** |

Route 2's cells *are* simpler, by 2.91× in mean facet count, exactly as the mechanism predicts. It
still loses, because it has **5.94× more of them** and the lifting penalty swamps the facet saving.

**Why, in one line: the $q$-cycle is a COMPLETE graph over its one-letter alphabet**
(`PCLF.is_complete` returns true — every node has an outgoing edge for the only mode). So
`prop:complete-redundant` applies and lifting yields $|S|$ near-isomorphic copies of the same
partition: $1\,694 \approx 5 \times 340$ against route 1's 285. The node records the past and
constrains nothing about the future, which is exactly the case where lifting is pure redundancy.

**The mechanism, and every measurement of it, is centralised in
[N1's plan](../N1_quotient_cost/plan.md)** — when route 2 wins, why, the per-node decomposition,
and the region-free case. N6 is one row of that table, and the one where the answer is negative.

### Two corrections this produced

- **"Route 1 cannot start" is false as stated.** It is true of *GLB's single node at a fixed
  four-facet template*, which is a statement about a template budget. The **induced common**
  ($V^\star(x)=\max_s\|R^{-s}x\|_\infty$, from `build_common_lyapunov`) starts fine and certifies at
  the same $\rho$ — with $4q$ facets in one node instead of 4 in each of $q$. Keep the two baselines
  distinct: the first is the feasibility leg, the second the cost leg.
- **The facet budgets match only at the TEMPLATE.** $5\times4$ against $1\times20$ is 20 either way,
  but the resulting *quotients* are 20 817 against 10 184 faces — route 2 pays 2.04× more. Total facet
  budget is not conserved through the construction, as N1 also found.

### What N6 is for, then

**The feasibility leg, not the cost leg.** Its value is that route 1's predecessor has no input at all
on a single linear map — the most accessible setting there is — and that is untouched by the cost
result. Quoting N6 as a speed result would be wrong, and a referee checking it would find route 1
faster.

## Why it works, and what the node means

The $q$-node cycle $s \to s{+}1$ over the one-letter alphabet is path-complete for $1^*$. Node $s$
carries $V_s(x) = \|R^{-s}x\|_\infty$, and

$$V_{s+1}(Ax) = \|R^{-(s+1)}\rho R x\|_\infty = \rho\|R^{-s}x\|_\infty = \rho V_s(x).$$

A one-step polytopic Lyapunov function needs $RB \subseteq (\gamma/\rho)B$; as $\gamma\to\rho$ the
unit ball must be invariant under the order-$q$ rotation group, and a polytope then needs $\ge q$
facets. **This is N5's orbit argument with the second mode deleted** — it never used it.

**With $M=1$ the graph node is a clock.** The memory is the phase $k \bmod q$, so the method becomes a
*periodic* (finite-step) Lyapunov function abstraction: the certificate decreases every $q$ steps
rather than every step, while the abstraction still respects per-step transitions and per-step
observations.

## What the problem becomes

Pure **verification** of an autonomous linear system against scLTL over polyhedral regions — there is
no switching signal to choose, so $\exists$ and $\forall$ coincide. This is not a degenerate problem:
questions about the orbit of a linear map relative to polyhedral sets are Skolem-flavoured, so a
finite bisimulation is a genuine object.

## The risk, and the defence

**The risk.** A referee may say that for one mode this is the classical finite-step / periodic
Lyapunov function idea (`Lazar_2011`, "On Infinity Norms as Lyapunov Functions: Alternative Necessary
and Sufficient Conditions"), and that using $A^q$ instead of $A$ is well known.

**The defence, which must be stated explicitly.** The contribution is the *abstraction*, not the
certificate. A periodic Lyapunov function is classical; using one to build a finite bisimulation that
respects per-step observations is not — and you cannot get it by lumping $q$ steps into one
transition, because the specification is checked at every step. The clock also lands in the abstract
state for free, so the resulting controller/monitor carries it without extra machinery.

**Decide deliberately how to use this.** Two options, and they are different papers:

- **As the opening example.** A single linear map is the most accessible possible setting, and the
  separation is exact and provable. It makes the whole thesis legible in half a page.
- **As a remark after the switched result.** Keeps switching as the headline and uses $M=1$ to show
  the mechanism is more general than it looks.

The first is more striking; the second is safer against the "that's just finite-step Lyapunov"
objection. Reading the `Lazar_2011` line of work before choosing is cheap and decides it.

## The better-motivated relative: periodic linear time-varying systems

$x_{k+1} = A_{k \bmod p}\,x$. Here the graph is **not a design choice — it is the period**, so the
memory is the problem's own, exactly as in [N2](../N2_constrained_switching/plan.md).

**The sharp case: every $A_i$ unstable while the monodromy $A_{p-1}\cdots A_0$ is stable.** Then no
common Lyapunov function exists at any complexity, since each mode alone diverges, so the predecessor
cannot start — while the $p$-cycle PCLF exists by construction. It is N2's argument with the most
constrained language possible (a single periodic word), in a far more classical application area:
sampled-data with periodic scheduling, rotating machinery, periodically time-varying plants.

**This is probably the most practically compelling non-switched configuration**, and it is easy to
instantiate: pick any $A_0, A_1$ with $\rho(A_i) > 1$ and $\rho(A_1A_0) < 1$.

## A softer relative: non-normal single systems

$\rho(A) < 1$ with large transient growth. A one-step polyhedral Lyapunov function exists but needs
many facets to wrap the transient; a clock-indexed $q$-step certificate needs fewer. Less clean than
the rotation argument, closer to real plants. Worth one row in a table, not its own experiment.

## What does *not* work, recorded so it is not attempted

**State-dependent switching / PWA**, $x_{k+1} = A_ix$ for $x \in X_i$. The framework treats modes as
free — chosen by the controller or an adversary — whereas here the mode is a function of the state.
That is a restriction by *state*, not by language, and the construction does not support it. A
plausible extension (Dionysos already has PWA machinery in
[`bemporad_morari.jl`](../../../../src/optim/bemporad_morari.jl)), but not a configuration where the
method helps today.

## Next steps

1. **Build the quotient**, not just the certificate — everything above is the Lyapunov side. Take
   $q = 5$, $\rho = 0.9$, a couple of polyhedral regions, and an scLTL formula, and confirm the
   construction terminates and the single-node arm has nothing to compare against.
2. **Instantiate the periodic LTV case** with unstable modes and a stable monodromy.
3. **Read `Lazar_2011`** and decide between the two framings above.
