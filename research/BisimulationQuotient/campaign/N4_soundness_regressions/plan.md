# N4 — Solver soundness: three open anomalies, no regression test

**Run this first.** Nothing else in the campaign matters if the object the paper calls a bisimulation
is not one. A referee who finds a deadend state in it stops reading.

## The three anomalies, all recorded in [`../../README.md`](../../README.md) and none chased

### 1. Deadend states in a bisimulation of a total system

Quotients contain states with minimum outgoing degree 0. The concrete system is total — every state
has a successor under every mode — so a bisimulation of it cannot have a state with no outgoing
transition. Either the quotient is not a bisimulation, or the deadends are an artefact of `atol`
coverage loss, in which case the *statement* of what is built needs qualifying.

**This is the serious one.** It is the only anomaly that, if real, invalidates the central claim.

### 2. `max_slices` does not cap the number of slices

Observed: `max_slices = 12` produced 20 slices. Whatever the attribute does, it is not what its name
says. This matters beyond tidiness because the invariance experiment (E1) explains a ±3 % spread in
covered volume by appealing to `max_slices` truncating at different levels — an explanation resting on
an attribute whose behaviour is not understood.

### 3. The JSR bound was non-monotone over nested conic orders

A certificate family ordered by inclusion must give a non-increasing bound. It did not, and that
artefact caused a retraction (`complexity-vs-memory.md` §5). The current workaround is to set
`tol_feas`/`tol_gap_*` to `1e-6`, under which the bound *is* monotone — but that is a setting, not a
guarantee, and every comparison across template complexities depends on it.

## What exists today

`test/optim/PCLFBisimulationQuotient/` is 778 lines and `test/utils/pclf.jl` another 445; the
statistics API is covered thoroughly, including a unit test for the `deadend_states` *function*. But
**no test asserts a semantic property of a quotient built from a real problem.** Counting deadends
correctly is not the same as there being none.

## MEASURED — the deadend anomaly is a numerical residue, not a soundness bug

Run over the four cached quotients, no rebuild needed:

| quotient | states | slices | deadends | deadend volume / total |
| :--- | --: | --: | --: | --: |
| E2 PCLF | 9 027 | 10 | **3** (0.03 %) | **3.0e-5 / 122.499 = 2.4e-7** |
| E2 induced CLF | 1 000 | 7 | **0** | — |
| E8 PCLF | 10 611 | 7 | **0** | — |
| E8 induced CLF | 2 680 | 6 | **0** | — |

**Three of four quotients are total.** The fourth has three deadends out of 9 027, carrying a combined
volume of $3\times10^{-5}$ against a domain of 122.5 — that is **two ten-millionths of the state
space**. All three sit in the unobserved region (obs 0); two are in inner slices (8 and 9), one in the
outermost, so it is not purely a `max_slices` truncation effect.

**Reading:** these are degenerate slivers left by the `atol = 1e-4` inset — hairline residues that have
no successors because they are numerically empty. The construction is not producing unreachable
behaviour; the arithmetic is leaving dust. That is a much better answer than the README's open
anomaly suggested, and it should be recorded as resolved-with-a-caveat rather than left hanging.

**What still has to be done about it.** A cell of volume $10^{-5}$ with no outgoing transition is
still a state in an object the paper calls a bisimulation. Two acceptable fixes, and one of them must
be chosen and stated:

1. **Prune** cells below a volume threshold tied to `atol` during construction, so the quotient is
   total by construction; or
2. **Document** that the quotient is a bisimulation of the system restricted to the covered set, up to
   the measure the inset discards — which the README already quantifies at ~0.4 % of the domain at
   `atol = 1e-4`, itself far larger than the deadend volume.

Option 2 is honest and costs nothing; option 1 is cleaner and makes the totality assertion of step 1
below unconditional. **Do not leave it as an unexplained anomaly.**

**Anomaly 2 did not reproduce.** `pclf_vs_clf.jl` requests `max_slices = 20` and its quotient has 10
slices; the Gol–Lazar–Belta runs request 50 and have 7 and 6. No cap is violated in any cached
quotient, so the observed 12 → 20 came from a configuration not represented here. Find that
configuration before calling it a bug — it may have been a different attribute or a stale run.

## MEASURED — anomaly 3 is a solver-tolerance artefact at conic order 3, and is now explained

Nested conic orders 1 ⊂ 2 ⊂ 3 give nested templates, so the certified rate must be non-increasing:

| | order 1 | order 2 | order 3 | monotone? |
| :--- | --: | --: | --: | :-- |
| **loosened `tol_* = 1e-6`** | | | | |
| gol_lazar_belta, 1 node | 0.938131 | 0.867381 | 0.859586 | yes |
| gol_lazar_belta, 2 nodes | 0.938131 | 0.867381 | 0.859586 | yes |
| observer_graph, 1 node | 0.913680 | 0.902140 | 0.891094 | yes |
| observer_graph, 2 nodes | 0.913685 | 0.902140 | 0.891089 | yes |
| two_mode, 1 node | 0.718616 | 0.713745 | 0.710852 | yes |
| two_mode, 2 nodes | 0.718616 | 0.713745 | 0.710852 | yes |
| **Clarabel defaults** | | | | |
| gol_lazar_belta, 1 node | 0.938138 | 0.867388 | **1.24874** | **NO** |
| gol_lazar_belta, 2 nodes | 0.938138 | 0.867388 | **0.957466** | **NO** |
| observer_graph, 1 node | 0.913702 | 0.902151 | **0.948782** | **NO** |
| observer_graph, 2 nodes | 0.913702 | 0.902151 | **1.053562** | **NO** |
| two_mode, 1 node | 0.718616 | 0.713751 | 0.710950 | yes |
| two_mode, 2 nodes | 0.718616 | 0.713751 | **Inf** | **NO** |

**Monotone in all six cases at `1e-6`; non-monotone in five of six at the defaults.** The failure is
sharply localized: orders 1 and 2 agree to the fifth decimal under both settings, and **order 3 is
where Clarabel's default `1e-8` feasibility tolerance stops converging**, returning a rate worse than
order 1 or `Inf` outright.

**So anomaly 3 is closed.** It is not a defect in the construction; it is a solver requirement that
the README prescribes (`tol_feas`/`tol_gap_* = 1e-6`) and that had never been demonstrated. Two
consequences:

1. **Any experiment at conic order $\ge 3$ must set the loosened tolerances**, and the paper should say
   so where it reports template complexity. The existing caches are safe — E8 uses order 2 and E2's
   certificate order 1 — so no published number is affected.
2. **This table is the regression test** that an earlier design note asked for, and it would have
   caught the retracted result immediately. Wire it into
   [`test/runtests.jl`](../../../../test/runtests.jl) asserting monotonicity at `1e-6`.

## RESOLVED — the construction is sound; every apparent violation is an escape from $X$

**Final classification of the E2 failures, 1 443 samples:**

| the image $A_m x$ … | count | verdict |
| :--- | --: | :--- |
| leaves the working set $X$ | **12** | legitimate — no successor is owed |
| is in $X$ but outside the covered slices | 0 | — |
| **is in $X$ and in a cell of the target node** | **0** | **no genuine missing transition** |

**Zero soundness violations.** The apparent ~1 % was entirely points escaping the working set, which
the abstraction is not required to follow.

This also explains both things that made it look real. The `atol` trend — smaller `atol` keeps more
thin cells alive near $\partial X$, so a larger fraction of sampled images escape, giving exactly the
observed 0.80 % → 1.21 % → 1.42 %. And why only E2 flagged: its $A_1$ has a row of norm $\approx 10.7$
on a domain of half-width 1.7, while Gol–Lazar–Belta's entries are below 1 on a half-width of 5.9, so
escapes are routine in one and rare in the other. The graph-nondeterminism lead was a coincidence.

**The half that remains unchecked.** This verifies the *simulation* direction: every concrete
transition is covered by an abstract one. A bisimulation also needs the converse — no *spurious*
abstract transitions, i.e. for each recorded $(q, m, q')$ there is some $x \in q$ with
$A_m x \in q'$. That check is the mirror image of this one and equally cheap; add it before claiming
the spot-check validates bisimulation rather than simulation.

## The investigation that got here, kept because the reasoning is reusable

The independent spot-check (sample interior points of a cell, step forward, ask whether **some**
recorded successor contains the image) found violations that survive every benign explanation tried.

| quotient | graph | checked | genuinely missing |
| :--- | :--- | --: | --: |
| E2 PCLF | observer graph, **nondeterministic** | 1 443 | **12 (0.83 %)** |
| E8 PCLF | De Bruijn order 1, deterministic | 1 902 | **0** |
| E8 induced CLF | single node | 5 022 | **0** |

**Three explanations tested and ruled out.**

1. *Boundary / first-match artefact.* The first version of the check located one cell containing the
   image and asked whether that exact cell was recorded. Cells are closed, so a point on a shared
   boundary lies in several. Fixed by testing against **all** recorded successors: 36 → 12, so
   two-thirds were this, and twelve were not.
2. *Floating point.* Violation distances are **0.058 – 0.221** (median 0.208) on a domain of
   half-width 1.7, and pulling sample points 95 % toward the cell centroid leaves the twelve failures
   and their distances **bit-identical**. Not proximity to a facet.
3. *The documented `atol` coverage loss.* Refuted by the sweep — the rate **increases** as `atol`
   shrinks, and coverage loss predicts the opposite:

| `atol` | cells | checked | missing | rate | worst violation |
| --: | --: | --: | --: | --: | --: |
| 1e-3 | 5 488 | 2 259 | 18 | 0.80 % | 1.081 |
| 1e-4 | 9 027 | 2 730 | 33 | 1.21 % | 3.722 |
| 1e-5 | 10 572 | 2 739 | 39 | 1.42 % | 3.252 |

Worst violations reach **3.7**, so the image is nowhere near the recorded successors. At `atol = 1e-5`
the quotient also reports 11 deadend states.

**The lead: graph nondeterminism.** E2's observer graph has node 2 carrying **three edges all labelled
mode 1**, to nodes 1, 3 and 4. E8's De Bruijn graph is deterministic. Only the nondeterministic one
fails. `refine_one_state!` and `group_edges_by_dest_mode`
([`pclf_bisimulation_quotient.jl:550`, `:629`](../../../../src/optim/hybrid_systems/PCLFBisimulationQuotient/pclf_bisimulation_quotient.jl#L550))
are where several same-labelled edges out of one node are handled, and that is where to look first.

**Two checks that would confirm or kill it, in order:**

1. **Is the image inside the covered region at all?** A point whose image leaves the outermost slice
   legitimately has no successor, and E2's $A_1$ has a row of norm $\approx 10.7$, so images escape
   easily. An earlier run reported all twelve images *inside* a cell of the target node, but with
   different sample placement — **re-run that classification against the current samples before
   anything else.** If they are outside the covered set, the construction is fine and only the
   theorem's scope needs stating.
2. **Run the same check on a deterministic incomplete graph** (the dual De Bruijn of N1). If it is
   clean, nondeterminism is implicated; if it fails too, incompleteness is.

*(Resolved above — the classification found all twelve to be escapes from $X$, and the
nondeterminism lead was a coincidence. The hold once placed on E2 results is lifted.)*

## What we will try

1. **A totality assertion.** On every benchmark in the folder, assert that the quotient has no deadend
   states. Where it fails, determine whether the deadends lie inside the working set or only in the
   `atol` inset — that distinction decides whether this is a bug or a documentation fix.
2. **A monotonicity regression test.** Nested conic partition orders 1 ⊂ 2 ⊂ 3 on a fixed system; assert
   the certified rate is non-increasing. An earlier design note already identified this as the first
   thing to add, and it is what would have caught the retracted result immediately.
3. **A `max_slices` contract test.** Either fix the cap or rename the attribute and document what it
   does. Then assert it.
4. **A soundness spot-check independent of the construction.** Sample concrete points, step them
   forward, and verify the abstract transition relation contains the concrete one — a direct check of
   the simulation property that does not trust the algorithm that produced it. This is the strongest
   possible answer to a referee, and it is a loop over samples.
5. **Wire all of it into [`test/runtests.jl`](../../../../test/runtests.jl)**, tagged `:slow` if
   needed, so it cannot regress silently.

## Why this outranks every cost experiment

The cost results say the method is cheaper. The soundness results say the method is *correct*. Only
one of those is load-bearing for the paper, and it is currently the untested one. It is also cheap:
items 1–3 are assertions over quotients that already exist in the caches.

## What would falsify the framing

Finding that the deadends are genuine — that is, inside the working set and not an `atol` artefact.
Then the object is not a bisimulation of a total system and the paper's statement has to change. That
is exactly why this is worth knowing now rather than at review time.
