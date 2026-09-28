# N1 — Building the bisimulation on the PCLF directly, instead of determinising it first

Working notes for the experiments section. Every number is either usable or marked **void** with the
reason. The one configuration issue that voids a row is a terminal set `D` overlapping an observation
region; §4b says exactly which rows that affects and which it does not. Nothing else is in the file.

## 1. The comparison

Given a path-complete Lyapunov function (PCLF) from an oracle, there are two ways to a bisimulation
quotient:

| | |
| :--- | :--- |
| **route 1** — the predecessor | induce a common Lyapunov function (`build_common_lyapunov`, the observer/subset construction), then run the single-node construction on it |
| **route 2** — ours | run the construction on the PCLF directly, one partition per node |

**Why the comparison is meaningful: graph invariance.** A PCLF and the common induced from it certify
the *same* contraction rate, and the bisimulations built from either certify the same property. The
two routes therefore differ only in what they compute on the way. Every script asserts the rates agree
and aborts otherwise.

**Each graph is its own experiment.** Route 1's input is whatever PCLF it is handed, so the common
induced from the *dual* PCLF and the one from the *primal* PCLF are different certificates — different
observer constructions, different geometry, different rates (0.717572 against 0.715897 on example A).
Every row below is a **paired** comparison: route 2 against the baseline *that graph* produces. Never
divide a route-1 number from one row by a route-2 number from another.

---

## 2. The mechanism — why the two graph orientations behave oppositely

With $M = 2$ at De Bruijn order 1 both graphs have two nodes and out-degree 2. What differs is what a
node *means*:

| | node $(m)$ means | outgoing (mode, dest) | |
| :--- | :--- | :--- | :--- |
| **primal** | "mode $m$ was just played" — the **past** | $(1,(1))$, $(2,(2))$ | two **different** modes |
| **dual** | "mode $m$ is what I play next" — the **future** | $(1,(1))$, $(1,(2))$ | **the same mode**, twice |

The primal graph is **complete** and deterministic; the dual is **incomplete** and non-deterministic.
That single change acts on the two routes through two separate channels.

### 2a. What determinisation costs — the dominant channel

Writing the edge condition as $V_d(A_mx)\le\gamma V_s(x)$ for every edge $(s,m,d)$:

- **complete** (every node an *outgoing* edge for every mode) $\Rightarrow$ $V_{\min}=\min_i V_i$ is a
  common Lyapunov function. If the minimum at $x$ sits at $i$, completeness supplies an edge
  $(i,m,d)$, so $\min_j V_j(A_mx)\le V_d(A_mx)\le\gamma V_i(x)$. Its sublevel set is
  $\bigcup_i\{V_i\le\gamma\}$ — the **union**: non-convex, but the larger region.
- **co-complete** (every node an *incoming* edge for every mode) $\Rightarrow$ $V_{\max}=\max_i V_i$
  is one. For each target $i$, co-completeness supplies $(s,m,i)$, so
  $V_i(A_mx)\le\gamma V_s(x)\le\gamma V_{\max}(x)$. Its sublevel set is
  $\bigcap_i\{V_i\le\gamma\}$ — the **intersection**: convex, but the smaller region.

**The two De Bruijn graphs sit at opposite extremes and each admits exactly one.** The primal records
the last mode: every node emits every mode but is entered by only one. The dual commits to the next:
every node emits one mode but is entered by all.

| | complete | co-complete | route 1's induced common |
| :--- | :--: | :--: | :--- |
| **primal** | ✓ | ✗ | $V_{\min}$ — the **union**, which fragments as the pieces separate |
| **dual** | ✗ | ✓ | $V_{\max}$ — the **intersection**, convex whatever the pieces |

**Determinising a dual graph costs coverage; determinising a primal one costs convexity.** Route 2
pays neither, because it never determinises.

Two consequences used throughout:

1. **On a co-complete graph a PCLF can never beat a convex common function in *rate*** — $V_{\max}$
   is convex and certifies the PCLF's own rate. The dual graph's advantage is therefore purely
   computational, which is what makes it the clean instance for an algorithmic claim.
2. **Piece diversity helps route 2 on the primal graph and hurts it on the dual.** On the primal it
   fragments route 1's union; on the dual route 1's common stays one convex set however different the
   pieces are, so diversity costs route 2 and costs route 1 nothing.

### 2b. What refinement costs — the second channel

**Redundancy of complete graphs.** If every node has an outgoing edge for every mode, lifting produces
$|S|$ near-isomorphic copies of one partition: the node records which mode was last played and
constrains nothing about what follows, so it must still serve the whole alphabet and does exactly the
single node's work. Route 2 then pays $|S|$ copies for nothing. Measured slice by slice on example A's
primal graph, the per-node counts are route 1's **exactly** — 1, 1, 4, 23, 116 against 1, 1, 4, 23, 115.

**Partition growth is a product.** Refining node $s$ superimposes, inside the new annulus, the
pre-images of its successors' cells — one family per outgoing edge. Superimposing families multiplies,
so $p_{i+1}^{(s)} \lesssim \prod_{(s,m,d)\in \mathrm{Out}(s)} p_i^{(d)}$, not a sum. The number of
modes a node must serve therefore matters superlinearly.

### 2c. The two cost laws, and the predictor

| | what its cost tracks |
| :--- | :--- |
| **route 2** | the **cell count** — its classes stay simple |
| **route 1** | the **fragmentation of the induced common**, not how many cells it makes |

On a **complete** graph route 2 builds *more* cells and wins anyway, on cell simplicity. On a
**co-complete** graph there is no fragmentation to avoid, so the only saving is the per-node
reduction

$$g \;=\; \frac{\text{route 1's cells}}{\text{route 2's cells per node}}$$

and route 2 pays $|S|$ copies of the pipeline, so **$g > |S|$ is break-even**. This is the one
quantity that predicts the *sign* of the gain before anything is built.

---

## 3. The benchmark set

Three examples, each placed where the competing explanation is structurally impossible. Piece
diversity is the largest relative gap between two node pieces' support functions over 64 directions —
exact for convex sets, no volume needed, reads as a percentage.

| | `A_identical_pieces` | `B_diverse_pieces` | `C_gol_lazar_belta` |
| :--- | :--- | :--- | :--- |
| system | `two_mode_problem` | the observer-study pair | Example 3.1 of arXiv:1208.5471 |
| certificate | shared 4-facet rotated template | shared conic partition, order 2 | shared conic partition, order 2 |
| piece gap, primal | 0.0278 | 0.4599 | **0.1778** |
| piece gap, dual | < 1e-4 | — | **0.1068** |
| induced common, primal | 3 parts | 13 parts | **17 parts** |
| induced common, dual | 1 part | — | 1 part |
| graphs | dual **and** primal | primal only | dual **and** primal |
| regions | two | none | three, and their co-safe formula |

**A is the control.** On its dual arm the two node pieces are the same set, so both routes carry the
*same function* and the entire difference is the graph. A functional gain is impossible there anyway
(§2a).

**B is primal-only on purpose** — the dual graph is where diversity works against route 2 (§2a).

**C is the external benchmark.** Their dynamics, their working set $\{\lVert Lx\rVert_\infty \le 10\}$,
their three observation regions, their co-safe formula and their initial point. What is *not* theirs
is the certificate: instead of their common Lyapunov function, an oracle supplies a path-complete one
and the two routes are compared on it. It is also the hardest case for the claim (§4).

---

## 4. Results

### 4a. Measured and usable

| example | graph | case | route 1 / route 2 cells | **time ratio** |
| :--- | :--- | :--- | --: | --: |
| A | dual | no regions | 2.27 | **2.59×** |
| A | primal | no regions | 0.57 | **2.63×** |
| B | primal | no regions | 0.81 | **7.61×** |
| **C (GLB)** | **primal** | **with regions** | **0.65** | **7.01×** |

The first three are `mid_D` (7 rungs, region-free, ladder pinned). The GLB row is their own working
set at 7 rungs with `D` fixed by the construction's rule: route 1 6 613 cells / 1 276 s against route
2 10 158 cells / 182 s. Route 2's 10 158 cells match the 10 611 of the published reference run, which
is the check that this is the real example.

Three things to read off it:

**Route 2 loses the cell count on every primal row (0.57–0.81) and wins the clock on all of them.**
That is §2b's redundancy plus §2a's fragmentation, and it is the clearest evidence in this folder that
**cell count is the wrong cost proxy**.

**The margin tracks fragmentation.** Induced common of 3 parts → 2.63×; 13 parts → 7.61×; 17 parts →
7.01×. Worst cell, route 1 against route 2: 399 parts against 42 on GLB, 128 against 16 on B.

**The cell ratio is invariant in problem size** — 2.26 then 2.27 on A's dual across quotients sixteen
times larger. That stability is worth more than any single margin.

### 4b. The dual graph, and the $g > \lvert S\rvert$ threshold

| example | graph | rungs used / needed | $g$ | time ratio | reading |
| :--- | :--- | :--- | --: | --: | :--- |
| A | dual | 7 / 5 — **valid** | **4.81** | 2.93× | $g > 2$ — route 2 wins |
| C (GLB) | dual | 9 / 13 — **void** | 1.55 | 0.72× | $g < 2$ — route 2 loses |

The two rows sit on opposite sides of the threshold and the outcome follows it in both. That is the
best support the predictor has, and it rests on **two points**, one of which still has to be re-run.

**On which rows are valid.** Forcing the rung count bypasses the construction's stopping rule — the
first level whose terminal set no longer meets an observation region — and a `D` that still meets a
region labels every point of the overlap with the terminal observation instead of its own. But the
clearing property is **monotone in the rung count**: more rungs means a smaller terminal set, so once
it clears it stays clear. Verified on example A, where it first holds at $k = 5$ and holds at every
$k \ge 5$. Forcing *more* rungs than needed is therefore harmless; forcing *fewer* is what breaks.

Example A was run at 7 rungs and needs 5, so **every example-A row in this file stands**, including
the "with regions" ones. GLB was run at 5 rungs where it needs 7 (primal) and at 9 where it needs 13
(dual), so those two are void. The GLB dual row has gone against route 2 twice (0.47× unwarmed, 0.72×
warmed), so the loss is likely real — but it is not established until it is re-run at 13.

### 4c. Soundness — the arms certify the same region

Every cost number is void unless both routes answer the specification identically. Cell counts cannot
show this (90 cells and 44 can cover the same area), so the check is on **volume**, on example A:

| | route 1 | route 2 | raw | normalised by each arm's covered area |
| :--- | --: | --: | --: | --: |
| dual | 2.1522 | 2.2011 | 2.27 % apart | 0.05500 vs 0.05497 — **0.05 %** |
| primal | 2.0857 | 2.2029 | 5.62 % apart | 0.05392 vs 0.05458 — **1.22 %** |

The raw primal gap exceeds 5 % for a measurable reason rather than a soundness one: the arms do not
tile the same region (§2a) and `atol` insets every cut, so the arm making more cuts loses more area.
Route 1 averages 9.12 parts per cell against 3.85 and ends with 4 % less covered area — very nearly
the whole 5.6 %.

On GLB both routes are additionally checked by solving the co-safe formula under **both quantifiers**
(∃ synthesis, ∀ verification) and comparing the certified sets. A cheaper quotient is only worth
having if it answers both questions the same way.

---

## 5. Protocol

**What is held fixed within a row.** Same system, same regions, same PCLF, same tolerance
(`atol = 1e-3`), same ladder. The ladder matters most: `ΓX` (the outer level) is computed once from
**route 1's** certificate and imposed on route 2, so the baseline keeps the geometry it would have
chosen and route 2 is made to conform. Taking it from route 2 would flatter route 2. Without pinning,
the two arms start from different $\tau_X$, tile different regions, and the counts measure coverage
rather than efficiency.

**The rung count is derived, not chosen**, from the construction's own stopping rule: the first level
whose terminal set no longer meets any observation region. This is also how GLB's $\Gamma_D = 5.063$ is
fixed. It differs between the two orientations — 7 on GLB's primal against 13 on its dual — because
the dual contracts more slowly ($\gamma = 0.924$ against $0.867$) and needs more levels to reach a
region-free terminal set. **Matching rung counts across orientations is wrong**; matching depth is
what is required.

**Timing.** Route 1's determinisation is inside its timed closure — handed a PCLF it must produce the
common before it can build anything. (Measured at **0.0 %** of route 1's time in every configuration:
the cost is not the determinisation step but what it returns.) Both arms run **back-to-back inside one
round**, so the ratio is paired; the median over rounds is reported with its spread. A warm-up call
precedes all rounds and every figure is drawn strictly after the last measurement.

**Three rules, all learned the hard way.**

1. **Only same-run ratios.** Two runs of identical code produced the same cell counts to the unit and
   absolute times differing by 3.6×–4.0×. Absolute timings from this folder are meaningless in
   isolation.
2. **A time ratio is quotable only from `mid_D` or larger.** Disjoint spreads are necessary but not
   sufficient: example B at `large_D` was measured twice, both with disjoint spreads, at 8.12× and
   4.18×. At builds of 0.3–2.7 s the machine state dominates.
3. **Cell, facet and part ratios are exempt** — they are deterministic and reproduce to the unit.

**One open measurement anomaly.** On the GLB primal row the *same* build took 666 s in the warm-up and
1 276 s in the timed round, producing identical output (6 613 cells) both times; route 2 moved only 8 %
(198 s → 182 s). The warm-up should be the slower of the two, since it compiles. Until this is
understood the 7.01× should be read as the interval **~3.7–7×**. Resolving it is the first thing to do
before quoting any single-round number.

---

## 6. Threats to validity, and what is *not* claimed

**No functional gain.** No instance here certifies a rate a common function cannot. The reason is
structural, not a failure to search: **polyhedral common Lyapunov functions are universal for the
JSR**, so with rich enough polyhedral templates a common function always matches — measured rates are
equal to four decimals. The classical separation lives in **restricted** classes, quadratics above
all, and that route is blocked in the library: `get_sublevel_set(::ObserverCLFPiece, …)` asserts
`Pi isa PolyhedralPiece`, so route 1's induced common cannot be built from ellipsoidal pieces. Lifting
that assertion is the prerequisite, and it is library work.

**Route 2 does not win everywhere.** On GLB's dual orientation it lost in every measurement so far.
The claim that survives all configurations is the conditional one: route 2 wins on a complete graph
through cell simplicity, and on a co-complete graph only when $g > |S|$.

**"Incomplete and mode-committed" is not sufficient.** A mode-committed graph sharing the dual's
properties — incomplete, identical pieces — has been measured *increasing* the quotient (305 against
93). The per-node gain must exceed $|S|$ and nothing guarantees it does.

**The baseline was corrected against us.** `build_common_lyapunov` always ran the full observer/subset
construction, enumerating subsets that cannot contribute: the sublevel set is a union over observer
states, so any state *containing* another is redundant. Dropping redundant supersets cut route 1's
primal cost by 40–60 % and **roughly halved our advantage there** (6.25× → ~3×). A reviewer will ask
whether the baseline is fair; this is the answer. Committed as `FIX utils: induced common shortcut and
bisection bracket`.

**One point in design space.** Everything is $M = 2$, De Bruijn order 1, in the plane. A $\theta$
sweep, order 2, $M = 3,4$ and dimension $> 2$ are untested.

**The regions' own contribution is unquantified.** Isolating it means varying region *geometry*, not
removing the regions, which only changes ladder depth.

---

## 7. Evidence that it is the graph, not the geometry

The natural alternative story — route 1 pays for node-dependent memory in facets — is false here, and
a controlled pair settles it. On two sheared modes $A_i = \rho\,[1\ {\pm}\alpha;\ 0\ 1]\,D$, ladder
pinned:

| $\alpha$ | piece gap | route 1 / route 2 cells |
| --: | --: | --: |
| **0.00** | **0.000** | **0.50** — route 2 loses |
| 0.20 | 0.023 | 0.58 |
| 0.40 | 0.076 | 0.68 |
| 0.60 | 0.076 | 0.86 |
| **0.80** | **0.000** | **1.84** — route 2 wins |

**Two points where the pieces are the same set, and opposite outcomes.** Geometry is held exactly
fixed while the result moves by 3.7×, so it cannot be the currency. What moves with it is $\alpha$,
the transversality of the modes.

**The underlying quantity is $B = A_2A_1^{-1}$.** Substituting $y = A_1x$ and using that a linear
bijection preserves a partition's combinatorics:

| arm | superposition | count |
| :--- | :--- | :--- |
| single node | $A_1^{-1}(P) \wedge A_2^{-1}(P)$ | $\lvert P \wedge B^{-1}P\rvert$ |
| primal node | $A_1^{-1}(P_1) \wedge A_2^{-1}(P_2)$ | $\lvert P_1 \wedge B^{-1}P_2\rvert$ |
| **dual node** | $A_1^{-1}(P_1) \wedge \mathbf{A_1^{-1}}(P_2)$ | $\lvert P_1 \wedge P_2\rvert$ |

The dual node's count is the only one free of $B$; the others refine a partition against a *linearly
twisted copy*, and the twist is $B$. So $A_1 = A_2 \Rightarrow B = I \Rightarrow$ no advantage, and
the design knob is the **relative** map, not $A_1$ and $A_2$ separately.

**Where it stands on the two examples, honestly.** A's $B$ has complex eigenvalues ($0.8407 \pm
0.227i$) — a genuine rotation with no real invariant direction. GLB's two modes differ by a **rank-one**
matrix, so $B = I + (\text{rank one})$ has eigenvalue exactly 1: a fixed direction, a shear. That
ordering matches the measured $g$ (4.81 against 1.55). But the scalar proxy tried for it,
$\lVert B - sI\rVert / \lVert B\rVert$, gives **0.955 for GLB against 0.368 for A** — the opposite
order — so it does **not** support the story and must not be quoted. What would settle it is measuring
$\lvert P \wedge B^{-1}P\rvert$ against $\lvert P_1 \wedge P_2\rvert$ directly.

**The advantage does not compound indefinitely.** Per-slice growth factors converge from the sixth
rung on, so the gain is a constant factor banked in the early slices, not a widening gap: within a
slice the cells form a planar arrangement bounded by Euler and by the number of distinct cut
directions. The reachable theorem is a bound with an explicit constant.

---

## 8. What is still to run

In priority order.

1. **GLB dual at 13 rungs.** The configuration that decides whether the conditional claim in §2c is
   stated with one supporting example or two. Estimated 4–6 h, ~74 000 cells per arm, memory risk.
2. **GLB primal at 7 rungs, 3 paired rounds** — the headline row, currently one round with the
   warm-up anomaly unresolved.
3. **Resolve the warm-up anomaly** in §5, without which no single-round number is quotable.
4. **Regenerate every figure.** None in this folder was produced under the corrected rule.
5. **`small_D` (8 rungs)** has not been run since the baseline correction.

---

## 9. Reproducing it

```
1_comparison.jl                    the matrix: every example × graph × case × size
2_piece_diversity.jl               the screen behind §7, and what selected example B
demo_primal_vs_dual_quotients.jl   [glb|A] — the 2×3 grid + paired timing (§4)
<example>/demo_route{1,2}_*.jl     one route per script, with the co-safe solve on GLB
```

Run from the repository root with `julia --project=test <path>`; `SAMPLES=n` sets the number of paired
rounds. One folder per example, each holding its own scripts and figures — see `README.md`.

| figure | what it shows |
| :--- | :--- |
| `fig_certificate_<graph>.png` | one panel per PCLF piece, then the common they induce — nested sublevel sets at the certificate's own rate (§2a, §3) |
| `fig_partitions_<graph>_<case>.png` | the abstraction: route 1's single plane, then route 2's $\lvert S\rvert$ lifted copies, on one shared window (§4a) |
| `fig_facets_<graph>_<case>.png` | facets-per-cell distribution, $\log_{10}$ on a linear axis — the tail is what costs |
| `fig_certified_<graph>.png` | the co-safe solve's certified region, both arms (§4c) |
| `primal_vs_dual_quotients.png` | the four quotients in one 2×3 grid, both orientations (§2) |
