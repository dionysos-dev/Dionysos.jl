# N1 — what a bisimulation quotient costs, and which of two routes builds it cheaper

Working notes behind the three experiments of `research/HSCC2027/`, which is where the finished runs
and the reported numbers live. What this file keeps is what that folder does not carry: why the two
De Bruijn graphs behave oppositely (§2, §7), the per-node threshold that predicts the sign of the
gain (§4a), the volume check that both routes certify the same region (§4b), and what is *not*
claimed (§6).

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

Three examples, and each is one of the three experiments of `research/HSCC2027/` — same dynamics,
same working set, same regions:

| here | there | system | placed where the competing explanation is impossible |
| :--- | :--- | :--- | :--- |
| `A_identical_pieces` | experiment 2 | `two_mode_problem` | the node pieces are the *same set*, so the whole difference is the graph |
| `B_diverse_pieces` | experiment 3 | the observer-study pair | A's working set, dynamics that force the pieces apart |
| `C_gol_lazar_belta` | experiment 1 | Example 3.1 of arXiv:1208.5471 | the external benchmark, and the hardest case for the claim |

Piece diversity is the largest relative gap between two node pieces' support functions over 64
directions — exact for convex sets, no volume needed, reads as a percentage.

| | A | B | C |
| :--- | --: | --: | --: |
| piece gap, primal | 0.0278 | 0.4599 | **0.1778** |
| piece gap, dual | < 1e-4 | — | **0.1068** |
| induced common, primal | 3 parts | 13 parts | **17 parts** |
| induced common, dual | 1 part | — | 1 part |

**A is the control.** On its dual arm the two node pieces are the same set, so both routes carry the
*same function*. A functional gain is impossible there anyway (§2a), so nothing can be confounded
with the algorithmic one.

**B is primal-only on purpose** — the dual graph is where diversity works against route 2 (§2a).

**C is the external benchmark.** Their dynamics, their working set $\{\lVert Lx\rVert_\infty \le
10\}$, their three observation regions, their co-safe formula and their initial point. What is *not*
theirs is the certificate: instead of their common Lyapunov function an oracle supplies a
path-complete one, and the two routes are compared on it.

---

## 4. Results

**The reported numbers are in `research/HSCC2027/README.md`**, measured under the construction's own
stopping rule with both arms alive for both timings. Nothing in this file supersedes them, and the
rows that were marked void here have been dropped rather than carried. Three findings live on
because that folder does not measure them.

### 4a. The per-node threshold $g > \lvert S\rvert$

On a **complete** graph route 2 builds more cells and wins anyway, on cell simplicity. On a
**co-complete** graph there is no fragmentation to avoid, so the only saving is the per-node
reduction $g$ of §2c — and route 2 pays $\lvert S\rvert$ copies of the pipeline, so
$g > \lvert S\rvert$ is break-even.

| example | graph | rungs used / needed | $g$ | time ratio | reading |
| :--- | :--- | :--- | --: | --: | :--- |
| A | dual | 7 / 5 — **valid** | **4.81** | 2.93× | $g > 2$ — route 2 wins |
| C (GLB) | dual | 9 / 13 — **void** | 1.55 | 0.72× | $g < 2$ — route 2 loses |

The two rows sit on opposite sides of the threshold and the outcome follows it in both. That is the
best support the predictor has, and it rests on **two points, one of which has never been re-run**.
GLB's dual needs 13 rungs and was measured at 9; it went against route 2 twice, 0.47× unwarmed and
0.72× warmed, so the loss is likely real, but it is not established.

**On which rows are valid.** Forcing the rung count bypasses the construction's stopping rule — the
first level whose terminal set no longer meets an observation region — and a `D` that still meets a
region labels every point of the overlap with the terminal observation instead of its own. But the
clearing property is **monotone in the rung count**: more rungs means a smaller terminal set, so once
it clears it stays clear. Forcing *more* rungs than needed is therefore harmless; forcing *fewer* is
what breaks. Example A was run at 7 and needs 5, so every example-A row here stands.

### 4b. Soundness — the arms certify the same region

Every cost number is void unless both routes answer the specification identically, and cell counts
cannot show it: 90 cells and 44 can cover the same area. The check is on **volume**, on example A.
**`research/HSCC2027/` does not carry this check**, so this is the only place the two routes are
measured to agree on what they certify.

| | route 1 | route 2 | raw | normalised by each arm's covered area |
| :--- | --: | --: | --: | --: |
| dual | 2.1522 | 2.2011 | 2.27 % apart | 0.05500 vs 0.05497 — **0.05 %** |
| primal | 2.0857 | 2.2029 | 5.62 % apart | 0.05392 vs 0.05458 — **1.22 %** |

The raw primal gap exceeds 5 % for a measurable reason rather than a soundness one: the arms do not
tile the same region (§2a) and `atol` insets every cut, so the arm making more cuts loses more area.
Route 1 averages 9.12 parts per cell against 3.85 and ends with 4 % less covered area — very nearly
the whole 5.6 %.

On GLB both routes are additionally checked by solving the co-safe formula under **both quantifiers**
and comparing the certified sets. A cheaper quotient is only worth having if it answers both
questions the same way.

### 4c. The cell ratio does not move with problem size

2.26 then 2.27 on A's dual arm, across quotients sixteen times larger. Each experiment in
`research/HSCC2027/` is reported at one depth, so this is the only evidence that the ratio is a
property of the pair rather than of the depth it happened to be measured at.

---

## 5. Protocol

**What is held fixed within a row.** Same system, same regions, same PCLF, same tolerance
(`atol = 1e-3`), same ladder. The ladder matters most: `ΓX`, the outer level, is computed once from
**route 1's** certificate and imposed on route 2, so the baseline keeps the geometry it would have
chosen and route 2 is made to conform. Taking it from route 2 would flatter route 2; without pinning,
the arms tile different regions and the counts measure coverage rather than efficiency.

**The rung count is derived, not chosen**, from the construction's own stopping rule. It differs
between the two orientations — 7 on GLB's primal against 13 on its dual — because the dual contracts
more slowly ($\gamma = 0.924$ against $0.867$). **Matching rung counts across orientations is wrong**;
matching depth is what is required.

**Timing.** Route 1's determinisation is inside its timed closure — handed a PCLF it must produce the
common before it can build anything. (Measured at **0.0 %** of route 1's time in every configuration:
the cost is not the determinisation step but what it returns.)

An anomaly recorded here was never resolved on its own configuration: on the GLB primal row the
*same* build took 666 s in the warm-up and 1 276 s in the timed round, for identical output, while
route 2 moved 8 %. It is consistent with what `research/HSCC2027/` later measured — a timed region
pays garbage collection proportional to the live heap, so an arm timed while more is alive is timed
slower — and that folder builds both arms before timing either for exactly this reason. It was not
confirmed here, so the 7.01× that row reported should still be read as an interval.

**Three rules, all learned the hard way.**

1. **Only same-run ratios.** Two runs of identical code produced the same cell counts to the unit and
   absolute times differing by 3.6×–4.0×. Absolute timings from this folder are meaningless in
   isolation.
2. **A time ratio is quotable only from `mid_D` or larger.** Disjoint spreads are necessary but not
   sufficient: example B at `large_D` was measured twice, both with disjoint spreads, at 8.12× and
   4.18×. At builds of 0.3–2.7 s the machine state dominates.
3. **Cell, facet and part ratios are exempt** — they are deterministic and reproduce to the unit.

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

## 8. What is still open

1. **GLB dual at 13 rungs.** The configuration that decides whether §4a's threshold is stated with
   one supporting example or two. Estimated 4–6 h, ~74 000 cells per arm, memory risk.
2. **Whether the two routes certify the same set on GLB**, not just on A. §4b checks volume on the
   control example only.
3. **More than one point in design space.** Everything is $M = 2$, De Bruijn order 1, in the plane. A
   $\theta$ sweep, order 2, $M = 3,4$ and dimension $> 2$ are untested.

---

## 9. Reproducing it

`comparison.jl` runs the three examples above, both graphs, both cases and two ladder depths, and
writes its figures next to itself, one folder per example and depth:

```
julia --project=test research/BisimulationQuotient/campaign/N1_quotient_cost/comparison.jl
```

`SAMPLES=n` sets the number of paired rounds. For the paper's numbers, run `research/HSCC2027/`
instead — three self-contained scripts, with the reported counts in its README.
