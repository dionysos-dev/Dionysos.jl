# Paper notes — the story, and the experiments that serve it

Design note, not a draft. Part one is the narrative chain for the introduction, with the two links
that do not currently hold. Part two decides *which* experiments go in and what each is there to
prove. Space is tight because the algorithm needs the room, so the rule is one claim per experiment
and nothing that needs a paragraph of theory to read.

> **References** are given by their key in `JCVD/sample-base.bib` wherever the entry already exists —
> which turned out to be every one that matters. Anything not in that file is flagged as such.

---

# Part one — the story

## The chain, repaired

Your draft is the right spine. Two links need work before a reviewer sees them, and I have marked
them ⚠.

**1. The problem.** Synthesis *and* verification of co-safe LTL specifications over regions of the
state space, for discrete-time switched linear systems whose modes are the control input.

**2. Why a finite abstraction.** Do not claim symbolic control is "the best" — that invites a reviewer
to name the alternatives (mixed-integer encodings, barrier certificates, funnel-based methods) and you
have to dispatch them. Claim *necessity* instead, which is unarguable: automata-theoretic synthesis
against an LTL specification needs a **finite** transition system, so the state space must be
partitioned. That is the definition of the setting, not a claimed superiority.
→ Tabuada, *Verification and Control of Hybrid Systems: A Symbolic Approach*, Springer 2009, is the
standard reference. Belta, Yordanov & Aydin Gol, *Formal Methods for Discrete-Time Dynamical Systems*,
Springer 2017, is closer to this exact setting — and it is the same group as the paper you build on,
which is a good citation to lean on.

### 3. "Simulation-based abstraction is conservative, so we want a bisimulation"

**The premise is right and it is the core of the argument** — I under-stated this earlier. An
abstraction related to the concrete system by a *simulation*, an *alternating simulation* or a
*feedback refinement relation* is **sound but not complete**, and the consequence is the one to lead
with: **when it fails, you cannot conclude anything.** You cannot tell "no controller exists" from
"my partition was too coarse". The synthesis question is left open, and refining once more is always
an available excuse.

In this setting the mechanism is concrete and worth one sentence: discretising a cell forces an
over-approximation of its reachable set, which introduces **spurious transitions**. Those spurious
transitions can make the abstract problem infeasible while the concrete one is perfectly feasible.
Conservatism is therefore not a small loss at the boundary — it can flip the answer.

References already in `sample-base.bib`:

| | |
| :--- | :--- |
| `tabuada2009verification` | Tabuada, *Verification and Control of Hybrid Systems: A Symbolic Approach*, Springer 2009 — the reference for simulation / alternating simulation and what each preserves |
| `reissig2016feedback` | Reissig, Weber & Rungger, *Feedback refinement relations for the synthesis of symbolic controllers*, IEEE TAC 2016 — the sound-but-incomplete synthesis relation, stated as such |
| `chutinan2001verification` | Chutinan & Krogh, *Verification of infinite-state dynamic systems using approximate quotient transition systems*, IEEE TAC 2001 — explicitly **over-approximating** quotients |
| `girard2009approximately` | Girard, Pola & Tabuada, *Approximately bisimilar symbolic models for incrementally stable switched systems*, IEEE TAC 2009 — the approximate route, and the closest neighbour to this paper's setting |
| `browne1988characterizing`, `baier2008principles` | what an **exact** bisimulation buys: preservation of everything expressible in CTL/LTL |

**Two further arguments, in this order:**

1. **A negative answer becomes a result.** Under a bisimulation, failure to find a controller is a
   **proof of infeasibility** for the concrete system. Refinement can never deliver that, because one
   can always refine once more. This is the strongest form of your own point and I would make it the
   headline sentence of the paragraph.
2. **Both quantifiers from one quotient.** A simulation relation transfers one direction. A
   bisimulation transfers both, so the same quotient answers "does *some* switching signal satisfy
   $\varphi$" and "does *every* one". Your experiments deliver exactly this on GLB's own formula, so
   the introduction should promise it.

**The one objection to pre-empt: "why not simply refine?"** Because refinement removes conservatism
only in the limit and carries **no stopping rule** — and the classical partition-refinement algorithm
is not guaranteed to terminate on continuous dynamics at all (`milner1989communication`,
`yordanov2010formal`). That is precisely the gap the stability certificate closes, and it is the pivot
into the next link: contraction bounds the number of slices, so the refinement *terminates by
construction*.

### ⚠ 4. "Bisimulation is infinite-time to construct" — true, but state it precisely

The honest statement is about **existence**, not just cost: for general dynamical systems a finite
bisimulation quotient **need not exist**, and classical partition refinement is not guaranteed to
terminate. Finiteness is known only for restricted classes. The draft introduction already says this
correctly and already cites the right entries — `alur1994theory` (timed automata, the canonical finite
bisimulation), `tabuada2006linear` (a class of linear systems), `milner1989communication` and
`yordanov2010formal` (non-termination). **Nothing to change here**; I raised it only because
"infinite-time to construct" understates it — the issue is that it may not exist at all.

**5. Gol, Ding, Lazar & Belta** (`gol2014finite`, IEEE TAC 2014, 59(12):3122–3134 — verified against
your `.bib`). They obtain a **finite** bisimulation for stable switched linear systems by using a
**stability certificate** — a common polyhedral Lyapunov function — to generate the partition: its
sublevel sets give a geometric ladder, the contraction rate bounds the number of rungs between the
working set and the terminal set, and the refinement therefore terminates *by construction*.

**Say explicitly that stability is the enabling hypothesis**, not a technicality. Contraction is what
makes the ladder finite, and it is why your paper inherits the same hypothesis and the same
restriction to scLTL (`kupferman2001model`).

### The pivot their own introduction hands you

Their introduction motivates the work as a bridge between control theory (complex dynamics, simple
specifications) and formal methods (simple systems, rich languages), then argues: finite bisimulations
exist only for restricted classes, the standard algorithm does not terminate, so restrict to a bounded
region and use a polyhedral Lyapunov function — which they note is *necessary* for stability under
arbitrary switching. And they present using **a single** Lyapunov function as an **advantage** over
neighbouring work that needs several.

*(That last point comes from a summary of the arXiv version, not the PDF — **verify it in the TAC
version before building a sentence on it**. If it holds, it is the best pivot available.)*

Because your paper inverts precisely that. One sentence:

> What their construction treats as a simplification — a single certificate for all modes — is
> where its geometric cost concentrates.

That is a sharper position than "we generalise their construction", and it costs one sentence.

### ⚠ 6. Our step — one claim here is false as usually phrased

The draft introduction currently says *"for many switched systems, simple common Lyapunov functions
may fail to exist altogether, even when the system is stable under arbitrary switching."* A reviewer
in this community will answer that **polyhedral common Lyapunov functions are universal** — for any
switched linear system with joint spectral radius below 1 one *exists* (`Jung09`), and GLB's own
introduction leans on exactly that when it calls polyhedral Lyapunov functions *necessary* for
stability under arbitrary switching. Someone reading both papers will see the tension. The campaign
also measured the consequence: with rich enough polyhedral templates a PCLF and a common function
certify the same rate to four decimals.

**The correct claim is a lower bound on complexity, and it is stronger than non-existence.** It is
also already in your bibliography, and it is the "papier limitant" of your own commented-out note:

| | |
| :--- | :--- |
| **`ahmadi2016lower`** | Ahmadi & Jungers, *Lower bounds on complexity of Lyapunov functions for switched linear systems*, Nonlinear Analysis: Hybrid Systems 2016 — **the reference that carries this link.** It is what licenses your abstract's phrasing, "that complexity admits no bound uniform over stable switched systems" |
| `AJPR:14` | Ahmadi, Jungers, Parrilo & Roozbehani, *Joint spectral radius and path-complete graph Lyapunov functions*, SICON 2014 — the foundational PCLF reference |
| `JunAhm_2017` | Jungers, Ahmadi, Parrilo & Roozbehani, *A characterization of Lyapunov inequalities for stability of switched systems*, IEEE TAC 2017 |
| `AAJP2017` | Angeli, Athanasopoulos, Jungers & Philippe, *Path-Complete Graphs and Common Lyapunov Functions*, HSCC 2017 — **the observer/subset construction** that determinises a PCLF; this is the object route 1 builds |
| `philippe2018path` | Philippe, Athanasopoulos, Angeli & Jungers, *On path-complete Lyapunov functions: geometry and comparison*, IEEE TAC 2018 |
| `Jung09` | Jungers, *The Joint Spectral Radius: Theory and Applications*, Springer 2009 — for universality, if you state it |

**Align the introduction with the abstract.** The abstract already has the right formulation; the
introduction is the loose one. Replace "may fail to exist" with the bounded-complexity statement and
cite `ahmadi2016lower`.

> **Housekeeping:** the bibliography holds each of two key papers twice — `AJPR:14` and
> `ahmadi2014joint` are the same SICON 2014 article, `AAJP2017` and `angeli2017path` the same HSCC
> 2017 paper. Citing both keys of a pair prints it twice in the reference list. Pick one of each.

So the argument becomes, cleanly:

> A PCLF is the cheaper certificate to obtain at a given template complexity. But the prior
> construction consumes a *common* function, so using a PCLF with it requires determinising it first —
> and determinisation is exactly what destroys the low complexity that made the PCLF attractive.

**And be precise about what determinisation costs**, because "exponential" alone is too loose and the
complete-graph case is not exponential:

| graph | induced common | cost |
| :--- | :--- | :--- |
| complete | $\min_i V_i$ | $\lvert S\rvert$ pieces, but the sublevel set is their **union** — non-convex, and every set operation must handle it as a union of convex parts |
| co-complete | $\max_i V_i$ | one convex piece — cheap, but it certifies a strictly smaller region |
| general | $\min_{\mathcal S}\max_{i\in\mathcal S} V_i$ over observer states | up to $2^{\lvert S\rvert}$ states — **this** is where the exponential lives |

The bisimulation algorithm's primitives — set difference, pre-image, emptiness — are superlinear in
the number of convex parts a cell carries. That is the entire cost argument, and your measurements
give it a number: on the Gol–Lazar–Belta example the worst cell holds **399** parts after
determinisation against **42** without it.

**7. The idea.** Do not build the induced common. Run the bisimulation construction **on the PCLF
directly**, one partition per graph node.

**8. The consequence, and sell it as a feature.** The quotient is augmented by the graph's node, so
the controller carries $\lvert S\rvert$ states of memory. Two things to say about that:

- the memory is **small and known in advance** — it is the certificate's own graph, not something
  discovered during synthesis;
- on the primal De Bruijn graph it is **free**: the node records the mode just played, which the
  controller already knows. No extra state has to be carried at run time.

**9. The extension.** Constrained switching — see the future-work section below.

## What I would add that is not in your draft

**A counterintuitive hook for the abstract.** *We show that a bisimulation quotient with 50 % more
states can be an order of magnitude cheaper to compute.* It is memorable, it is what your data says,
and it reframes the contribution from "an optimisation" to "the cost model was wrong".

**A one-line thesis that covers both the result and the extension.** *The structure of the certificate
should be preserved, not flattened.* Determinisation flattens a graph-indexed certificate into a
single function; everything downstream then pays for the geometry that flattening creates. This one
sentence carries the speed result *and* the constrained-switching extension, which is why I would put
it in the abstract.

**Constrained switching is not merely an extension — it is where the argument becomes structural.**
Under a switching constraint the admissible mode sequences are given by an automaton, so the problem
*is* graph-structured. A common Lyapunov function has no place to put that structure; a path-complete
certificate does, because its graph can be composed with the constraint automaton. Determinising then
throws away precisely the object the problem hands you. Framed this way the contribution stops being
"we are faster" and becomes "this is the right certificate for this class of systems" — a much better
position to defend.
→ Philippe, Essick, Dullerud & Jungers, *Stability of discrete-time switching systems with
constrained switching sequences*, Automatica 2016, is the reference for the constrained-switching
setting with multiple Lyapunov functions.

**What you cannot claim, and should pre-empt.** No instance in the campaign certifies a rate that a
common function cannot (universality again). So the contribution is **computational, not
functional** — say it yourself, in one clause, before a reviewer says it for you. It costs nothing
because the computational claim is the one you are making anyway.

---

# Part two — the experiments

## What the two experiments have to establish

| | claim | why it needs its own experiment |
| :--- | :--- | :--- |
| **1** | Building directly on the PCLF is **several times faster** than determinising it first, on a benchmark the reader already knows, and it answers the same specification under both quantifiers. | credibility — external problem, external specification |
| **2** | The speed-up comes from **two different mechanisms**, and **cell count is the wrong cost proxy**: route 2 is faster whether it builds more cells or fewer. | otherwise a reader assumes "faster = smaller quotient" and the contribution looks like an optimisation |

Everything else the campaign measured — the $g > \lvert S\rvert$ threshold, the piece-diversity sweep,
the $B = A_2A_1^{-1}$ analysis — is deferred. It is the material that needs theory to read.

---

## Experiment 1 — Gol–Lazar–Belta Example 3.1, primal De Bruijn

**The setup.** Their dynamics, their working set $\{\lVert Lx\rVert_\infty \le 10\}$, their three
observation regions, their co-safe formula

$$\varphi = (\lnot R_2 \,\mathcal{U}\, D) \wedge \mathbf{F} R_1 \wedge \big((R_3 \Rightarrow \mathbf{X}\lnot R_1) \,\mathcal{U}\, D\big)$$

and their initial point $a = (-4, -7)$. The certificate is **ours**: an oracle supplies a
path-complete Lyapunov function on the order-1 primal De Bruijn graph, two nodes, conic order-2
template. The terminal set $D$ is fixed by the construction's own rule — the first level clearing
$R_1, R_2, R_3$ — which is also how their $\Gamma_D = 5.063$ is fixed.

**The two arms.** Route 1 determinises the PCLF into a common Lyapunov function and runs the
single-node construction on it. Route 2 runs the construction on the PCLF directly, one partition per
node.

**The numbers** (7 rungs, measured once — see *Still to run*):

| | cells | worst cell | build |
| :--- | --: | --: | --: |
| route 1 — determinise first | 6 613 | **399 parts** | 1 276 s |
| route 2 — directly on the PCLF | **10 158** | **42 parts** | **182 s** |

The sentence the experiment exists to produce: **route 2 builds 1.5× more cells, with cells an order
of magnitude simpler, and finishes 7× sooner.** Route 2's 10 158 cells sit next to the 10 611 of the
published reference run, which is the check that this is genuinely their example and not a shrunken
version of it.

**Figures — three panels, matching the layout already used in the thesis/paper:**

1. **the problem** — the three regions, $D$, $x_0$, and the synthesised closed-loop trajectory
2. **synthesis (∃)** — the certified region for $\varphi$ when the modes are the controller's
3. **verification (∀)** — the certified region when the modes are the environment's

Both quantifiers, because their example poses both and because *a cheaper quotient is only worth
having if it answers both questions identically*. The two green sets are therefore the soundness
check, not decoration — if route 1's and route 2's disagree, every timing above is void.

If a fourth panel fits, the strongest one is route 1's fragmented plane beside route 2's two clean
per-node planes: it shows the 399-against-42 in a way no table does.

### One framing to get right

**Do not write "faster than Gol–Lazar–Belta's algorithm."** Their algorithm is the $\lvert S\rvert = 1$
case of ours, and route 1 *is* their construction — but applied to the common induced from our PCLF,
not to their published $\lVert Lx\rVert_\infty$. Comparing against their own certificate would compare
*certificates*, not algorithms: a certificate determines its own slice family and terminal set, so the
two quotients would partition different regions and the cell counts would answer different questions.

The defensible framing, which is also the more useful one:

> Given a path-complete certificate, the prior construction can only be used by determinising it
> first. That is route 1. We compare against it on the same certificate.

Stated that way the baseline is unimpeachable and the comparison is purely algorithmic.

---

## Experiment 2 — the same system on both graph orientations

**Your instinct is right, the example choice needs one correction.** You proposed using a
non-aligned-pieces example to show the dual case better. That is exactly backwards for the dual, and
here is why in one line:

- on a **complete** (primal) graph route 1's induced common is $\min_i V_i$, whose sublevel set is the
  **union** of the node pieces — it fragments as the pieces separate, so **piece diversity is what
  makes the primal case shine**;
- on a **co-complete** (dual) graph it is $\max_i V_i$, the **intersection** — convex however
  different the pieces are. Diversity there costs route 2 and costs route 1 nothing.

Measured: forcing the pieces apart on the dual graph (gap 0.34–0.46) puts route 2's cell ratio at
0.29–0.41, i.e. route 2 *loses*. A diverse-piece example on the dual would demonstrate the opposite of
what you want.

**So: diverse pieces for the primal claim (that is GLB, gap 0.178, 17-part common — experiment 1
already has it), interchangeable pieces for the dual claim.**

**The proposal.** `A_identical_pieces` — `two_mode_problem`, two observation regions, shared rotated
4-facet template — run on **both** orientations, changing nothing but the graph:

| | route 1 cells | route 2 cells | | speed-up |
| :--- | --: | --: | :--- | --: |
| **primal** (complete) | 3 814 | 6 833 | route 2 has **1.8× more**, and simpler | 2.11× |
| **dual** (co-complete) | 3 615 | 1 502 | route 2 has **2.4× fewer** | 2.93× |

**These two rows are valid.** They were built at 7 rungs where the stopping rule needs 5, and the
clearing property is monotone in the rung count — more rungs, smaller terminal set, so once `D`
clears the observation regions it stays clear. Verified on this example: it first holds at $k = 5$
and holds at every $k \ge 5$. (Only GLB was affected by the forced-ladder defect, because there the
forced count was *below* what the rule requires.)

This is the cleanest control the campaign has: **one system, one certificate family, one
specification; only the meaning of a node changes.** The cell count swings by a factor of 4.5 between
the two rows and the time verdict does not move. That is the "cell count is the wrong proxy" claim,
demonstrated rather than asserted.

**Figure:** one 2×3 grid — route 1's single plane against route 2's two per-node planes, primal above,
dual below. The script that produces it is `research/HSCC2027/experiment2_primal_and_dual.jl`. At a
glance the primal row shows route 2 finer but clean, the dual row shows it visibly coarser.

### Show it, do not explain it

**Report the difference and stop there.** Do not open the union/intersection dichotomy, the two
channels or the $g$ threshold. The unexplained observation is worth more here than a compressed
explanation would be: it is what *motivates* the future work, and an explanation given in half a
paragraph would only invite the reviewer to ask for the other half.

What the reader needs is the minimum to interpret the figure — one sentence of description, no
mechanism:

> In the primal graph a node records the mode just played; in the dual it commits to the mode played
> next. Nothing else changes: same system, same certificate family, same specification.

Then the observation, flatly:

> The quotient is 1.8× *larger* than route 1's in the first case and 2.4× *smaller* in the second,
> while route 2 is the faster arm in both. The orientation of the graph therefore changes the size of
> the quotient by a factor of 4.5 without changing which construction is cheaper.

And the hand-off, one clause:

> What governs this — and how to predict it from the graph — is left to future work.

Three short blocks, no theory. The reader is left with a concrete, quantified puzzle rather than a
half-argument, and the future-work paragraph inherits it.

---

## What I would leave out

| | why |
| :--- | :--- |
| `B_diverse_pieces` | region-free, so no specification to pose; its 7.61× is impressive but it duplicates GLB's message with a synthetic system, and GLB is more credible |
| the $g > \lvert S\rvert$ threshold | needs the two-channel decomposition to state → **future work** |
| the piece-diversity sweep (§7 of `plan.md`) | it is a *control* against an alternative explanation. Worth a single sentence if a reviewer pushes, not a figure → **future work** |
| GLB on the dual graph | route 2 **loses** there ($g = 1.55 < 2$). Not hidden — it is one of the two points the future-work threshold rests on |

**Announce the scope explicitly**: the experiments are on **complete** path-complete graphs, the
standard De Bruijn construction, plus one controlled co-complete instance. Saying so costs one clause
and removes the only line of attack.

---

## Future work

Announcing *"the impact of the graph structure, such as the primal/dual property, on the bisimulation
quotient, and the constrained-switching case"* is the right future-work paragraph, and experiment 2
sets it up rather than competing with it: the experiment leaves a **quantified, unexplained
observation** on the table — a factor of 4.5 in quotient size from flipping the orientation — and this
paragraph is where it gets picked up. Showing without explaining is what makes the sequence work.

This is also where the $g > \lvert S\rvert$ threshold, the piece-diversity control and the
$B = A_2A_1^{-1}$ analysis belong.

**Constrained switching is the stronger half of the paragraph**, and worth leading with. There, the
graph is **imposed by the switching constraint** rather than chosen: the automaton comes from the
application. So "which orientation do I get, and what does it cost me" stops being a free design
choice and becomes a question one is forced to answer — which is precisely what makes the study worth
a paper rather than a remark.

**What you already hold in reserve for it**, and which I would *not* reveal here: the two mechanisms
decomposition, the break-even $g > \lvert S\rvert$ with two measured points on either side of it
(4.81 → route 2 wins; 1.55 → route 2 loses), the union/intersection dichotomy behind it, and the
counter-example showing "incomplete and mode-committed" is not sufficient. That is a follow-up paper's
worth of content, and it is already measured. Saying "future work" while holding it is honest; the
paragraph promises a systematic study, not that nothing is known.

---

## Still to run — nothing here is blocked on design

1. **GLB primal, 3 paired rounds.** The 7× is currently a single round, and on that round the same
   build took 666 s in the warm-up and 1 276 s in the timed round — identical output both times, while
   route 2 moved 8 %. Until that is understood the headline is the interval **~3.7–7×**, not a number.
   About 45 min.
2. **Example A, both orientations, with regions, under the corrected terminal-set rule.** The two rows
   in experiment 2 were measured with the rung count forced, which leaves $D$ overlapping the
   observation regions. The *ordering* is not in doubt (it reproduced across runs) but the numbers are
   not quotable.
3. **Regenerate every figure.** None in the folder was produced under the corrected rule.

Rough total: a few hours, none of it the 4–6 h GLB-dual run, which this design does not need.

## Space budget

| | |
| :--- | :--- |
| experiment 1 | 3 panels + a 3-row table + ~6 sentences |
| experiment 2 | 1 grid figure + a 2-row table + ~4 sentences |
| scope + protocol | ~4 sentences (what is held fixed; paired timing; the terminal-set rule) |

Two figures and two small tables. If only one figure survives, keep experiment 1's synthesis/
verification pair — it carries both the credibility and the soundness check.
