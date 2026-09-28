# Path-Complete Abstractions of Switched Systems

Code and data for the numerical experiments of the paper, submitted to HSCC 2027.

Finite-bisimulation constructions for switched systems partition the state space with the sublevel
sets of a *common* Lyapunov function [1]. The paper replaces that certificate by a *path-complete*
one and lifts the abstraction over the nodes of the certificate's graph. Given a path-complete
certificate there are then two ways to obtain a finite bisimulation quotient, and the experiments
compare them:

* **Approach 1.** Determinise the path-complete function into a common one, by the observer
  construction of [2], then run the single-node construction of [1] on the result. This is what the
  earlier method requires when it is handed a path-complete certificate.
* **Approach 2.** Run the construction on the path-complete function directly, with one partition per
  graph node. This is the paper's generalisation of [1] to the path-complete setting.

Both arms start from the same stability certificate and both produce a bisimulation. They therefore
answer the specification exactly, and give the same answer; what separates them is cost. The
experiments report it as three quantities:

* **the number of abstract states**: the cells of the quotient, and the size of the graph the fixed
  point behind a specification iterates over;
* **the complexity of a cell**: a cell is a semilinear set, a finite union of polytopes, and it is
  its polytopes and their facets that every geometric primitive of the construction is charged for,
  namely set difference, pre-image and emptiness;
* **the time** to build the quotient and to solve the specification on it.

The cases where the three disagree are the informative ones, and each experiment is built around one.

## Installation and use

```
julia --project=. -e 'using Pkg; Pkg.develop(path="../.."); Pkg.instantiate()'
julia --project=. run_all.jl
```

The first command points the environment at the copy of Dionysos in this repository and installs the
remaining dependencies. It is needed once and takes a few minutes of precompilation. `run_all.jl`
then runs the three experiments and writes the figures to `figures/`.

The experiments can also be run separately:

```
julia --project=. experiment1_gol_lazar_belta.jl
julia --project=. experiment2_primal_and_dual.jl
julia --project=. experiment3_diverse_pieces.jl
```

`SAMPLES=n` reports the median of `n` timed rounds instead of one. Each round rebuilds everything, so
`n = 3` roughly triples the runtime. `FIGURES=0` skips the figures. In experiment 1, `CACHE=0` forces
a rebuild instead of loading the quotient from `cache/`, which is what the reported build time
requires.

After precompilation, experiment 1 takes about five minutes, experiment 2 about one and experiment 3
about nine. Experiment 1 is the only one that caches its quotients, so with `CACHE=0` it spends six
minutes more rebuilding them. Every arm is built twice in any case, once for the figures and the
counts and once to be timed.

Every experiment builds both quotients first and times them afterwards, with both alive. This is not
incidental: a timed region pays garbage collection proportional to the live heap, so an arm timed
before the other exists is timed on a lighter one, and building and timing in a single pass favours
whichever arm goes first by enough to reverse a verdict. Each timed section is also preceded by an
untimed warm-up call, so Julia's compilation is excluded, and `SAMPLES=3` takes the median of three.

Even so, absolute times on a loaded machine are not reliable: identical work has differed by a factor
of two between runs here. The counts are deterministic and reproduce exactly, and it is those that
carry the mechanism. Where a time was of the same order as its own noise, it was replaced by the
count of geometric operations underneath it, which is what decides the time and can be checked.

## Contents

| file | |
| :--- | :--- |
| `experiment1_gol_lazar_belta.jl` | experiment 1 |
| `experiment2_primal_and_dual.jl` | experiment 2 |
| `experiment3_diverse_pieces.jl` | experiment 3 |
| `run_all.jl` | all three, in order |
| `Project.toml`, `Manifest.toml` | the pinned environment |
| `cache/` | experiment 1's quotients, once built |
| `figures/` | output |

Each experiment is a single self-contained file: it defines its own system, certificate and figures,
and shares no code with the others.

## Experiment 1: the Gol–Lazar–Belta benchmark

The system, the three observation regions, the co-safe LTL formula and the initial point are those of
Example 3.1 of [1].

The system is `x⁺ = A_σ x` with two modes:

```
A₁ = ⎡-0.65   0.32⎤      A₂ = ⎡ 0.65   0.32⎤
     ⎣-0.42  -0.92⎦           ⎣-0.42  -0.92⎦
```

The Lyapunov certificate is not theirs. Theirs is a common polyhedral function, built by hand, and it
certifies a contraction rate of 0.94. The path-complete framework searches instead, over the Lyapunov
functions a graph admits, and on the same system returns a markedly tighter certificate: rate
0.864529 on the order-1 primal De Bruijn graph, over a shared template of 7 placed cones.

The certificate itself. A polyhedral piece is the gauge `V_s(x) = maxᵢ |(G_s x)ᵢ| / wᵢ`, one weight
per row, and the solver is free to return any weights it likes. Here it returned 1 on every row of
both pieces, so the weights drop out and each piece is the infinity norm `V_s(x) = ‖G_s x‖_∞`, whose
Γ-sublevel set is the symmetric polytope `{x : |G_s x| ≤ Γ}`. Only `G` is shown for that reason,
`G₁` for node (1,) and `G₂` for node (2,):

```
       ⎡ 0.6116   0.3108⎤          ⎡ 0.6116   0.4224⎤
       ⎢ 0.5114   0.4674⎥          ⎢ 0.6116   0.4224⎥
       ⎢ 0.3236   0.5518⎥          ⎢ 0.3236   0.5518⎥
G₁ =   ⎢ 0.0866   0.5456⎥   G₂ =   ⎢ 0.0865   0.5456⎥
       ⎢-0.2170   0.4192⎥          ⎢-0.2169   0.4193⎥
       ⎢-0.4706   0.1956⎥          ⎢-0.5402   0.1343⎥
       ⎣-0.6116  -0.0695⎦          ⎣-0.6116   0.0000⎦
```

Both approaches then build a quotient and solve the specification under both quantifiers: synthesis,
where the modes belong to the controller, and verification, where they belong to the environment.

| | approach 1 | approach 2 | ratio 1/2 |
| :--- | --: | --: | --: |
| cells | 2434 | 4625 | 0.53 |
| Σ polytopes | 22438 | 16847 | **1.33** |
| mean polytopes / cell | 9.22 | 3.64 | 2.53 |
| max polytopes / cell | 156 | 51 | **3.06** |
| mean facets / cell | 36.70 | 14.52 | 2.53 |
| max facets / cell | 622 | 205 | 3.03 |
| build (s) | 255.50 | 130.54 | **1.96** |
| solve ∃, a graph fixed point (s) | 0.06 | 0.09 | 0.66 |
| solve ∀, a graph fixed point (s) | 0.04 | 0.08 | 0.52 |

The table answers three different questions, and they do not all favour the same arm.

**Building the quotient is geometry**, charged per polytope of a cell: set difference, pre-image and
emptiness. Approach 2 has 4625 cells against 2434 but 16847 polytopes against 22438, and its worst
cell is three times simpler. It builds the larger quotient in half the time.

**Solving the specification is a graph fixed point** on the quotient, and nothing else: the automaton
is assembled from transitions the build already stored, so no set operation runs. It therefore
follows the number of cells, and approach 1 wins with half as many.

**Abstracting a concrete set is geometry again.** Asking which abstract states a concrete set meets
tests it against every cell, and a cell is a union of polytopes, so the answer costs one feasibility
LP per polytope. That is the `Σ polytopes` row, and approach 2 carries a third fewer. Such a request
is therefore answered faster on its quotient, although that quotient has nearly twice the cells.

The complexity of a cell decides the geometry, which is building the quotient and abstracting a
concrete problem onto it, the specification's own sets included. Both are priced in polytopes. The
number of cells decides what comes after, once the quotient is built and the specification
abstracted: the fixed point, which is graph work and nothing else. Building the quotient dominates by
far, so approach 2 is much faster overall.

## Experiment 2: the primal and the dual De Bruijn graph

One system, one certificate family and one specification; only the meaning of a graph node changes.
In the primal graph a node records the mode just played, so every node emits every mode. In the dual
graph a node commits to the mode played next, so every node emits one.

The system is `x⁺ = A_σ x`, with two observation regions:

```
A₁ = ⎡0.70  0.10⎤      A₂ = ⎡0.60  -0.15⎤
     ⎣0.00  0.65⎦           ⎣0.10   0.55⎦
```

The two certificates: primal, complete, certified rate 0.715897; dual, co-complete, rate 0.717572.

All four pieces, both graphs and both nodes, share the same `G`: the rotation by π/6 that the
template fixes. The whole certificate therefore sits in the per-row weights of
`V_s(x) = maxᵢ |(G x)ᵢ| / wᵢ`.

```
G = ⎡0.8660  -0.5000⎤
    ⎣0.5000   0.8660⎦
```

| | w₁ | w₂ |
| :--- | --: | --: |
| primal, node (1,) | 2.9622 | 2.1975 |
| primal, node (2,) | 2.8774 | 2.2152 |
| dual, node (1,) | 2.8516 | 2.1648 |
| dual, node (2,) | 2.8517 | 2.1647 |

The two dual rows agree to four decimals while the two primal rows differ in the second. Pieces this
close are what this experiment needs, since it isolates the graph and nothing else.

| | primal (1) | primal (2) | dual (1) | dual (2) |
| :--- | --: | --: | --: | --: |
| cells | 400 | **769** | 378 | **172** |
| Σ polytopes | 3619 | **2930** | 1401 | **881** |
| mean polytopes / cell | 9.05 | 3.81 | 3.71 | 5.12 |
| max polytopes / cell | 174 | **68** | 38 | 34 |
| mean facets / cell | 36.09 | **15.19** | 14.74 | 20.24 |
| max facets / cell | 699 | **277** | 154 | 135 |
| build (s) | 18.67 | **6.72** | 3.95 | **2.99** |

The quotient of approach 2 is larger than approach 1's on the primal graph and smaller on the dual,
and approach 2 builds faster on both. Which of the two De Bruijn graphs is used therefore changes the
size of the quotient substantially without changing which construction is cheaper, so the number of
cells does not predict the build. What it does predict is the total number of polytopes, which falls
in both cases.

## Experiment 3: the same two graphs with pieces that do not coincide

Experiments 1 and 2 use certificates whose node pieces are close to one another. This one changes the
dynamics so that they are not. The diversity is discovered rather than imposed: both nodes are given
the same conic template and the solver returns two markedly different pieces.

The system is `x⁺ = A_σ x`:

```
A₁ = 1/10 · ⎡1.5519  0.4474⎤      A₂ = 1/10 · ⎡0.4750  9.1755⎤
            ⎣7.6412  7.4716⎦                  ⎣1.8955  0.1850⎦
```

There are no observation regions here. With regions, part of the refinement serves to respect them,
and that part is work both approaches do identically; removing them leaves the certificate's own
geometry as the only thing driving the partition.

The primal and the dual order-1 De Bruijn graphs are both run, as in experiment 2, and read beside it
they separate the two channels behind the result. The two certificates: primal, complete, rate
0.902151; dual, co-complete, rate 0.869447.

Here the template is a conic partition of order 2, eight rows per piece, and every weight came back
at 1, so each piece is again an infinity norm. The two pieces of the primal certificate:

```
       ⎡ 0.7519   0.0527⎤          ⎡ 0.4061   0.3627⎤
       ⎢ 0.6981   0.1605⎥          ⎢ 0.3023   0.5704⎥
       ⎢ 0.5856   0.2729⎥          ⎢ 0.1995   0.6731⎥
G₁ =   ⎢ 0.2905   0.4205⎥   G₂ =   ⎢ 0.0915   0.7272⎥
       ⎢-0.3370   0.4205⎥          ⎢-0.0258   0.7272⎥
       ⎢-0.5713   0.3034⎥          ⎢-0.1481   0.6660⎥
       ⎢-0.6926   0.1820⎥          ⎢-0.2781   0.5360⎥
       ⎣-0.7519   0.0634⎦          ⎣-0.4061   0.2799⎦
```

Set these beside experiment 2, where one rotation served all four pieces. Nothing was done to force
them apart: the same template went in for both nodes and the solver came back with two markedly
different gauges, `G₁` reaching far along x₁ and `G₂` along x₂. The dual certificate's two pieces
stay much closer to one another, and its node (1,) even repeats one row three times, leaving six
distinct facet normals out of eight; the script prints both.

| | primal (1) | primal (2) | dual (1) | dual (2) |
| :--- | --: | --: | --: | --: |
| cells | 1429 | 1706 | 3815 | **2463** |
| Σ polytopes | 6490 | **4391** | 13871 | **11274** |
| mean polytopes / cell | 4.54 | **2.57** | 3.64 | 4.58 |
| max polytopes / cell | 128 | **16** | 62 | **42** |
| mean facets / cell | 18.14 | **10.31** | 14.50 | 18.24 |
| max facets / cell | 510 | **64** | 244 | **170** |
| build (s) | 13.86 | **4.61** | 31.81 | **14.10** |

On the primal graph approach 2 builds more cells, and much simpler ones: 16 polytopes in its worst
cell against 128. On the dual graph it builds fewer cells of roughly the same complexity. Approach 2
is therefore the faster of the two on both De Bruijn graphs, the primal and the dual, but for a
different reason on each.

## References

[1] E. Aydin Gol, X. Ding, M. Lazar and C. Belta. Finite bisimulations for switched linear systems.
*IEEE Transactions on Automatic Control*, 59(12):3122–3134, 2014.

[2] D. Angeli, N. Athanasopoulos, R. M. Jungers and M. Philippe. Path-complete graphs and common
Lyapunov functions. In *Proceedings of the 20th International Conference on Hybrid Systems:
Computation and Control (HSCC)*, pages 81–90, 2017.
