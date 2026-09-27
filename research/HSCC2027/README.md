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
julia --project=. experiment2_graph_orientation.jl
julia --project=. experiment3_diverse_pieces.jl
```

`SAMPLES=n` reports the median of `n` timed rounds instead of one. Each round rebuilds everything, so
`n = 3` roughly triples the runtime. `FIGURES=0` skips the figures. In experiment 1, `CACHE=0` forces
a rebuild instead of loading the quotient from `cache/`, which is what the reported build time
requires.

After precompilation and with `CACHE=0`, experiment 1 spends about six minutes on its two builds,
experiment 2 about half a minute and experiment 3 about twenty seconds. Figures add several more
minutes to experiment 1, whose quotient figure colours every cell separately.

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
| `experiment2_graph_orientation.jl` | experiment 2 |
| `experiment3_diverse_pieces.jl` | experiment 3 |
| `run_all.jl` | all three, in order |
| `figures/` | output |

Each experiment is a single self-contained file: it defines its own system, certificate, ladder and
figures, and shares no code with the others.

## Experiment 1: the Gol–Lazar–Belta benchmark

The system, three observation regions, co-safe LTL formula and initial point are those
of Example 3.1 of [1].

The system is `x⁺ = A_σ x` with two modes:

```
A₁ = [-0.65  0.32]      A₂ = [ 0.65  0.32]
     [-0.42 -0.92]           [-0.42 -0.92]
```

The Lyapunov certificate is not theirs. Theirs is a common polyhedral function, hand-crafted, and it
certifies a contraction rate of 0.94; reproducing their construction from it takes 11 slices. The
path-complete framework replaces that hand by a systematic search over the Lyapunov functions a graph
admits, and on the same system it returns a markedly tighter certificate: rate 0.864529 on the
order-1 primal De Bruijn graph, over a shared template of 7 placed cones. A tighter rate contracts
the ladder faster, and this one reaches a region-free terminal level in 6 slices.

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

The complexity of a cell decides the geometry: building the quotient, and abstracting a concrete
problem onto it, the sets a specification refers to included. Both are priced in polytopes. The
number of cells decides what comes after, once the quotient is built and the specification
abstracted: the fixed point, which is graph work and nothing else. Building the quotient dominates by
far, so approach 2 is much faster overall.

## Experiment 2: the two orientations of the graph

One system, one certificate family and one specification; only the meaning of a graph node changes.
In the primal graph a node records the mode just played, so every node emits every mode. In the dual
graph a node commits to the mode played next, so every node emits one.

The system is `x⁺ = A_σ x`, with two observation regions:

```
A₁ = [0.70  0.10]      A₂ = [0.60 -0.15]
     [0.00  0.65]           [0.10  0.55]
```

The two certificates: primal, complete, certified rate 0.715897, induced common in 3 polytopes,
ladder ΓX = 1.2950 over 5 slices; dual, co-complete, rate 0.717572, induced common convex, ladder
ΓX = 1.3252 over 5 slices.

| | primal (1) | primal (2) | dual (1) | dual (2) |
| :--- | --: | --: | --: | --: |
| cells | 400 | **769** | 378 | **172** |
| Σ polytopes | 3619 | **2930** | 1401 | **881** |
| mean polytopes / cell | 9.05 | 3.81 | 3.71 | 5.12 |
| max polytopes / cell | 174 | **68** | 38 | 34 |
| mean facets / cell | 36.09 | **15.19** | 14.74 | 20.24 |
| max facets / cell | 699 | **277** | 154 | 135 |
| build (s) | 18.67 | **6.72** | 3.95 | **2.99** |

The quotient of approach 2 is larger than approach 1's in the first case and smaller in the second,
and approach 2 builds faster in both. The orientation therefore changes the size of the quotient
substantially without changing which construction is cheaper, so the number of cells does not predict
the build. What it does predict is the total number of polytopes, which falls in both cases.

## Experiment 3: the same flip with pieces that do not coincide

Experiments 1 and 2 use certificates whose node pieces are close to one another. This one changes the
dynamics so that they are not. The diversity is discovered rather than imposed: both nodes are given
the same conic template and the solver returns two markedly different pieces, a support-function gap
of 0.46 against 0.03.

```
A₁ = 1/10 · [1.5519  0.4474]      A₂ = 1/10 · [0.4750  9.1755]
            [7.6412  7.4716]                  [1.8955  0.1850]
```

There are no observation regions here. With regions, part of the refinement serves to respect them,
and that part is work both approaches do identically; removing them leaves the certificate's own
geometry as the only thing driving the partition. The ladder then has to be given explicitly, since
the construction's stopping rule is satisfied immediately when there is nothing to clear, and the
outer level is taken from a probe that lets approach 1 derive its own.

Both orientations are run, and read beside experiment 2 they separate the two channels behind the
result. The two certificates: primal, complete, rate 0.902151, induced common in 13 polytopes,
ladder ΓX = 1.7628 over 7 rungs; dual, co-complete, rate 0.869447, induced common convex, ladder
ΓX = 2.2139 over 7 rungs.

| | primal (1) | primal (2) | dual (1) | dual (2) |
| :--- | --: | --: | --: | --: |
| piece gap | 0.460 | | 0.201 | |
| cells | 308 | 380 | 600 | **458** |
| Σ polytopes | 2348 | **1265** | 3128 | 3179 |
| mean polytopes / cell | 7.62 | **3.33** | 5.21 | 6.94 |
| max polytopes / cell | 128 | **16** | 34 | 35 |
| mean facets / cell | 30.56 | **13.36** | 20.85 | 27.80 |
| max facets / cell | 510 | **64** | 134 | 142 |
| build (s) | 11.89 | **2.46** | 7.21 | 6.92 |

Read beside experiment 2, the two orientations separate the two channels, and the polytope row says
which one is open. On the primal graph approach 2 builds 23 % more cells out of 46 % fewer polytopes,
and is 4.8× faster. On the dual graph it builds 24 % fewer cells, but each is more complex, and the
two totals land within 2 % of one another: 3128 polytopes against 3179. The build times land there
too, 7.21 s against 6.92 s. Nothing is gained because there is nothing to gain, which is a cleaner
statement than the one earlier runs of this arm supported, when the measurement protocol let the
verdict wander between 1.9× slower and 1.3× faster.

The piece gap is not a property of the system alone. The same dynamics and the same template give
0.46 on the primal graph and 0.20 on the dual, because the solver is answering a different question
on each. Imposing diversity instead, by giving each node its own template, does carry across and
removes the dual's advantage entirely: the cell ratio of approach 2 then goes from 0.51 to between
2.25 and 3.29. That is worth knowing before reading a gap as a property of the problem.

## References

[1] E. Aydin Gol, X. Ding, M. Lazar and C. Belta. Finite bisimulations for switched linear systems.
*IEEE Transactions on Automatic Control*, 59(12):3122–3134, 2014.

[2] D. Angeli, N. Athanasopoulos, R. M. Jungers and M. Philippe. Path-complete graphs and common
Lyapunov functions. In *Proceedings of the 20th International Conference on Hybrid Systems:
Computation and Control (HSCC)*, pages 81–90, 2017.
