# N2 — At a bounded per-node budget, the common construction does not exist

**The safest result in the folder.** Analytic, 2-D, no optimizer enters the claim. It is already
implemented in [`../../experiments/memory_vs_geometry.jl`](../../experiments/memory_vs_geometry.jl);
what remains is writing it up as the headline rather than as one experiment among eight.

## What we want to show

> Fix the per-node geometric budget at four facets. There are planar, two-mode, asymptotically stable
> switched systems for which **no** common polyhedral Lyapunov function of that complexity certifies,
> while a path-complete family at the *same* per-node budget certifies at the optimal rate. The
> predecessor returns no abstraction; we return a quotient.

This is `rem:two-regimes` (ii) instantiated, and `cor:bounded-local-forces-memory` made concrete.

## The construction

$$A_1 = \rho\,R(2\pi/q), \qquad A_2 = \rho\,R(-2\pi/q)$$

on the cycle graph $\mathbb{Z}_q$ with edges $(s, s{+}1, 1)$ and $(s, s{-}1, 2)$, carrying

$$V_s(x) = \|R^{-s}x\|_\infty \qquad \text{— four facets per node, for every } q.$$

Along $(s, s{+}1, 1)$, $V_{s+1}(A_1x) = \|R^{-(s+1)}\rho R x\|_\infty = \rho\|R^{-s}x\|_\infty =
\rho V_s(x)$, and symmetrically for mode 2. So the cycle certifies at rate exactly $\rho$ at four
facets per node, for every $q$.

A *common* function with unit ball $B$ must satisfy $R^{\pm1}B \subseteq (\gamma/\rho)B$; as
$\gamma \to \rho$ this forces $B$ invariant under the order-$q$ rotation group, and a polytopic $B$
then needs at least $q$ facets. Restricted to four, the achievable rate is bounded away from $\rho$:
for $q = 5$ the gauge of $R_{72°}(1,1)$ is 1.260, so the best rate is $\approx 1.26\rho$ and any
$\rho > 1/1.26 \approx 0.794$ defeats it.

Measured: at 4 facets/node, one node gives rate $\ge 1$ and **no quotient**; the $\mathbb{Z}_q$ cycle
gives rate $= \rho$ and **1 797 cells**.

## Why this is the result to lead with

- **It is a proof, not a search.** Nothing here can be defeated by a referee giving the baseline a
  larger search budget — which is what killed two candidate separations in
  [`../../notes/paper-experiments.md`](../../notes/paper-experiments.md) §E5.
- **It is planar and two-mode**, so it costs nothing to reproduce.
- **It is a separation in kind**, not a cost ratio that has to survive a timing argument.

## What we will try

1. **Write it up as the primary experiment**, with the analytic argument in the body and
   `verify_pclf_rate` (sampling, no solver) as the numerical check.
2. **Report the single-node side honestly**: it is searched over many template orientations, and the
   search budget must be stated, because under-searching it would manufacture the separation.
3. **Tabulate over $q$**: for each $q$, the four-facet common's best achievable rate, the cycle's rate
   $\rho$, and the threshold $\rho$ above which the separation holds. A table over $q = 4 \ldots 8$
   makes the mechanism visible rather than anecdotal.
4. **State the limitation in the same breath** — see below. A paper that marks its own boundary is
   much harder to attack.

## The limitation that must be stated

This is a **construction, not a property of typical systems**. A separation was searched for over 24
random planar systems and 20 random systems in dimensions 3 and 4, and **none was found**: for
generic low-dimensional systems the common certificate already attains the JSR to five decimals, so
$\mathfrak F(\Sigma)$ is small and there is nothing to redistribute. Gol–Lazar–Belta's own system is
an instance — one node, two nodes and four nodes all certify at identical rates at every template
budget tried.

Pair this with **N3**, which turns that limitation into a usable rule: the length of the
spectrum-maximizing product says in advance whether a system is in the regime where memory can help.

## What would falsify it

`verify_pclf_rate` returning anything above $\rho$ would mean the analytic certificate is wrong. The
single-node side certifying at rate $< 1$ at four facets, once searched properly, would mean the
geometric lower-bound argument is wrong.
