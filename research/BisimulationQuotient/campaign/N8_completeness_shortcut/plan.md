# N8 — Stop handicapping route 1

Not an experiment so much as a fairness fix, but it must land before any timing table is published,
because route 1 currently wins E8 *while handicapped* and a referee will notice.

## The problem

`th:induced_common` gives three cases. In general
$V^\star(x) = \min_{Q\in S_{\mathrm{obs}}}\max_{s\in Q} V_s(x)$ over the observer graph. But:

- if $\mathcal G$ is **complete**, $V^\star = \min_s V_s$;
- if $\mathcal G$ is **co-complete**, $V^\star = \max_s V_s$.

`build_common_lyapunov` ([`pclf.jl:914`](../../../../src/utils/pclf.jl#L914)) always takes the general
path: it runs `build_observer_graph` and wraps the result in an `ObserverCLFPiece`, whose
`get_sublevel_set` builds $\bigcup_Q \bigcap_{s\in Q}\mathcal P^{(s)}_\Gamma$ by repeated
`set_difference_decompose`.

**On E8 that is strictly wasteful.** The graph is De Bruijn order 1 over two modes, hence complete, so
$V^\star = \min_s V_s$ applies. And the two pieces are *nested* — measured
$\mathrm{vol}(P_1 \cap P_2) = \mathrm{vol}(P_1)$ exactly, so $P_1 \subseteq P_2$ — hence
$\min_s V_s = V_2$: a **single 8-row polytope, 16 facets**. What the code builds instead:

| $\gamma$ | general observer path | what the theorem permits |
| --: | :--- | :--- |
| 1.0 | 21 pieces, 94 facets, 24.0 ms | 1 piece, 16 facets |
| 0.5 | 17 pieces, 80 facets, 14.9 ms | 1 piece, 16 facets |

## What we will try

1. **Dispatch on completeness in `build_common_lyapunov`.** `PCLF.is_complete` already exists (it is
   used in [`../../examples/incomplete_showcase.jl`](../../examples/incomplete_showcase.jl)); add the
   co-complete case too.
2. **Re-run the E8 arm** and report route 1's new build time. Expect it to drop, making route 2 look
   worse — publish that.
3. **Check it changes nothing else.** The certified winning set must be identical; if it is not,
   either `th:induced_common` or its implementation is wrong, and that is a far more important
   finding than any timing.
4. **Note the interaction with N1.** The dual De Bruijn graph of N1 is *incomplete*, so the shortcut
   does not apply there and N1's 1.74× is unaffected. Confirm this explicitly so the two results do
   not appear to contradict.

## Expected size of the effect

Per sublevel set the saving is large in ratio (24 ms against a direct H-polytope construction), but
sublevel sets are built once per level in `build_sublevel_sequence`, so at 20–50 levels the total is
0.3–1.2 s against route 1's 168.7 s — **under 1 %**. So this will not overturn any conclusion. It is
worth doing because it removes an objection, not because it moves a number.

## What would falsify the reasoning

`is_complete` returning false for the De Bruijn order-1 graph, or the nesting $P_1 \subseteq P_2$ not
holding at other $\gamma$ — both are one-line checks and should be done before the code change.
