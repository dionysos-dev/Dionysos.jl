# Full-benchmark SAT-filtered conjunction curves

deepseek-v4-pro, T=0.7, 4736 prompts, 47360 saved responses. Existing outputs only; no model/API calls.

Parse → allowed APs → individual SAT → discard failures → semantic equivalence clusters → rank by frequency (earliest query breaks ties) → conjoin min(k,m) representatives. All rates use all prompts; zero clusters means abstention.

| k | P_good | P_not-guaranteed | P_abstain | P_unsat | Top-1 rescue | Semantic rescue | Representatives used |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | 72.61% | 27.32% | 0.06% | 0.00% | 0.00% | 12.73% | 4733 |
| 2 | 78.06% | 21.58% | 0.06% | 0.30% | 5.60% | 12.73% | 6666 |
| 3 | 80.36% | 19.17% | 0.06% | 0.40% | 7.92% | 13.22% | 7725 |
| 4 | 81.08% | 18.37% | 0.06% | 0.49% | 8.68% | 13.43% | 8384 |
| 5 | 81.78% | 17.63% | 0.06% | 0.53% | 9.40% | 13.96% | 8804 |
| 6 | 82.05% | 17.34% | 0.06% | 0.55% | 9.69% | 14.15% | 9082 |
| 7 | 82.14% | 17.25% | 0.06% | 0.55% | 9.78% | 14.21% | 9239 |
| 8 | 82.18% | 17.21% | 0.06% | 0.55% | 9.82% | 14.25% | 9326 |
| 9 | 82.20% | 17.19% | 0.06% | 0.55% | 9.84% | 14.27% | 9362 |
| 10 | 82.22% | 17.17% | 0.06% | 0.55% | 9.86% | 14.29% | 9374 |

Good means SAT(Φ_k) and Φ_k implies the reference, over unconstrained infinite traces. Unsatisfiable conjunctions do not count as good.

Top-1 rescue means Good(Φ_k) and the top representative is not sufficient. Semantic rescue means Good(Φ_k) and no included representative is equivalent to the reference. Sufficient means SAT(candidate) and candidate implies reference. Representatives used is the sum of min(k,m) across all prompts.

Excluded candidates: OrderedCollections.OrderedDict{String, Any}("invalid_ltl" => 498, "disallowed_propositions" => 129, "unsatisfiable" => 24). Retained samples: 46709; semantic clusters summed over prompts: 9374.

These are logical metrics; controller synthesis is not part of this analysis. All counts, cluster members, formulas, per-prompt outcomes, source SHA-256, and policy details are saved in summary.json.
