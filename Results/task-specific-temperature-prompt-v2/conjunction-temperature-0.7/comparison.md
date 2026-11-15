Frequency-ranked cluster conjunctions — T=0.7

40 prompts; 400 saved samples. No new LLM calls. Clusters contain valid outputs only; ties use the earliest sampled member.

| k | Eligible prompts | Satisfiable | Implies reference | Logical good | Not guaranteed | Unsatisfiable | Semantic rescue |
|---|---:|---:|---:|---:|---:|---:|---:|
| 1 | 40 | 100.0% | 62.5% | 62.5% | 37.5% | 0.0% | 15.0% |
| 2 | 40 | 100.0% | 75.0% | 75.0% | 25.0% | 0.0% | 20.0% |
| 3 | 40 | 97.5% | 80.0% | 77.5% | 20.0% | 2.5% | 22.5% |
| 4 | 40 | 97.5% | 80.0% | 77.5% | 20.0% | 2.5% | 22.5% |
| 5 | 40 | 97.5% | 82.5% | 80.0% | 17.5% | 2.5% | 25.0% |
| 6 | 40 | 97.5% | 82.5% | 80.0% | 17.5% | 2.5% | 25.0% |
| 7 | 40 | 97.5% | 82.5% | 80.0% | 17.5% | 2.5% | 25.0% |
| 8 | 40 | 97.5% | 82.5% | 80.0% | 17.5% | 2.5% | 25.0% |
| 9 | 40 | 97.5% | 82.5% | 80.0% | 17.5% | 2.5% | 25.0% |
| 10 | 40 | 97.5% | 82.5% | 80.0% | 17.5% | 2.5% | 25.0% |

Missing-cluster policy: saturate. Strict rows estimate probabilities conditional on at least k distinct clusters; their populations change with k. Saturated rows conjoin min(k,m) clusters. Prompts with zero valid clusters are unavailable under both policies.

Semantic rescue definition: no_individual_equivalent. The selected default is SAT(Φ_k) and Φ_k implies the reference, while no included cluster is individually equivalent to the reference. This includes stronger sufficient formulas, so rescue can be nonzero at k=1.

Raw implication includes vacuous truth for unsatisfiable conjunctions. Logical good explicitly requires satisfiability. Logical good + not guaranteed + unsatisfiable = 100% for every nonempty row.

All checks use standard LTL over unconstrained infinite traces. Every ranked representative, frequency, and Φ_k is saved in summary.json. A second table with the alternative missing-cluster policy is in summary-alternative-policy.json.
