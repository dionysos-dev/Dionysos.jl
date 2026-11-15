Frequency-ranked cluster conjunctions — deepseek-v4-pro, T=0.3

40 prompts; 400 saved samples. No new LLM calls. Candidates must be valid and individually satisfiable before semantic clustering; ties use the earliest sampled member.

Filtered 6 invalid candidates, 0 candidates with disallowed propositions, and 0 unsatisfiable candidates before clustering.

| k | Eligible prompts | Satisfiable | Implies reference | Logical good | Not guaranteed | Unsatisfiable | Semantic rescue | Top-1 rescue | Synergistic rescue |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | 40 | 100.0% | 57.5% | 57.5% | 42.5% | 0.0% | 12.5% | 0.0% | 0.0% |
| 2 | 40 | 100.0% | 65.0% | 65.0% | 35.0% | 0.0% | 12.5% | 7.5% | 0.0% |
| 3 | 40 | 100.0% | 67.5% | 67.5% | 32.5% | 0.0% | 15.0% | 10.0% | 0.0% |
| 4 | 40 | 100.0% | 70.0% | 70.0% | 30.0% | 0.0% | 17.5% | 12.5% | 2.5% |
| 5 | 40 | 100.0% | 70.0% | 70.0% | 30.0% | 0.0% | 17.5% | 12.5% | 2.5% |
| 6 | 40 | 100.0% | 70.0% | 70.0% | 30.0% | 0.0% | 17.5% | 12.5% | 2.5% |
| 7 | 40 | 100.0% | 70.0% | 70.0% | 30.0% | 0.0% | 17.5% | 12.5% | 0.0% |
| 8 | 40 | 100.0% | 70.0% | 70.0% | 30.0% | 0.0% | 17.5% | 12.5% | 0.0% |
| 9 | 40 | 100.0% | 70.0% | 70.0% | 30.0% | 0.0% | 17.5% | 12.5% | 0.0% |
| 10 | 40 | 100.0% | 70.0% | 70.0% | 30.0% | 0.0% | 17.5% | 12.5% | 0.0% |

Missing-cluster policy: saturate. Strict rows estimate probabilities conditional on at least k distinct clusters; their populations change with k. Saturated rows conjoin min(k,m) clusters. Prompts with zero eligible clusters are unavailable under both policies; eligible means valid and individually satisfiable.

Semantic rescue definition: no_individual_equivalent. The selected default is SAT(Φ_k) and Φ_k implies the reference, while no included cluster is individually equivalent to the reference. This includes stronger sufficient formulas, so rescue can be nonzero at k=1.

Top-1 rescue: Good(Φ_k) and the top-ranked representative is not individually sufficient. Synergistic rescue: Good(Φ_k) and none of the included representatives is individually sufficient. Here sufficient means SAT(ψ) and ψ implies the reference. Both rescue metrics are zero at k=1. Synergistic rescue can disappear when an individually sufficient cluster is added, even if the conjunction stays logically good.

Raw implication includes vacuous truth for unsatisfiable conjunctions. Logical good explicitly requires satisfiability. Logical good + not guaranteed + unsatisfiable = 100% for every nonempty row.

All checks use standard LTL over unconstrained infinite traces. Every ranked representative, frequency, and Φ_k is saved in summary.json. A second table with the alternative missing-cluster policy is in summary-alternative-policy.json.
