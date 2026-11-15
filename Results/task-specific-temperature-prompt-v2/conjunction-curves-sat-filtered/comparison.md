SAT-filtered conjunction curves — deepseek-v4-pro

Existing outputs only; no model/API calls. Parse, allowed APs, and individual SAT precede semantic clustering. Rank by frequency, breaking ties by earliest sampled member. Conjoin min(k,m) distinct representatives. All probabilities use all selected prompts. With m=0, abstain; no empty conjunction is evaluated.

T=0.7: 40 prompts, 400 saved outputs, 393 retained candidate samples, 91 distinct semantic clusters summed across prompts.
Excluded before clustering: 6 invalid/incomplete, 0 disallowed-AP, 1 individually UNSAT.

| k | P(good) | P(not guaranteed) | P(abstain) | P(unsat) | P(top-1 rescue) | Semantic rescue | Included representatives (total) | Prompts with ≥k clusters |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | 62.5% | 37.5% | 0.0% | 0.0% | 0.0% | 15.0% | 40 | 40 |
| 2 | 75.0% | 25.0% | 0.0% | 0.0% | 12.5% | 20.0% | 59 | 19 |
| 3 | 77.5% | 22.5% | 0.0% | 0.0% | 15.0% | 22.5% | 72 | 13 |
| 4 | 77.5% | 22.5% | 0.0% | 0.0% | 15.0% | 22.5% | 79 | 7 |
| 5 | 80.0% | 20.0% | 0.0% | 0.0% | 17.5% | 25.0% | 83 | 4 |

Semantic rescue: Good(Φ_k) and none of its included representatives is equivalent to the reference. Top-1 rescue: Good(Φ_k) and the top representative is not sufficient. Sufficient means SAT(candidate) and candidate implies reference, over unconstrained infinite traces.

Available candidate samples include duplicate generations after filtering; available semantic clusters are distinct formulas modulo semantic equivalence, counted separately for each prompt. The included-representative count sums min(k,m) over all prompts. Good + not guaranteed + abstain + unsat = 100% in every row.
