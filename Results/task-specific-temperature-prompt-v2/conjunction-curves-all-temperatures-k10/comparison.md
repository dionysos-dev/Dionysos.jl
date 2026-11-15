SAT-filtered conjunction curves — deepseek-v4-pro

Existing outputs only; no model/API calls. Parse, allowed APs, and individual SAT precede semantic clustering. Rank by frequency, breaking ties by earliest sampled member. Conjoin min(k,m) distinct representatives. All probabilities use all selected prompts. With m=0, abstain; no empty conjunction is evaluated.

T=0.3: 40 prompts, 400 saved outputs, 394 retained candidate samples, 69 distinct semantic clusters summed across prompts.
Excluded before clustering: 6 invalid/incomplete, 0 disallowed-AP, 0 individually UNSAT.

| k | P(good) | P(not guaranteed) | P(abstain) | P(unsat) | P(top-1 rescue) | Semantic rescue | Included representatives (total) | Prompts with ≥k clusters |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | 57.5% | 42.5% | 0.0% | 0.0% | 0.0% | 12.5% | 40 | 40 |
| 2 | 65.0% | 35.0% | 0.0% | 0.0% | 7.5% | 12.5% | 56 | 16 |
| 3 | 67.5% | 32.5% | 0.0% | 0.0% | 10.0% | 15.0% | 62 | 6 |
| 4 | 70.0% | 30.0% | 0.0% | 0.0% | 12.5% | 17.5% | 65 | 3 |
| 5 | 70.0% | 30.0% | 0.0% | 0.0% | 12.5% | 17.5% | 67 | 2 |
| 6 | 70.0% | 30.0% | 0.0% | 0.0% | 12.5% | 17.5% | 68 | 1 |
| 7 | 70.0% | 30.0% | 0.0% | 0.0% | 12.5% | 17.5% | 69 | 1 |
| 8 | 70.0% | 30.0% | 0.0% | 0.0% | 12.5% | 17.5% | 69 | 0 |
| 9 | 70.0% | 30.0% | 0.0% | 0.0% | 12.5% | 17.5% | 69 | 0 |
| 10 | 70.0% | 30.0% | 0.0% | 0.0% | 12.5% | 17.5% | 69 | 0 |

T=0.7: 40 prompts, 400 saved outputs, 393 retained candidate samples, 91 distinct semantic clusters summed across prompts.
Excluded before clustering: 6 invalid/incomplete, 0 disallowed-AP, 1 individually UNSAT.

| k | P(good) | P(not guaranteed) | P(abstain) | P(unsat) | P(top-1 rescue) | Semantic rescue | Included representatives (total) | Prompts with ≥k clusters |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | 62.5% | 37.5% | 0.0% | 0.0% | 0.0% | 15.0% | 40 | 40 |
| 2 | 75.0% | 25.0% | 0.0% | 0.0% | 12.5% | 20.0% | 59 | 19 |
| 3 | 77.5% | 22.5% | 0.0% | 0.0% | 15.0% | 22.5% | 72 | 13 |
| 4 | 77.5% | 22.5% | 0.0% | 0.0% | 15.0% | 22.5% | 79 | 7 |
| 5 | 80.0% | 20.0% | 0.0% | 0.0% | 17.5% | 25.0% | 83 | 4 |
| 6 | 80.0% | 20.0% | 0.0% | 0.0% | 17.5% | 25.0% | 87 | 4 |
| 7 | 80.0% | 20.0% | 0.0% | 0.0% | 17.5% | 25.0% | 89 | 2 |
| 8 | 80.0% | 20.0% | 0.0% | 0.0% | 17.5% | 25.0% | 90 | 1 |
| 9 | 80.0% | 20.0% | 0.0% | 0.0% | 17.5% | 25.0% | 91 | 1 |
| 10 | 80.0% | 20.0% | 0.0% | 0.0% | 17.5% | 25.0% | 91 | 0 |

T=1.0: 40 prompts, 400 saved outputs, 386 retained candidate samples, 101 distinct semantic clusters summed across prompts.
Excluded before clustering: 13 invalid/incomplete, 1 disallowed-AP, 0 individually UNSAT.

| k | P(good) | P(not guaranteed) | P(abstain) | P(unsat) | P(top-1 rescue) | Semantic rescue | Included representatives (total) | Prompts with ≥k clusters |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | 57.5% | 42.5% | 0.0% | 0.0% | 0.0% | 10.0% | 40 | 40 |
| 2 | 70.0% | 30.0% | 0.0% | 0.0% | 12.5% | 17.5% | 63 | 23 |
| 3 | 75.0% | 25.0% | 0.0% | 0.0% | 17.5% | 17.5% | 79 | 16 |
| 4 | 75.0% | 25.0% | 0.0% | 0.0% | 17.5% | 17.5% | 85 | 6 |
| 5 | 77.5% | 22.5% | 0.0% | 0.0% | 20.0% | 20.0% | 90 | 5 |
| 6 | 80.0% | 20.0% | 0.0% | 0.0% | 22.5% | 22.5% | 94 | 4 |
| 7 | 80.0% | 20.0% | 0.0% | 0.0% | 22.5% | 22.5% | 97 | 3 |
| 8 | 80.0% | 20.0% | 0.0% | 0.0% | 22.5% | 22.5% | 99 | 2 |
| 9 | 80.0% | 20.0% | 0.0% | 0.0% | 22.5% | 22.5% | 100 | 1 |
| 10 | 80.0% | 20.0% | 0.0% | 0.0% | 22.5% | 22.5% | 101 | 1 |

Semantic rescue: Good(Φ_k) and none of its included representatives is equivalent to the reference. Top-1 rescue: Good(Φ_k) and the top representative is not sufficient. Sufficient means SAT(candidate) and candidate implies reference, over unconstrained infinite traces.

Available candidate samples include duplicate generations after filtering; available semantic clusters are distinct formulas modulo semantic equivalence, counted separately for each prompt. The included-representative count sums min(k,m) over all prompts. Good + not guaranteed + abstain + unsat = 100% in every row.
