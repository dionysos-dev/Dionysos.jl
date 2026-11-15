# DeepSeek controller-aware risk and issuance

All percentages use uniform semantic-class means, matching calibration.
300 calibration classes; 826 test classes. No new LLM calls.

| Target alpha | Selected k | Test risk (%) | Issuance (%) | Guaranteed issuance (%) | Status |
|---:|---:|---:|---:|---:|---|
| 0.05 | infinity | 0.00 | 0.00 | 0.00 | complete_always_abstain |
| 0.10 | infinity | 0.00 | 0.00 | 0.00 | complete_always_abstain |
| 0.15 | 3 | 13.49 | 84.55 | 71.06 | complete |
| 0.20 | 2 | 16.40 | 86.33 | 69.94 | complete |
| 0.25 | 1 | 22.89 | 89.54 | 66.65 | complete |

Issuance uses the successful prefix product, not merely synthesis at selected k.
Test risk + guaranteed issuance = issuance. Computational failure causes abstention.
The historical 16.32% used the earlier solver/resources and is not inserted into this new evaluation.
This resplit has prior benchmark evaluation history; it is not a pristine held-out test.
