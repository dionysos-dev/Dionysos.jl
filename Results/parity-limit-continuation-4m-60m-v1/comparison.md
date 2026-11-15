# Ground-truth audit after larger parity limits

Production algorithms unchanged. Six resolved pilot cases reused; all 159 remaining parity cases attempted sequentially with cooldown breaks.

| Outcome | Semantic classes | % of classes | Records | % of records |
|---|---:|---:|---:|---:|
| controller_issued | 992 | 88.10% | 1049 | 88.60% |
| plant_unrealizable | 38 | 3.37% | 38 | 3.21% |
| logically_inconsistent | 0 | 0.00% | 0 | 0.00% |
| computational_failure | 96 | 8.53% | 97 | 8.19% |

Computational failures remain unknown controller existence. All successful parity controllers passed both native and separate Julia infinite-horizon verification. Original audit reports are preserved.
