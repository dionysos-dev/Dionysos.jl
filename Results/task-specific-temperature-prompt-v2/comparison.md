Temperature pilot — deepseek-v4-pro

40 fixed development prompts × 10 samples × 3 temperature(s). Complete: true.

| T | Complete prompts | Valid LTL | Mean clusters | Sample accuracy | Modal accuracy | Correct appears | Sufficient appears | Mean output tokens |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| 0.3 | 40 | 98.5% | 1.73 | 47.2% | 45.0% | 52.5% | 70.0% | 31.54 |
| 0.7 | 40 | 98.5% | 2.27 | 46.0% | 47.5% | 55.0% | 80.0% | 31.79 |
| 1.0 | 40 | 96.5% | 2.52 | 44.8% | 47.5% | 57.5% | 75.0% | 32.29 |

Sufficiency means satisfiable candidate ⇒ reference under unconstrained infinite-trace LTL. Ties use the cluster containing the earliest sample; JSON/CSV also report tie rate and lower/upper modal accuracy. All-invalid groups count as incorrect and have zero clusters.

Candidates are filtered for individual satisfiability before semantic clustering and modal selection. Output-token cost still includes excluded candidates; per-prompt JSON reports their counts.

Token cost is measured in API-reported output tokens, including invalid generations. Missing usage is excluded from that average and its coverage is reported. Sample accuracy includes every scheduled sample. API failures are recorded separately. The pilot is for development selection; it does not estimate held-out test performance.

Excluded 0 benchmark records with malformed references or propositions outside the five allowed colors. Their IDs and reasons are saved in manifest.json; the source benchmark is unchanged.
