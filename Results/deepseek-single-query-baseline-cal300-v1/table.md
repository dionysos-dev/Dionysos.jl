All percentages use uniform semantic-class means; 300 calibration classes and 826 test classes.
The baseline uses query 1 exactly. No new LLM calls.

| Method | Statistical guarantee | α | k | Test risk (%) | Issuance (%) |
|---|---|---:|---:|---:|---:|
| Single LLM query (equivalence) | No | -- | -- | 34.52 | 88.54 |
| Single LLM query (language inclusion) | No | -- | -- | 22.68 | 88.54 |
| NL2RCtrl | Yes | 0.10 | infinity | 0.00 | 0.00 |
| NL2RCtrl | Yes | 0.15 | 3 | 13.49 | 84.55 |
| NL2RCtrl | Yes | 0.20 | 2 | 16.40 | 86.33 |
| NL2RCtrl | Yes | 0.25 | 1 | 22.89 | 89.54 |
