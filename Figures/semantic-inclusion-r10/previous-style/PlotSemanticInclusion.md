# Non-trivial semantic inclusion figure

Run from Terminal; Python 3.11+ and Matplotlib are required. On this machine:

```sh
/Users/abayat/miniforge3/bin/python -m venv /tmp/ltl-plot-env
/tmp/ltl-plot-env/bin/pip install -r /Users/abayat/Documents/PhD/Code/Dio/Dionysos.jl_1/PlotSemanticInclusion_requirements.txt
/tmp/ltl-plot-env/bin/python /Users/abayat/Documents/PhD/Code/Dio/Dionysos.jl_1/PlotSemanticInclusion.py
```

The script reads the completed R=10 full-benchmark conjunction analysis, recounts
`satisfiable && implies_reference` from every prompt, checks all ten summary
rows, and plots the exact ratios without intermediate percentage rounding.
No model requests or controller synthesis occur.

Outputs in `Figures/semantic-inclusion-r10/`:

- `semantic_inclusion.pdf`, `.svg`, `.png`: detailed 65–85% vertical scale.
- `semantic_inclusion_full_scale.pdf`, `.svg`, `.png`: 0–100% scale.
- `semantic_inclusion_data.json`: all exact counts, rates, metadata and hashes.
- `caption.txt`: caption with model, temperature, repetitions and denominator.

Source: `Results/deepseek-v4-pro-task-specific-full-T0.7-K10-conjunctions/summary.json`.
Override paths with `--input /path/to/summary.json --output /path/to/figures`.

Interpretation: DeepSeek V4 Pro, T=0.7, R=10, all 4,736 active NL prompts in the
task-specific benchmark. The denominator includes the three no-candidate cases
and all unsatisfiable conjunctions. Results are prompt-weighted, rather than
averaged uniformly over the 1,126 ground-truth semantic classes.

The event is SAT(Phi_k) AND Phi_k implies the reference. Its language inclusion
need not be strict and does not require equivalence. Individually invalid and
unsatisfiable candidates were removed before clustering. Frequency ranks the
remaining semantic classes, and Phi_k contains min(k,m(p)) representatives.
Consequently k is a cluster budget, not a claim that every prompt has k classes.
This is a logical metric, independent of controller issuance.

Verified endpoints: k=1: 3439/4736 = 72.6140202703%;
k=10: 3894/4736 = 82.2212837838%. Difference: 9.6072635135 percentage points.
