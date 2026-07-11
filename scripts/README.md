# Case-study scripts

Runnable driver scripts that exercise Dionysos on concrete examples — used to reproduce the
figures of the slides and papers presenting Dionysos and its algorithms.

These are **not** library code (that lives in [`../src/`](../src/), the `Dionysos` package) and
they are **not** the problem definitions (those live in [`../problems/`](../problems/), which each
script `include`s). This folder was previously named `utils/`; it was renamed to `scripts/` to
avoid confusion with the `src/utils/` library module (`Dionysos.Utils`).

## Organization

One folder per example (`PathPlanning/`, `Pendulum/`, `Thermostat/`, …), mirroring the layout of
[`../problems/`](../problems/). Thematic collections (`CDC2024/`, `PCLF/`,
`BisimulationQuotient/`, `Settings/`, `TrajectoryCertificationOptimizer/`) group scripts by topic
rather than by a single problem. Superseded experiments are kept under `Deprecated/`.

See the comments in the first few lines of each file for what it does.
