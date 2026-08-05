# Research scripts

Driver scripts reproducing the figures and experiments of our papers and internal studies — our own
simulations. These are **not** user-facing examples (those live in [`../examples/`](../examples/))
and **not** library code ([`../src/`](../src/)).

## Contents

| Folder | Associated work |
| :--- | :--- |
| [`CDC2024/`](CDC2024/) | Lazy ellipsoidal abstraction / local transitions — CDC 2024 paper. |
| [`BisimulationQuotient/`](BisimulationQuotient/) | PCLF-based bisimulation quotients — HSCC 2027 (see the folder README). |
| [`PCLF/`](PCLF/) | Path-complete Lyapunov function experiments (PCLF ↔ CLF, polyhedral/quadratic pieces). |
| [`TrajectoryCertificationOptimizer/`](TrajectoryCertificationOptimizer/) | Trajectory-certification studies. |
| [`Settings/`](Settings/) | Feature / configuration demos (distributed build, periodic mapping, hybrid abstraction). |

Superseded experiments are kept under [`Deprecated/`](Deprecated/).

Scripts load their problem definitions from [`../problems/`](../problems/) via `pathof(Dionysos)`, so
they run from any working directory.
