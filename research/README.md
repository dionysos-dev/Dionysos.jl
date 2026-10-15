# Research scripts

Driver scripts reproducing the figures and experiments of our papers and internal studies: our own
simulations. These are **not** user-facing examples (those live in [`../examples/`](../examples/))
and **not** library code ([`../src/`](../src/)).

## Contents

| Folder | Associated work |
| :--- | :--- |
| [`CDC2024/`](CDC2024/) | Lazy ellipsoidal abstraction and local transitions, CDC 2024 paper. |
| [`BisimulationQuotient/`](BisimulationQuotient/) | PCLF-based bisimulation quotients, HSCC 2027 (see the folder README). |
| [`HSCC2027/`](HSCC2027/) | Path-complete abstractions of switched systems: code and data of the experiments of the HSCC 2027 paper (see its README). |
| [`ECC2027/`](ECC2027/) | Velocity-controlled mechanical systems in symbolic control: the CodeOcean capsule of the ECC 2027 paper (see its README). |
| [`TrajectoryCertificationOptimizer/`](TrajectoryCertificationOptimizer/) | Trajectory-certification studies. |
| [`Settings/`](Settings/) | Feature and configuration demos (distributed build, periodic mapping, hybrid abstraction). |

Scripts load their problem definitions from [`../problems/`](../problems/) via `pathof(Dionysos)`, so
they run from any working directory.
