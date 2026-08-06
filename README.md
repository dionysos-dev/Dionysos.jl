# Dionysos
<picture>
  <source srcset="assets/images/logo-dark.png" media="(prefers-color-scheme: dark)">
  <img src="assets/images/logo.png"  height="240">
</picture>

| **Documentation** | **Build Status** | **Citation** | **Community** |
|:-----------------:|:----------------:|:------------:|:-------------:|
| [![][docs-latest-img]][docs-latest-url] [![][docs-stable-img]][docs-stable-url] | [![Build Status][build-img]][build-url] [![Codecov][codecov-img]][codecov-url] [![Aqua QA][aqua-img]][aqua-url] [![PkgEval][pkgeval-img]][pkgeval-url] | [![Paper DOI][paper-img]][paper-url] | [![Slack][slack-img]][slack-url] |


[docs-stable-img]: https://img.shields.io/badge/docs-stable-blue.svg
[docs-latest-img]: https://img.shields.io/badge/docs-latest-blue.svg
[docs-stable-url]: https://dionysos-dev.github.io/Dionysos.jl/stable
[docs-latest-url]: https://dionysos-dev.github.io/Dionysos.jl/dev

[build-img]: https://github.com/dionysos-dev/Dionysos.jl/actions/workflows/ci.yml/badge.svg?branch=master
[build-url]: https://github.com/dionysos-dev/Dionysos.jl/actions?query=workflow%3ACI
[codecov-img]: https://codecov.io/github/dionysos-dev/Dionysos.jl/coverage.svg
[codecov-url]: https://app.codecov.io/github/dionysos-dev/Dionysos.jl
[aqua-img]: https://juliatesting.github.io/Aqua.jl/dev/assets/badge.svg
[aqua-url]: https://github.com/JuliaTesting/Aqua.jl
[pkgeval-img]: https://juliaci.github.io/NanosoldierReports/pkgeval_badges/D/Dionysos.svg
[pkgeval-url]: https://juliaci.github.io/NanosoldierReports/pkgeval_badges/report.html

[paper-img]: https://img.shields.io/badge/Paper-10.21105%2Fjcon.00160-blue.svg
[paper-url]: https://doi.org/10.21105/jcon.00160

<!-- Software (Zenodo) DOI badge. Dionysos has no Zenodo record yet: create one by enabling the
     GitHub–Zenodo integration for this repository and making a release, then replace ZENODO_ID
     below with the *concept* DOI (the one that always resolves to the latest release) and add
     `[![Software DOI][zenodo-img]][zenodo-url]` to the Citation column above.
[zenodo-img]: https://img.shields.io/badge/Software-10.5281%2Fzenodo.ZENODO__ID-blue.svg
[zenodo-url]: https://doi.org/10.5281/zenodo.ZENODO_ID
-->

[slack-img]: https://img.shields.io/badge/chat-%23control-4A154B?logo=slack&logoColor=white
[slack-url]: https://julialang.org/slack/
[slack-control-url]: https://julialang.slack.com/archives/CKH1UTZT9

**Dionysos** is a [Julia](https://julialang.org/) framework for **correct-by-construction controller
synthesis** through symbolic (abstraction-based) control. Its guiding vision is
**Control as a Service (CaaS)**: making certified controller design an automated, on-demand capability
rather than a bespoke, months-long expert effort. It is the software of the ERC project
[Learning to Control](https://perso.uclouvain.be/raphael.jungers/content/erc-consolidator-grant) (L2C).

The following article showcases the basic functionality, highlighting some of the key design choices:

> [**Dionysos.jl: a Modular Platform for Smart Symbolic Control**][paper-url]
> Julien Calbert, Adrien Banse, Benoît Legat, Raphaël M. Jungers.
> *JuliaCon Proceedings*, 6(66), 160, 2024.

See [How to cite](#how-to-cite) below.

## Why Dionysos — Control as a Service (CaaS)

Controlling a complex system traditionally means a team of experts hand-crafting an ad hoc controller
over months, at significant cost. Dionysos pursues a different paradigm — **Control as a Service** —
turning controller design into an automated, on-demand pipeline that is accessible even to teams
without a dedicated control or IT department:

> **describe the system → select the specification → pick a solver → obtain a controller together with
> a formal certificate.**

You bring the model and the goal; Dionysos returns the controller and its certificate. Under the hood,
the system is *abstracted* into a finite-state automaton by discretizing its variables; a controller
is synthesized on that finite object with graph algorithms and then *concretized* back to the original
system with a formal guarantee. Dionysos is an **ecosystem, not a single algorithm**: every solver is
a [MathOptInterface](https://jump.dev/MathOptInterface.jl) optimizer driven through
[JuMP](https://jump.dev/JuMP.jl), so a control task can be re-solved, compared, and benchmarked by
*swapping the solver* rather than rewriting the model.

It solves reach-avoid, safety, reach-and-stay, and co-safe LTL specifications with a catalog of
solvers: uniform grid abstraction (SCOTS-style), uniform and lazy ellipsoidal abstractions,
hybrid-system abstraction, a PCLF bisimulation quotient, discrete-automaton synthesis, and
optimization-based solvers (Bemporad–Morari, branch and bound).

## Quick start

Describe a system and a specification in a JuMP model, pick the solver by choosing the optimizer, and
read back the controller:

```julia
using Dionysos, JuMP, StaticArrays

model = Model(Dionysos.Optimizer)

# State and input variables, with the initial state as `start`.
@variable(model, x_low[i] <= x[i = 1:3] <= x_upp[i], start = x0[i])
@variable(model, -1 <= u[1:2] <= 1)

# Dynamics, written with the `∂` (derivative) operator ...
@constraint(model, ∂(x[1]) == u[1] * cos(α + x[3]) * sec(α))
@constraint(model, ∂(x[2]) == u[1] * sin(α + x[3]) * sec(α))
@constraint(model, ∂(x[3]) == u[1] * tan(u[2]))

# ... and the target set, written with `final`.
@constraint(model, final(x[1]) in MOI.Interval(3.0, 3.6))
@constraint(model, final(x[2]) in MOI.Interval(0.3, 0.9))

# Solver parameters: discretization grids and time step.
set_attribute(model, "time_step", 0.3)
set_attribute(model, "state_grid", Dionysos.Mapping.GridFree(x0, hx))
set_attribute(model, "input_grid", Dionysos.Mapping.GridFree(u0, hu))

optimize!(model)
concrete_controller = get_attribute(model, "concrete_controller")
```

[Getting started](https://dionysos-dev.github.io/Dionysos.jl/dev/generated/getting_started/) walks
through a complete run step by step, and the
[Path planning example](https://dionysos-dev.github.io/Dionysos.jl/dev/generated/path_planning/)
adds obstacles and a 3-D state space.

## Repository layout

| Path | Description |
| :--- | :--- |
| [`src/`](src/) | The `Dionysos` library (`Utils`, `System`, `Problem`, `Mapping`, `Symbolic`, `Optim`) and the JuMP front-end ([`Wrapper`](src/wrapper/README.md)). |
| [`ext/`](ext/) | Package extensions behind optional dependencies (Plots, Symbolics, Spot, CSV, …). |
| [`problems/`](problems/) | Reusable benchmark problem library (path planning, DC-DC, pendulum, …). |
| [`examples/`](examples/) | Runnable example drivers (user-facing), one folder per problem. |
| [`research/`](research/) | Our paper / experiment simulations (CDC 2024, HSCC 2027, …). |
| [`docs/`](docs/) | Documentation site and examples. |
| [`test/`](test/) | Test suite, mirroring the `src/` layout. |

## Installation

Install [Julia](https://julialang.org/downloads/) (≥ 1.10), then add Dionysos from the Julia REPL:

```julia
import Pkg
Pkg.add("Dionysos")
```

## How to cite

If you use this package in your work, please cite the paper describing Dionysos:

```bibtex
@article{Calbert2024,
  author       = {Julien Calbert and
                  Adrien Banse and
                  Beno{\^\i}t Legat and
                  Rapha{\"e}l M. Jungers},
  title        = {{Dionysos.jl}: a Modular Platform for Smart Symbolic Control},
  journal      = {JuliaCon Proceedings},
  volume       = {6},
  number       = {66},
  pages        = {160},
  year         = {2024},
  month        = dec,
  publisher    = {The Open Journal},
  issn         = {2642-4029},
  url          = {https://doi.org/10.21105/jcon.00160},
  doi          = {10.21105/jcon.00160}
}
```

To cite the *software release* you actually ran — which is what makes an experiment reproducible —
cite the archived version of the package as well, alongside the paper above.

## Documentation

The full documentation, including the manual, examples, and API reference, is available for the
[released version][docs-stable-url] and the [development version][docs-latest-url]. If Dionysos is
useful in your research, please see [How to cite](#how-to-cite).

## Community & support

Questions, ideas, and discussion are welcome in the **`#control`** channel on the
[JuliaLang Slack](https://julialang.org/slack/) — join the workspace, then open
[`#control`][slack-control-url]. For bug reports and feature requests, please open an
[issue](https://github.com/dionysos-dev/Dionysos.jl/issues).

## Contributing

Contributions are welcome. Please open an [issue](https://github.com/dionysos-dev/Dionysos.jl/issues)
to report a bug or discuss a feature, and see the Developer Docs in the documentation for the setup,
conventions, and Git workflow.

## Acknowledgements

This project has received funding from the European Research Council (ERC) under the European Union's
Horizon 2020 research and innovation programme under grant agreement No 864017 - L2C.
