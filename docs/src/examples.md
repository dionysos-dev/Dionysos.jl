# Examples

Each example is one complete run: a model, an abstraction, a synthesized controller, and an animated
simulation of the closed loop. They are ordered by what they add, so reading them in order introduces
one new idea at a time.

Start with [Getting started](generated/getting_started.md), which walks the same five steps in detail
on a pendulum.

## What each example covers

The **System** column classifies the plant — that is usually what you are scanning for. *Switched*
means the mode is chosen by the control input; *hybrid* means it changes when the state enters a
guard.

| Example | System | Specification | Front-end features |
| :--- | :--- | :--- | :--- |
| [Getting started](generated/getting_started.md) | 2-D continuous, nonlinear | reach-avoid | `∂`, `start`, `final`, `∉` |
| [Path planning](generated/path_planning.md) | 3-D continuous, nonlinear | reach-avoid | `∉` over several coordinates, `@expression` |
| [Unicycle robot](generated/unicycle_robot.md) | 3-D discrete-time, nonlinear | reach-avoid | `Δ`, custom `growthbound_map` |
| [DC-DC converter](generated/dcdc_converter.md) | 2-D continuous, switched by the input | safety | `Always`, dynamics as a Julia function, `time_domain`, two solvers on one problem |
| [Thermostat](generated/thermostat.md) | 1-D hybrid, 2 modes with guards | reach | `@mode`, `add_transition!`, guards, per-mode bounds |

## Not covered yet

The front-end does more than these four pages show. The gaps are listed rather than hidden: each row
is a page worth writing, and the feature itself works today.

| Missing example | Specification | Front-end features |
| :--- | :--- | :--- |
| a timed hybrid system | reach within a time window | clocks, reset maps, timed specifications |
| an LTL task | co-safe LTL | `Label`, `@specification` |
| a stabilisation task | reach-and-stay | `EventuallyAlways` |
| a deadline | reach within a horizon | the `horizon` attribute |
| a set-valued initial condition | any | `Start` over a region rather than a point |
| a mode controlled only by switching | any | a hybrid mode with no continuous input — needs zero-dimensional input spaces in `Mapping`/`Symbolic` |

## Beyond the front-end

The [Solver families](generated/dcdc_converter.md) examples drive one optimizer each through
MathOptInterface directly. They exist both to document those algorithms and to reach the ones the
JuMP front-end cannot express — piecewise-affine systems, ellipsoidal cells, SDP-based local
controllers.

To add an example, see [Adding an example](@ref "Adding an example").
