"""
    Dionysos.Wrapper

The JuMP/MOI front-end: it compiles a JuMP model into a `MathematicalSystems` /
`HybridSystems` object plus a [`Dionysos.Problem.ProblemType`](@ref), then drives an existing
solver with them. It owns no control semantics of its own.

The pipeline is one-directional:

```
JuMP model ──parse──▶ ModelIR ──lower──▶ (system, problem) ──▶ solver ──▶ controller
```

[`ModelIR`](@ref) is a plain, dependency-free description of what the user wrote, which is why
this module lives in `src/` rather than in an extension: parsing needs no Symbolics, only
*compiling a dynamics expression into a callable* does (see [`AbstractDynamicsBackend`](@ref)).

See `src/wrapper/README.md` for the user guide.
"""
module Wrapper

import StaticArrays: SVector, SMatrix
import MathematicalSystems as MS
import MathOptInterface as MOI
import HybridSystems
import JuMP
import LazySets

using ..Utils
const UT = Utils

using ..System
const ST = System

using ..Problem
const PR = Problem

using ..Mapping
const MP = Mapping

using ..Optim
const OP = Optim
const OPDS = OP.DiscreteSystems

include("operators.jl")
include("variables.jl")
include("specification.jl")
include("modes.jl")
include("model_ir.jl")
include("dynamics_backend.jl")
include("solver_selection.jl")
# `optimizer.jl` defines the `Optimizer` type, so it precedes the files that add methods on it.
include("optimizer.jl")
include("parse.jl")
include("parse_scoped.jl")
include("lower_system.jl")
include("lower_problem.jl")
include("clock.jl")
include("lower_hybrid.jl")
include("results.jl")

end # module Wrapper
