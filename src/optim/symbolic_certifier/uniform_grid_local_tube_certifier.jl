module UniformGridLocalTubeCertifier

using StaticArrays
import MathOptInterface as MOI

import Dionysos
import MathematicalSystems as MS
const UT = Dionysos.Utils
const DO = Dionysos.Domain
const SY = Dionysos.Symbolic



# A very simple certifier:
# - Build X_local = union of hyperrectangles around each x_k in traj.X with radius r
# - Set system.X = X_local
# - Run uniform-grid abstraction + reachability solver on that reduced domain
mutable struct CertifierUniformGridLocalTube{T} <: AbstractSymbolicCertifier
    # Inputs
    concrete_problem::Any
    traj::Any

    # Tube parameters
    radius::Any                # Number or SVector{N,T}
    margin::T
    incl_mode::DO.INCL_MODE    # e.g. DO.INNER/OUTER/CENTER depending on your conventions

    # Uniform-grid abstraction settings (forwarded)
    state_grid::Union{Nothing, DO.Grid}
    h::Any
    input_grid::Union{Nothing, DO.Grid}
    approx_mode::Any
    time_step::Union{Nothing, T}
    threaded::Bool
    print_level::Int

    # Outputs
    abstract_system::Any
    abstract_controller::Any
    controllable_set::Any
    concrete_value_function::Any
    success::Bool
    solve_time_sec::T

    function CertifierUniformGridLocalTube{T}() where {T}
        return new{T}(
            nothing, nothing,
            0.1, 0.0, DO.OUTER,
            nothing, nothing, nothing,
            Dionysos.Optim.Abstraction.UniformGridAbstraction.GROWTH,
            nothing,
            false,
            1,
            nothing, nothing, nothing, nothing,
            false, 0.0
        )
    end
end

CertifierUniformGridLocalTube() = CertifierUniformGridLocalTube{Float64}()

# -------------------------
# Minimal required API
# -------------------------
function set_problem!(c::CertifierUniformGridLocalTube, p)
    c.concrete_problem = p
    return c
end

function set_trajectory!(c::CertifierUniformGridLocalTube, traj)
    c.traj = traj
    return c
end

get_controller(c::CertifierUniformGridLocalTube) = c.abstract_controller
get_success(c::CertifierUniformGridLocalTube) = c.success
get_solve_time(c::CertifierUniformGridLocalTube) = c.solve_time_sec

# -------------------------
# Local tube -> union of rectangles
# -------------------------
# radius can be scalar or per-dim vector
@inline _rad_i(r, i) = r isa Number ? r : r[i]


function build_tube(trajX, radius; margin=0.0)
    # Build vector of HyperRectangle{SVector} and wrap in LazyUnionSetArray
    rects = Vector{UT.HyperRectangle}()
    for x in trajX
        N = length(x)
        lb = ntuple(i -> x[i] - (_rad_i(radius, i) + margin), N)
        ub = ntuple(i -> x[i] + (_rad_i(radius, i) + margin), N)
        push!(rects, UT.HyperRectangle(SVector(lb), SVector(ub)))
    end
    return UT.LazyUnionSetArray(rects)
end

# -------------------------
# "Restrict system.X" helper
# -------------------------
"""
Return a system identical to `sys` but with state constraint X replaced by X_local.

You MUST adapt this depending on how your system is represented.
Common cases:
- sys is mutable and has field X -> `sys2 = deepcopy(sys); sys2.X = X_local`
- sys is immutable -> reconstruct via its constructor
"""
function restrict_system_X(sys, X_local)
    sys2 = deepcopy(sys)
    sys2.X = X_local
    return sys2
end

# -------------------------
# Main certification
# -------------------------
function certify!(c::CertifierUniformGridLocalTube)
    t0 = time()
    @assert c.concrete_problem !== nothing
    @assert c.traj !== nothing
    @assert c.input_grid !== nothing
    @assert (c.state_grid !== nothing) || (c.h !== nothing)

    # 1) Build local union set around trajectory
    trajX = getproperty(c.traj, :X)  # expects traj.X
    X_local = build_tube(trajX, c.radius; margin=c.margin)

    # 2) Make reduced system/problem (only X changes)
    sys = c.concrete_problem.system
    sys_local = restrict_system_X(sys, X_local)

    empty_problem = Dionysos.Problem.EmptyProblem(sys_local)  # adapt if your ctor differs

    # 3) Abstraction on reduced domain (reuse your optimizer)
    abs_opt = Dionysos.Optim.Abstraction.UniformGridAbstraction.OptimizerEmptyProblem()
    MOI.set(abs_opt, MOI.RawOptimizerAttribute("empty_problem"), empty_problem)
    if c.state_grid !== nothing
        MOI.set(abs_opt, MOI.RawOptimizerAttribute("state_grid"), c.state_grid)
    else
        MOI.set(abs_opt, MOI.RawOptimizerAttribute("h"), c.h)
    end
    MOI.set(abs_opt, MOI.RawOptimizerAttribute("input_grid"), c.input_grid)
    c.time_step !== nothing && MOI.set(abs_opt, MOI.RawOptimizerAttribute("time_step"), c.time_step)
    MOI.set(abs_opt, MOI.RawOptimizerAttribute("approx_mode"), c.approx_mode)
    MOI.set(abs_opt, MOI.RawOptimizerAttribute("threaded"), c.threaded)
    MOI.set(abs_opt, MOI.RawOptimizerAttribute("print_level"), c.print_level)

    MOI.optimize!(abs_opt)
    c.abstract_system = MOI.get(abs_opt, MOI.RawOptimizerAttribute("abstract_system"))

    # 4) Solve abstract control problem (reuse your optimizer)
    ocp_opt = Dionysos.Optim.OptimizerOptimalControlProblem()
    MOI.set(ocp_opt, MOI.RawOptimizerAttribute("concrete_problem"), c.concrete_problem)
    MOI.set(ocp_opt, MOI.RawOptimizerAttribute("abstract_system"), c.abstract_system)
    MOI.set(ocp_opt, MOI.RawOptimizerAttribute("print_level"), c.print_level)

    MOI.optimize!(ocp_opt)

    c.abstract_controller = MOI.get(ocp_opt, MOI.RawOptimizerAttribute("abstract_controller"))
    c.controllable_set    = MOI.get(ocp_opt, MOI.RawOptimizerAttribute("controllable_set"))
    c.concrete_value_function = MOI.get(ocp_opt, MOI.RawOptimizerAttribute("concrete_value_function"))
    c.success             = MOI.get(ocp_opt, MOI.RawOptimizerAttribute("success"))

    c.solve_time_sec = time() - t0
    return c
end

end # module
