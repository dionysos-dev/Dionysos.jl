using StaticArrays
import MathOptInterface as MOI
import Dionysos
const DI = Dionysos
const UT = DI.Utils
const ST = DI.System
const DO = DI.Domain
const OP = DI.Optim
const AB = OP.Abstraction

# Your API
import ..SymbolicCertifier: AbstractSymbolicCertifier, set_trajectory!, certify!,
                              get_success, get_solve_time

"""
Minimal certifier:
- user provides a configured AB.UniformGridAbstraction.Optimizer
- certifier provides traj + radius (+ margin + incl_mode)
- no deepcopy: we temporarily modify the system.X in the *problem* given to the optimizer,
  run optimize!, then restore it.
"""
mutable struct UniformGridLocalTubeCertifier{T} <: AbstractSymbolicCertifier
    # Required
    traj::Any

    # Tube params
    radius::Any
    margin::T
    incl_mode::DO.INCL_MODE

    # Injected solver (already configured by user)
    optimizer::Union{Nothing, AB.UniformGridAbstraction.Optimizer}

    # Outputs
    solve_time_sec::T

    function UniformGridLocalTubeCertifier{T}() where {T}
        return new{T}(nothing, 0.1, 0.0, DO.INNER, nothing, 0.0)
    end
end

UniformGridLocalTubeCertifier() = UniformGridLocalTubeCertifier{Float64}()

set_optimizer!(c::UniformGridLocalTubeCertifier, optimizer::AB.UniformGridAbstraction.Optimizer) = (c.optimizer = optimizer; c)
set_trajectory!(c::UniformGridLocalTubeCertifier, x_traj) = (c.traj = x_traj; c)
get_success(c::UniformGridLocalTubeCertifier) = c.optimizer.control_solver.success
get_solve_time(c::UniformGridLocalTubeCertifier) = c.solve_time_sec
get_controller(c::UniformGridLocalTubeCertifier) = c.optimizer.concrete_controller


function certify!(c::UniformGridLocalTubeCertifier)
    t0 = time()
    @assert c.optimizer !== nothing "Set cert.optimizer = your configured AB.UniformGridAbstraction.Optimizer"
    @assert c.traj !== nothing "Call set_trajectory! first."
    
    x_traj = c.traj
    # Build tube set
    X_local = build_tube(x_traj, c.radius; margin=c.margin)

    # Build abstracion & solve problem
    optimizer = c.optimizer
    concrete_problem = MOI.get(optimizer, MOI.RawOptimizerAttribute("concrete_problem"))
    @assert concrete_problem !== nothing "Your uniformGridOptimizer must already have concrete_problem set."
    MOI.set(optimizer, MOI.RawOptimizerAttribute("abstraction_region"), X_local)
    MOI.set(optimizer, MOI.RawOptimizerAttribute("incl_mode"), c.incl_mode)
    MOI.optimize!(optimizer)

    c.solve_time_sec = time() - t0
    return c
end


# radius can be scalar or per-dim vector
@inline _rad_i(r, i) = r isa Number ? r : r[i]

function build_tube(x_traj::ST.Trajectory, radius; margin=0.0)
    rects = Vector{UT.HyperRectangle}()
    for x in ST.enum_elems(x_traj)
        N = length(x)
        lb = ntuple(i -> x[i] - (_rad_i(radius, i) + margin), N)
        ub = ntuple(i -> x[i] + (_rad_i(radius, i) + margin), N)
        push!(rects, UT.HyperRectangle(SVector(lb), SVector(ub)))
    end
    return UT.LazyUnionSetArray(rects)
end
