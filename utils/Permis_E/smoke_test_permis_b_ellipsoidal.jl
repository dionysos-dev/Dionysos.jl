using Dionysos
using JuMP
import Clarabel
import IntervalArithmetic as IA
import MathOptInterface as MOI
import StaticArrays: SVector

import MosekTools

const DI = Dionysos
const UT = DI.Utils
const OP = DI.Optim
const SC = OP.SymbolicCertifier
const ST = DI.System

include(joinpath(@__DIR__, "get_nominal_traje.jl"))

nominal = get_nominal_traje()
@assert nominal.candidate !== nothing

function print_state_list(state_list)
    println("State list (" * string(length(state_list)) * " states):")
    for (i, x) in enumerate(state_list)
        println("[" * string(i) * "] " * string(x))
    end
end

function print_input_list(input_list)
    println("Input list (" * string(length(input_list)) * " controls):")
    for (i, u) in enumerate(input_list)
        println("[" * string(i) * "] " * string(u))
    end
end

# je devrais mettre ça dans /Users/florentincanart/Dionysos.jl-updated/src/domain/grid_domain/periodic_domain.jl
function unwrap_periodic_state_list(state_list, periodic_dims, periodic_periods)
    isempty(state_list) && return state_list
    length(periodic_dims) == length(periodic_periods) ||
        error("periodic_dims et periodic_periods doivent avoir la meme longueur.")

    nx = length(state_list[1])
    xs = [collect(Float64, x) for x in state_list]

    for i in eachindex(periodic_dims)
        d = Int(periodic_dims[i])
        p = Float64(periodic_periods[i])
        1 <= d <= nx || error("Dimension periodique invalide: $d.")
        p > 0 || error("Periode invalide (<= 0): $p.")

        for k in 2:length(xs)
            Δ_raw = xs[k][d] - xs[k - 1][d]
            Δ_min = Δ_raw - round(Δ_raw / p) * p
            xs[k][d] = xs[k - 1][d] + Δ_min
        end
    end

    return [SVector{nx, Float64}(x) for x in xs]
end

Wdom = UT.HyperRectangle(SVector(0.0), SVector(0.0))
# backend = optimizer_with_attributes(Clarabel.Optimizer, MOI.Silent() => false)

backend = optimizer_with_attributes(
    MosekTools.Optimizer,
    MOI.Silent() => true,
    MOI.RawOptimizerAttribute("MSK_IPAR_NUM_THREADS") => 1,
    MOI.RawOptimizerAttribute("MSK_IPAR_LOG") => 10,
)

opts = (
    maxδx = 100.0,
    maxδu = 200.0,
    λ = 0.001,
    rayon_terminal = 0.45,
    ΔX = IA.IntervalBox(
        IA.interval(-1.0, 1.0),
        IA.interval(-1.0, 1.0),
        IA.interval(-1.0, 1.0),
        IA.interval(-1.0, 1.0),
    ),
    ΔU = IA.IntervalBox(IA.interval(-1.0, 1.0), IA.interval(-1.2, 1.2)),
    ΔW = IA.IntervalBox(IA.interval(0.5, 0.5), 1),
    symbolic_rk4_substeps = 1,
)

cert_cfg = SC.EllipsoidalBackwardConfig(
    nominal.problem.system.X,
    nominal.problem.system.U,
    Wdom,
    backend,
    opts,
)

params = AV.Params(; L1 = 1.0, L2 = 1.0, Lc = 0.5)

symbolic_builder = function (prob, candidate, cfg)
    o = cfg.options
    return AV.symbolic_system(
        prob.system.X;
        _U_ = prob.system.U,
        params = params,
        Ts = candidate.Ts,
        ΔX = o.ΔX,
        ΔU = o.ΔU,
        ΔW = o.ΔW,
        rk4_num_substeps = o.symbolic_rk4_substeps,
    )
end

cand_raw = nominal.candidate
xs_raw = collect(ST.enum_elems(cand_raw.x_traj))
us_raw = collect(ST.enum_elems(cand_raw.u_traj))
print_state_list(xs_raw)
print_input_list(us_raw)
xs_unwrapped = unwrap_periodic_state_list(xs_raw, SVector(3, 4), SVector(2pi, 2pi))
cand = OP.CandidateTrajectory(
    ST.Trajectory(xs_unwrapped),
    cand_raw.u_traj;
    Ts = cand_raw.Ts,
    source = cand_raw.source,
    metadata = cand_raw.metadata,
)

cert = SC.EllipsoidalBackwardCertifier{
    typeof(nominal.problem),
    typeof(cand),
    typeof(cert_cfg),
    Any,
    typeof(symbolic_builder),
}(
    nothing,
    nothing,
    cert_cfg,
    nothing,
    false,
    0.0,
    symbolic_builder,
)

SC.set_problem!(cert, nominal.problem)
SC.set_trajectory!(cert, cand)
SC.certify!(cert)

res = SC.get_result(cert)
@assert res !== nothing
@assert length(res.steps) >= 1

println("generator_success = ", nominal.success)
println("certifier_success = ", SC.get_success(cert))
println("n_steps = ", length(res.steps))
