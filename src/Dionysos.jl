module Dionysos

include("utils/utils.jl")
include("system/system.jl")
include("problem/problem.jl")
include("mapping/mapping.jl")
include("symbolic/symbolic.jl")
include("optim/optim.jl")

# ----- Wrapper functions for optional dependencies ---------
import JuMP
function Optimizer end
export ∂, Δ, final, start

function _diff end
const ∂ = JuMP.NonlinearOperator(_diff, :∂)

function _delta end
const Δ = JuMP.NonlinearOperator(_delta, :Δ)

function _final end
const final = JuMP.NonlinearOperator(_final, :final)

function _start end
const start = JuMP.NonlinearOperator(_start, :start)

# ----- CSV functions for optional dependencies ---------

function export_controller_csv(args...; kwargs...)
    return error("export_controller_csv requires CSV.jl and DataFrames.jl.")
end

function import_controller_csv(args...; kwargs...)
    return error("import_controller_csv requires CSV.jl and DataFrames.jl.")
end

# ----- Spot functions for optional dependencies ---------

function spot_stepper(args...; kwargs...)
    return error(
        "spot_stepper is available only when Spot.jl is installed and loaded with `using Spot`.",
    )
end

end
