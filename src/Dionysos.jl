module Dionysos

include("utils/utils.jl")
include("system/system.jl")
include("problem/problem.jl")
include("mapping/mapping.jl")
include("symbolic/symbolic.jl")
include("optim/optim.jl")
include("MOI_wrapper.jl")

function export_controller_csv(args...; kwargs...)
    return error("export_controller_csv requires CSV.jl and DataFrames.jl.")
end

function import_controller_csv(args...; kwargs...)
    return error("import_controller_csv requires CSV.jl and DataFrames.jl.")
end

function spot_stepper(args...; kwargs...)
    return error(
        "spot_stepper is available only when Spot.jl is installed and loaded with `using Spot`.",
    )
end

end
