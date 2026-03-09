import Dionysos

include(joinpath(@__DIR__, "..", "helpers", "helpers.jl"))
include(joinpath(@__DIR__, "..", "helpers", "plot.jl"))

"""
Placeholder benchmark for cercle trajectory scenario.
Keep this script as the future entry point for circle-tracking experiments.
"""
function main()
    outputs = joinpath(@__DIR__, "outputs")
    mkpath(outputs)
    println("TODO: benchmark cercle a implementer.")
    println("Outputs directory prepared at: ", outputs)
    return (; outputs)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
