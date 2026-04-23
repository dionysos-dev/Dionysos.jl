module ControlServerDeployment

using Dionysos
using ..Controller4Server

function to_server_controller(ctrl::Dionysos.AbstractDeployableController)
    x0 = Vector{Float64}(Dionysos.initial_state(ctrl))
    f = (x, y) -> Vector{Float64}(Dionysos.update_state(ctrl, x, y))
    g = (x, y) -> Vector{Float64}(Dionysos.output_control(ctrl, x, y))
    return Controller4Server.Controller(x0, f, g)
end

function load_server_controller(filename::AbstractString)
    ctrl = Dionysos.import_controller_jld2(filename)
    return to_server_controller(ctrl)
end

end
