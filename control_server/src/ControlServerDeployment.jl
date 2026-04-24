module ControlServerDeployment

using Dionysos
const ST = Dionysos.System

using JLD2
using ..Controller4Server

export to_server_controller, load_server_controller

function to_server_controller(ctrl::ST.AbstractController)
    x0 = ST.initial_state(ctrl)
    f = (x, y) -> ST.update_state(ctrl, x, y)
    g = (x, y) -> ST.output_control(ctrl, x, y)
    return Controller4Server.Controller(x0, f, g)
end

function load_server_controller(filename::AbstractString)
    ctrl = JLD2.load(filename, "controller")
    return to_server_controller(ctrl)
end

end # module
