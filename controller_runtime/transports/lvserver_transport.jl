module LVServerTransport

using LVServer
const CS = Main.ControllerRuntime.Server.ControllerServer

export start_lvserver

function start_lvserver(; initOK::Bool = true)
    load_controller(; path) = CS.load_controller!(; path = path)
    unload_controller(;) = CS.unload_controller!()
    reset_controller(; kwargs...) = CS.reset_controller!(; kwargs...)
    step_controller(; x, t = 0.0, kwargs...) = CS.step_controller(; x = x, t = t, kwargs...)

    return LVServer.server_0mq4lv(
        (; load_controller, unload_controller, reset_controller, step_controller);
        initOK = initOK,
    )
end

end
