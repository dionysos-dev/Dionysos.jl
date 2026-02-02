module LVServerTransport

# ---------------------------
# Shared controller entry point
# ---------------------------
const CS = Main.ControllerRuntime.Server.ControllerServer

# =====================================================
# Approach A: Use LVServer.server_0mq4lv (LVServer protocol)
# =====================================================
module ApproachA_LVServer

    using ..LVServerTransport: CS
    using LVServer

    export start_lvserver

    # Helper: normalize x into Vector{Float64} if it is a vector-like
    _to_vecfloat(x) = x
    _to_vecfloat(x::AbstractVector) = Float64.(x)

    """
    Start an LVServer-protocol ZMQ server (whatever LVServer.server_0mq4lv expects on wire).

    LabVIEW must use a matching LVServer client/protocol.
    """
    function start_lvserver(; initOK::Bool = true)
        load_controller(; path) = CS.load_controller!(; path = path)
        unload_controller(;) = CS.unload_controller!()
        reset_controller(; kwargs...) = CS.reset_controller!(; kwargs...)

        function step_controller(; x, t = 0.0, kwargs...)
            u = CS.step_controller(; x = x, t = t, kwargs...)

            # If u is an array, return it under bigarrs.
            if u isa AbstractArray
                return (; bigarrs = (; u,))
            else
                # scalar / small output: OK to return directly
                return (; u = u)
            end
        end

        return LVServer.server_0mq4lv(
            (; load_controller, unload_controller, reset_controller, step_controller);
            initOK = initOK,
        )
    end

end # module ApproachA_LVServer

# =====================================================
# Approach B: Simple JSON over ZMQ (explicit protocol)
# =====================================================
module ApproachB_JSON

    using ..LVServerTransport: CS
    using ZMQ
    using JSON3

    export start_json_server

    """
    Simple JSON protocol over ZMQ REQ/REP.

    Client sends a JSON string:

      {"cmd":"load_controller","args":{"path":"..."}}
      {"cmd":"reset_controller","args":{"seed":1}}
      {"cmd":"step_controller","args":{"x":[...], "t":0.01}}
      {"cmd":"unload_controller","args":{}}
      {"cmd":"shutdown","args":{}}

    Server replies:

      {"ok":true, "result": ...}
    or
      {"ok":false, "error":"...", "stack":"..."}
    """
    function start_json_server(;
        endpoint::AbstractString = "tcp://*:5556",
        verbose::Bool = true,
    )
        ctx = Context()
        sock = Socket(ctx, REP)
        bind(sock, endpoint)

        verbose && println("[ApproachB_JSON] Listening on $endpoint")

        # --- handlers ---
        function load_controller(args)
            path = get(args, :path, nothing)
            path === nothing && error("Missing required argument: path")
            CS.load_controller!(; path = String(path))
            return true
        end

        function unload_controller(args)
            CS.unload_controller!()
            return true
        end

        function reset_controller(args)
            kwargs = Pair{Symbol, Any}[]
            for (k, v) in pairs(args)
                push!(kwargs, Symbol(k) => v)
            end
            CS.reset_controller!(; kwargs...)
            return true
        end

        function step_controller(args)
            x = get(args, :x, nothing)
            x === nothing && error("Missing required argument: x")

            t = Float64(get(args, :t, 0.0))

            xval = x isa AbstractVector ? Float64.(collect(x)) : Float64(x)

            kwargs = Pair{Symbol, Any}[]
            for (k, v) in pairs(args)
                if k == :x || k == :t
                    continue
                end
                push!(kwargs, Symbol(k) => v)
            end

            u = CS.step_controller(; x = xval, t = t, kwargs...)
            return u
        end

        handlers = Dict{String, Function}(
            "load_controller" => load_controller,
            "unload_controller" => unload_controller,
            "reset_controller" => reset_controller,
            "step_controller" => step_controller,
        )

        # --- main loop ---
        running = true
        while running
            msg = String(recv(sock))
            verbose && println("[ApproachB_JSON] Received: $msg")

            reply = Dict{String, Any}()

            try
                req = JSON3.read(msg)
                cmd = String(get(req, :cmd, ""))
                args = get(req, :args, JSON3.Object())

                cmd == "" && error("Missing field: cmd")

                if cmd == "shutdown"
                    running = false
                    reply["ok"] = true
                    reply["result"] = "shutting down"
                else
                    f = get(handlers, cmd, nothing)
                    f === nothing && error("Unknown cmd: $cmd")

                    result = f(args)
                    reply["ok"] = true
                    reply["result"] = result
                end
            catch err
                reply["ok"] = false
                reply["error"] = sprint(showerror, err)
                reply["stack"] = sprint(showerror, err, catch_backtrace())
            end

            reply_json = JSON3.write(reply)
            send(sock, reply_json)
            verbose && println("[ApproachB_JSON] Replied: $reply_json")
        end

        close(sock)
        close(ctx)

        verbose && println("[ApproachB_JSON] Server stopped.")
        return nothing
    end

end # module ApproachB_JSON

# =====================================================
# Convenience: run either or both
# =====================================================
export run_lvserver, run_json_server, run_both

"""
Run Approach A (LVServer protocol). Blocks (server runs in foreground).
"""
run_lvserver(; initOK::Bool = true) = ApproachA_LVServer.start_lvserver(; initOK = initOK)

"""
Run Approach B (JSON protocol). Blocks (server runs in foreground).
"""
run_json_server(; endpoint::AbstractString = "tcp://*:5556", verbose::Bool = true) =
    ApproachB_JSON.start_json_server(; endpoint = endpoint, verbose = verbose)

"""
Run BOTH servers concurrently:
- LVServer protocol (Approach A) in one task
- JSON protocol (Approach B) in another task (default port 5556)

Keeps Julia alive until you Ctrl+C.
"""
function run_both(;
    initOK::Bool = true,
    json_endpoint::AbstractString = "tcp://*:5556",
    verbose::Bool = true,
)
    tA = @async begin
        try
            ApproachA_LVServer.start_lvserver(; initOK = initOK)
        catch err
            @error "Approach A crashed" err
            rethrow()
        end
    end

    tB = @async begin
        try
            ApproachB_JSON.start_json_server(; endpoint = json_endpoint, verbose = verbose)
        catch err
            @error "Approach B crashed" err
            rethrow()
        end
    end

    println("[LVServerTransport] Running BOTH servers.")
    println("  - Approach A: LVServer protocol (LVServer.server_0mq4lv)")
    println("  - Approach B: JSON protocol on $json_endpoint")

    wait(tA)  # If A exits, we exit. (Adjust if you want different behavior.)
    wait(tB)
    return nothing
end

end # module LVServerTransport
