module ServerRuntime
using Sockets

import ..Controller4Server as C

export start_control_server

const HEADER_BUF = Vector{UInt8}(undef, 4)

"""
    start_control_server(controller; port=5000, max_packets=100000,
                         log_data=false, received_data_size=1,
                         state_to_vector=nothing)

Starts a server that listens for incoming measurement vectors.

- `controller.x` may be any controller-memory object.
- `controller.f(x, y)` returns the next controller state.
- `controller.g(x, y)` returns the control output to send back as numeric data.
- If `log_data=true`, the function logs time and measurements.
- If `log_data=true` and `state_to_vector !== nothing`, controller states are also logged
  through `state_to_vector(controller.x)`.
- Control outputs are logged whenever `log_data=true`.

Returns:
- `nothing` if `log_data=false`
- `(t, measurements, controls, states)` if `log_data=true`
  where `states === nothing` when `state_to_vector === nothing`.
"""
function start_control_server(
    controller::C.Controller;
    port = 5000,
    max_packets = 100000,
    log_data = false,
    received_data_size = 1,
    state_to_vector = nothing,
)
    server = listen(port)
    println("Server listening on port $port...")
    keep_running = true

    x0 = deepcopy(controller.x)

    measurement_history =
        log_data ? Matrix{Float64}(undef, received_data_size, max_packets) : nothing
    time_history = log_data ? Vector{Float64}(undef, max_packets) : nothing

    state_history = nothing
    if log_data && state_to_vector !== nothing
        state_dim = length(state_to_vector(controller.x))
        state_history = Matrix{Float64}(undef, state_dim, max_packets)
    end

    control_history = nothing
    idx = 1

    try
        while keep_running
            sock = accept(server)
            println("Client connected!")
            controller.x = deepcopy(x0)
            idx = 1
            start_time = time()

            try
                while isopen(sock)
                    # 1. Read header (4 bytes)
                    read!(sock, HEADER_BUF)
                    msg_len = ntoh(reinterpret(UInt32, HEADER_BUF)[1])

                    # Special reset message
                    if msg_len == UInt32(0xFFFFFFFF)
                        controller.x = deepcopy(x0)
                        continue
                    end

                    # 2. Read payload as raw bytes
                    nbytes = Int(msg_len)
                    nbytes % 8 == 0 ||
                        error("Payload size $nbytes is not a multiple of 8 bytes.")

                    payload_bytes = Vector{UInt8}(undef, nbytes)
                    read!(sock, payload_bytes)

                    # 3. Convert network-order bytes to Float64 vector
                    payload_u64 = reinterpret(UInt64, payload_bytes)
                    measurements = reinterpret(Float64, ntoh.(payload_u64))

                    # 4. Update the controller state and compute control output
                    x_plus = controller.f(controller.x, measurements)
                    controller.x = x_plus
                    control = controller.g(controller.x, measurements)

                    # Ensure control is a Float64 vector for transmission/logging
                    control_vec = Float64[control...]

                    # 5. Data logging
                    if log_data && idx <= max_packets
                        time_history[idx] = time() - start_time
                        measurement_history[:, idx] = measurements

                        if state_history !== nothing
                            state_history[:, idx] = Float64[controller.x...]
                        end

                        if control_history === nothing
                            control_history =
                                Matrix{Float64}(undef, length(control_vec), max_packets)
                        end
                        control_history[:, idx] = control_vec

                        idx += 1
                    end

                    # 6. Send response (Length + Data) in network byte order
                    resp_u64 = hton.(reinterpret(UInt64, control_vec))
                    resp_bytes = reinterpret(UInt8, resp_u64)
                    write(sock, hton(UInt32(length(resp_bytes))))
                    write(sock, resp_bytes)
                end
            catch e
                if !(e isa EOFError)
                    @error "Session error" exception = e
                end
            finally
                close(sock)
                println("Client disconnected.")
                # I think we need here:
                # keep_running = false
                # Otherwise, the server keeps running
            end
        end
    catch e
        if e isa InterruptException
            println("\nStopping server...")
        else
            rethrow(e)
        end
    finally
        println("Closing server...")
        close(server)
    end

    # 7. Return logged data
    if log_data
        n = idx - 1
        t = time_history[1:n]
        measurements = measurement_history[:, 1:n]
        controls = control_history === nothing ? nothing : control_history[:, 1:n]
        states = state_history === nothing ? nothing : state_history[:, 1:n]
        println("Logged $n packets.")
        return (t, measurements, controls, states)
    end

    return nothing
end

end # module
