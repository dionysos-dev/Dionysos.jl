module ServerRuntime
using Sockets

import ..Controller4Server as C

export start_control_server

# Use a pre-allocated buffer for the length header (4 bytes)
const HEADER_BUF = Vector{UInt8}(undef, 4)

"""
    start_labview_server(controller; port=5000)

Starts a server that listens for LabVIEW data. 
'controller' is a structure containing a Vector{Float64} and two functions (f and g) that take two Vector{Float64} as arguments and returns a Vector{Float64}.
"""
function start_control_server(controller::C.Controller; port=5000, max_packets=100000, log_data=false, received_data_size=1, )
    server = listen(port)
    println("Server listening on port $port...")
    keep_running = true
    
    x0 = controller.x

    # Pre-allocate buffers if logging is enabled
    history = log_data ? (
        measurements = Matrix{Float64}(undef, received_data_size, max_packets), 
        states = Matrix{Float64}(undef, length(controller.x), max_packets),
        t = Vector{Float64}(undef, max_packets)
    ) : nothing

    idx = 1
    
    while keep_running
        sock = accept(server)
        println("Client connected!")
        controller.x = x0
        idx = 1 # Packet counter
        start_time = time()
        
        try
            while isopen(sock)
                # 1. Read header (4 bytes)
                read!(sock, HEADER_BUF) 
                msg_len = ntoh(reinterpret(UInt32, HEADER_BUF)[1])

                if msg_len == UInt32(4294967295)
                    # 2a. Reset controller if "special" input is received
                    controller.x = x0
                else
                    # 2b. Read payload into a Float64 vector otherwise
                    payload = Vector{Float64}(undef, Int(msg_len / 8))
                    read!(sock, payload)
                    
                    # 3. Fix endianness
                    measurements = ntoh.(payload)

                    # 4. Update the controller state and send back the controller output
                    x_plus = controller.f(controller.x, measurements)
                    control = controller.g(controller.x, measurements)
                    controller.x = x_plus

                    # 5. Data Logging (if enabled)
                    if log_data && idx <= max_packets
                        history.t[idx] = time() - start_time
                        history.measurements[:, idx] = measurements
                        history.states[:, idx] = controller.x
                        idx += 1
                    end

                    # 6. Send response (Length + Data)
                    resp_data = hton.(control)
                    write(sock, hton(UInt32(length(resp_data) * 8)))
                    write(sock, resp_data)
                end
            end
        catch e
            if !(e isa EOFError)
                @error "Session error" exception=e
            end
        finally
            close(sock)
            println("Client disconnected.")
            keep_running = false
        end

    end

    println("Closing server...")
    close(server)

    # 6. Return or process the matrices
    if log_data
        # Trim the matrices to the actual number of packets received
        t = history.t[1:idx-1]
        measurements = history.measurements[:, 1:idx-1]
        states = history.states[:, 1:idx-1]
        println("Logged $(idx-1) packets.")
    end

    return log_data ? (t, measurements, states) : nothing
end

end #module