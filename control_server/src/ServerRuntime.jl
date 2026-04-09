module ServerRuntime
using Sockets

import ..Controller4Server as C

export start_control_server

# Use a pre-allocated buffer for the length header (4 bytes)
const HEADER_BUF = Vector{UInt8}(undef, 4)

"""
    start_labview_server(process_fn; port=5000)

Starts a server that listens for LabVIEW data. 
'process_fn' should be a function that takes a Vector{Float64} and returns a Vector{Float64}.
"""
function start_control_server(controller::C.Controller; port=5000)
    server = listen(port)
    println("Server listening on port $port...")
    x0 = controller.x
    
    while true
        sock = accept(server)
        println("Client connected!")
        controller.x = x0
        
        @async begin # Use @async to handle multiple clients simultaneously
            try
                while isopen(sock)
                    # 1. Read header (4 bytes)
                    read!(sock, HEADER_BUF) 
                    msg_len = ntoh(reinterpret(UInt32, HEADER_BUF)[1])

                    # 2. Read payload into a Float64 vector
                    payload = Vector{Float64}(undef, Int(msg_len / 8))
                    read!(sock, payload)
                    
                    # 3. Fix endianness
                    values = ntoh.(payload)

                    # 4. Update the controller state and send back the controller output
                    x_plus = controller.f(controller.x, values)
                    controller.x = x_plus
                    result = controller.g(controller.x)

                    # 5. Send response (Length + Data)
                    resp_data = hton.(result)
                    write(sock, hton(UInt32(length(resp_data) * 8)))
                    write(sock, resp_data)
                end
            catch e
                if !(e isa EOFError)
                    @error "Session error" exception=e
                end
            finally
                close(sock)
                println("Client disconnected.")
            end
        end
    end
end

end