# Stand in for the robot against `deploy_biped_footstep.jl`.
#
# Speaks the control server's wire protocol -- a 4-byte big-endian length, then
# Float64s in network byte order, both ways -- and plays the plant the 4-D biped
# model describes: velocity control, so one step of the loop is
#
#     x⁺ = x + tstep * u        (tstep = 0.1 s)
#
# which is exact for this model, the same property that makes the abstraction a
# bisimulation. Use it to check a deployment end to end before the real robot:
# the certified footstep is 30 steps for the slew controller (27 for plain),
# ending with the swing foot on the foothold.
#
#   julia --project=control_server control_server/scripts/deploy_biped_footstep.jl   # terminal 1
#   julia --project=control_server control_server/scripts/biped_footstep_client.jl   # terminal 2
#
# Environment: BIPED_PORT (default 5000), BIPED_STEPS (default 400 max).

using Sockets

const PORT = parse(Int, get(ENV, "BIPED_PORT", "5000"))
const MAXSTEPS = parse(Int, get(ENV, "BIPED_STEPS", "400"))
const TSTEP = 0.1
const X0 = [0.2, 0.0, -0.2, 0.0]      # θ1..θ4, the example's initial posture

"Send `v` as length-prefixed network-order Float64s."
function send_vec(sock, v::Vector{Float64})
    payload = reinterpret(UInt8, hton.(reinterpret(UInt64, v)))
    write(sock, hton(UInt32(length(payload))))
    write(sock, payload)
    return flush(sock)
end

"Read one length-prefixed network-order Float64 vector."
function recv_vec(sock)
    header = Vector{UInt8}(undef, 4)
    read!(sock, header)
    nbytes = Int(ntoh(reinterpret(UInt32, header)[1]))
    payload = Vector{UInt8}(undef, nbytes)
    read!(sock, payload)
    return collect(reinterpret(Float64, ntoh.(reinterpret(UInt64, payload))))
end

sock = connect(ip"127.0.0.1", PORT)
println("connected to the control server on port $PORT")

x = copy(X0)
us = Vector{Vector{Float64}}()

for k in 1:MAXSTEPS
    send_vec(sock, x)
    local u
    try
        u = recv_vec(sock)
    catch err
        # The server closes the session when the controller is undefined at the
        # measured state -- the certificate says nothing beyond its domain.
        println("step $k: server closed the session ($(typeof(err)))")
        break
    end
    push!(us, u)
    global x = x + TSTEP * u
end
close(sock)

println("steps commanded : ", length(us))
println("final state     : ", round.(x; digits = 4))
if length(us) > 1
    slew = maximum(maximum(abs.(us[k + 1] - us[k])) for k in 1:(length(us) - 1))
    println("max |Δu|        : ", round(slew; digits = 3), " rad/s")
end
println("max |u|         : ", isempty(us) ? 0.0 : maximum(maximum(abs.(u)) for u in us))
