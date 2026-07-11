module TestMain

import Dionysos
include(joinpath(dirname(dirname(pathof(Dionysos))), "test", "testsetup.jl"))

const PID = ST.PIDControllers

@testset "PID controller" begin
    # Proportional-only controller, reference (1, 0), plain error e = r - x.
    pid = PID.PIDControllerVector(;
        Kp = SVector(2.0, 0.0),
        Ki = SVector(0.0, 0.0),
        Kd = SVector(0.0, 0.0),
        ref = PID.ConstantSignal(SVector(1.0, 0.0)),
        error = (x, r, t) -> r .- x,
        dt = 0.1,
        e0 = SVector(0.0, 0.0),
    )

    @test ST.controller_kind(pid) isa ST.DynamicKind

    mem0 = ST.initial_state(pid)
    @test mem0 isa PID.PIDMemory
    @test mem0.initialized == false

    # Kp · (r - x) = dot([2,0], [1,0]) = 2 at x = 0.
    u = ST.output_control(pid, mem0, SVector(0.0, 0.0))
    @test u ≈ 2.0

    # State advances and marks itself initialized.
    mem1 = ST.update_state(pid, mem0, SVector(0.0, 0.0))
    @test mem1 isa PID.PIDMemory
    @test mem1.initialized == true

    # Saturation clamps the output.
    pid_sat = PID.PIDControllerVector(;
        Kp = SVector(100.0, 0.0),
        Ki = SVector(0.0, 0.0),
        Kd = SVector(0.0, 0.0),
        ref = PID.ConstantSignal(SVector(1.0, 0.0)),
        error = (x, r, t) -> r .- x,
        dt = 0.1,
        umin = -1.0,
        umax = 1.0,
        e0 = SVector(0.0, 0.0),
    )
    u_sat = ST.output_control(pid_sat, ST.initial_state(pid_sat), SVector(0.0, 0.0))
    @test all(u_sat .<= 1.0) && all(u_sat .>= -1.0)
end

end # module TestMain
