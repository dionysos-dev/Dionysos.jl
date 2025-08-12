module TestMain

using Test, StaticArrays, Plots
using Dionysos

const DI = Dionysos
const UT = DI.Utils
const DO = DI.Domain
const ST = DI.System
const SY = DI.Symbolic

function F_sys(x::SVector{N, T}, u) where {N, T}
    return u
end

function get_traj(symmodel, sys, x0, nstep)
    u = SVector(1.0, 1.0)
    x = copy(x0)
    s = SY.get_abstract_state(symmodel, x)
    xpos = DO.get_pos_by_coord(symmodel.Xdom, x)
    traj = [(x, s, xpos)]
    for k in 1:nstep
        x = sys.sys_map(x, u, sys.tstep)
        s = SY.get_abstract_state(symmodel, x)
        xpos = DO.get_pos_by_coord(symmodel.Xdom, x)
        push!(traj, (x, s, xpos))
    end
    return traj
end

function test()
    #system
    tstep = 0.8
    measnoise = SVector(0.0, 0.0)
    sys = ST.NewSimpleSystem(tstep, F_sys, measnoise, 4)

    #input space
    Udom = DO.CustomList([SVector(1.0, 1.0)])

    #state space
    X = UT.HyperRectangle(SVector(0.0, 0.0), SVector(30.0, 30.0))
    obstacle = UT.HyperRectangle(SVector(10.0, 10.0), SVector(15.0, 15.0))

    # Grid settings
    hx = SVector(3.0, 0.3)

    # Create a LazySetMinus: free space = X minus obstacle
    free_space = UT.LazySetMinus(X, UT.LazyUnionSetArray([obstacle]))

    x0 = SVector(0.0, 0.0)
    nstep = 10
    @testset "LazySymbolic" begin
        Xdom = DO.DomainList(hx)
        symmodel = SY.LazySymbolicModelList(Xdom, Udom, free_space)
        traj = get_traj(symmodel, sys, x0, nstep)
        (x, s, xpos) = traj[1]
        @test (s, xpos) == (1, (0, 0))

        (x, s, xpos) = traj[2]
        @test (s, xpos) == (2, (0, 3))

        (x, s, xpos) = traj[3]
        @test (s, xpos) == (3, (1, 5))

        (x, s, xpos) = traj[4]
        @test (s, xpos) == (4, (1, 8))

        (x, s, xpos) = traj[5]
        @test (s, xpos) == (5, (1, 11))

        (x, s, xpos) = traj[6]
        @test (s, xpos) == (6, (1, 13))

        (x, s, xpos) = traj[7]
        @test (s, xpos) == (7, (2, 16))

        (x, s, xpos) = traj[8]
        @test (s, xpos) == (8, (2, 19))

        (x, s, xpos) = traj[9]
        @test (s, xpos) == (9, (2, 21))

        (x, s, xpos) = traj[10]
        @test (s, xpos) == (10, (2, 24))

        (x, s, xpos) = traj[11]
        @test (s, xpos) == (11, (3, 27))

        @test SY.get_n_state(symmodel) == 11

        @test all(SY.enum_states(symmodel) .== 1:11)

        translist = [(1, 2, 1), (1, 8, 1)]
        SY.add_transitions!(symmodel.autom, translist)

        fig = plot(; aspect_ratio = :equal)
        lyap_fun = Dict(state => 2.0 * state for state in SY.enum_states(symmodel))
        plot!(fig, symmodel; arrowsB = true, cost = true, lyap_fun = lyap_fun)
        @test isa(fig, Plots.Plot{Plots.GRBackend})

        fig = plot(; aspect_ratio = :equal)
        plot!(fig, symmodel; arrowsB = true, cost = false)
        @test isa(fig, Plots.Plot{Plots.GRBackend})
    end
end

test()

end
