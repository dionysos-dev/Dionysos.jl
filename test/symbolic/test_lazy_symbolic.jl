module TestMain

using Test, StaticArrays, Plots
using Dionysos

const DI = Dionysos
const UT = DI.Utils
const DO = DI.Domain
const ST = DI.System
const SY = DI.Symbolic

function f1(x)
    return x
end
function fi1(x)
    return x
end
function f2(x)
    return SVector(x[2] + sin(x[1]), x[1])
end
function fi2(x)
    return SVector(x[2], x[1] - sin(x[2]))
end

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
    hx = [3.0, 0.3]
    d = DO.RectangularObstacles(X, [obstacle])

    x0 = SVector(0.0, 0.0)
    nstep = 10

    @testset "LazySymbolic" begin
        Xdom = DO.GeneralDomainList(hx; elems = d, f = f1, fi = fi1, fit = true)
        symmodel = SY.LazySymbolicModelList(Xdom, Udom)
        traj = get_traj(symmodel, sys, x0, nstep)
        (x, s, xpos) = traj[1]
        @test (s, xpos) == (1, (0, 0))

        (x, s, xpos) = traj[2]
        @test (s, xpos) == (2, (0, 2))

        (x, s, xpos) = traj[3]
        @test (s, xpos) == (3, (0, 5))

        (x, s, xpos) = traj[4]
        @test (s, xpos) == (4, (0, 8))

        (x, s, xpos) = traj[5]
        @test (s, xpos) == (5, (1, 10))

        (x, s, xpos) = traj[6]
        @test (s, xpos) == (6, (1, 13))

        (x, s, xpos) = traj[7]
        @test (s, xpos) == (7, (1, 16))

        (x, s, xpos) = traj[8]
        @test (s, xpos) == (8, (1, 18))

        (x, s, xpos) = traj[9]
        @test (s, xpos) == (9, (2, 21))

        (x, s, xpos) = traj[10]
        @test (s, xpos) == (10, (2, 24))

        (x, s, xpos) = traj[11]
        @test (s, xpos) == (11, (2, 26))

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

    @testset "Deformed grid and LazySymbolic" begin
        Xdom = DO.GeneralDomainList(hx; elems = d, f = f2, fi = fi2, fit = true)
        symmodel = SY.LazySymbolicModelList(Xdom, Udom)
        traj = get_traj(symmodel, sys, x0, nstep)
        (x, s, xpos) = traj[1]
        @test (s, xpos) == (1, (0, 0))

        (x, s, xpos) = traj[2]
        @test (s, xpos) == (1, (0, 0))

        (x, s, xpos) = traj[3]
        @test (s, xpos) == (2, (0, 2))

        (x, s, xpos) = traj[4]
        @test (s, xpos) == (3, (0, 5))

        (x, s, xpos) = traj[5]
        @test (s, xpos) == (4, (1, 10))

        (x, s, xpos) = traj[6]
        @test (s, xpos) == (5, (1, 15))

        (x, s, xpos) = traj[7]
        @test (s, xpos) == (6, (1, 19))

        (x, s, xpos) = traj[8]
        @test (s, xpos) == (7, (1, 20))

        (x, s, xpos) = traj[9]
        @test (s, xpos) == (8, (2, 20))

        (x, s, xpos) = traj[10]
        @test (s, xpos) == (9, (2, 21))

        (x, s, xpos) = traj[11]
        @test (s, xpos) == (10, (2, 23))

        @test SY.get_n_state(symmodel) == 10
        @test all(SY.enum_states(symmodel) .== 1:10)
    end
end

test()

end
