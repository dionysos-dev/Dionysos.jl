module TestMain

using Test
using StaticArrays, Plots
using Dionysos
const DI = Dionysos
const UT = DI.Utils
const DO = DI.Domain
const ST = DI.System
const SY = DI.Symbolic

sleep(0.1) # used for good printing
println("Started test")

@testset "SymbolicModel" begin
    x0 = SVector(0.0, 0.0)
    h = SVector(1.0, 2.0)
    Xgrid = DO.GridFree(x0, h)
    Xfull = DO.DomainList(Xgrid)
    DO.add_pos!(Xfull, (1, 1))
    DO.add_pos!(Xfull, (2, 2))

    u0 = SVector(0.0)
    h = SVector(0.5)
    Ugrid = DO.GridFree(u0, h)
    Ufull = DO.DomainList(Ugrid)
    DO.add_pos!(Ufull, (0,))

    symmodel = SY.NewSymbolicModelListList(Xfull, Ufull)

    @test SY.is_xpos(symmodel, (1, 1)) == true
    @test SY.is_xpos(symmodel, (2, 2)) == true
    @test SY.is_xpos(symmodel, (3, 4)) == false
    @test SY.is_xpos(symmodel, (0, 0)) == false

    stateslist = Int[]
    push!(stateslist, SY.get_state_by_xpos(symmodel, (1, 1)))
    push!(stateslist, SY.get_state_by_xpos(symmodel, (2, 2)))
    sort!(stateslist)
    @test all(stateslist .== [1, 2])

    xposlist = Tuple{Int, Int}[]
    push!(xposlist, SY.get_xpos_by_state(symmodel, 1))
    push!(xposlist, SY.get_xpos_by_state(symmodel, 2))
    sort!(xposlist)
    @test all(xposlist .== [(1, 1), (2, 2)])

    @test SY.get_concrete_input(symmodel, 1) == [0.0]

    @test SY.get_abstract_state(symmodel, SVector(1.0, 2.5)) == 1
    @test SY.get_abstract_state(symmodel, SVector(1.5, 3.0)) == 2

    subDomain = SY.get_domain_from_states(symmodel, [1])
    positions = [pos for pos in DO.enum_pos(subDomain)]
    @test all(positions .== [(1, 1)])

    subDomain = SY.get_domain_from_states(symmodel, [1, 2])
    positions = [pos for pos in DO.enum_pos(subDomain)]
    @test all(positions .== [(1, 1), (2, 2)])

    translist = [(1, 2, 1), (2, 1, 1)]
    SY.add_transitions!(symmodel.autom, translist)

    fig = plot(; aspect_ratio = :equal)
    lyap_fun = Dict(state => 2.0 * state for state in SY.enum_states(symmodel))
    plot!(fig, symmodel; arrowsB = true, cost = true, lyap_fun = lyap_fun)
    @test isa(fig, Plots.Plot{Plots.GRBackend})

    fig = plot(; aspect_ratio = :equal)
    plot!(fig, symmodel; arrowsB = true, cost = false)
    @test isa(fig, Plots.Plot{Plots.GRBackend})
end

sleep(0.1) # used for good printing
println("End test")

end  # module TestMain
