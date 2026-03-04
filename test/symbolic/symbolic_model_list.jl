module TestSymbolicModelList

using Test
using StaticArrays
using Plots
using Dionysos
const DI = Dionysos
const MP = DI.Mapping
const SY = DI.Symbolic

sleep(0.1)
println("Started test")

@testset "SymbolicModelList (Mapping-based)" begin
    # -----------------------
    # Build X mapping (finite explicit grid)
    # -----------------------
    x0 = SVector(0.0, 0.0)
    hx = SVector(1.0, 2.0)
    Xgrid = MP.GridFree(x0, hx)
    Xmap = MP.ExplicitGridMapping(Xgrid)

    q11 = MP.add_pos!(Xmap, (1, 1))
    q22 = MP.add_pos!(Xmap, (2, 2))
    @test MP.get_n_state(Xmap) == 2
    @test Set(MP.enum_states(Xmap)) == Set([1, 2])

    # -----------------------
    # Build U mapping (finite explicit grid)
    # -----------------------
    u0 = SVector(0.0)
    hu = SVector(0.5)
    Ugrid = MP.GridFree(u0, hu)
    Umap = MP.ExplicitGridMapping(Ugrid)

    u0id = MP.add_pos!(Umap, (0,))
    @test MP.get_n_state(Umap) == 1
    @test u0id == 1

    # -----------------------
    # Construct SymbolicModelList with defaults
    # (Xset=all states, Uset=all inputs, Rset=Xset)
    # -----------------------
    sym = SY.SymbolicModelList(Xmap, Umap)

    @test SY.get_state_mapping(sym) === Xmap
    @test SY.get_input_mapping(sym) isa MP.AbstractMapping

    @test SY.get_n_state(sym) == 2
    @test SY.get_n_input(sym) == 1
    @test Set(SY.enum_states(sym)) == Set([1, 2])
    @test Set(SY.enum_inputs(sym)) == Set([1])

    # concrete ↔ abstract conversions
    @test SY.get_concrete_state(sym, 1) == MP.get_coord_by_state(Xmap, 1)
    @test SY.get_concrete_state(sym, 2) == MP.get_coord_by_state(Xmap, 2)
    @test SY.get_concrete_input(sym, 1) ==
          MP.get_coord_by_state(SY.get_input_mapping(sym), 1)

    @test SY.get_abstract_state(sym, SVector(1.0, 2.5)) == 1
    @test SY.get_abstract_state(sym, SVector(1.5, 3.0)) == 2

    # -----------------------
    # Build a state-set from states
    # -----------------------
    S1 = SY.get_state_set_from_states(sym, [1])
    @test MP.contains_state(S1, Xmap, 1)
    @test !MP.contains_state(S1, Xmap, 2)
    @test Set(MP.enum_states(S1, Xmap)) == Set([1])

    S12 = SY.get_state_set_from_states(sym, [1, 2])
    @test Set(MP.enum_states(S12, Xmap)) == Set([1, 2])

    # -----------------------
    # Restrict Xset / Uset explicitly
    # -----------------------
    Xset_only1 = MP.stateset_from_states(Xmap, [1])
    Uset_only1 = MP.stateset_from_states(SY.get_input_mapping(sym), [1])

    sym_restricted = SY.SymbolicModelList(
        Xmap,
        Umap;
        Xset = Xset_only1,
        Uset = Uset_only1,
        Rset = Xset_only1,
    )

    @test SY.get_n_state(sym_restricted) == 1
    @test Set(SY.enum_states(sym_restricted)) == Set([1])
    @test SY.get_n_input(sym_restricted) == 1

    # -----------------------
    # Transitions + plotting recipe
    # -----------------------
    translist = [(1, 2, 1), (2, 1, 1)]  # (q′, q, u)
    SY.add_transitions!(sym.autom, translist)
    @test SY.get_n_transitions(sym) == 2

    fig = plot(; aspect_ratio = :equal)
    lyap_fun = Dict(q => 2.0 * q for q in SY.enum_states(sym))
    plot!(fig, sym; with_arrows = true, value_function = lyap_fun)  # uses your recipe kwargs
    @test isa(fig, Plots.Plot{Plots.GRBackend})

    fig2 = plot(; aspect_ratio = :equal)
    plot!(fig2, sym; with_arrows = true)
    @test isa(fig2, Plots.Plot{Plots.GRBackend})
end

sleep(0.1)
println("End test")

end # module
