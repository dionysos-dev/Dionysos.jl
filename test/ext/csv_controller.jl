module TestMain

import Dionysos
include(joinpath(dirname(dirname(pathof(Dionysos))), "test", "testsetup.jl"))

import CSV
import DataFrames   # importing both loads the DionysosCSVExt extension
import LazySets

# The DionysosCSVExt controller export/import: a full export → import round-trip against a
# real symbolic model + `DiscreteStaticController`, plus a hand-built import case exercising
# the grid table and the parsing path in isolation.
const CSVExt = Base.get_extension(Dionysos, :DionysosCSVExt)

@testset "CSV controller import round-trip" begin
    @test CSVExt !== nothing   # extension is loaded once CSV + DataFrames are

    tmp = mktempdir()
    base = joinpath(tmp, "ctrl")

    # Grid table via the extension's own builder, then write it as export would.
    grid = MP.GridFree(SVector(0.0, 0.0), SVector(1.0, 2.0))
    CSV.write(base * "_Grid.csv", CSVExt.build_grid_df(grid); delim = ';', decimal = ',')

    # The remaining three tables in the on-disk export format.
    CSV.write(
        base * "_StateMap.csv",
        DataFrames.DataFrame(; abstract_state = [1, 2], x1 = [0.0, 1.0], x2 = [0.0, 2.0]);
        delim = ';',
        decimal = ',',
    )
    CSV.write(
        base * "_ControllerMap.csv",
        DataFrames.DataFrame(; abstract_state = [1, 2], abstract_input = [10, 20]);
        delim = ';',
        decimal = ',',
    )
    CSV.write(
        base * "_InputMap.csv",
        DataFrames.DataFrame(; abstract_input = [10, 20], u1 = [-1.0, 1.0]);
        delim = ';',
        decimal = ',',
    )

    origin, h, coord2state, state2input, input2u =
        Dionysos.import_controller_csv(base; delim = ';', decimal = ',')

    @test origin == [0.0, 0.0]
    @test h == [1.0, 2.0]
    @test coord2state[(0.0, 0.0)] == 1
    @test coord2state[(1.0, 2.0)] == 2
    @test state2input == Dict(1 => 10, 2 => 20)
    @test input2u[10] == (-1.0,)
    @test input2u[20] == (1.0,)
end

@testset "CSV controller export/import round-trip" begin
    # A small grid-backed symbolic model with three inputs.
    Xgrid = MP.GridFree(SVector(0.0, 0.0), SVector(1.0, 1.0))
    Xmap = MP.ImplicitGridMapping(
        Xgrid,
        LazySets.Hyperrectangle(; low = SVector(0.0, 0.0), high = SVector(2.0, 2.0)),
    )
    Umap = MP.ListMapping([SVector(-1.0), SVector(0.0), SVector(1.0)])
    symmodel = SY.SymbolicModelList(
        Xmap,
        Umap;
        Xset = MP.FullStateSet{2}(),
        Rset = nothing,
        Uset = MP.FullStateSet{1}(),
        convert_U_to_list = false,
    )

    states = collect(SY.enum_states(symmodel))
    inputs = collect(SY.enum_inputs(symmodel))
    @test length(states) >= 2 && length(inputs) == 3

    # A static controller that acts on exactly two states.
    table = ST.ControlTable(MP.get_n_state(Xmap))
    ST.set_control!(table, states[1], inputs[1])
    ST.set_control!(table, states[2], inputs[3])
    controller = ST.DiscreteStaticController(Set([states[1], states[2]]), table, false)

    base = joinpath(mktempdir(), "ctrl")
    Dionysos.export_controller_csv(symmodel, controller, base)

    @test isfile(base * "_Grid.csv")
    @test isfile(base * "_StateMap.csv")
    @test isfile(base * "_ControllerMap.csv")
    @test isfile(base * "_InputMap.csv")

    origin, h, coord2state, state2input, input2u = Dionysos.import_controller_csv(base)

    @test origin == [0.0, 0.0]
    @test h == [1.0, 1.0]
    # The protocol-based `get_input_symbol` exports the real symbols of controlled states and
    # `-1` elsewhere (the old `controller.X`/`.h` access would have crashed the export).
    @test state2input[states[1]] == inputs[1]
    @test state2input[states[2]] == inputs[3]
    @test all(state2input[q] == -1 for q in states if !(q in (states[1], states[2])))
    # Each exported symbol resolves back to its concrete input.
    @test input2u[inputs[1]] == (-1.0,)
    @test input2u[inputs[3]] == (1.0,)
end

end # module TestMain
