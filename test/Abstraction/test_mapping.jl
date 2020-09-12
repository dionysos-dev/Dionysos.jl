include(joinpath(@__DIR__, "../../src/Abstraction/abstraction.jl"))

module TestMain

using Test
using Main.Abstraction
ABS = Main.Abstraction
using StaticArrays

println("")

@testset "Abstraction/mapping" begin
map = ABS.ListGridMap(String)
append!(map.coll, ["a", "b", "c", "d"])
get_symbol = ABS.cell2symbol(map)
get_symbol_try = ABS.cell2symbol_try(map)
get_cell = ABS.symbol2cell(map)
ST = ABS._symboltype(typeof(map))
symb = get_symbol("a")
# Main.@code_warntype hash(symb, UInt(0))
@test symb === ST(1)
@test get_symbol_try("a") === (ST(1), true)
@test !get_symbol_try("e")[2]
@test get_cell(ST(2)) === "b"
@test ABS._nsymbols(map) === 4
print("")
end

end  # module TestMain
