module DionysosCSVExt

import Dionysos
import CSV
import DataFrames
import Dionysos: export_controller_csv, import_controller_csv

const DI = Dionysos
const MP = DI.Mapping
const ST = DI.System
const SY = DI.Symbolic
const AB = DI.Optim.Abstraction

using DataFrames: DataFrame, eachrow, ncol, Not

# --------------------------------------- #
# ----- CSV Controller Export/Import ---- #
# --------------------------------------- #

# --------------- Export ---------------- #

function export_controller_csv(
    optimizer::AB.UniformGridAbstraction.Optimizer,
    basename::String;
    delim = ';',
    decimal = ',',
)
    abs_sys = optimizer.abstraction_solver.abstract_system
    ctrl = optimizer.control_solver.abstract_controller
    ctrl === nothing && error("Controller not available")
    return export_controller_csv(abs_sys, ctrl, basename; delim = delim, decimal = decimal)
end

function export_controller_csv(
    sym::SY.SymbolicModel,
    controller,
    basename::String;
    delim = ';',
    decimal = ',',
)
    Xmap = SY.get_state_mapping(sym)

    if Xmap isa MP.GridMapping
        grid = MP.get_grid(Xmap)
        CSV.write(
            basename * "_Grid.csv",
            build_grid_df(grid);
            delim = delim,
            decimal = decimal,
        )
    end

    CSV.write(
        basename * "_StateMap.csv",
        build_state_map_df(sym);
        delim = delim,
        decimal = decimal,
    )
    CSV.write(
        basename * "_ControllerMap.csv",
        build_controller_map_df(sym, controller);
        delim = delim,
        decimal = decimal,
    )
    CSV.write(
        basename * "_InputMap.csv",
        build_input_map_df(sym);
        delim = delim,
        decimal = decimal,
    )

    return nothing
end

function build_grid_df(grid)
    origin = MP.get_origin(grid)
    h = MP.get_h(grid)
    ndims = length(origin)

    header = ["key"; ["x$(j)" for j in 1:ndims]]
    rows = [["origin"; origin], ["h"; h]]

    df = DataFrame()
    for j in eachindex(header)
        df[!, Symbol(header[j])] = getindex.(rows, j)
    end
    return df
end

function build_state_map_df(sym::SY.SymbolicModel)
    states = SY.enum_states(sym)
    isempty(states) && error("Symbolic model has no states.")

    x1 = SY.get_concrete_state(sym, first(states))
    ndims = length(x1)

    headers = ["abstract_state"; ["x$(j)" for j in 1:ndims]]
    rows = [(q, SY.get_concrete_state(sym, q)...) for q in states]

    return DataFrame([headers[i] => getindex.(rows, i) for i in eachindex(headers)])
end

function build_controller_map_df(sym::SY.SymbolicModel, controller)
    states = SY.enum_states(sym)
    rows = [(q, get_input_symbol(controller, q)) for q in states]
    return DataFrame([
        "abstract_state" => getindex.(rows, 1),
        "abstract_input" => getindex.(rows, 2),
    ])
end

# The abstract input symbol the controller assigns to abstract state `state`, or `-1` when
# the controller is undefined there. Routes through the controller protocol
# (`output_control`) so it works for any `AbstractController` (static or dynamic); a static
# controller ignores the `nothing` memory argument.
function get_input_symbol(controller, state)
    u = ST.output_control(controller, nothing, state)
    return u === nothing ? -1 : u
end

function build_input_map_df(sym::SY.SymbolicModel)
    inputs = SY.enum_inputs(sym)
    isempty(inputs) && error("Symbolic model has no inputs.")

    u1 = SY.get_concrete_input(sym, first(inputs))
    ndims_u = length(u1)

    headers = ["abstract_input"; ["u$(j)" for j in 1:ndims_u]]
    rows = [(i, SY.get_concrete_input(sym, i)...) for i in inputs]

    return DataFrame([headers[i] => getindex.(rows, i) for i in eachindex(headers)])
end

# --------------- Import ---------------- #

function import_controller_csv(basename::String; delim = ';', decimal = ',')
    grid_path = basename * "_Grid.csv"
    state_path = basename * "_StateMap.csv"
    ctrl_path = basename * "_ControllerMap.csv"
    input_path = basename * "_InputMap.csv"

    grid_df =
        isfile(grid_path) ?
        CSV.read(grid_path, DataFrame; delim = delim, decimal = decimal) : nothing
    state_df = CSV.read(state_path, DataFrame; delim = delim, decimal = decimal)
    ctrl_df = CSV.read(ctrl_path, DataFrame; delim = delim, decimal = decimal)
    input_df = CSV.read(input_path, DataFrame; delim = delim, decimal = decimal)

    return parse_controller_tables(grid_df, state_df, ctrl_df, input_df)
end

function parse_controller_tables(grid_df, state_df, ctrl_df, input_df)
    origin = nothing
    h = nothing

    if grid_df !== nothing
        origin = Vector{Float64}(grid_df[grid_df.key .== "origin", Not(:key)][1, :])
        h = Vector{Float64}(grid_df[grid_df.key .== "h", Not(:key)][1, :])
    end

    ndims_x = ncol(state_df) - 1
    coord2state = Dict{NTuple{ndims_x, Float64}, Int}()

    for row in eachrow(state_df)
        x = ntuple(i -> Float64(row[Symbol("x$i")]), ndims_x)
        coord2state[x] = Int(row.abstract_state)
    end

    state2input = Dict(Int.(ctrl_df.abstract_state) .=> Int.(ctrl_df.abstract_input))

    ndims_u = ncol(input_df) - 1
    input2u = Dict{Int, NTuple{ndims_u, Float64}}()
    for row in eachrow(input_df)
        u = ntuple(i -> Float64(row[Symbol("u$i")]), ndims_u)
        input2u[Int(row.abstract_input)] = u
    end

    return origin, h, coord2state, state2input, input2u
end

end
