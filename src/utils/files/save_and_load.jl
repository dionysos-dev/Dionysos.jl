using JLD2

function save_problem(filename, abstract_system, abstract_problem, abstract_controller, value_function)
    jldsave(filename; abstract_system, abstract_problem, abstract_controller, value_function)
    return 
end

function rebuild_concrete_controller(abstract_system, abstract_controller)
    function concrete_controller(x; param = false)
        xpos = DO.get_pos_by_coord(abstract_system.Xdom.grid, x)
        if !(xpos ∈ abstract_system.Xdom)
            @warn("State out of domain")
            return nothing
        end
        source = SY.get_state_by_xpos(abstract_system, xpos)
        symbollist = UT.fix_and_eliminate_first(abstract_controller, source)
        if isempty(symbollist)
            @warn("Uncontrollable state")
            return nothing
        end
        if param
            temp = rand(collect(symbollist))
            symbol, tstep = temp[1], temp[2]
        else
            temp = first(symbollist)
            symbol, tstep = temp[1], temp[2] # symbol, timestep
        end
        upos = SY.get_upos_by_symbol(abstract_system, symbol) 
        u = DO.get_coord_by_pos(abstract_system.Udom.grid, upos)
        return u, tstep
    end
end

function load_problem(filename) 
    file = jldopen(filename, "r")
    abstract_system = file["abstract_system"]
    abstract_problem = file["abstract_problem"]
    abstract_controller = file["abstract_controller"]
    value_function = file["value_function"]
    concrete_controller = rebuild_concrete_controller(abstract_system, abstract_controller)
    return abstract_system, abstract_problem, abstract_controller, value_function, concrete_controller
end