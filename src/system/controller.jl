abstract type Controller end

# Returns all admissible controls at `state`, if multiple exist.
function get_all_controls(controller::Controller, state) end
# Returns a control input for the given state.
function get_control(controller::Controller, state) end
# Returns `true` if the controller is defined at the given state.
function is_defined(controller::Controller, state) end
# Returns the set of all states for which the controller is defined.
function domain(controller::Controller) end

abstract type SymbolicController <: Controller end

function add_control!(controller::SymbolicController, state::Int, symbol::Int) end
# Delete previous control law for 'state' and 'symbol'
function set_control!(controller::SymbolicController, state::Int, symbol::Int) end

#################################################
############ Symbolic implementations ###########
#################################################
abstract type SymbolicController <: Controller end
function add_control!(controller::SymbolicController, state::Int, symbol::Int) end
get_control(controller::SymbolicController, state::Int) =
    first(get_all_controls(controller, state))

struct SymbolicControllerList <: SymbolicController
    transitions::UT.SortedTupleSet{2, NTuple{2, Int}}  # (state, symbol)
end

SymbolicControllerList() = SymbolicControllerList(UT.SortedTupleSet{2, NTuple{2, Int}}())

get_all_controls(controller::SymbolicControllerList, state::Int) = collect(
    UT.fix_and_eliminate_first(controller.transitions, state; drop = UT.drop_first_int),
)
is_defined(controller::SymbolicControllerList, state::Int) =
    !isempty(get_all_controls(controller, state))
domain(controller::SymbolicControllerList) =
    unique(map(first, UT.get_data(controller.transitions)))
add_control!(controller::SymbolicControllerList, state::Int, symbol::Int) =
    UT.push_new!(controller.transitions, (state, symbol))

function set_control!(controller::SymbolicControllerList, state::Int, symbol::Int)
    delete!(controller.transitions, (state,), (a, b) -> a[1] == b[1])
    return add_control!(controller, state, symbol)
end

struct SymbolicControllerDict <: SymbolicController
    control_map::Dict{Int, Vector{Int}}  # state => list of symbols
end

SymbolicControllerDict() = SymbolicControllerDict(Dict{Int, Vector{Int}}())

get_all_controls(controller::SymbolicControllerDict, state::Int) =
    controller.control_map[state]
is_defined(controller::SymbolicControllerDict, state::Int) =
    haskey(controller.control_map, state)
domain(controller::SymbolicControllerDict) = keys(controller.control_map)
add_control!(controller::SymbolicControllerDict, state::Int, symbol::Int) =
    push!(get!(controller.control_map, state, Int[]), symbol)

function set_control!(controller::SymbolicControllerDict, state::Int, symbol::Int)
    return controller.control_map[state] = [symbol]
end
##########################################################
########## Continuous Concrete implementations ###########
##########################################################
abstract type ContinuousController <: Controller end

struct BlackBoxContinuousController{F, G} <: ContinuousController
    f::F
    is_defined::G
end
BlackBoxContinuousController(f::F) where {F} = BlackBoxContinuousController(f, x -> true)

get_control(controller::BlackBoxContinuousController, x) = controller.f(x)
get_all_controls(controller::BlackBoxContinuousController, x) = [get_control(controller, x)]
is_defined(controller::BlackBoxContinuousController, x) = controller.is_defined(x)
domain(controller::BlackBoxContinuousController) = nothing

struct ConstantController{VT} <: ContinuousController
    c::VT
end
get_control(controller::ConstantController, x) = controller.c
get_all_controls(controller::ConstantController, x) = [get_control(controller, x)]
is_defined(controller::ConstantController, x) = true
domain(controller::ConstantController) = nothing

struct AffineController{MT, VT1, VT2} <: ContinuousController
    K::MT
    c::VT1
    ℓ::VT2
end
get_control(controller::AffineController, x::AbstractVector) =
    controller.K * (x - controller.c) + controller.ℓ
get_all_controls(controller::AffineController, x::AbstractVector) =
    [get_control(controller, x)]
is_defined(controller::AffineController, x::AbstractVector) = true
domain(controller::AffineController) = nothing

##########################################################
############## Sub-routine for debugging #################
##########################################################

struct ControllerAnalysis{C <: Controller, F, W, S1, S2, CE}
    controller::C
    f_eval::F
    Wset::W
    domain_set::S1
    target_set::Union{Nothing, S2}
    cost_eval::Union{Nothing, CE}
    dims::Vector{Int}
    N::Int
end

function ControllerAnalysis(
    controller::C,
    f_eval,
    Wset,
    domain_set;
    target_set = nothing,
    cost_eval = nothing,
    dims = [1, 2],
    N = 100,
) where {C <: Controller}
    return ControllerAnalysis(
        controller,
        f_eval,
        Wset,
        domain_set,
        target_set,
        cost_eval,
        dims,
        N,
    )
end

@recipe function f(
    analysis::ControllerAnalysis;
    arrowsB = true,
    cost = true,
    sphere_radius = 0.01,
)
    dims = analysis.dims
    samples = UT.sample(analysis.domain_set; N = analysis.N)

    controller = analysis.controller
    nx = UT.get_dims(analysis.domain_set)
    # Show domain and target sets
    @series begin
        color := :green
        return analysis.domain_set
    end

    if analysis.target_set !== nothing
        @series begin
            color := :red
            return analysis.target_set
        end
    end

    # Compute cost (if enabled)
    if cost && analysis.cost_eval !== nothing
        costs = Float64[]
        for x in samples
            u = get_control(controller, x)
            push!(costs, analysis.cost_eval(x, u))
        end
        vmin, vmax = extrema(costs)
        mycolorMap = UT.Colormap([vmin, vmax], Colors.colormap("Blues"))
        P = (1.0 / sphere_radius) * Matrix{Float64}(I(2))  # small ellipsoid

        for (i, x) in enumerate(samples)
            color = UT.get_color(mycolorMap, costs[i])
            ell = UT.Ellipsoid(P, x[dims])
            @series begin
                color := color
                lw := 0
                return ell
            end
        end
        @series begin
            mycolorMap
        end
    end

    if arrowsB
        for x in samples
            u = get_control(controller, x)
            w = UT.sample(analysis.Wset)
            x2 = analysis.f_eval(x, u, w)

            @series begin
                color := :black
                return UT.DrawArrow(SVector{nx}(x[dims]), SVector{nx}(x2[dims]))
            end
        end
    end
end

# data-driven check
function check_feasibility(
    ell1,
    ell2,
    f_eval,
    controller::Controller,
    Uset,
    Wset;
    N = 500,
    input_check = true,
    noise_check = true,
)
    samples = UT.sample(ell1; N = N)
    nw = UT.get_dims(Wset)
    for x in samples
        unew = get_control(controller, x)
        if input_check && !(unew ∈ Uset)
            println("Not feasible input")
            return false
        end
        noise_check ? wnew = UT.sample(Wset) : wnew = zeros(nw)
        xnew = f_eval(x, unew, wnew)
        if !(xnew ∈ ell2)
            println("Not in the target ellipsoid")
            return false
        end
    end
    return true
end
