include(joinpath("..", "..", "Abstraction", "abstraction.jl"))
include("../general_domain.jl")
include("../utils.jl")
include("../alternating_simulation.jl")
include("../partition.jl")
include("../lazy_abstraction.jl")
include("../branch_and_bound.jl")
include("../optimal_control.jl")

module Test
"""
    Simple two-dimensional reachability problem
"""

using Plots, StaticArrays, JuMP

using ..Abstraction
const AB = Abstraction

using ..DomainList
const D = DomainList

using ..Utils
const U = Utils

using ..BranchAndBound
const BB = BranchAndBound

using ..OptimalControl
const OC = OptimalControl

## problem specific-functions required by the algo
function compute_reachable_set(rect::AB.HyperRectangle{SVector{N,T}},contsys,Udom) where {N,T}
    tstep = contsys.tstep
    r = (rect.ub-rect.lb)/2.0 + contsys.measnoise
    Fr = r
    x = U.center(rect)
    n =  U.dims(rect)
    lb = SVector(ntuple(i -> Inf, Val(N)))
    ub = SVector(ntuple(i -> -Inf, Val(N)))
    for upos in AB.enum_pos(Udom)
        u = AB.get_coord_by_pos(Udom.grid,upos)
        Fx = contsys.sys_map(x, u, tstep)
        lb = min.(lb,Fx .- Fr)::SVector{N,T}
        ub = max.(ub,Fx .+ Fr)::SVector{N,T}
    end
    return AB.HyperRectangle(lb,ub)
end
function minimum_transition_cost(symmodel,contsys,source,target)
    return 1.0
end
function post_image(symmodel,contsys,xpos,u)
    Xdom = symmodel.Xdom
    x = AB.get_coord_by_pos(Xdom.grid, xpos)
    tstep = contsys.tstep
    Fx = contsys.sys_map(x, u, tstep)
    r = Xdom.grid.h/2.0 + contsys.measnoise
    Fr = r

    rectI = AB.get_pos_lims_outer(Xdom.grid, AB.HyperRectangle(Fx .- Fr, Fx .+ Fr))
    ypos_iter = Iterators.product(AB._ranges(rectI)...)
    over_approx = Int[]
    for ypos in ypos_iter
        ypos = D.set_in_period_pos(Xdom,ypos)
        if !(ypos in Xdom)
            empty!(over_approx)
            break
        end
        target = AB.get_state_by_xpos(symmodel, ypos)
        push!(over_approx, target)
    end
    return over_approx
end

function pre_image(symmodel,contsys,xpos,u)
    grid = symmodel.Xdom.grid
    x = AB.get_coord_by_pos(grid, xpos)
    tstep = contsys.tstep
    potential = Int[]
    x_prev = x-tstep*u
    xpos_cell = AB.get_pos_by_coord(grid,x_prev)
    n = 2
    for i=-n:n
        for j=-n:n
            x_n = (xpos_cell[1]+i,xpos_cell[2]+j)
            x_n = D.set_in_period_pos(symmodel.Xdom,x_n)
            if x_n in symmodel.Xdom
                cell = AB.get_state_by_xpos(symmodel, x_n)
                if !(cell in potential)
                    push!(potential,cell)
                end
            end
        end
    end
    return potential
end
## problem data

struct SimpleSystem{N,T,F<:Function} <: AB.ControlSystem{N,T}
    tstep::Float64
    measnoise::SVector{N,T}
    sys_map::F
end

function NewSimpleControlSystem(tstep,measnoise)
    function sys_map(x, u, tstep)
        return x+tstep*u
    end
    return SimpleSystem(tstep, measnoise, sys_map)
end

function build_system()
    tstep = 0.8
    measnoise = SVector(0.0, 0.0)
    return NewSimpleControlSystem(tstep,measnoise)
end

function build_input()
    U = AB.HyperRectangle(SVector(-2.0, -2.0), SVector(2.0, 2.0))
    x0 = SVector(0.0, 0.0); hu = SVector(0.5,0.5)
    Ugrid = AB.GridFree(x0,hu)
    Udom = AB.DomainList(Ugrid)
    AB.add_set!(Udom, U, AB.OUTER)
    box = AB.HyperRectangle(SVector(-0.5, -0.5), SVector(0.5, 0.5))
    AB.remove_set!(Udom, box, AB.OUTER)
    return Udom
end

function transition_cost(x,u)
    return 1.0
end

##

function test()
    X = AB.HyperRectangle(SVector(0.0,0.0),SVector(60.0,60.0))
    contsys = build_system()
    Udom = build_input()

    # control problem
    _T_ = AB.HyperRectangle(SVector(44.0, 38.0), SVector(49.0, 41.0))
    #_T_ = AB.HyperRectangle(SVector(5.0, 10.0), SVector(8.0, 12.0))
    x0 = SVector(7.0,7.0)
    _I_ = AB.HyperRectangle(SVector(6.5, 6.5), SVector(7.5, 7.5))
    # problem-specific functions

    hx_coarse = SVector(10.0, 10.0)*1.0
    hx_medium = SVector(1.0, 1.0)*0.5*2.0
    hx_fine = SVector(0.5, 0.5)*1.0
    periodic = [1,2]
    periods = Float64[60.0,60.0]
    T0 = Float64[0.0,0.0]
    fig = plot(aspect_ratio = 1,legend = false)
    plot!(U.rectangle(_I_.lb,_I_.ub), opacity=.8,color=:green)
    plot!(U.rectangle(_T_.lb,_T_.ub), opacity=.8,color=:red)
    functions = (compute_reachable_set,minimum_transition_cost,post_image,pre_image)

    optimal_control_prob = OC.OptimalControlProblem(x0,_I_,_T_,contsys,periodic,periods,T0,Udom,transition_cost,(X,hx_coarse,[]),hx_medium,hx_fine,functions)
    max_iter = 5
    max_time = 1000
    optimizer = BB.Optimizer(optimal_control_prob,max_iter,max_time,log_level=2)
    @time MOI.optimize!(optimizer)
    println(optimizer.status)
    println(optimizer.best_sol)
    display(fig)
    (traj,cost,sucess) = OC.simulate_trajectory(optimal_control_prob, optimizer.best_sol)
    #println(traj)
end


test()
end # module




#=
function build_partition(α)
    H1 = AB.HyperRectangle(SVector(0.0,0.0),SVector(8.0,4.0))
    H2 = AB.HyperRectangle(SVector(8.0,0.0),SVector(12.0,4.0))
    H3 = AB.HyperRectangle(SVector(12.0,0.0),SVector(16.0,4.0))
    H4 = AB.HyperRectangle(SVector(0.0,4.0),SVector(3.0,9.0))
    H5 = AB.HyperRectangle(SVector(3.0,4.0),SVector(10.0,9.0))
    H6 = AB.HyperRectangle(SVector(10.0,4.0),SVector(16.0,7.0))
    H7 = AB.HyperRectangle(SVector(0.0,9.0),SVector(5.0,13.0))
    H8 = AB.HyperRectangle(SVector(5.0,9.0),SVector(10.0,16.0))
    H9 = AB.HyperRectangle(SVector(12.0,9.0),SVector(16.0,14.0))
    H10 = AB.HyperRectangle(SVector(0.0,13.0),SVector(5.0,16.0))
    H11 = AB.HyperRectangle(SVector(10.0,14.0),SVector(16.0,16.0))
    H12 = AB.HyperRectangle(SVector(7.0,16.0),SVector(14.0,20.0))
    H = [H1,H2,H3,H4,H5,H6,H7,H8,H9,H10,H11,H12]

    O1 = AB.HyperRectangle(SVector(10.0,7.0),SVector(16.0,9.0))
    O2 = AB.HyperRectangle(SVector(10.0,9.0),SVector(12.0,14.0))
    O = [O1,O2]
    H = [U.scale(h,α) for h in H]
    O = [U.scale(o,α) for o in O]
    fig = plot(aspect_ratio = 1,legend = false)
    for h in H
        plot!(U.rectangle(h.lb,h.ub), opacity=.4,color=:blue)
    end
    for o in O
        plot!(U.rectangle(o.lb,o.ub), opacity=.4,color=:black)
    end
    display(fig)
    return H,O
end

=#
