include(joinpath("..", "Abstraction", "abstraction.jl"))
#include(joinpath("..", "CDC2021","CDC2021.jl"))
include(joinpath("..", "CDC2021","general_domain.jl"))
include(joinpath("..", "CDC2021","utils.jl"))
include(joinpath("..", "CDC2021","lazy_abstraction.jl"))
include(joinpath("..", "CDC2021","alternating_simulation.jl"))
module TestM2
using Test, StaticArrays,Plots, LinearAlgebra, LazySets
using Distributions

using ..Abstraction
const AB = Abstraction

using ..DomainList
const D = DomainList

using ..Utils
const UT = Utils

using ..Lazy_abstraction
const LA = Lazy_abstraction

using ..AlternatingSimulation
const AS = AlternatingSimulation

using ForwardDiff, IntervalArithmetic
using HybridSystems
include("sorted_vector_set.jl")
include("automaton.jl")
include("domain.jl")
include("markov_chain.jl")
include("symbolic.jl")
include("algo.jl")
include("system.jl")


function get_pre(symmodel,sys,source,u,tstep)
    N = 500
    (l,xpos) = get_xpos_by_state(symmodel, source)
    points = sample_elem(symmodel.Xdom.domains[l].grid, xpos, N)
    targetlist = Int[]
    for x in points
        Fx = sys.reverse_sys_map(x, u, tstep)
        target = get_cell(symmodel,Fx)
        push!(targetlist,target)
    end
    unique!(targetlist)
    return targetlist
end

function pre_image(symmodel,contsys,source,u,tstep)
    # (l,xpos) = get_xpos_by_state(symmodel, source)
    # Xdom = symmodel.Xdom.domains[l]
    # x = AB.get_coord_by_pos(Xdom.grid, xpos)
    # #tstep = contsys.tstep #############
    # Fx = contsys.reverse_sys_map(x, u, tstep)
    # r = Xdom.grid.h/2.0 + contsys.measnoise
    # Fr = contsys.reverse_growthbound_map(r, u, tstep, x, nstep=10)
    # println(Fr)
    # post_rec = AB.HyperRectangle(Fx .- Fr, Fx .+ Fr)
    # rectI = AB.get_pos_lims_outer(Xdom.grid, AB.HyperRectangle(Fx .- Fr, Fx .+ Fr))
    # ypos_iter = Iterators.product(AB._ranges(rectI)...)
    # over_approx = Int[]
    # for ypos in ypos_iter
    #     ypos = D.set_in_period_pos(Xdom,ypos)
    #     if !(ypos in Xdom)
    #         empty!(over_approx)
    #         break
    #     end
    #     target = get_state_by_xpos(symmodel, ypos, l)
    #     push!(over_approx, target)
    # end


    over_approx = get_pre(symmodel,contsys,source,u,tstep)
    return over_approx
end



function build_dom(f,fi,hx,r)
    X = AB.HyperRectangle(SVector(-r, -r), SVector(r, r))
    obstacle = AB.HyperRectangle(SVector(-100.0, -100.0), SVector(-99.0, -99.0))
    d = D.RectanglularObstacles(X, [obstacle])
    Xdom = D.GeneralDomainList(hx,d;fit=false)
    return X,Xdom
end

function build_system(X)
    tstep = 0.08
    measnoise = SVector(0.0,0.0)
    sysnoise = SVector(0.0,0.0)
    nsys = 5
    ngrowthbound = 5

    μ = 2.0
    k = 1.0
    f = system_with_boundary(Van_der_pol_oscillator(μ;k=k),10.0,-0.8)
    f_reverse = reverse_Van_der_pol_oscillator(μ;k=k)
    # f2 = Van_der_pol_oscillator(2.0)
    jacobian = Jacobian_with_boundary(Van_der_pol_oscillator_Jacobian(2.0),12.0,-0.8)
    contsys = NewControlSystemGrowthLx(tstep, f, jacobian, f_reverse, measnoise, sysnoise, nsys, ngrowthbound, X)
    return contsys
end


function eliminate!(prob::AnalysisProblem)
    symmodel = prob.symmodel
    Xdom = symmodel.Xdom
    mc = symmodel.mc
    indexes = sortperm(mc.steady_state)
    d = 0
    bool2 = true
    for (i,index) in enumerate(indexes)
        s = mc.symmodel_from_mc[index]
        if get_steady_state(mc, s) > 0.0
            u = SVector(0.0,0.0)
            tstep = 0.08*1.8 #1.8 #4.0
            cells = prob.pre_image(symmodel,prob.sys,s,u,tstep)
            println(cells)
            #plot_some(symmodel,s,cells)

            # println(cells)
            bool2 = true
            for cell in cells
                #println(get_steady_state(mc, cell))
                if cell != 0 && get_steady_state(mc, cell) > 0.0
                    bool2 = false
                    break
                end
            end
            if bool2
                d = s
                break
            end
            # if i == 10
            #     break
            # end
        end
    end
    if !bool2
        println("Nothing found")
        return false
    else
        i = mc.mc_from_symmodel[d]
        println(d)
        # println(get_steady_state(mc, d))
        mc.steady_state[i] = 0.0

        # fig = plot(aspect_ratio = 1,legend = false,title="steady-state probability")
        # plot_steady_state!(symmodel,fact=0.0)
        # (l,xpos) = get_xpos_by_state(symmodel, d)
        # dom = Xdom.domains[l]
        # AB.plot_elem!(dom.grid, xpos,opacity=1.0,color=:red)
        # display(fig)
        return true
    end

end

function test()
    f_grid = nothing
    f_grid_i = nothing
    hx = [1.0, 1.0]*0.3#0.3
    X,Xdom = build_dom(f_grid,f_grid_i,hx,12.0) #10.0
    Udom = UT.CustomList([SVector(0.0,1.0)])
    sys =  build_system(X)

    symmodel = NewMultiSymbolicModel(Xdom, Udom)
    prob = NewAnalysisProblem(symmodel,sys,nothing,pre_image,false)
    time = @elapsed begin
        symmodel = symmodel
        println("Symbolic model")
        compute_symbolic_full_domain!(prob,symmodel,sys)
        println("Markov chain")
        update_MC!(symmodel)
        s = get_state_by_coord(symmodel,SVector(0.0,0.0))
        println(get_steady_state(symmodel.mc,s))
    end
    println("total time: lazy abstraction + controller: ", time)
    println("Eliminate")
    fig = plot(aspect_ratio = 1,legend = false,title="steady-state probability")
    plot_steady_state!(prob.symmodel,fact=0.0)
    for k in 1:1000
        println(k)
        if !eliminate!(prob)
            println("cells deleted: ",k)
            break
        end
    end
    plot_steady_state!(prob.symmodel,fact=0.0,color=:blue)
    x0 = SVector(2.0,2.0)
    plot_trajectory!(sys,x0,200)
    display(fig)
end

function test_2()
    tstep = 0.08
    measnoise = SVector(0.0,0.0)
    sysnoise = SVector(0.0,0.0)
    nsys = 5
    μ = 2.0
    k = 1.0

    f = system_with_boundary(Van_der_pol_oscillator(μ;k=k),10.0,-0.8)
    f_reverse = reverse_Van_der_pol_oscillator(μ;k=k)
    sys1 = NewControlSystemGrowthLx(tstep, f, nothing, nothing, measnoise, sysnoise, nsys, nothing, nothing)
    sys2 = NewControlSystemGrowthLx(tstep, f_reverse, nothing, nothing, measnoise, sysnoise, nsys, nothing, nothing)
    fig = plot(aspect_ratio = 1,legend = false,title="Simu")
    plot_trajectory!(sys1,SVector(2.0,2.0),300)
    # plot_trajectory!(sys2,SVector(2.0,2.0),13;color=:blue)
    # plot_trajectory!(sys2,SVector(2.1,2.1),10;color=:blue)
    # plot_trajectory!(sys2,SVector(1.0,1.0),100;color=:blue)
    # plot_trajectory!(sys2,SVector(1.0,-1.0),14;color=:blue)
    # plot_trajectory!(sys2,SVector(-1.0,-1.0),100;color=:blue)
    plot_trajectory!(sys1,SVector(0.01,0.01),100;color=:blue)
    plot_trajectory!(sys1,SVector(0.01,0.0),100;color=:blue)
    plot_trajectory!(sys1,SVector(0.0,0.01),100;color=:blue)
    display(fig)
end

test()
end








# function Van_der_pol_oscillator(μ)
#     function f(x::SVector{N,T}, u) where {N,T}
#         y1 = x[2]
#         y2 = μ*(1-x[1]*x[1])*x[2]-x[1]
#         y = SVector(y1,y2)
#         return y
#     end
#     return f
# end
#
# function system_with_boundary(f1,radius,a)
#     function f(x::SVector{N,T}, u) where {N,T}
#         y = f1(x,u)
#         if norm(x,2)>radius
#             A = @SMatrix [ a 0.0 ;
#                            0.0 a ]
#             y = A*x
#         end
#         return y
#     end
#     return f
# end
#
#
# function f1(x)
#     return x
# end
# function fi1(x)
#     return x
# end
#
# #counter-clockwise
# function rotate(x,θ)
#     R = @SMatrix [ cos(θ) -sin(θ) ;
#                    sin(θ)  cos(θ)]
#     return R*x
# end
# function f2(x)
#     c =  SVector(0.0,0.0) #SVector(2.0,3.0)
#     θ = π/3.0
#     return rotate(x-c,θ)+c
# end
# function fi2(x)
#     c =  SVector(0.0,0.0) #SVector(2.0,3.0)
#     θ = π/3.0
#     return rotate(x-c,-θ)+c
# end
#
# function f3(x)
#     return SVector(x[2]+sin(x[1]),x[1])
# end
# function fi3(x)
#     return SVector(x[2],x[1]-sin(x[2]))
# end
#
# function f4(x)
#     return SVector(x[1]*cos(x[2]),x[1]*sin(x[2]))
# end
# function fi4(x)
#     return SVector(sqrt(x[1]*x[1] + x[2]*x[2]),atan(x[2],x[1]))
# end
#
# function build_f_rotation(θ,c =  SVector(0.0,0.0))
#     function f(x)
#         return rotate(x-c,θ)+c
#     end
#     function fi(x)
#         return rotate(x-c,-θ)+c
#     end
#     return f,fi
# end
#
#
#
#
# struct SimpleSystem{N,T,F<:Function,F2<:Function} <: AB.ControlSystem{N,T}
#     tstep::Float64
#     measnoise::SVector{N,T}
#     sys_map::F
#     f::F2
# end
#
# function RungeKutta4(F, x, u, tstep, nsub::Int)
#     τ = tstep/nsub
#     for i in 1:nsub
#         Fx1 = F(x, u)
#         xrk = x + Fx1*(τ/2.0)
#         Fx2 = F(xrk, u)
#         xrk = x + Fx2*(τ/2.0)
#         Fx3 = F(xrk, u)
#         xrk = x + Fx3*τ
#         Fx4 = F(xrk, u)
#         x = x + (Fx1 + Fx2*2.0 + Fx3*2.0 + Fx4)*(τ/6.0)
#     end
#     return x
# end
#
# function NewControlSystemGrowthRK4(tstep, F_sys, measnoise::SVector{N,T}, nsys) where {N,T}
#     sys_map = let nsys = nsys
#         (x::SVector{N,T}, u, tstep) ->
#             RungeKutta4(F_sys, x, u, tstep, nsys)::SVector{N,T}
#     end
#     return SimpleSystem(tstep, measnoise, sys_map, F_sys)
# end
#
# function build_system_Van_der_pol_oscillator(μ = 2.0,tstep = 0.08,measnoise = SVector(0.0, 0.0),nsys = 4)
#     function F_sys(x::SVector{N,T}, u) where {N,T}
#         y1 = x[2]
#         y2 = μ*(1-x[1]*x[1])*x[2]-x[1]
#         y = SVector(y1,y2)
#         if norm(x,2)>8.0
#             a = -0.8
#             b = 0.0
#             A = @SMatrix [ a -b ;
#                            b a ]
#             y = A*x
#         end
#         return y
#     end
#     return NewControlSystemGrowthRK4(tstep,F_sys,measnoise,nsys)
# end
