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
# include("algo.jl")
include("multigrid.jl")
include("system.jl")


function build_dom(f,fi,hx,r)
    X = AB.HyperRectangle(SVector(-r, -r), SVector(r, r))
    obstacle = AB.HyperRectangle(SVector(-100.0, -100.0), SVector(-99.0, -99.0))
    d = D.RectanglularObstacles(X, [obstacle])
    Xdom = D.GeneralDomainList(hx,d;fit=true)  #il faut fit = true !!!
    return X,Xdom
end

function build_system(X,μ)
    tstep = 0.08
    measnoise = SVector(0.0,0.0)
    sysnoise = SVector(0.0,0.0)
    nsys = 5
    ngrowthbound = 5

    #μ = 2.0 #2.0
    k = 1.0
    f = system_with_boundary(Van_der_pol_oscillator(μ;k=k),10.0,-0.8)
    f_reverse = reverse_Van_der_pol_oscillator(μ;k=k)
    # f2 = Van_der_pol_oscillator(2.0)
    jacobian = Jacobian_with_boundary(Van_der_pol_oscillator_Jacobian(2.0),12.0,-0.8)
    contsys = NewControlSystemGrowthLx(tstep, f, jacobian, f_reverse, measnoise, sysnoise, nsys, ngrowthbound, X)
    return contsys
end

function test_limit_cycle()
    μ = 2.0
    hx = [1.0, 1.0]*0.3  #0.3
    X,Xdom = build_dom(nothing,nothing,hx,12.0) #10.0
    Udom = UT.CustomList([SVector(0.0,1.0)])
    sys =  build_system(X,μ)
    symmodel = NewMultiSymbolicModel(Xdom, Udom)
    compute_symbolic_full_domain!(symmodel,sys)
    update_MC!(symmodel)
    plot_steady_state(sys,symmodel,xlims=[-3.0,3.0],ylims=[-4.0,4.0])
    # plot_steady_state(sys,symmodel,fact=0.0015,xlims=[-5.0,5.0])
    plot_steady_state(sys,symmodel,fact=0.0015,xlims=[-5.0,5.0])
    plot_steady_state(sys,symmodel,fact=0.001,xlims=[-5.0,5.0])
    # plot_steady_state(sys,symmodel,fact=0.0005,xlims=[-5.0,5.0])
    # plot_steady_state(sys,symmodel,fact=0.0001,xlims=[-5.0,5.0])
    # plot_steady_state(sys,symmodel,fact=0.0,xlims=[-5.0,5.0])
    #plot_shannon_entropy(sys,symmodel,xlims=[-8.0,8.0], ylims=[-5.0,10.0])

    # region = get_steady_state_region_threshold(symmodel, 0.0005)#0.0008
    # plot_region(symmodel,sys,region)
    # fixed_point!(sys,symmodel,region)
    # println("LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL")
    # println(region)
    # plot_region(symmodel,sys,region)
end
# test_limit_cycle()

function update!(symmodel,sys,max)
    Xdom = symmodel.Xdom
    mc = symmodel.mc
    indexes = sortperm(mc.steady_state)#sortperm(mc.entropy_bis,rev = true) # sortperm(mc.steady_state) #sortperm(mc.entropy_bis,rev=true)
    count = 0
    for (i,index) in enumerate(indexes)
        s = mc.symmodel_from_mc[index]
        if get_steady_state(mc, s) > 0.0
            count += 1
            split_node!(symmodel, sys, s)
        end
        if count == max
            break
        end
    end

end


function algo(symmodel,sys,tab)
    compute_symbolic_full_domain!(symmodel,sys)
    update_MC!(symmodel)
    plot_limit_cycle(symmodel,sys)
    for steps in tab
        update!(symmodel,sys,steps)
        update_MC!(symmodel)
        plot_limit_cycle(symmodel,sys)
    end
end

function algo!(symmodel,sys,tab)
    compute_symbolic_full_domain!(symmodel,sys)
    update_MC!(symmodel)
    plot_steady_state!(symmodel,fact=0.0)
    colors = [:blue,:green,:purple,:black]
    for (i,steps) in enumerate(tab)
        update!(symmodel,sys,steps)
        update_MC!(symmodel)
        plot_steady_state!(symmodel,fact=0.0,color=colors[i])
    end
end



function update_finale!(symmodel,sys)
    Xdom = symmodel.Xdom
    mc = symmodel.mc
    indexes = sortperm(mc.entropy_SS_finale, rev = true)#sortperm(mc.entropy_bis,rev = true) # sortperm(mc.steady_state) #sortperm(mc.entropy_bis,rev=true)
    s = mc.symmodel_from_mc[indexes[1]]
    split_node!(symmodel, sys, s)
end
function test_finale()
    μ = 2.0
    hx = [1.0, 1.0]*1.0  #0.3
    X,Xdom = build_dom(nothing,nothing,hx,12.0) #10.0
    Udom = UT.CustomList([SVector(0.0,1.0)])
    sys =  build_system(X,μ)

    symmodel = NewMultiSymbolicModel(Xdom, Udom)
    compute_symbolic_full_domain!(symmodel,sys)
    for i in 1:150
        println()
        println(i)
        update_MC!(symmodel)
        println("number of cells : ", get_active_ncells(symmodel))
        println("entropy: ", sum(symmodel.mc.entropy_SS_finale))
        #plot_steady_state(sys,symmodel)
        update_finale!(symmodel,sys)
        #plot_limit_cycle(symmodel,sys)
    end
    update_MC!(symmodel)
    plot_steady_state(sys,symmodel)
end

function test_finale_2()
    μ = 2.0
    hx = [1.0, 1.0]*0.8  #0.3
    X,Xdom = build_dom(nothing,nothing,hx,12.0) #10.0
    Udom = UT.CustomList([SVector(0.0,1.0)])
    sys =  build_system(X,μ)

    symmodel = NewMultiSymbolicModel(Xdom, Udom)
    compute_symbolic_full_domain!(symmodel,sys)
    update_MC!(symmodel)
    println("number of cells : ", get_active_ncells(symmodel))
    println("entropy: ", sum(symmodel.mc.entropy_SS_finale))
    plot_steady_state(sys,symmodel)
end


test_finale()
test_finale_2()

function test1()
    μ = 2.0
    hx = [1.0, 1.0]*1.0  #0.3
    X,Xdom = build_dom(nothing,nothing,hx,12.0) #10.0
    Udom = UT.CustomList([SVector(0.0,1.0)])
    sys =  build_system(X,μ)

    symmodel = NewMultiSymbolicModel(Xdom, Udom)
    tab = [100, 200, 200, 100, 100]#[100, 200, 200, 100, 100]
    algo(symmodel,sys,tab)
end

function test2()
    μ = 2.0
    hx = [1.0, 1.0]*2.0  #0.3
    X,Xdom = build_dom(nothing,nothing,hx,12.0) #10.0
    Udom = UT.CustomList([SVector(0.0,1.0)])
    sys =  build_system(X,μ)

    symmodel = NewMultiSymbolicModel(Xdom, Udom)

    fig = plot(aspect_ratio = 1,legend = false,title="steady-state probability")
    tab = [80,300,300]
    algo!(symmodel,sys,tab)

    x0 = SVector(2.0,2.0)
    plot_trajectory!(sys,x0,200)
    display(fig)
end


function test3()
    μ = 2.0
    hx = [1.0, 1.0]*2.0  #0.3
    X,Xdom = build_dom(nothing,nothing,hx,12.0) #10.0
    Udom = UT.CustomList([SVector(0.0,1.0)])
    sys =  build_system(X,μ)

    symmodel = NewMultiSymbolicModel(Xdom, Udom)

    fig = plot(aspect_ratio = 1,legend = false,title="steady-state probability")
    tab = [20,80,80]
    algo!(symmodel,sys,tab)

    x0 = SVector(2.0,2.0)
    plot_trajectory!(sys,x0,200)
    display(fig)
end

function test4()
    μ = 1.0
    hx = [1.0, 1.0]*1.0  #0.3
    X,Xdom = build_dom(nothing,nothing,hx,12.0) #10.0
    Udom = UT.CustomList([SVector(0.0,1.0)])
    sys =  build_system(X,μ)

    symmodel = NewMultiSymbolicModel(Xdom, Udom)

    fig = plot(aspect_ratio = 1,legend = false,title="steady-state probability")
    tab = [20,120,200]
    algo!(symmodel,sys,tab)

    x0 = SVector(2.0,2.0)
    plot_trajectory!(sys,x0,200)
    display(fig)
end

function test5()
    μ = 1.0
    hx = [1.0, 1.0]*1.0  #0.3
    X,Xdom = build_dom(nothing,nothing,hx,12.0) #10.0
    Udom = UT.CustomList([SVector(0.0,1.0)])
    sys =  build_system(X,μ)

    symmodel = NewMultiSymbolicModel(Xdom, Udom)

    tab = [40,200,300,400] # [20,400,400]
    algo(symmodel,sys,tab)
end

function test()
    f_grid = nothing
    f_grid_i = nothing
    hx = [1.0, 1.0]*2.0  #0.3
    X,Xdom = build_dom(f_grid,f_grid_i,hx,12.0) #10.0
    Udom = UT.CustomList([SVector(0.0,1.0)])
    sys =  build_system(X)

    symmodel = NewMultiSymbolicModel(Xdom, Udom)

    #test1(symmodel, sys)
    test3(symmodel, sys)

end

function test6()
    proba = [0.1,0.1,0.8]
    proba2 = [0.4,0.3,0.3]
    e = 0.0
    for p in proba2
        e = e - p*log2(p)
    end
    println(e)
    println(log(4,2))
    println(log2(4))
end
# test1()
# test2()
# test3()
# test4()
# test1()
end
