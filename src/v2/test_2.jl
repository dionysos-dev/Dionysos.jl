include(joinpath("..", "Abstraction", "abstraction.jl"))
#include(joinpath("..", "CDC2021","CDC2021.jl"))
include(joinpath("..", "CDC2021","general_domain.jl"))
include(joinpath("..", "CDC2021","utils.jl"))
include(joinpath("..", "CDC2021","lazy_abstraction.jl"))
include(joinpath("..", "CDC2021","alternating_simulation.jl"))


# test lazy abstraction with the simple 2D example (should be deleted as particular case of periodic (general) domain)
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


using HybridSystems
include("sorted_vector_set.jl")
include("automaton.jl")
include("domain.jl")
include("markov_chain.jl")
include("symbolic.jl")
# include("algo.jl")
include("multigrid.jl")
include("system.jl")

function build_system(X,μ)
    tstep = 0.08
    measnoise = SVector(0.0,0.0)
    sysnoise = SVector(0.0,0.0)
    nsys = 5
    ngrowthbound = 5

    #μ = 2.0 #2.0
    k = 1.0
    f = system_with_boundary(Van_der_pol_oscillator(μ;k=k),10.0,-1.0)
    f_reverse = reverse_Van_der_pol_oscillator(μ;k=k)
    # f2 = Van_der_pol_oscillator(2.0)
    jacobian = Jacobian_with_boundary(Van_der_pol_oscillator_Jacobian(2.0),12.0,-0.8)
    contsys = NewControlSystemGrowthLx(tstep, f, jacobian, f_reverse, measnoise, sysnoise, nsys, ngrowthbound, X)
    return contsys
end

function build_dom(f,fi,hx,r)
    X = AB.HyperRectangle(SVector(-r, -r), SVector(r, r))
    obstacle = AB.HyperRectangle(SVector(-100.0, -100.0), SVector(-99.0, -99.0))
    d = D.RectanglularObstacles(X, [obstacle])
    Xdom = D.GeneralDomainList(hx,d;fit=true,f=f,fi=fi)
    return X,Xdom
end

function test()
    f = nothing
    fi = nothing
    hx = [1.0, 1.0]*0.5
    X,Xdom = build_dom(f,fi,hx,10.0) #10.0
    Udom = UT.CustomList([SVector(0.0,1.0)])
    F_sys = system_with_boundary(Van_der_pol_oscillator(2.0),8.0,-0.8)
    sys = NewControlSystemGrowthRK4(0.08,F_sys,SVector(0.0, 0.0),4)

    #sys = build_system_Van_der_pol_oscillator()
    symmodel = NewMultiSymbolicModel(Xdom, Udom)
    prob = NewAnalysisProblem(symmodel,sys,nothing,false)
    time = @elapsed begin
    solve!(prob)
    end
    println("total time: lazy abstraction + controller: ", time)
    plot_steady_state(sys,symmodel)
    # plot_steady_state(sys,symmodel,fact=0.004)
    # plot_steady_state(sys,symmodel,fact=0.003)
    # plot_steady_state(sys,symmodel,fact=0.002)
    # plot_steady_state(sys,symmodel,fact=0.0015)
    # plot_steady_state(sys,symmodel,fact=0.001)
    # plot_steady_state(sys,symmodel,fact=0.0005)
    # plot_steady_state(sys,symmodel,fact=0.00001)
    # plot_steady_state(sys,symmodel,fact=0.000005)
    plot_steady_state(sys,symmodel,fact=0.0)
    h = get_entropy_chain(symmodel.mc)
    V = AB.get_volume(Xdom.grid)
    ϵ = 0.001
    l = log2(V/(4*ϵ*ϵ))
    println("h: ",h)
    println(V)
    println(log2(V/(4*ϵ*ϵ)))
    println("hV: ",h+l)
    #detect(symmodel.mc.steady_state)
    #plot_automaton(prob)
end

# test()


###############################################################################

function test_size()
    f = nothing
    fi = nothing
    #f,fi = build_f_rotation((15/180)*π)
    sizes = reverse([h for h in 0.2:0.05:5.0])

    F_sys = system_with_boundary(Van_der_pol_oscillator(2.0),8.0,-0.8)
    sys = NewControlSystemGrowthRK4(0.08,F_sys,SVector(0.0, 0.0),4)
    #sys = build_system_Van_der_pol_oscillator()
    entropy_SS = []
    entropy_V = []
    hx_vec = []
    for (i,c) in enumerate(sizes)
        hx = [1.0, 1.0]*c
        push!(hx_vec,hx[1])
        X,Xdom = build_dom(f,fi,hx,10.0)
        Udom = UT.CustomList([SVector(0.0,1.0)])

        symmodel = NewMultiSymbolicModel(Xdom, Udom)
        compute_symbolic_full_domain!(symmodel,sys)
        update_MC!(symmodel)

        # plot_steady_state(sys,symmodel)
        # plot_steady_state(sys,symmodel,bool=true)
        # plot_shannon_entropy(sys,symmodel)
        h = get_entropy_chain(symmodel.mc)
        V = AB.get_volume(Xdom.grid)
        ϵ = 0.001
        l = log2(V/(4*ϵ*ϵ))
        hV = h+l
        push!(entropy_SS,h)
        push!(entropy_V,hV)
        println(hx, " , entropy steady-state , ",h)
        println(hx, " , with volume entropy , ",hV)
    end
    fig = plot(legend = false)
    plot!(hx_vec,entropy_SS)
    display(fig)
    fig = plot(legend = false)
    plot!(hx_vec,entropy_V)
    display(fig)
    #plot_automaton(prob)
end

# test_size()

function detect(ps)
    threshold = [t for t in 0.0:0.001:0.1]
    countSet = []
    for t in threshold
        count = 0
        for p in ps
            if p>t
                count+=1
            end
        end
        push!(countSet,count)
    end
    fig = plot(legend = false)
    plot!(threshold,countSet)
    display(fig)
end

function test_volume_entropy_ϵ()
    V = 1.0
    ϵSet = [ϵ for ϵ in 0.00001:0.0001:0.1]
    entropy = []
    for ϵ in ϵSet
        h = log2(V/(4*ϵ*ϵ))
        push!(entropy,h)
    end
    fig = plot(legend = false)
    plot!(ϵSet,entropy)
    display(fig)
end
function test_volume_entropy_V()
    Vset = [v for v in 0.1:0.1:10.0]
    ϵ = 0.0001
    entropy = []
    for V in Vset
        h = log2(V/(4*ϵ*ϵ))
        push!(entropy,h)
    end
    fig = plot(legend = false)
    plot!(Vset,entropy)
    display(fig)
end
# test_volume_entropy_ϵ()
# test_volume_entropy_V()

function test_compare()
    X = AB.HyperRectangle(SVector(-10.0, -10.0), SVector(10.0, 10.0))
    obstacle = AB.HyperRectangle(SVector(-100.0, -100.0), SVector(-99.0, -99.0))
    hx = [1.0, 1.0]*0.5
    d = D.RectanglularObstacles(X, [obstacle])

    Udom = UT.CustomList([SVector(0.0,1.0)])
    sys = build_system_Van_der_pol_oscillator()
    angles = [θ for θ in 0:((15/180)*π):π/2.0]
    entropy = Float64[]
    for θ in angles
        f,fi = build_f_rotation(θ)
        Xdom = D.GeneralDomainList(hx,d;fit=false,f=f,fi=fi)
        symmodel = NewMultiSymbolicModel(Xdom, Udom)
        compute_symbolic_full_domain!(symmodel,sys)
        update_MC!(symmodel)
        h = get_entropy_chain(symmodel.mc)
        println()
        println(θ)
        println(h)
        push!(entropy,h)
        # plot_steady_state(sys, symmodel)
        # plot_steady_state(sys, symmodel, fact=0.0)
        # plot_shannon_entropy(sys, symmodel)
    end
    fig = plot(aspect_ratio = 1,legend = false)
    println(angles)
    println(entropy)
    plot!(angles,entropy)
    display(fig)
    #plot_automaton(prob)
end

# test_compare()


function compare_angles()
    μ = 2.0
    hx = [1.0, 1.0]*0.5  #0.3

    angles = [θ for θ in 0:((15/180)*π):π/2.0]
    entropy_SS = Float64[]
    for θ in angles
        f,fi = build_f_rotation(θ)
        X,Xdom = build_dom(f,fi,hx,12.0) #10.0
        Udom = UT.CustomList([SVector(0.0,1.0)])
        sys =  build_system(X,μ)

        # Xdom = D.GeneralDomainList(hx,d;fit=false,f=f,fi=fi)
        symmodel = NewMultiSymbolicModel(Xdom, Udom)
        compute_symbolic_full_domain!(symmodel,sys)
        update_MC!(symmodel)

        # plot_steady_state(sys,symmodel)
        # plot_steady_state(sys,symmodel,fact=0.0015)
        # plot_steady_state(sys,symmodel,fact=0.001)
        # plot_steady_state(sys,symmodel,fact=0.0005)
        # plot_steady_state(sys,symmodel,fact=0.0001)
        # plot_steady_state(sys,symmodel,fact=0.0)
        push!(entropy_SS,get_entropy_SS(symmodel.mc))
        println("entropy SS: ",θ," is ", get_entropy_SS(symmodel.mc))
    end
    fig = plot(aspect_ratio = 1,legend = false)
    println(angles)
    println(entropy_SS)
    plot!(angles,entropy_SS)
    display(fig)
end
# compare_angles()

function test_rotated()
    μ = 2.0
    hx = [1.0, 1.0]*0.3  #0.3
    θ = 0.7853981633974483#0.2617993877991494#1.3#1.8
    f,fi = build_f_rotation(θ)
    X,Xdom = build_dom(f,fi,hx,12.0) #10.0
    Udom = UT.CustomList([SVector(0.0,1.0)])
    sys =  build_system(X,μ)

    # Xdom = D.GeneralDomainList(hx,d;fit=false,f=f,fi=fi)
    symmodel = NewMultiSymbolicModel(Xdom, Udom)
    compute_symbolic_full_domain!(symmodel,sys)
    update_MC!(symmodel)
    a = 5.0
    b = 5.0
    h = get_entropy_chain(symmodel.mc)
    println(h)
    #plot_shannon_entropy(sys,symmodel)
    plot_steady_state(sys,symmodel,xlims=[-b,b],ylims=[-a,a])
    plot_steady_state(sys,symmodel)
    plot_steady_state(sys,symmodel,fact=0.0015)
    plot_steady_state(sys,symmodel,fact=0.001)
    plot_steady_state(sys,symmodel,fact=0.0005)
    plot_steady_state(sys,symmodel,fact=0.0001)
    plot_steady_state(sys,symmodel,fact=0.0)
    println("entropy SS: ", get_entropy_SS(symmodel.mc))

    # plot_steady_state(sys,symmodel,fact=0.0015,xlims=[-5.0,5.0])
    #plot_steady_state(sys,symmodel,fact=0.001,xlims=[-5.0,5.0])
    # plot_steady_state(sys,symmodel,fact=0.0005,xlims=[-5.0,5.0])
    # plot_steady_state(sys,symmodel,fact=0.0001,xlims=[-5.0,5.0])
    # plot_steady_state(sys,symmodel,fact=0.0,xlims=[-5.0,5.0])
end
test_rotated()
 # 0.2617993877991494
 # 0.7853981633974483

###############################################################################
function plot_cell_sample!(symmodel, sys, coord)
    grid = symmodel.Xdom.domains[1].grid
    N = 1000
    xpos = AB.get_pos_by_coord(grid, coord)
    points = sample_elem(grid, xpos, N)
    targetlist = Int[]
    u = SVector(0.0,0.0)
    symbol = 1
    Fpoints = []
    source = get_state_by_coord(symmodel,coord)
    for x in points
        Fx = sys.sys_map(x, u, sys.tstep)
        push!(Fpoints,Fx)
        target = get_state_by_coord(symmodel,Fx)
        push!(targetlist,target)
    end
    occurences = count_occurences(targetlist)
    probaList = [(target,occ/N) for (target,occ) in occurences]
    transitions = []
    for (target,proba) in probaList
        push!(transitions, (target, source, symbol, proba))
    end
    plot!(symmodel,annotate=true)
    for x in points
        scatter!([x[1]],[x[2]],color =:green,markersize=1)
    end
    for x in Fpoints
        scatter!([x[1]],[x[2]],color =:red,markersize=1)
    end
    println()
    println("Transitions:  ", transitions)
    add_transitions!(symmodel.autom, transitions)
    h = get_entropy(transitions)
    println("Shannon entropy:  ", h)
end

function test_plot(coord)
    fig = plot(aspect_ratio = 1,legend = false)
    f = nothing
    fi = nothing
    hx = [1.0, 1.0]*0.5
    X,Xdom = build_dom(f,fi,hx,14.0)
    Udom = UT.CustomList([SVector(0.2,1.0)])
    F_sys = system_with_boundary(Van_der_pol_oscillator(2.0),12.0,-0.8)
    sys = NewControlSystemGrowthRK4(0.08,F_sys,SVector(0.0, 0.0),4)
    symmodel = NewMultiSymbolicModel(Xdom, Udom)
    plot_cell_sample!(symmodel, sys, coord)
    display(fig)
end
# test_plot(SVector(0.0,-5.0))
# test_plot(SVector(5.0,0.01))
# test_plot(SVector(0.0,-0.1))
# test_plot(SVector(8.0,8.0))

#############################################################################
function plot_circular_grid()
    h = SVector(1.0,π/6.0)
    grid = build_circular_grid(SVector(0.0,0.0), h)
    fig = plot(aspect_ratio = 1,legend = false)
    for i=0:10
        for j=0:12
            pos = (i,j)
            AB.plot_elem!(grid, pos, opacity=1.0,color=:yellow,N=12)
        end
    end
    display(fig)
end

# plot_circular_grid()

function build_dom2()
    X = AB.HyperRectangle(SVector(-10.0, -10.0)*2.0, SVector(10.0, 10.0)*2.0)
    obstacle = AB.HyperRectangle(SVector(-100.0, -100.0), SVector(-99.0, -99.0))
    h = SVector(1.0,π/6.0)
    origin = SVector(0.0,0.0)
    grid = build_circular_grid(origin, h)
    d = D.RectanglularObstacles(X, [obstacle])
    Xdom = D.GeneralDomainList(grid,d)
    return X,Xdom
end

function test_2()
    X,Xdom = build_dom2()
    Udom = UT.CustomList([SVector(0.0,1.0)])
    F_sys = system_with_boundary(Van_der_pol_oscillator(2.0),8.0,-0.8)
    #sys = NewControlSystemGrowthRK4(0.08,F_sys,SVector(0.0, 0.0),4)
    μ = 2.0
    sys =  build_system(X,μ)

    symmodel = NewMultiSymbolicModel(Xdom, Udom)
    compute_symbolic_full_domain!(symmodel,sys)
    update_MC!(symmodel)
    a = 5.0
    b = 5.0
    plot_steady_state(sys,symmodel,xlims=[-b,b],ylims=[-a,a])
    plot_steady_state(sys,symmodel,fact=0.001,xlims=[-b,b],ylims=[-a,a])
    plot_steady_state(sys,symmodel,fact=0.0001,xlims=[-b,b],ylims=[-a,a])
    plot_steady_state(sys,symmodel,fact=0.0,xlims=[-b,b],ylims=[-a,a])
end
# test_2()

end # end module
