include("../partition.jl")

using StaticArrays,LazySets
using Main.PartitionAbstraction
PA = PartitionAbstraction

a = 3
b = 1
function build_zonotopes()
    centers = [[-2.5*a,-12*b],[-5*a,8*b],[3*a,15*b],[5*a,4*b]]./3
    G1 = PA.build_generators(1,centers,[2,3,4])
    Z1 = Zonotope(centers[1],G1)
    G2 = PA.build_generators(2,centers,[1,3,4])
    Z2 = Zonotope(centers[2],G2)
    G3 = PA.build_generators(3,centers,[1,2,4])
    Z3 = Zonotope(centers[3],G3)
    G4 = PA.build_generators(4,centers,[1,2,3])
    Z4 = Zonotope(centers[4],G4)
    Z = [Z1,Z2,Z3,Z4]
    return Z
end

function build_constrainedZonotope(vertices)
    P = VPolytope(vertices)
    return PA.convert_HorV_CG(P)
end

function build_constrainedZonotopes(vertices)
    v = [1,2,5]
    CZ1 = build_constrainedZonotope(vertices[v])
    v = [3,4,8]
    CZ2 = build_constrainedZonotope(vertices[v])
    v = [1,5,6,9]
    CZ3 = build_constrainedZonotope(vertices[v])
    v = [1,6,7]
    CZ4 = build_constrainedZonotope(vertices[v])
    v = [4,8,11,13]
    CZ5 = build_constrainedZonotope(vertices[v])
    v = [10,12,21]
    CZ6 = build_constrainedZonotope(vertices[v])
    v = [12,14,15]
    CZ7 = build_constrainedZonotope(vertices[v])
    v = [12,15,21]
    CZ8 = build_constrainedZonotope(vertices[v])
    v = [15,20,21]
    CZ9 = build_constrainedZonotope(vertices[v])
    v = [16,18,19]
    CZ10 = build_constrainedZonotope(vertices[v])
    v = [17,18,22]
    CZ11 = build_constrainedZonotope(vertices[v])
    v = [18,19,20,22]
    CZ12 = build_constrainedZonotope(vertices[v])
    v = [20,21,22]
    CZ13 = build_constrainedZonotope(vertices[v])

    CZ = [CZ1,CZ2,CZ3,CZ4,CZ5,CZ6,CZ7,CZ8,CZ9,CZ10,CZ11,CZ12,CZ13]

    return CZ
end

function build_obstacles()
    c = 5
    d = -7.5
    O1 = VPolytope([[-1.0*a+c,5*b+d],[0*a+c,-1*b+d],[2*a+c,0*b+d],[1.5*a+c,4.5*b+d]])
    O2 = VPolytope([[-2.0*a,-2*b],[-1*a,-2*b],[-2*a,-6*b],[-1*a,-6*b]])
    O3 = VPolytope([[-5.2*a,-3*b],[-5.2*a,5.8*b],[-4*a,-3*b],[-4*a,5.8*b]])
    O4 = VPolytope([[-5.2*a,-3*b],[-2*a,-3*b],[-5.2*a,-0.5*b],[-2*a,-0.5*b]])
    O5 = VPolytope([[1.0*a,0.0*b],[-3.0*a,0.0*b],[1.0*a,0.5*b],[-3.0*a,0.5*b]])
    return [O1,O2,O3,O5]
end

function build_system()
    tstep = 0.5
    nsys = 5
    function F_sys(x, u)
        return SVector{2}(u[1],u[2])
    end
    measnoise = SVector(0.0, 0.0)
    Lip = 0.0
    return PA.NewControlSystemLipschitzRK4(tstep, F_sys, Lip,measnoise, nsys)
end

function build_domain_partition()
    c = [0.0, 0.0]; r = [5.0*a, 10.0*b]
    X = Hyperrectangle(c, r)
    Z = build_zonotopes()
    v = PA.get_vertices(X,Z)
    CZ = build_constrainedZonotopes(v)
    L = [Z...,CZ...]
    D = PA.Partition{AbstractPolytope}(X,L)
    return D_e = PA.expansion(1.2,D) #1.1
end

function print_ntransitions(control_list)
    ntransisitions = [length(control_data[1].autom.transitions) for control_data in control_list]
    println("Number of transisitions       : ",ntransisitions)
    println("Total number of transisitions : ",sum(ntransisitions))
end

function test1()
    # build the domain partition
    D = build_domain_partition()
    PA.plot_partition(D)

    # build the obstacles
    O = build_obstacles()

    #build the toplogical graph
    A = PA.topological_graph(D,O)
    PA.plot_topological_graph(D,A,O)

    # definition of the system
    contsys = build_system()

    # abstraction parameters
    _U_ = Hyperrectangle(low = [-2.0,-2.0],high = [2.0,2.0])
    hu = SVector(0.5,0.5)
    hx_l = [[0.4,0.4],[0.8,0.8],[0.5,0.5],[0.4,0.4]]

    # control problem
    X0 = Hyperrectangle([-12.7, -7.2],[0.6,0.6])
    Xf = Hyperrectangle([13.5, 8.0],[0.8,0.8])

    control_list = PA.control(X0,Xf,D,A,contsys,O,_U_,hu,hx_l)

    x0 = SVector(-8.4-4, -7.5)
    trajectory = PA.trajectory_reach_full(contsys,control_list,x0; randchoose = false)

    print_ntransitions(control_list)
    PA.print_trajectory_full(control_list,trajectory;full=false,partition=D,obstacles=O)
    PA.print_trajectory_full(control_list,trajectory;full=true,obstacles=O,X0=X0,Xf=Xf)

end


function test2()
    # build the domain
    c = [0.0, 0.0]; r = [5.0*a, 10.0*b]
    X = Hyperrectangle(c, r)
    D = PA.Partition{AbstractPolytope}(X,[X])
    # build the obstacles
    O = build_obstacles()

    #build the toplogical graph
    A = PA.topological_graph(D,O)
    PA.plot_topological_graph(D,A,O)

    # definition of the system
    contsys = build_system()

    # abstraction parameters
    _U_ = Hyperrectangle(low = [-2.0,-2.0],high = [2.0,2.0])
    hu = SVector(0.5,0.5)
    hx_l = [[0.5,0.5]]


    # control problem
    X0 = Hyperrectangle([-12.7, -7.2],[0.6,0.6])#Singleton([-8.4, -7.5])
    Xf = Hyperrectangle([13.5, 8.0],[0.8,0.8])#Hyperrectangle([-4.0,-7.0],[1.0,0.8])#Hyperrectangle([-7.5,-5.0],[1.0,0.8])#Hyperrectangle([-8.0,-6.4],[1.4,0.9])#PA.expansion(0.9,D.L[7])#Singleton([-7.0, -4.0]) # SVector(8.0, 6.0) #SVector(-7.0, -4.0)

    control_list = PA.control(X0,Xf,D,A,contsys,O,_U_,hu,hx_l)

    x0 = SVector(-12.4, -7.5)
    trajectory = PA.trajectory_reach_full(contsys,control_list,x0; randchoose = false)

    print_ntransitions(control_list)
    PA.print_trajectory_full(control_list,trajectory;full=false,partition=D,obstacles=O)
    PA.print_trajectory_full(control_list,trajectory;full=true,obstacles=O,X0=X0,Xf=Xf)

end

println()
println()
println("START TEST")

test1()
#test2()
