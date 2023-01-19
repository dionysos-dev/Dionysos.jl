#  Copyright 2022, Lucas N. Egidio, Thiago Alves Lima, and contributors
#############################################################################
# Dionysos
# See https://github.com/dionysos-dev/Dionysos.jl
#############################################################################

using JuMP
using LinearAlgebra
using Polyhedra
using HybridSystems
using SpecialFunctions
using ProgressMeter
using ..Domain
using ..Utils


"""
    _has_transition(subsys::Union{HybridSystems.ConstrainedAffineControlDiscreteSystem,HybridSystems.HybridSystems.ConstrainedAffineControlMap},#
    P,c,Pp,cp,W,L,U, optimizer)

Verifies whether a controller u(x)=K(x-c)+ell exists for `subsys` satisfying input requirements
defined by `U` ans performs a sucessful transitions from a starting set Bs = {(x-c)'P(x-c) ≤ 1}
to the final set Bs = {(x-c)'P(x-c) ≤ 1}. A tight upper bound on the transition cost c(x,u) is 
minimized where c(x,u) = |L*[x; u; 1]|^2, from the parameter `L` and each columm of the matrix  
`W` defines a vertex of the polytope from which additive disturbance are drawn. 

The input restrictions are defined by the list U as |U[i]*u| ≤ 1. `optimizer` must be a JuMP 
SDP optimizer (e.g., Mosek, SDPA, COSMO, ...).


## Note
This implements the optimization problem presented in Corollary 1 of the following paper 
https://arxiv.org/pdf/2204.00315.pdf

"""
function _has_transition(subsys::Union{HybridSystems.ConstrainedAffineControlDiscreteSystem,HybridSystems.HybridSystems.ConstrainedAffineControlMap},#
    P,c,Pp,cp,W,L,U, optimizer)
    
    eye(n) = diagm(ones(n))
    A = subsys.A
    B = subsys.B
    g = subsys.c
    n = length(c);
    m = size(U[1],2);
    N = size(W,2);
    p = size(W,1);
    Nu = length(U);


    model = Model(optimizer)
    @variable(model, K[i=1:m,j=1:n]) 
    @variable(model, ell[i=1:m,j=1:1])
    @variable(model, bta[i=1:N] >= 0)
    @variable(model, tau[i=1:Nu] >= 0)
    @variable(model, gamma >= 0)
    @variable(model, J >= 0)

    @expressions(model, begin
        At, A+B*K
        gt, g+B*ell
    end)


    t(x) = transpose(x);

    z = zeros(n,1);
    
    for i in 1:N
        w = W[:,i];
        aux = A*hcat(c)+hcat(gt)-hcat(cp)+hcat(w)
        @constraint(model,
            [bta[i]*P        z         t(At)
            t(z)        1-bta[i]        t(aux)
             At       aux      inv(Pp)           ] >= eye(2*n+1)*1e-4, PSDCone())

    end
    
    for i=1:Nu
        n_ui = size(U[i],1);
        @constraint(model,
        [tau[i]*P        z          t(U[i]*K)
        t(z)        1-tau[i]        t(U[i]*ell)
        U[i]*K      U[i]*ell        eye(n_ui)   ] >= eye(n+n_ui+1)*1e-4, PSDCone())
        
    end
    n_L = size(L,1);
    @constraint(model,
    [gamma*P                     z           [I t(K) z]*t(L)
     t(z)                   J-gamma          [t(c) t(ell) 1]*t(L)
     L*t([I t(K) z])    L*t([t(c) t(ell) 1])        eye(n_L)       ] >= eye(n+n_L+1)*1e-4, PSDCone())
    
    @objective(model, Min, J)

    #print(model)
    optimize!(model)

    kappa = [value.(K) value.(ell)];
    cost = value(J);
    ans = solution_summary(model).termination_status == MOI.OPTIMAL
    #println("$(solution_summary(model).solve_time) s")
    return ans, cost, kappa
end

"""
    _compute_base_cell(r::SVector{S})

Computes a polyhedron containing the base hyperrectangular cell, centered at the origin
and with the i-th side lenght given by `2*r[i]`. 

"""
function _compute_base_cell(r::SVector{S}) where S
    baseCellList = []
    for i in 1:S
        vec = SVector{S}(1:S .==i)
        append!(baseCellList, [ HalfSpace(-vec, r[i]) ∩ HalfSpace(vec, r[i])  ])
    end 
    return polyhedron(intersect(baseCellList...))
end



"""
    compute_symmodel_from_hybridcontrolsystem!(symmodel::SymbolicModel{N}, transitionCost::AbstractDict, transitionKappa::AbstractDict,
    hybridsys::AbstractHybridSystem, W, L, U, opt_sdp, opt_qp)

Builds an abstraction `symmodel` where the transitions have costs given in `transitionCost`
and are parameterized by affine-feedback controllers in `transitionKappa`. The concrete system 
is `hybridsys` and `W`, `L` and `U` are defined as in `_has_transition`. An SDP optimizer `opt_sdp`
and a QP optimizer `opt_qp` must be provided as JuMP optimizers.

"""
function compute_symmodel_from_hybridcontrolsystem!(symmodel::SymbolicModel{N}, transitionCost::AbstractDict, transitionKappa::AbstractDict,
    hybridsys::AbstractHybridSystem, W, L, U, opt_sdp, opt_qp) where N
    println("compute_symmodel_from_hybridcontrolsystem! started")
    Xdom = symmodel.Xdom
    

    r = Xdom.grid.h/2.0

    n_sys = length(r)
    if Xdom.grid isa Domain.GridEllipsoidalRectangular
        Pm = Xdom.grid.P
        P = Pm
        R = _get_min_bounding_box(P, opt_qp)

    else
        Pm = (1/n_sys) * diagm(inv.(r.^2))
        P =  Pm
        R = r
    end
 
    # get affine mode number for a point x
    get_mode(x) = findfirst(m -> (x ∈ m.X), hybridsys.resetmaps)

    vec_list = collect(Iterators.product(eachcol(repeat(hcat([-1,1]),1,n_sys))...))[:] # list of vertices of a hypersquare centered at the origin and length 2

    bds2rectverts(lb,ub) = hcat([v.*((ub-lb)/2)+(ub+lb)/2  for v in vec_list ]...) #generate a matrix containing vertices of a hyperrectangle with lower vertice lb and upper one ub
    function _compute_xpost(A,x,B,U,c,R)
        Axcell = A*bds2rectverts(x-R,x+R)
        
        Bu = B*hcat(points(U)...)
        

        return [min(eachcol(Axcell)...) + min(eachcol(Bu)...) + c-R, #
                max(eachcol(Axcell)...) + max(eachcol(Bu)...) + c+R]
    end

    trans_count = 0
    @showprogress 1 "Computing symbolic control system: " (
    for xpos in Xdom.elems
        
        source = get_state_by_xpos(symmodel, xpos)
        x = Domain.get_coord_by_pos(Xdom.grid, xpos)
        m = get_mode(x)

        A = hybridsys.resetmaps[m].A
        B = hybridsys.resetmaps[m].B
        c = hybridsys.resetmaps[m].c
        Upoly = hybridsys.resetmaps[m].U
        
        xpost = _compute_xpost(A,x,B,Upoly,c,R)
                
        rectI = Domain.get_pos_lims_outer(Xdom.grid, Xdom.grid.rect ∩ Utils.HyperRectangle(xpost[1],xpost[2]))

        xmpos_iter = Iterators.product(Domain._ranges(rectI)...)
        for xmpos in xmpos_iter
            xm = Domain.get_coord_by_pos(Xdom.grid, xmpos) 
            ans, cost, kappa =_has_transition(hybridsys.resetmaps[m],P,x,Pm,xm,W,L,U,opt_sdp)
        
            if(ans)
                trans_count += 1
                target = get_state_by_xpos(symmodel, xmpos)
                symbol = get_symbol_by_upos(symmodel, xmpos);
                #println("->Added $(trans_count+1)\nfrom\t $(source)\n","to\t $(target)\n\n")
                add_transition!(symmodel.autom, source, target, symbol)
                transitionCost[(source,target)] = cost
                transitionKappa[(source,target)] = kappa
            end
        end
        
        
    end
    )
    println("compute_symmodel_from_controlsystem! terminated with success: ",
    "$(trans_count) transitions created")    
end


"""
    _provide_P(subsys::HybridSystems.ConstrainedAffineControlDiscreteSystem, optimizer)

If `subsys` is a stabilizable system, finds the matrix `P` and the state-feedback gain `K`
that satisfy the discrete-time Lyapunov inequality (A+BK)'P(A+BK)-P < 0. The condition number
of `P` is minimized. `optimizer` must be a JuMP SDP optimizer.
"""
function _provide_P(subsys::HybridSystems.ConstrainedAffineControlDiscreteSystem, optimizer)
    
    eye(n) = diagm(ones(n))
    A = subsys.A
    B = subsys.B
    n = size(A,1);
    m = size(B,2);


    model = Model(optimizer)
    @variable(model, L[i=1:m,j=1:n]) 
    @variable(model, S[i=1:n,j=1:n], PSD) 
    @variable(model, gamma >= 0)



    t(x) = transpose(x);

    
    @constraint(model, [S      t(A*S+B*L);
                        A*S+B*L    S]        >= 0, PSDCone())
    @constraint(model, eye(n) >= S, PSDCone())
    @constraint(model, S >= -gamma*eye(n), PSDCone())

    
    @objective(model, Min, gamma)

    #print(model)
    optimize!(model)

    P = inv(value.(S));
    K = value.(L)*P;
    gamma = value(gamma);
    ans = solution_summary(model).termination_status == MOI.OPTIMAL
    return ans, K, P, gamma
end


"""
    ellipsoid_vol(P,r)
    
Calculates the n-volume of the n-ellipsoid defined as {x'Px < r}.
"""
function ellipsoid_vol(P,r) 
    N = size(P,1)
    return pi^(N/2)/(gamma(N/2+1))*det(P/r)^(-1/2)
end

"""
    _get_min_bounding_box(P, optimizer) 

Finds the minimum bounding box containing the ellipsoid {x'Px < 1}. 
"""
function _get_min_bounding_box(P, optimizer) 
    n = size(P,1)
    R = zeros(n)
    
    model = Model(optimizer)
    @variable(model, x[i=1:n])

    @constraint(model, x'P*x  <= 1) 
    
    for i in 1:n
        new_model, reference_map = copy_model(model)
        set_optimizer(new_model,optimizer)
        @objective(new_model, Max, reference_map[x[i]])
        optimize!(new_model)
        R[i] = abs(value(reference_map[x[i]]))
    end
    return R
end
