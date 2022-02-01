using JuMP, GLPK, SCS, Gurobi, Mosek, MosekTools
using LinearAlgebra
using Polyhedra
using HybridSystems
using SpecialFunctions
using ..Abstraction


function _has_transition(subsys::HybridSystems.ConstrainedAffineControlDiscreteSystem,P,c,Pp,cp,W,L,U)
    
    eye(n) = diagm(ones(n))
    A = subsys.A
    B = subsys.B
    g = subsys.c
    n = length(c);
    m = size(U[1],2);
    N = size(W,2);
    p = size(W,1);
    Nu = length(U);
    optimizer = optimizer_with_attributes(Mosek.Optimizer, MOI.Silent() => true)


    model = Model(optimizer)
    @variable(model, K[i=1:m,j=1:n]) 
    @variable(model, ell[i=1:m,j=1])
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
        @SDconstraint(model,
            [bta[i]*P        z         t(At)
            t(z)        1-bta[i]        t(aux)
             At       aux      inv(Pp)           ] ⪰ eye(2*n+1)*1e-4)

    end
    
    for i=1:Nu
        n_ui = size(U[i],1);
        @SDconstraint(model,
        [tau[i]*P        z          t(U[i]*K)
        t(z)        1-tau[i]        t(U[i]*ell)
        U[i]*K      U[i]*ell        eye(n_ui)   ]⪰ eye(n+n_ui+1)*1e-4)
        
    end
    n_L = size(L,1);
    @SDconstraint(model,
    [gamma*P                     z           [I t(K) z]*t(L)
     t(z)                   J-gamma          [t(c) t(ell) 1]*t(L)
     L*t([I t(K) z])    L*t([t(c) t(ell) 1])        eye(n_L)       ]>=eye(n+n_L+1)*1e-4)
    
    @objective(model, Min, J)

    #print(model)
    optimize!(model)

    kappa = [value.(K) value.(ell)];
    cost = value(J);
    ans = solution_summary(model).termination_status == MOI.OPTIMAL
    #println("$(solution_summary(model).solve_time) s")
    return ans, cost, kappa
end

# Computs a basic cell from 
function _compute_base_cell(r::SVector{S}) where S
    baseCellList = []
    for i in 1:S
        vec = SVector{S}(1:S .==i)
        append!(baseCellList, [ HalfSpace(-vec, r[i]) ∩ HalfSpace(vec, r[i])  ])
    end 
    return polyhedron(intersect(baseCellList...))
end

function compute_symmodel_from_hybridcontrolsystem!(symmodel::SymbolicModel{N}, transitionCost::AbstractDict, transitionKappa::AbstractDict,
    hybridsys::AbstractHybridSystem, W, L, Uc) where N
    println("compute_symmodel_from_hybridcontrolsystem! started")
    Xdom = symmodel.Xdom
    Udom = symmodel.Udom

    r = Xdom.grid.h/2.0
    n_sys = length(r)


    Pm = (1/n_sys) * diagm(inv.(r.^2))
    P =  Pm

    function get_mode(x)
        o = map(m-> (x ∈ m.X), hybridsys.modes).*(1:length(hybridsys.modes))
        return o[o.>0][1]
    end

    vec_list = collect(Iterators.product(eachcol(repeat(hcat([-1,1]),1,n_sys))...))[:]

    bds2rectverts(lb,ub) = hcat([v.*((ub-lb)/2)+(ub+lb)/2  for v in vec_list ]...)
    function _compute_xpost(A,x,B,U,c,r)
        Axcell = A*bds2rectverts(x-r,x+r)
        
        Bu = B*hcat(points(U)...)
        

        return [min(eachcol(Axcell)...) + min(eachcol(Bu)...) + c, #
                max(eachcol(Axcell)...) + max(eachcol(Bu)...) + c]
    end

    trans_count = 0
    @showprogress 1 "Computing symbolic control system: " (
    for xpos in Xdom.elems
        
        source = get_state_by_xpos(symmodel, xpos)
        x = get_coord_by_pos(Xdom.grid, xpos)
        m = get_mode(x)

        A = hybridsys.modes[m].A
        B = hybridsys.modes[m].B
        c = hybridsys.modes[m].c
        Upoly = hybridsys.modes[m].U
        
        xpost = _compute_xpost(A,x,B,Upoly,c,r)
                
        rectI = get_pos_lims_outer(Xdom.grid, Xdom.grid.rect ∩ HyperRectangle(xpost[1],xpost[2]))

        xmpos_iter = Iterators.product(_ranges(rectI)...)
        for xmpos in xmpos_iter
            xm = get_coord_by_pos(Xdom.grid, xmpos) 
            ans, cost, kappa =_has_transition(hybridsys.modes[m],P,x,Pm,xm,W,L,Uc)
        
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



function _provide_P(subsys::HybridSystems.ConstrainedAffineControlDiscreteSystem)
    
    eye(n) = diagm(ones(n))
    A = subsys.A
    B = subsys.B
    n = size(A,1);
    m = size(B,2);
    optimizer = optimizer_with_attributes(Mosek.Optimizer, MOI.Silent() => true)


    model = Model(optimizer)
    @variable(model, L[i=1:m,j=1:n]) 
    @variable(model, S[i=1:n,j=1:n], PSD) 
    @variable(model, gamma >= 0)



    t(x) = transpose(x);

    
    @SDconstraint(model, [S      t(A*S+B*L);
                        A*S+B*L    S]        ⪰ 0)
    @SDconstraint(model, eye(n) ⪰ S)
    @SDconstraint(model, S ⪰ -gamma*eye(n))

    
    @objective(model, Min, gamma)

    #print(model)
    optimize!(model)

    P = inv(value.(S));
    K = value.(L)*P;
    gamma = value(gamma);
    ans = solution_summary(model).termination_status == MOI.OPTIMAL
    return ans, K, P, gamma
end


function ellipsoid_vol(P,r) 
    N = size(P,1)
    return pi^(N/2)/(gamma(N/2+1))*det(P/r)^(-1/2)
end