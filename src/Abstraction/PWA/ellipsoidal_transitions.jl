using JuMP, GLPK, SCS, Gurobi, Mosek, MosekTools
using LinearAlgebra
using Polyhedra
using HybridSystems
using ..Abstraction


function _has_transition(subsys::HybridSystems.ConstrainedAffineControlDiscreteSystem,P,c,Pp,cp,W,L,U)
    
    eye(n) = diagm(ones(n))
    A = subsys.A
    Bc = subsys.B
    g= subsys.c
    n = length(c);
    m = size(U[1],2);
    N = size(W,2);
    p = size(W,1);
    Nu = length(U);
    B = Bc[:,1:m]
    H = Bc[:,(m+1):(m+p)]

    model = Model(Mosek.Optimizer)
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
        w = H*W[:,i];
        @SDconstraint(model,
            [bta[i]*P        z         t(At)
            t(z)        1-bta[i]    t(A*c+gt-cp+w)
             At       A*c+gt-cp+w      inv(Pp)           ] ⪰ eye(2*n+1)*1e-4)

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

function compute_symmodel_from_hybridcontrolsystem!(symmodel::SymbolicModel{N},
    hybridsys::AbstractHybridSystem, W, L, Uc) where N
    println("compute_symmodel_from_hybridcontrolsystem! started")
    Xdom = symmodel.Xdom
    Udom = symmodel.Udom

    r = Xdom.grid.h/2.0

    baseCell = _compute_base_cell(r)


    P = Diagonal(inv.(r.^2))
    Pm = 1/(r'P*r) * P

    function get_mode(x)
        return sum(map(m-> (x ∈ m.X), hybridsys.modes).*(1:length(hybridsys.modes)))
    end



    
    for xpos in Xdom.elems

        source = get_state_by_xpos(symmodel, xpos)
        x = get_coord_by_pos(Xdom.grid, xpos)
        m = get_mode(x)
        xcell = vrep([x]) + baseCell

        A = hybridsys.modes[m].A
        B = hybridsys.modes[m].B
        c = vrep([Vector(hybridsys.modes[m].c)])
        U = hybridsys.modes[m].U
        
        transitionCost = Dict()
        transitionKappa = Dict()
        xpost = A*xcell + Matrix(B)*Diagonal([1.0; 0; 0])*U + c

        v_xpost = hcat(points(xpost)...)
           
        
        rectI = get_pos_lims_outer(Xdom.grid, Xdom.grid.rect ∩ HyperRectangle(min(eachcol(v_xpost)...), max(eachcol(v_xpost)...)))

        xmpos_iter = Iterators.product(_ranges(rectI)...)
        for xmpos in xmpos_iter
            xm = get_coord_by_pos(Xdom.grid, xmpos) 
            ans, cost, kappa =_has_transition(hybridsys.modes[m],P,x,Pm,xm,W,L,Uc)
            if(ans)
                target = get_state_by_xpos(symmodel, xmpos)
                add_transition!(symmodel.autom, source, target, target)
                transitionCost[(source,c)] = cost
                transitionKappa[(source,c)] = kappa
            end
        end

        # get over approximation of Axcell?
        # get points in this over approximation (create domain)
        # for each target in the domain, solve LMI, add transition if feasible
        
    end    
end
