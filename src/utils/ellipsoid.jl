using SpecialFunctions
using Plots
using LazySets
using LinearAlgebra
using Polyhedra
using IntervalArithmetic
using JuMP, Ipopt

struct testJ
    v1::Int
    v2::Int
    testJ(v1,v2) = new(v1+v2,v2)
end

struct Ellipsoid{T<:Real,MT<:AbstractMatrix{T},VT<:AbstractVector{T}}
    P::MT
    c::VT
    function Ellipsoid(P::MT,c::VT) where {T<:Real,MT<:AbstractMatrix{T},VT<:AbstractVector{T}}
        M = (P+P')./2
        if !isposdef(M) # see if delete this test
            error("matrix must be positive definite")
        end
        return new{T,MT,VT}(M,c)
    end
end


function Base.in(elli1::Ellipsoid, elli2::Ellipsoid; bisectionB=true)
    e_max = eigmax(elli1.P-elli2.P)
    if e_max<0
        return false
    elseif  elli1.c==elli2.c
        return e_max>=0
    elseif !(elli1.c ∈ elli2)
        return false
    else 
        L = cholesky((elli2.P+elli2.P')/2).L
        P = L\elli1.P/L';
        c = L'*(elli1.c -elli2.c)
        specDecomp = eigen(P)
        vp = specDecomp.values
        ct = specDecomp.vectors'*c
        lb = 1/min(vp...)
        ub = 1-ct'*ct
        if(lb > ub)
            return false
        end
        f(β) = -(1-β + β*sum((vp./(1 .- β*vp)).*(ct.^2)))
        if bisectionB
            (val, _) = bisection(f, interval=[lb+1e-15, ub], verbose=false, stopIfNegative=true)
        else
            df(β) = -(1 - sum((vp./(1 .- β*vp).^2).*(ct.^2)))
            ddf(β) = -(-2*sum(((vp.^2)./(1 .- β*vp).^3).*(ct.^2)))
            (val, _) = newton_method(f, df, ddf; x0=4.0, interval=[lb,ub], verbose=false, stopIfNegative=true)
        end
        return val<=0
    end
end

function Base.in(x::AbstractVecOrMat, elli::Ellipsoid)
    return (x-elli.c)'elli.P*(x-elli.c) ≤ 1
end

function centerDistance(elli1::Ellipsoid,elli2::Ellipsoid)
    return norm(get_center(elli1)-get_center(elli2))
end

function pointCenterDistance(elli::Ellipsoid, x)
    return norm(get_center(elli)-x)
end

function volume(elli::Ellipsoid)
    N = size(elli.P,1)
    return pi^(N/2)/(gamma(N/2+1))*det(elli.P)^(-1/2)
end

function get_center(elli::Ellipsoid)
    return elli.c
end

function get_shape(elli::Ellipsoid)
    return elli.P
end

function get_dims(elli::Ellipsoid)
    return length(elli.c)
end

function scale(elli::Ellipsoid, α)
    return Ellipsoid(elli.P*(1/α),elli.c*α)
end

function expand(elli::Ellipsoid, α)
    return Ellipsoid(elli.P*(1/α),elli.c)
end

function plotE!(elli::Ellipsoid; color=:blue, opacity=1.0)
    P = get_shape(elli)
    Q = inv(P)
    Q = (Q+Q')./2
    E = LazySets.Ellipsoid(get_center(elli), Q) #not optimal, require to inverse
    Plots.plot!(E; color=color, opacity=opacity)
end

"""
    ellipsoid_vol(P,r)
    
Calculates the n-volume of the n-ellipsoid defined as {x'Px < r}.
"""
function ellipsoid_vol(elli::Ellipsoid) 
    N = size(P,1)
    return pi^(N/2)/(gamma(N/2+1))*det(elli.P)^(-1/2)
end


function get_farthest_point(elli::Ellipsoid, d)
    d = d/norm(d)
    Q = inv(elli.P)
    a = Q*d
    return a/sqrt(d'a)
end

"""
    _get_min_bounding_box(P, optimizer) 

Finds the minimum bounding box containing the ellipsoid {x'Px < 1}. 
"""
function get_min_bounding_box(elli::Ellipsoid; optim=false, optimizer=optimizer_with_attributes(Ipopt.Optimizer, MOI.Silent() => true)) 
    P = elli.P
    n = size(P,1)
    R = zeros(n)
    
    if optim
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
    else
        for i in 1:n
            ei = zeros(n)
            ei[i] = 1
            R[i] = get_farthest_point(elli, ei)[i]
        end
    end
    box = IntervalBox(elli.c .- R,elli.c .+ R)
    return box
end


