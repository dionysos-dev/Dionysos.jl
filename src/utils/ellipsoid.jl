using SpecialFunctions
using Plots
using LazySets
using LinearAlgebra
using Polyhedra
using IntervalArithmetic
using JuMP

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

function get_center(elli::Ellipsoid)
    return elli.c
end

function get_shape(elli::Ellipsoid)
    return elli.P
end

function get_dims(elli::Ellipsoid)
    return length(elli.c)
end

function centerDistance(elli1::Ellipsoid,elli2::Ellipsoid)
    return norm(get_center(elli1)-get_center(elli2))
end

function pointCenterDistance(elli::Ellipsoid, x)
    return norm(get_center(elli)-x)
end

function get_volume(elli::Ellipsoid)
    N = size(elli.P,1)
    return pi^(N/2)/(gamma(N/2+1))*det(elli.P)^(-1/2)
end

function Base.:*(elli::Ellipsoid, r::Real)
    Ellipsoid(elli.P*(1/r), elli.c)
end

function Base.:*(r::Real, elli::Ellipsoid)
    elli*r
end

function Base.:/(elli::Ellipsoid, r::Real)
    elli*(1/r)
end

function scale(elli::Ellipsoid, α)
    return Ellipsoid(elli.P*(1/α),elli.c*α)
end

function expand(elli::Ellipsoid, α)
    return Ellipsoid(elli.P*(1/α),elli.c)
end

function plotE!(elli::Ellipsoid; color=:blue, opacity=1.0, label="")
    P = get_shape(elli)
    Q = inv(P)
    Q = (Q+Q')./2
    E = LazySets.Ellipsoid(get_center(elli), Q) #not optimal, require to inverse
    Plots.plot!(E; color=color, opacity=opacity, label=label)
end


# get the farthest point of the ellipsoid in direction d
function get_farthest_point(elli::Ellipsoid, d)
    d = d/norm(d)
    Q = inv(elli.P)
    a = Q*d
    return a/sqrt(d'a)
end

"""
    get_min_bounding_box(elli, optimizer) 

Finds the minimum bounding box containing the ellipsoid {(x-c)'P(x-c) < 1}. 
"""
function get_min_bounding_box(elli::Ellipsoid; optim=false, optimizer) 
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


include("ellipsoid_inclusion.jl")
include("ellipsoid_intersection.jl")


# compress E1 if E1∩E2≠∅
# return nothing if impossible
function compress_if_intersection(E1::Ellipsoid, E2::Ellipsoid)
    if E1∩E2
        return scale_for_noninclusion_contact_point(E1, E2)
    else
        return E1
    end
end
