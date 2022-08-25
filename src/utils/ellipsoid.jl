using SpecialFunctions
struct Ellipsoid{T<:Real,MT<:AbstractMatrix{T},VT<:AbstractVector{T}}
    P::MT
    c::VT
    # function Ellipsoid{T,MT,VT}(P::MT,c::VT) where {T<:Real,MT<:AbstractMatrix{T},VT<:AbstractVector{T}}
    #     new((P+P')./2,c)
    # end
end



function Base.in(x, elli::Ellipsoid)
    return (x-elli.c)'elli.P*(x-elli.c) ≤ 1
end



function Base.in(elli1::Ellipsoid, elli2::Ellipsoid)
    e_max = eigmax(elli1.P-elli2.P)
    if e_max<0
        return false
    elseif  elli1.c==elli2.c
        return e_max>=0
    elseif !(elli1.c ∈ elli2)
        return false
    else 
        L = cholesky((elli2.P+elli2.P')/2).U #TODO remove ' when Ellipsoid constructor fixed
        P = L'\elli1.P/L;
        specDecomp = eigen(P)
        lb = specDecomp.values
        ct = specDecomp.vectors'*L*(elli1.c -elli2.c)
        α = 1/min(lb...)

        polPos(β) = -(1-β + β*sum((lb./(1 .- β*lb)).*(ct.^2)))
        if(α==1)
            return polPos(1)<=0
        end
        (val, _) = bisection(polPos, interval=[α+1e-15, 1-norm(ct)], verbose=false, stopIfNegative=true)

        return val<=0
    end

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

function get_dims(elli::Ellipsoid)
    return length(elli.c)
end

function scale(elli::Ellipsoid, α)
    return Ellipsoid(elli.P*(1/α),elli.c*α)
end
