struct Ellipsoid{MT,VT} 
    P::MT
    c::VT
        Ellipsoid(P,c)= new((P+P')./2,c)
    
end



function Base.in(x, elli::Ellipsoid)
    return (x-elli.c)'elli.P*(x-elli.c) ≤ 1
end



function Base.in(elli1::Ellipsoid, elli2::Ellipsoid)
    e_max = eigmax(elli1.P-elli2.P)
    if e_max<0
        return false
    elseif  elli1.c==elli2.c
        return e_max>0
    elseif ~(elli1.c ∈ elli2)
        return false
    else 
        L = cholesky(elli2.P).U
        P = L'\P/L;
        specDecomp = eigen(P)
        lb = specDecomp.values
        ct = specDecomp.vectors'*L*(elli1.c -elli2.c)
        α = max(lb...)
        polPos(β) = -(1-β + β*sum((lb./(1-β*lb)).*(ct.^2)))
        (val, _) = bissection(polPos, interval=[α+1e-5, 1], stopIfNegative=true)

        return val<=0
    end

end

function volume(elli::Ellipsoid)
    N = size(elli.P,1)
    return pi^(N/2)/(gamma(N/2+1))*det(P)^(-1/2)
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
