# Polyhedron defined by {x : |[Ax]_i| ≦ b_i ∀ i}
struct CenteredPolyhedron{MT,VT}
    A::MT
    b::VT
end

# Move each face so that it contains H + x where |x_i| ≦ e_i ∀ i.
# function inflate(H::CenteredPolyhedron, e) where {T,S}
#     b = H.b + abs.(H.A)*e
#     return CenteredPolyhedron(H.A, b)
# end

function Base.in(x, H::CenteredPolyhedron)
    return all(abs.(H.A*x) .<= H.b)
end
# all(x .<= y) is (surprisingly) faster than all(i -> x[i] <= y[i], eachindex(x))
# See also test_performances
