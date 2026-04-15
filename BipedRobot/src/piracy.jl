module Piracy

# To avoid ambiguities
import LinearAlgebra, Quaternions, RigidBodyDynamics, Symbolics

function Base.:/(q::Quaternions.Quaternion, x::Symbolics.Num)
    return Quaternions.Quaternion(q.s / x.val, q.v1 / x.val, q.v2 / x.val, q.v3 / x.val)
end

# Needed for `show`ing a `Transform3D`.
# See https://github.com/JuliaSymbolics/Symbolics.jl/issues/97#issuecomment-798844158
Symbolics.@register_symbolic Base.rem2pi(x, y::RoundingMode)

end
