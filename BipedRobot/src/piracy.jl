# This file include type piracy that is needed in order to make Symbolics.jl
# types interact will with Quaternions, etc...
# https://docs.julialang.org/en/v1/manual/style-guide/#Avoid-type-piracy
# In the long term, we should try to merge these fixes upstream so we should
# keep references here to open issues and/or pull requests to keep track of when we can remove these type piracies

module Piracy

# To avoid ambiguities
import LinearAlgebra, Quaternions, RigidBodyDynamics, Symbolics
function LinearAlgebra.normalize(v::RigidBodyDynamics.FreeVector3D, p::Real)
    return RigidBodyDynamics.FreeVector3D(v.frame, LinearAlgebra.normalize(v.v, p))
end
function Base.:/(q::Quaternions.Quaternion, x::Symbolics.Num)
    return Quaternions.Quaternion(q.s / x.val, q.v1 / x.val, q.v2 / x.val, q.v3 / x.val)
end

# Needed for `show`ing a `Transform3D`.
# See https://github.com/JuliaSymbolics/Symbolics.jl/issues/97#issuecomment-798844158
Symbolics.@register_symbolic Base.rem2pi(x, y::RoundingMode)

end
