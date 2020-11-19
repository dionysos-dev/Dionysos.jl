export reachable_any_to_any

using MathematicalSystems, Polyhedra

"""
    reachable_any_to_any(set_from, set_to, map)

Returns a `Bool` indicating whether any element of `set_from` can reach any
element of `set_to` using map `map`.
"""
function reachable_any_to_any end

function reachable_any_to_any(set_from::HRep, set_to::HRep, m::ConstrainedLinearControlMap)
    return !isempty(polyhedron((hrep(set_from) * m.U) ∩ ([m.A m.B] \ hrep(set_to)), Polyhedra.library(set_from)))
end

#function reachable_any_to_any(sets, maps::AbstractVector{ConstrainedLinearControlMap})
#end
