# Temporary fix to push upstream
using JuMP, Polyhedra
Polyhedra.convexhull!(v::Polyhedra.Hull, r::Polyhedra.VStruct) = Polyhedra.convexhull!(v.rays, r)
Polyhedra.convexhull!(v::Polyhedra.RaysHull, r::Ray) = push!(v.rays, r)
function _LPHRep(model::MOI.ModelLike, T::Type = Float64)
    _model = Polyhedra._MOIModel{T}()
    bridged = MOI.Bridges.LazyBridgeOptimizer(_model)
    # Only enable constraint bridges that don't create variables and don't add
    # any variable bridge so that there is an identity mapping betwenen
    # variables of `model` and polyhedra dimensions.
    MOI.Bridges.add_bridge(bridged, MOI.Bridges.Constraint.GreaterToLessBridge{T})
    MOI.Bridges.add_bridge(bridged, MOI.Bridges.Constraint.LessToGreaterBridge{T})
    MOI.Bridges.add_bridge(bridged, MOI.Bridges.Constraint.NonnegToNonposBridge{T})
    MOI.Bridges.add_bridge(bridged, MOI.Bridges.Constraint.NonposToNonnegBridge{T})
    MOI.Bridges.add_bridge(bridged, MOI.Bridges.Constraint.ScalarizeBridge{T})
    MOI.Bridges.add_bridge(bridged, MOI.Bridges.Constraint.VectorizeBridge{T})
    MOI.Bridges.add_bridge(bridged, MOI.Bridges.Constraint.ScalarFunctionizeBridge{T})
    MOI.Bridges.add_bridge(bridged, MOI.Bridges.Constraint.VectorFunctionizeBridge{T})
    MOI.Bridges.add_bridge(bridged, MOI.Bridges.Constraint.SplitIntervalBridge{T})
    MOI.Bridges.add_bridge(bridged, MOI.Bridges.Constraint.NormInfinityBridge{T})
    # This one creates variables so the user need to consider the polyhedra as the
    # feasible set of the extended formulation.
    MOI.Bridges.add_bridge(bridged, MOI.Bridges.Constraint.NormOneBridge{T})
    MOI.Bridges.add_bridge(bridged, Polyhedra.PolyhedraToLPBridge{T})
    MOI.copy_to(bridged, model)
    return LPHRep(_model)
end
