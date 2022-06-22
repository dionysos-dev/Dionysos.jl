module PathPlanning

using StaticArrays
using MathematicalSystems, HybridSystems
using Dionysos
const UT = Dionysos.Utils
const DO = Dionysos.Domain
const CO = Dionysos.Control

function _initTargetSets()
    _I_ = UT.HyperRectangle(SVector(0.4, 0.4, 0.0), SVector(0.4, 0.4, 0.0))
    _T_ = UT.HyperRectangle(SVector(3., 0.5, -100.0), SVector(3.6, 0.8, 100.0))
    return _I_, _T_
end

function system(;
    X1_lb = [1.0, 2.2,  2.2, 3.4,  4.6, 5.8,  5.8,  7.0, 8.2, 8.4,  9.3, 8.4,  9.3, 8.4,  9.3],
    X1_ub = [1.2, 2.4,  2.4, 3.6,  4.8, 6.0,  6.0,  7.2, 8.4, 9.3, 10.0, 9.3, 10.0, 9.3, 10.0],
    X2_lb = [0.0, 0.0,  6.0, 0.0,  1.0, 0.0,  7.0,  1.0, 0.0, 8.2,  7.0, 5.8,  4.6, 3.4,  2.2],
    X2_ub = [9.0, 5.0, 10.0, 9.0, 10.0, 6.0, 10.0, 10.0, 8.5, 8.6,  7.4, 6.2,  5.0, 3.8,  2.6]
)
    # Domains   
    _X_ = UT.HyperRectangle(SVector(0.0, 0.0, -pi-0.4), SVector(4., 10.0, pi+0.4))
    _U_ = UT.HyperRectangle(SVector(-1.0, -1.0), SVector(1.0, 1.0))

    # Remove obstacles from _X_
    _I_, _T_ = _initTargetSets() # Maybe remove this if not useful to check ∩ below
    obstacles = typeof(_I_)[]
    for (x1lb, x2lb, x1ub, x2ub) in zip(X1_lb, X2_lb, X1_ub, X2_ub)
        ob = UT.HyperRectangle(SVector(x1lb, x2lb, _X_.lb[3]), SVector(x1ub, x2ub, _X_.ub[3]))
        if ob ⊆ _X_ && isempty(ob ∩ _I_) && isempty(ob ∩ _T_)
            push!(obstacles, ob)
        end
    end
    obstacles_LU = UT.LazyUnionSetArray(obstacles)
    _X_ = UT.LazySetMinus(_X_, obstacles_LU)

    # System eq x' = F_sys(x, u)
    function F_sys(x, u)
        α = atan(tan(u[2])/2)
        return SVector{3}(
            u[1]*cos(α + x[3])/cos(α),
            u[1]*sin(α + x[3])/cos(α),
            u[1]*tan(u[2]))
    end
    system = ConstrainedBlackBoxControlContinuousSystem(F_sys, 3, 2, _X_, _U_)
    return system
end

""""
    problem()

This function create the system with `PathPlanning.system`.

Then, we define initial and target domains for the state of the system.

Finally, we instantiate our Reachability Problem with the system, the 
initial and target domains.
"""
function problem()
    sys = system()
    _I_, _T_ = _initTargetSets()
    problem = CO.ReachabilityProblem(sys, _I_, _T_)
    return problem
end

end