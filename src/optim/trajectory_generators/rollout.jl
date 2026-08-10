# Typed rollout kernels (plan.md ★S1): simulate an input sequence and score it with
# no trajectory object and no per-sample allocation. Shared by every sampling
# generator (MPPI, CEM).
#
# Hard-constraint semantics (plan.md §3.1-E2): a rollout that leaves the state set is
# *frozen* — the state stops advancing, the step count keeps going — and the caller
# charges `violation_penalty` per frozen step. Rollouts are never truncated: every
# sample has a full-horizon, comparable cost, and crashing early is never cheap.

# Score `us` from `x0` with a CompositeCost, online. Returns (cost, n_frozen).
function rollout_cost(
    f,
    x0,
    us,
    wrap,
    cost::CompositeCost,
    X,
    hard::Bool,
    violation_penalty::Float64,
)
    x = wrap(x0)
    accs = _composite_init(cost)
    frozen = 0
    k = 0
    for u in us
        k += 1
        accs = _composite_step(cost, accs, x, u, k)
        if frozen == 0
            xn = wrap(f(x, u))
            if hard && !(xn ∈ X)
                frozen = 1
            else
                x = xn
            end
        else
            frozen += 1
        end
    end
    return _composite_final(cost, accs, x) + violation_penalty * frozen, frozen
end

# Full state rollout for callers that need the trajectory (the elitist candidate, the
# closure-cost fallback). Frozen steps repeat the last valid state so the channel
# invariant (states = inputs + 1) always holds. Returns (trajectory, n_frozen).
function rollout_trajectory(f, x0, us, wrap, X, hard::Bool)
    x = wrap(x0)
    xs = [x]
    frozen = 0
    for u in us
        if frozen == 0
            xn = wrap(f(x, u))
            if hard && !(xn ∈ X)
                frozen = 1
            else
                x = xn
            end
        else
            frozen += 1
        end
        push!(xs, x)
    end
    return ST.Trajectory(xs; inputs = collect(us)), frozen
end
