module AdaptiveCruiseControl

import LazySets
using StaticArrays
using MathematicalSystems
using Plots

using Dionysos

const DI = Dionysos
const UT = DI.Utils
const PR = DI.Problem

"""
    Params

Longitudinal car-following parameters, from the adaptive-cruise-control benchmark of Ames,
Xu, Grizzle and Tabuada (IEEE TAC 2017). `c0`, `c1`, `c2` are the coefficients of the
resistance force `F_r(v) = c0 + c1 v + c2 v²`; `τ_h` is the time headway, also the ISO 15622
default.

`v_desired` deliberately exceeds `v_lead`: the set-point is faster than the car in front,
which is what makes the headway constraint bind rather than decorate the problem.

The operating envelope lives here rather than in a keyword of [`system`](@ref) because
[`jacobian_bound`](@ref) depends on it — the resistance term is bounded over the declared
speed range.
"""
Base.@kwdef struct Params{T}
    m::T = 1650.0            # vehicle mass [kg]
    c0::T = 0.1              # [N]
    c1::T = 5.0              # [N⋅s/m]
    c2::T = 0.25             # [N⋅s²/m²]
    v_lead::T = 13.89        # lead-vehicle speed, held constant [m/s]
    v_desired::T = 24.0      # cruising set-point [m/s]
    τ_h::T = 1.8             # time headway [s]
    a_max::T = 0.3 * 9.81    # wheel-friction acceleration limit [m/s²]
    z_min::T = 0.0           # gap bounds [m]
    z_max::T = 100.0
    v_min::T = 0.0           # ego speed bounds [m/s]
    v_max::T = 30.0
end

"Deceleration [m/s²] contributed by rolling resistance and drag at speed `v`."
rolling_resistance(p::Params, v) = (p.c0 + p.c1 * v + p.c2 * v^2) / p.m

"""
    dynamic(p)

`ż = v_lead - v`, `v̇ = a - F_r(v)/m`, with the state `x = (z, v)`: gap to the lead vehicle
and ego speed.

The input is an **acceleration**, not a wheel force, so it shares units with `a_max` and an
input-grid step reads directly as m/s².
"""
function dynamic(p::Params = Params())
    return (x, u) -> SVector{2}(p.v_lead - x[2], u[1] - rolling_resistance(p, x[2]))
end

"Exact Jacobian `∂f/∂x`, for the `LINEARIZED` approximation mode."
function jacobian(p::Params = Params())
    return (x, u) -> SMatrix{2, 2}(0.0, 0.0, -1.0, -(p.c1 + 2 * p.c2 * x[2]) / p.m)
end

"""
    jacobian_bound(p)

Growth bound `L` with `L[i,j] ≥ |J[i,j]|` off the diagonal and `L[i,i] ≥ J[i,i]` on it.

`∂v̇/∂v = -(c1 + 2 c2 v)/m` is a contraction that only strengthens with speed, so the diagonal
entry is taken at `v_min`, where it is weakest. The off-diagonal `|∂ż/∂v| = 1` is the only
coupling in the system: after a step, the gap uncertainty is the gap uncertainty plus the time
step times the speed uncertainty.
"""
function jacobian_bound(p::Params = Params())
    L22 = -(p.c1 + 2 * p.c2 * p.v_min) / p.m
    return u -> SMatrix{2, 2}(0.0, 0.0, 1.0, L22)
end

state_set(p::Params = Params()) = LazySets.Hyperrectangle(;
    low = SVector(p.z_min, p.v_min),
    high = SVector(p.z_max, p.v_max),
)

input_set(p::Params = Params()) =
    LazySets.Hyperrectangle(; low = SVector(-p.a_max), high = SVector(p.a_max))

function system(;
    params::Params = Params(),
    _X_ = state_set(params),
    _U_ = input_set(params),
)
    return MathematicalSystems.ConstrainedBlackBoxControlContinuousSystem(
        dynamic(params),
        LazySets.dim(_X_),
        LazySets.dim(_U_),
        _X_,
        _U_,
    )
end

"""
    safe_set(p)

`{(z, v) ∈ X : z ≥ τ_h v}` — the operating envelope cut by the time-headway constraint. This
is the `□` of every specification below; leaving it is a tailgating event.
"""
function safe_set(p::Params = Params())
    return LazySets.HPolytope([
        LazySets.HalfSpace(SVector(1.0, 0.0), p.z_max),
        LazySets.HalfSpace(SVector(-1.0, 0.0), -p.z_min),
        LazySets.HalfSpace(SVector(0.0, 1.0), p.v_max),
        LazySets.HalfSpace(SVector(0.0, -1.0), -p.v_min),
        LazySets.HalfSpace(SVector(-1.0, p.τ_h), 0.0),
    ])
end

# A speed band, restricted to the gaps that satisfy the headway at its fastest point (plus
# `margin`), so that every target set is a subset of `safe_set`.
function _speed_band(p::Params, v_center, ε, margin, gap_high)
    return LazySets.Hyperrectangle(;
        low = SVector(p.τ_h * (v_center + ε) + margin, v_center - ε),
        high = SVector(gap_high, v_center + ε),
    )
end

"""
    cruise_set(p; ε = 0.5, margin = 5.0, gap_high = p.z_max)

Band around the set-point `v_desired`.

It cannot be *held*: with `v_desired > v_lead` the gap closes at `v_desired - v_lead` and no
gap in the envelope lasts forever. Handing it alone to [`reach_and_stay_problem`](@ref) is
expected to be infeasible, which is the benchmark proving something true rather than failing.
"""
cruise_set(p::Params = Params(); ε = 0.5, margin = 5.0, gap_high = p.z_max) =
    _speed_band(p, p.v_desired, ε, margin, gap_high)

"""
    follow_set(p; ε = 0.5, margin = 5.0, gap_high = p.z_max)

Band around the lead's speed, with `margin` metres of gap in hand beyond the headway
requirement. This one *is* holdable — at `v = v_lead` the gap is stationary — and it is what
an adaptive cruise controller settles into when the set-point is out of reach.

`gap_high` closes the band from above. The default leaves it open, which asks only "match the
lead's speed at a safe distance" — from a large gap that is satisfied without ever closing it.
Setting it makes the specification a *following distance* rather than a floor, so the ego has
to overspeed, catch up, and then settle: the behaviour an adaptive cruise controller is bought
for.
"""
follow_set(p::Params = Params(); ε = 0.5, margin = 5.0, gap_high = p.z_max) =
    _speed_band(p, p.v_lead, ε, margin, gap_high)

"""
    safety_problem(; params, _I_)

`□ (z ≥ τ_h v)`. The certificate is the maximal controlled-invariant set, which
[`min_safe_gap`](@ref) brackets in closed form.
"""
function safety_problem(;
    params::Params = Params(),
    _I_ = LazySets.Hyperrectangle(; low = SVector(59.0, 19.8), high = SVector(61.0, 20.2)),
)
    return PR.SafetyProblem(system(; params = params), _I_, safe_set(params))
end

"""
    reach_and_stay_problem(; params, _I_, _T_, stay_on_first_entry = true)

`□ (z ≥ τ_h v) ∧ ◇□ _T_`. The default target is "cruise at the set-point **or** follow the
lead", which is what an adaptive cruise controller is actually asked to do. Pass
`_T_ = cruise_set(params)` for the set-point alone — see [`cruise_set`](@ref).

`stay_on_first_entry` defaults to `true` here, unlike
[`PR.ReachAndStayProblem`](@ref Dionysos.Problem.ReachAndStayProblem). Once the ego has settled
into a band it should not drop out of it and come back, so the stronger reading is the one that
matches the task; it is also the cheap one, being a safety solve followed by a reachability
solve rather than the nested fixed point ◇□ needs.
"""
function reach_and_stay_problem(;
    params::Params = Params(),
    _I_ = LazySets.Hyperrectangle(; low = SVector(69.0, 7.8), high = SVector(71.0, 8.2)),
    _T_ = UT.set_union([cruise_set(params), follow_set(params)]),
    stay_on_first_entry::Bool = true,
)
    return PR.ReachAndStayProblem(
        system(; params = params),
        _I_,
        _T_,
        safe_set(params);
        stay_on_first_entry = stay_on_first_entry,
    )
end

"""
    optimal_control_problem(; params, _I_, _T_, state_cost, transition_cost)

`(z ≥ τ_h v) U (matched speed)` — settle behind the lead from a slow, distant start without
ever violating the headway. The value function is the worst-case time to settle.
"""
function optimal_control_problem(;
    params::Params = Params(),
    _I_ = LazySets.Hyperrectangle(; low = SVector(69.0, 7.8), high = SVector(71.0, 8.2)),
    _T_ = follow_set(params),
    state_cost = nothing,
    transition_cost = nothing,
)
    return PR.OptimalControlProblem(
        system(; params = params),
        _I_,
        _T_,
        state_cost,
        transition_cost,
        PR.Infinity();
        safe_set = safe_set(params),
    )
end

"Strongest deceleration available at speed `v`: the friction limit plus the resistance that assists it."
max_deceleration(p::Params, v) = p.a_max + rolling_resistance(p, v)

"""
    min_safe_gap(p, v; a_brake = p.a_max)

Smallest gap from which `z ≥ τ_h v` can be held forever, at speed `v`, when the strongest
available deceleration is `a_brake`.

The barrier `h = z - τ_h v` obeys `ḣ = (v_lead - v) - τ_h v̇`, which depends on `(v, u)` alone
— the gap itself drops out — so it admits an exact one-dimensional analysis. Braking at
`a_brake` gives `ḣ = v⋆ - v` with `v⋆ = v_lead + τ_h a_brake`, so `h` can only fall while
`v > v⋆`, and riding the transient down to `v⋆` costs `(v - v⋆)²/(2 a_brake)` of gap.

`a_brake = p.a_max` ignores the resistance force, which assists braking, and so returns a gap
that is larger than necessary: an **inner** bound on the controlled-invariant set. Passing
`a_brake = max_deceleration(p, p.v_max)` returns an **outer** one. Together they bracket the
set an abstraction should compute — the one check on this benchmark that does not depend on
the abstraction being right.
"""
function min_safe_gap(p::Params, v; a_brake = p.a_max)
    v_star = p.v_lead + p.τ_h * a_brake
    excess = max(v - v_star, zero(v))
    return p.τ_h * v + excess^2 / (2 * a_brake)
end

function _rectangle!(fig, xlo, xhi, ylo, yhi; kwargs...)
    return plot!(
        fig,
        [xlo, xhi, xhi, xlo, xlo],
        [ylo, ylo, yhi, yhi, ylo];
        label = "",
        kwargs...,
    )
end

"""
    system_plot!(; kwargs...)

Closure `(fig, x, u) -> fig` drawing the two vehicles on the road, for
`Dionysos.animate_trajectory_dashboard`. The bar ahead of the ego is the *required* headway
`τ_h v`: the specification holds exactly while it stops short of the lead vehicle, and the bar
turns red when it does not.

The frame is **pinned to the ego**, and the state is the *gap*, not distance travelled — so
neither car slides along the road, and what moves is the lead, closing in or pulling away. Both
speeds are drawn so the absolute motion the frame hides stays readable.
"""
function system_plot!(;
    params::Params = Params(),
    car_length = 4.5,
    car_width = 2.2,
    xlims = (-12.0, 112.0),
    ylims = (-9.0, 9.0),
    show_state = true,
)
    return function (fig, x, u)
        z = Float64(x[1])
        v = Float64(x[2])
        a = Float64(u[1])

        required = params.τ_h * v

        _rectangle!(fig, xlims[1], xlims[2], -4.0, 4.0; color = :gray, fill = (true, 0.15))
        plot!(
            fig,
            [xlims[1], xlims[2]],
            [0.0, 0.0];
            lw = 1,
            ls = :dot,
            color = :gray,
            label = "",
        )

        _rectangle!(
            fig,
            -car_length,
            0.0,
            -car_width / 2,
            car_width / 2;
            lw = 2,
            color = :blue,
            fill = (true, 0.5),
        )
        _rectangle!(
            fig,
            z,
            z + car_length,
            -car_width / 2,
            car_width / 2;
            lw = 2,
            color = :black,
            fill = (true, 0.35),
        )

        # The frame hides that both cars are moving; the labels put it back.
        annotate!(
            fig,
            -car_length / 2,
            2.6,
            text("ego $(round(v; digits = 1)) m/s", 8, :center, :blue),
        )
        annotate!(
            fig,
            z + car_length / 2,
            2.6,
            text("lead $(round(params.v_lead; digits = 1)) m/s", 8, :center),
        )
        annotate!(
            fig,
            (xlims[1] + xlims[2]) / 2,
            8.0,
            text("view pinned to the ego — the road is not scrolling", 7, :center, :gray),
        )

        # The headway requirement drawn to scale, so the constraint is read off the picture.
        plot!(
            fig,
            [0.0, required],
            [-5.5, -5.5];
            lw = 4,
            color = required > z ? :red : :green,
            label = "",
        )
        annotate!(
            fig,
            required / 2,
            -7.0,
            text("τ_h⋅v = $(round(required; digits = 1)) m", 8),
        )

        if show_state
            annotate!(
                fig,
                xlims[1] + 0.03 * (xlims[2] - xlims[1]),
                6.5,
                text(
                    "v = $(round(v; digits = 1)) m/s\n" *
                    "gap = $(round(z; digits = 1)) m\n" *
                    "a = $(round(a; digits = 2)) m/s²",
                    9,
                    :left,
                ),
            )
        end

        xlims!(fig, xlims...)
        ylims!(fig, ylims...)

        return fig
    end
end

end # module
