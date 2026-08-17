# Empirical inflation stress test (PR #569's methodology, hardened): replay the
# certified feedback chain from the entry ellipsoid inflated by α ≥ 1 and
# report, per α, WHERE the rollouts fail — not just how many. Results at α > 1
# are NOT certificates: they measure the empirical success basin beyond the
# formal guarantee.

# Uniform sample in the α-inflated ellipsoid shell α ∈ [alpha_inner, alpha_outer]
# (alpha_inner = 0 gives the full ellipsoid): uniform direction, radial law
# ρ = (aᵢⁿ + U·(aₒⁿ − aᵢⁿ))^(1/n), mapped through an eigen square root of the
# shape matrix — exact for degenerate (needle) shapes, no rejection. Shells keep
# the per-α statistics from being dominated by points near the center.
function _sample_shell(E::LazySets.Ellipsoid, alpha_inner, alpha_outer, rng)
    c = collect(Float64, LazySets.center(E))
    n = length(c)
    F = LA.eigen(LA.Symmetric(Matrix{Float64}(LazySets.shape_matrix(E))))
    A = F.vectors * LA.Diagonal(sqrt.(max.(F.values, 0.0)))
    u = rng === nothing ? randn(n) : randn(rng, n)
    nu = LA.norm(u)
    nu == 0.0 && return c
    r01 = rng === nothing ? rand() : rand(rng)
    ρ = (alpha_inner^n + r01 * (alpha_outer^n - alpha_inner^n))^(1 / n)
    return c .+ A * ((ρ / nu) .* u)
end

_in_ellipsoid(x, E) =
    LA.dot(
        x .- collect(LazySets.center(E)),
        Matrix{Float64}(UT.get_quadratic_form(E)) * (x .- collect(LazySets.center(E))),
    ) <= 1.0 + 1e-9

"""
    inflation_stress(f, kappas, ellipsoids, target; alphas, n_samples, rng,
                     input_set = nothing, project_input = identity,
                     domain = nothing) -> Vector{NamedTuple}

Replay the certified feedback chain on the plant map `f((x, u) -> x⁺)` from
samples of the entry ellipsoid `ellipsoids[1]` inflated by each `α ∈ alphas`,
and report the outcome decomposed by FAILURE MODE. All arguments live in the
frame the chain was certified in (pass the z-frame map and sets when the chain
is normalized).

- `kappas` — the chain's affine feedbacks, forward-ordered (length `K`);
- `ellipsoids` — the funnel ellipsoids, forward-ordered (length `K + 1`,
  `ellipsoids[k]` is the validity region of `kappas[k]`);
- `target` — endpoint membership set;
- `alphas` — inflation factors (1.0 is forced in and sampled over the full
  ellipsoid; each α > 1 is sampled over the shell between it and its
  predecessor, so the rates are per-annulus, not cumulative);
- `input_set` — when given, the RAW feedback output `κ(x)` is checked against
  it before `project_input` is applied (an input violation means the
  certificate's feedback demanded more authority than the plant has);
- `domain` — when given, every visited state is checked against it.

Returns one row per α:
`(; alpha, n, success_rate, certified_rate, tube_exit_rate,
   domain_violation_rate, input_violation_rate, target_miss_rate,
   first_exit_median)` where `success` = endpoint in `target` with no domain
or input violation, `certified` additionally requires never leaving the funnel
chain, `tube_exit` counts rollouts that left their ellipsoid at some step
(`first_exit_median` locates where), and the two violation rates split the
constraint failures. At α = 1 the certificate guarantees
`certified_rate == 1.0` (up to solver/sampling roundoff); every α > 1 row is
an empirical measurement, NOT a certificate.
"""
function inflation_stress(
    f,
    kappas,
    ellipsoids,
    target;
    alphas = [1.0, 1.25, 1.5, 2.0, 3.0],
    n_samples = 200,
    rng = nothing,
    input_set = nothing,
    project_input = identity,
    domain = nothing,
)
    K = length(kappas)
    length(ellipsoids) == K + 1 ||
        error("expected length(ellipsoids) == length(kappas) + 1 (forward order).")
    E1 = ellipsoids[1]

    αs = sort(unique(vcat(1.0, collect(Float64, alphas))))
    rows = NamedTuple[]
    for (i, α) in enumerate(αs)
        α_inner = i == 1 ? 0.0 : αs[i - 1]
        n_success = 0
        n_certified = 0
        n_tube_exit = 0
        n_domain = 0
        n_input = 0
        n_miss = 0
        first_exits = Int[]
        for _ in 1:n_samples
            x = collect(Float64, _sample_shell(E1, α_inner, α, rng))
            first_exit = nothing
            domain_viol = false
            input_viol = false
            for k in 1:K
                first_exit === nothing &&
                    !_in_ellipsoid(x, ellipsoids[k]) &&
                    (first_exit = k)
                u_raw = Matrix(kappas[k].A) * x .+ collect(kappas[k].c)
                input_set !== nothing && u_raw ∉ input_set && (input_viol = true)
                x = collect(Float64, f(x, project_input(u_raw)))
                domain !== nothing && x ∉ domain && (domain_viol = true)
            end
            hit = x ∈ target
            success = hit && !domain_viol && !input_viol
            n_success += success
            n_certified += success && first_exit === nothing
            first_exit === nothing || (n_tube_exit += 1; push!(first_exits, first_exit))
            n_domain += domain_viol
            n_input += input_viol
            n_miss += !hit
        end
        push!(
            rows,
            (;
                alpha = α,
                n = n_samples,
                success_rate = n_success / n_samples,
                certified_rate = n_certified / n_samples,
                tube_exit_rate = n_tube_exit / n_samples,
                domain_violation_rate = n_domain / n_samples,
                input_violation_rate = n_input / n_samples,
                target_miss_rate = n_miss / n_samples,
                first_exit_median = isempty(first_exits) ? nothing :
                                    first_exits[div(length(first_exits) + 1, 2)],
            ),
        )
    end
    return rows
end

"""
    print_inflation_stress(rows; io = stdout)

Print the [`inflation_stress`](@ref) rows as a fixed-width table, one line per
inflation factor.
"""
function print_inflation_stress(rows; io = stdout)
    println(io, "  α      success  certified  tube-exit  domain  input  miss   first-exit")
    for r in rows
        println(
            io,
            "  ",
            rpad(round(r.alpha; digits = 2), 6),
            " ",
            rpad(round(r.success_rate; digits = 3), 8),
            " ",
            rpad(round(r.certified_rate; digits = 3), 10),
            " ",
            rpad(round(r.tube_exit_rate; digits = 3), 10),
            " ",
            rpad(round(r.domain_violation_rate; digits = 3), 7),
            " ",
            rpad(round(r.input_violation_rate; digits = 3), 6),
            " ",
            rpad(round(r.target_miss_rate; digits = 3), 6),
            " ",
            r.first_exit_median === nothing ? "—" : "k=$(r.first_exit_median)",
        )
    end
    return nothing
end
