# Prefix re-planning (plan.md §6-2): when the backward chain fails at `failed_k`,
# the suffix funnel (states failed_k+1 … K+1) is already certified and stays valid.
# Re-plan ONLY the prefix — the generator is retargeted at the suffix's entry
# ellipsoid — then certify the prefix with its terminal pinned inside that entry,
# and splice the two funnels into one controller. Every round reuses all certified
# work; the splice is gated on exact ellipsoid inclusion, so it is sound by the
# same containment argument as the bidirectional handoff.

"""
    prefix_replan_certify!(gen, cert, failed_result;
                           gen_problem, seed = nothing,
                           prepare = identity, backmap = identity,
                           retarget_cost! = E -> nothing,
                           margin = 5, terminal_shrink = 0.5, max_rounds = 3)

Recover a failed backward certification by re-planning its prefix.

- `failed_result` — the failed [`CertificationResult`](@ref) (its `lmi_data` holds
  the certified suffix, entry ellipsoid first);
- `gen_problem` — the generator-frame problem; the prefix problem reuses its
  initial set but targets the suffix entry (mapped through `backmap`, e.g. a
  de-normalization);
- `seed` — warm-start trajectory for the prefix (typically the failed trajectory;
  the generator trims it to the prefix horizon);
- `prepare` — generator-frame → certifier-frame trajectory transform (the same
  hook the driver uses);
- `retarget_cost!` — caller hook to re-aim the generator's shaping terms at the
  new (generator-frame) target ellipsoid, e.g. mutating a `TerminalEllipsoidCost`'s
  arrays in place;
- `margin` — extra steps granted to the prefix beyond `failed_k`;
- `terminal_shrink` — the prefix chain's terminal is the entry ellipsoid scaled by
  this factor, centered at the prefix endpoint; the splice gate then requires it
  to lie inside the entry ellipsoid (`UT.is_included`, exact kernel).

Returns `(; success, controller, k_prefix, rounds, prefix_result)` — on success
`controller` is the spliced [`ST.FunnelController`](@ref) covering the whole
horizon.
"""
function prefix_replan_certify!(
    gen,
    cert,
    failed_result::CertificationResult;
    gen_problem,
    seed = nothing,
    prepare = identity,
    backmap = identity,
    retarget_cost! = E -> nothing,
    margin::Int = 5,
    terminal_shrink::Float64 = 0.5,
    max_rounds::Int = 3,
)
    k_f = failed_result.failed_k
    k_f === nothing &&
        error("prefix_replan_certify! needs a *failed* result with a certified suffix")
    suffix_ell = failed_result.lmi_data.ellipsoids
    suffix_kap = failed_result.lmi_data.kappas
    isempty(suffix_ell) && error("the failed result has no certified suffix to splice on")

    E_entry = suffix_ell[1]                       # certifier frame, at state k_f + 1
    E_entry_gen = backmap(E_entry)                # generator frame

    prefix_problem = PR.OptimalControlProblem(
        gen_problem.system,
        gen_problem.initial_set,
        E_entry_gen,
        gen_problem.state_cost,
        gen_problem.transition_cost,
        gen_problem.time,
        gen_problem.safe_set,
    )
    retarget_cost!(E_entry_gen)

    set_problem!(gen, prefix_problem)
    seed === nothing || set_seed_trajectory!(gen, seed)
    gen.nstep = k_f + margin

    # The warm-start seed passes through the entry's center, so it already
    # "succeeds" against the prefix target — the generator must keep optimizing
    # (smoothness, terminal centering) past success or no exploration happens.
    restore_stop = hasproperty(gen, :stop_on_success) ? gen.stop_on_success : nothing
    restore_stop === nothing || (gen.stop_on_success = false)

    # The prefix chain's own terminal: a shrunk copy of the entry ellipsoid — the
    # chain centers it at the prefix endpoint, so containment in E_entry is what
    # the splice gate verifies afterwards.
    cert.options.terminal_shape =
        terminal_shrink^2 .* Matrix{Float64}(LazySets.shape_matrix(E_entry))
    cert.options.check_terminal = false           # the suffix already closes the spec

    restore!() = restore_stop === nothing || (gen.stop_on_success = restore_stop)

    for round in 1:max_rounds
        generate!(gen)
        get_success(gen) || continue

        prefix_traj = prepare(get_trajectory(gen))
        set_trajectory!(cert, prefix_traj)
        certify!(cert)
        pres = get_result(cert)
        (pres === nothing || !pres.success) && continue

        prefix_ell = pres.lmi_data.ellipsoids
        prefix_kap = pres.lmi_data.kappas

        # Splice gate: the prefix's terminal ellipsoid must sit inside the
        # certified suffix's entry.
        UT.is_included(prefix_ell[end], E_entry) || continue

        kappas = vcat(collect(prefix_kap), collect(suffix_kap))
        ellipsoids = vcat(collect(prefix_ell[1:(end - 1)]), collect(suffix_ell))
        controller =
            length(ellipsoids) == length(kappas) + 1 ?
            ST.FunnelController(kappas, ellipsoids) : nothing

        restore!()
        return (;
            success = controller !== nothing,
            controller,
            k_prefix = length(prefix_kap),
            rounds = round,
            prefix_result = pres,
        )
    end

    restore!()
    return (;
        success = false,
        controller = nothing,
        k_prefix = 0,
        rounds = max_rounds,
        prefix_result = get_result(cert),
    )
end
