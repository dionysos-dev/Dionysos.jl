# Bidirectional handoff (plan.md §4.4/§6): run the forward chain from the initial
# set and the backward chain from the target on the same trajectory; the
# specification is certified as soon as the forward tube at some state s is
# contained in the backward funnel at s — neither direction has to close the whole
# horizon alone. The combined controller splices the forward κ's before the
# handoff with the backward κ's after it; soundness of the splice is exactly the
# containment: x ∈ E^f_s ⊆ E^b_s, so the backward guarantee takes over at s.

"""
    bidirectional_certify!(forward_cert::ForwardCertifier,
                           backward_cert::BackwardCertifier,
                           problem, traj)
        -> (; success, k_handoff, controller, forward_result, backward_result)

Certify `traj` with both directions and hand off at the first state `s` where the
forward tube is contained in the backward funnel (`UT.is_included`, exact
ellipsoid kernel). On success `controller` is the spliced
[`ST.FunnelController`](@ref); `k_handoff` is the 1-based state index of the
handoff (`nothing` on failure). Either chain may have failed partway — only the
overlap of their certified ranges is searched.
"""
function bidirectional_certify!(
    forward_cert::ForwardCertifier,
    backward_cert::BackwardCertifier,
    problem,
    traj,
)
    set_problem!(forward_cert, problem)
    set_trajectory!(forward_cert, traj)
    certify!(forward_cert)

    set_problem!(backward_cert, problem)
    set_trajectory!(backward_cert, traj)
    certify!(backward_cert)

    fres = get_result(forward_cert)
    bres = get_result(backward_cert)
    failed = (;
        success = false,
        k_handoff = nothing,
        controller = nothing,
        forward_result = fres,
        backward_result = bres,
    )
    (fres === nothing || bres === nothing) && return failed

    f_ell = fres.lmi_data.ellipsoids           # state s ↦ f_ell[s], s = 1..len_f
    b_ell = bres.lmi_data.ellipsoids           # state s ↦ b_ell[s − b_offset]
    b_offset = bres.failed_k === nothing ? 0 : bres.failed_k
    (isempty(f_ell) || isempty(b_ell)) && return failed

    # The backward funnel must actually reach the target for the handoff to close
    # the specification (its terminal gate); a backward chain that failed its
    # terminal containment cannot anchor a certificate.
    bres.terminal_contained === false && return failed

    for s in 1:length(f_ell)
        j = s - b_offset
        1 <= j <= length(b_ell) || continue
        UT.is_included(f_ell[s], b_ell[j]) || continue

        f_k = fres.lmi_data.kappas             # κ_i at step i, i = 1..len_f−1
        b_k = bres.lmi_data.kappas             # κ at steps b_offset+1 .. K
        (s - 1 <= length(f_k) && j <= length(b_k) + 1) || continue

        kappas = vcat(collect(f_k[1:(s - 1)]), collect(b_k[j:end]))
        ellipsoids = vcat(collect(f_ell[1:s]), collect(b_ell[(j + 1):end]))
        controller =
            length(ellipsoids) == length(kappas) + 1 ?
            ST.FunnelController(kappas, ellipsoids) : nothing

        return (;
            success = true,
            k_handoff = s,
            controller,
            forward_result = fres,
            backward_result = bres,
        )
    end

    return failed
end
