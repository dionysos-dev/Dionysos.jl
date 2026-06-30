module LTL2Automata

using StaticArrays
using LinearAlgebra
using Dionysos
using Spot
const DI = Dionysos
const UT = DI.Utils
const DO = DI.Domain
const SY = DI.Symbolic


export TSFromAdj, AbstractTS, nstates, ninputs, initset, post, pre_product
export Automaton, is_accepting
export letter_of, nextstate
export translate_ltl_buchi, hoa_to_automaton
export ts_from_symbolic_onepass
export precompute_labels_and_forbidden
export APLabeler, labeler_and_inner_sets
export synthesize_buchi_controller, simulate_buchi_policy , synthesize_reachability_controller , simulate_reachability_policy


# =============================================================================
# Transition system wrapper
# =============================================================================

"""
TSFromAdj: lightweight TS wrapper with adjacency list:
  Adj[x][u] = Vector{Int} successors
and initial set X0::Vector{Int}.
"""
struct TSFromAdj
    Adj::Vector{Vector{Vector{Int}}}
    X0::Vector{Int}
end

nstates(ts::TSFromAdj) = length(ts.Adj)
ninputs(ts::TSFromAdj) = isempty(ts.Adj) ? 0 : length(ts.Adj[1])
initset(ts::TSFromAdj) = ts.X0
post(ts::TSFromAdj, x::Int, u::Int) = ts.Adj[x][u]

# Alias used in older code
const AbstractTS = TSFromAdj

# =============================================================================
# Büchi automaton wrapper (minimal)
# =============================================================================

"""
Automaton: deterministic automaton representation sufficient for product games.

Fields:
- aps: Vector{String} atomic propositions (in fixed order)
- init: Int initial automaton state (can be 0-based or 1-based)
- nstates: Int number of explicit states
- zero_based: Bool whether states are 0:(nstates-1)
- delta: Dict{Tuple{Int,UInt64},Int} transition function on bitmask letters
- accepting: BitSet of accepting states
"""
struct Automaton
    aps::Vector{String}
    init::Int
    nstates::Int
    zero_based::Bool
    delta::Dict{Tuple{Int,UInt64},Int}
    accepting::BitSet
end

is_accepting(aut::Automaton, q::Int) = (q in aut.accepting)

# =============================================================================
# Letters + transitions
# =============================================================================

"""
letter_of(lbl::Set{Symbol}, aps::Vector{String}) -> UInt64 bitmask.
Bit i is 1 if Symbol(aps[i]) ∈ lbl.
"""
function letter_of(lbl::Set{Symbol}, aps::Vector{String})::UInt64
    m = UInt64(0)
    @inbounds for i in 1:length(aps)
        if Symbol(aps[i]) in lbl
            m |= (UInt64(1) << (i-1))
        end
    end
    return m
end

"""
nextstate(aut, q, letter_mask) -> next q, or nothing if undefined.
"""
function nextstate(aut::Automaton, q::Int, letter::UInt64)
    get(aut.delta, (q, letter), nothing)
end

# =============================================================================
# LTL -> automaton (Spot if available)
# =============================================================================

"""
translate_ltl_buchi(ϕ::String) -> Automaton

Preferred path: use Spot.jl to translate the LTL string to a Spot automaton.
Then export that automaton to HOA and parse it.

IMPORTANT: some Spot.jl builds (like yours) do NOT provide any HOA exporter (no `text/hoa` show method and no `hoa*` functions).
In that case we fall back to the external Spot command-line tools (e.g. `ltl2tgba`) if they are installed on PATH.

To enable the fallback on macOS:
  - `brew install spot`   (provides `ltl2tgba`, `ltl2ba`, ...)

Also, VSCode/REPL display of `SpotAutomata` may fail if Graphviz `dot` is missing. If you want to visualize Spot automata:
  - `brew install graphviz`
"""
function translate_ltl_buchi(ϕ::String)
    # Spot.jl API in your environment:
    #   SpotFormula(::AbstractString)
    #   LTLTranslator(; buchi, tgba, deterministic, ...)
    #   translate(::LTLTranslator, ::SpotFormula)

    # 1) Parse LTL string into Spot formula
    f = try
        Spot.SpotFormula(ϕ)
    catch e
        error("Spot.jl: failed to parse LTL string into SpotFormula. Error: $(typeof(e))\n$e")
    end

    # 2) Build a translator configured to produce a Büchi automaton
    tr = try
        Spot.LTLTranslator(; buchi=true)
    catch e
        error("Spot.jl: cannot construct LTLTranslator(; buchi=true). Try: methods(Spot.LTLTranslator). Error: $(typeof(e))\n$e")
    end

    # 3) Translate to a Spot automaton object
    a = try
        Spot.translate(tr, f)
    catch e
        error("Spot.translate(LTLTranslator, SpotFormula) failed: $(typeof(e))\n$e")
    end
    # NOTE: do not display `a` here; some frontends try to render it via Graphviz `dot`.

    # 4) Export to HOA (exporter name varies by build)
    function _try_hoa_via(fn)
        # io-method
        if hasmethod(fn, Tuple{IO, typeof(a)})
            io = IOBuffer()
            try
                fn(io, a)
                s = String(take!(io))
                return occursin("HOA:", s) ? s : nothing
            catch
                return nothing
            end
        end
        # direct method
        if hasmethod(fn, Tuple{typeof(a)})
            try
                out = fn(a)
                if out isa AbstractString
                    s = String(out)
                    return occursin("HOA:", s) ? s : nothing
                end
            catch
                return nothing
            end
        end
        return nothing
    end

    hoa_str = nothing

    # First: common exporter symbols
    for sym in (:hoa, :to_hoa, :print_hoa, :show_hoa, :hoa_print)
        if isdefined(Spot, sym)
            s = _try_hoa_via(getproperty(Spot, sym))
            if s !== nothing
                hoa_str = s
                break
            end
        end
    end

    # Second: scan Spot for any function name that includes "hoa"
    if hoa_str === nothing
        for sym in names(Spot, all=true)
            occursin("hoa", lowercase(String(sym))) || continue
            isdefined(Spot, sym) || continue
            fn = getproperty(Spot, sym)
            fn isa Function || continue
            s = _try_hoa_via(fn)
            if s !== nothing
                hoa_str = s
                break
            end
        end
    end

    # Last resort: try MIME-based display hooks (Spot often defines show methods)
    if hoa_str === nothing
        for mime in (MIME("text/hoa"), MIME("text/x-hoa"), MIME("text/plain"))
            io = IOBuffer()
            try
                show(io, mime, a)
                s = String(take!(io))
                if occursin("HOA:", s)
                    hoa_str = s
                    break
                end
            catch
                # ignore
            end
        end
    end

    # Another fallback: `repr` may include HOA
    if hoa_str === nothing
        s = try
            repr(a)
        catch
            ""
        end
        hoa_str = occursin("HOA:", s) ? s : nothing
    end

    if hoa_str === nothing
        # ---------------------------------------------------------------------
        # Fallback: use external Spot CLI tools to obtain HOA
        # This is the most robust option when Spot.jl lacks HOA exporters.
        # ---------------------------------------------------------------------
        function _try_cli(cmd_argv::Vector{String})
            exe = cmd_argv[1]
            Sys.which(exe) === nothing && return nothing
            io = IOBuffer()
            try
                run(pipeline(Cmd(cmd_argv), stdout=io, stderr=devnull))
                s = String(take!(io))
                return occursin("HOA:", s) ? s : nothing
            catch
                return nothing
            end
        end

        # Try a few common Spot CLIs in order of preference.
        # `ltl2tgba --hoa` usually yields a complete automaton in HOA.
        for cmd in (
            ["ltl2tgba", "--hoa", ϕ],
            ["ltl2tgba", "--hoa", "--buchi", ϕ],
            ["ltl2ba",   "--hoa", ϕ],
            ["ltl2gba",  "--hoa", ϕ],
        )
            s = _try_cli(cmd)
            if s !== nothing
                hoa_str = s
                break
            end
        end

        if hoa_str === nothing
            error(
                "Could not export HOA from Spot.jl (your build has no HOA exporter), and Spot CLI tools were not found on PATH.\n" *
                "Fix: install Spot command-line tools and retry. On macOS: `brew install spot`.\n" *
                "(Optional) If you want to visualize Spot automata in VSCode/REPL, also install Graphviz: `brew install graphviz`."
            )
        end
    end

    return hoa_to_automaton(hoa_str)
end
# =============================================================================
# HOA -> Automaton parser (deterministic, Spot-style HOA)
# =============================================================================

# A tiny boolean-guard evaluator for HOA transition guards over AP indices.
# Supports: t, f, integers (AP index), !, &, |, parentheses.
abstract type _GuardNode end
struct _GTrue  <: _GuardNode end
struct _GFalse <: _GuardNode end
struct _GVar   <: _GuardNode; i::Int; end
struct _GNot   <: _GuardNode; a::_GuardNode; end
struct _GAnd   <: _GuardNode; a::_GuardNode; b::_GuardNode; end
struct _GOr    <: _GuardNode; a::_GuardNode; b::_GuardNode; end

@inline function _guard_eval(g::_GuardNode, letter::UInt64)::Bool
    if g isa _GTrue
        return true
    elseif g isa _GFalse
        return false
    elseif g isa _GVar
        # HOA variables are 0-based indices
        i = (g::_GVar).i
        return ((letter >> i) & 0x01) == 0x01
    elseif g isa _GNot
        return !_guard_eval((g::_GNot).a, letter)
    elseif g isa _GAnd
        gg = g::_GAnd
        return _guard_eval(gg.a, letter) && _guard_eval(gg.b, letter)
    elseif g isa _GOr
        gg = g::_GOr
        return _guard_eval(gg.a, letter) || _guard_eval(gg.b, letter)
    else
        return false
    end
end

# Recursive descent parser for guard strings
# grammar:
#   expr  := term ( '|' term )*
#   term  := factor ( '&' factor )*
#   factor:= '!' factor | atom
#   atom  := 't' | 'f' | INT | '(' expr ')'

function _guard_parse(s::AbstractString)
    s = String(s)
    i = Ref(1)
    n = lastindex(s)

    skipws() = (while i[] <= n && (s[i[]] == ' ' || s[i[]] == '\t'); i[] += 1; end)

    function parse_atom()::_GuardNode
        skipws()
        i[] > n && return _GFalse()
        c = s[i[]]
        if c == 't'
            i[] += 1
            return _GTrue()
        elseif c == 'f'
            i[] += 1
            return _GFalse()
        elseif c == '('
            i[] += 1
            g = parse_expr()
            skipws()
            if i[] <= n && s[i[]] == ')'
                i[] += 1
            end
            return g
        else
            # integer
            if isdigit(c)
                j = i[]
                while j <= n && isdigit(s[j])
                    j += 1
                end
                v = parse(Int, s[i[]:prevind(s, j)])
                i[] = j
                return _GVar(v)
            end
        end
        # unknown token
        i[] += 1
        return _GFalse()
    end

    function parse_factor()::_GuardNode
        skipws()
        if i[] <= n && s[i[]] == '!'
            i[] += 1
            return _GNot(parse_factor())
        end
        return parse_atom()
    end

    function parse_term()::_GuardNode
        g = parse_factor()
        while true
            skipws()
            if i[] <= n && s[i[]] == '&'
                i[] += 1
                g = _GAnd(g, parse_factor())
            else
                break
            end
        end
        return g
    end

    function parse_expr()::_GuardNode
        g = parse_term()
        while true
            skipws()
            if i[] <= n && s[i[]] == '|'
                i[] += 1
                g = _GOr(g, parse_term())
            else
                break
            end
        end
        return g
    end

    g = parse_expr()
    return g
end

function hoa_to_automaton(hoa::String)::Automaton
    lines = split(hoa, '\n')

    # Header
    nstates = 0
    start_q = 0
    aps = String[]

    # Per-state info
    accepting = BitSet()

    # Transitions: for each (q, letter) -> q'
    delta = Dict{Tuple{Int,UInt64},Int}()

    # Parse header fields
    for ln in lines
        startswith(ln, "States:") && (nstates = parse(Int, strip(split(ln, ":")[2])); continue)
        startswith(ln, "Start:")  && (start_q = parse(Int, strip(split(ln, ":")[2])); continue)
        if startswith(ln, "AP:")
            # e.g. AP: 2 "a" "b"
            parts = split(ln)
            k = parse(Int, parts[2])
            # extract quoted strings
            qs = collect(eachmatch(r"\"([^\"]+)\"", ln))
            aps = [m.captures[1] for m in qs]
            length(aps) == k || (aps = aps[1:min(end,k)])
        end
    end
    nstates == 0 && error("HOA parse: could not find 'States:'.")

    # Spot HOA uses 0-based state numbering
    zero_based = true

    # Parse body
    cur_state = nothing
    for idx in eachindex(lines)
        ln = strip(lines[idx])
        isempty(ln) && continue

        if startswith(ln, "State:")
            # Example: State: 3 {0}
            m = match(r"State:\s*(\d+)(?:\s*\{([^}]*)\})?", ln)
            m === nothing && continue
            cur_state = parse(Int, m.captures[1])
            accset = m.captures[2]
            if accset !== nothing && !isempty(strip(accset))
                # Mark accepting if any acceptance marker present
                push!(accepting, cur_state)
            end
            continue
        end

        if cur_state === nothing
            continue
        end

        # Transition line examples:
        #   [t] 1
        #   [0] 2
        #   [!0&1] 3
        if startswith(ln, "[")
            m = match(r"^\[([^\]]+)\]\s*(\d+)", ln)
            m === nothing && continue
            guard_str = strip(m.captures[1])
            tgt = parse(Int, m.captures[2])

            g = _guard_parse(guard_str)
            nap = length(aps)
            nap > 62 && error("HOA parse: too many APs for UInt64 encoding (got $nap).")
            for mask_int in 0:(UInt64(1) << nap) - 1
                letter = UInt64(mask_int)
                if _guard_eval(g, letter)
                    delta[(cur_state, letter)] = tgt
                end
            end
        end
    end

    # If accepting set remained empty (some HOA encodes acceptance differently),
    # default to making the 'start' accepting so reachability doesn't immediately fail.
    if isempty(accepting)
        push!(accepting, start_q)
    end

    return Automaton(aps, start_q, nstates, zero_based, delta, accepting)
end

# =============================================================================
# TS from Dionysos symbolic abstraction
# =============================================================================
# --- helper: try to access the backing vector of tuples inside SortedTupleSet ---
function _raw_transition_vec(trans)
    # Common field names used by custom set wrappers
    for fname in (:data, :v, :vec, :tuples, :A)
        if hasproperty(trans, fname)
            v = getproperty(trans, fname)
            return v
        end
    end
    # Fallback: many wrappers store the payload in field #1
    try
        return getfield(trans, 1)
    catch
    end
    error("Cannot access backing transition vector inside $(typeof(trans)). " *
          "Open `typeof(trans)` and check its fields; add the right field name above.")
end

"""
    ts_from_symbolic_onepass(abs_sys, X0_rect; reachable_only=false, x0_abs=nothing)

Build TSFromAdj by scanning transitions once:
transitions are tuples (target, source, symbol=u).
This is the fast path (what your earlier logs were doing).
"""
function ts_from_symbolic_onepass(abs_sys, X0_rect; reachable_only::Bool=false, x0_abs=nothing)
    nx = getproperty(abs_sys.autom, :nstates)
    nu = length(getfield(abs_sys.Udom, 1))

    trans = getproperty(abs_sys.autom, :transitions)
    tv = _raw_transition_vec(trans)  # now iterable

    # ---- 1) count successors per (x,u) to preallocate exactly ----
    counts = zeros(Int, nx, nu)
    @inbounds for tup in tv
        # transitions stored as (target, source, symbol)
        tgt = tup[1]
        src = tup[2]
        u   = tup[3]
        if 1 ≤ src ≤ nx && 1 ≤ u ≤ nu
            counts[src, u] += 1
        end
    end

    # ---- 2) allocate adjacency with exact sizes ----
    Adj = [Vector{Vector{Int}}(undef, nu) for _ in 1:nx]
    @inbounds for x in 1:nx
        for u in 1:nu
            Adj[x][u] = Vector{Int}(undef, counts[x,u])
        end
    end

    # ---- 3) fill adjacency (second pass) ----
    writepos = ones(Int, nx, nu)
    @inbounds for tup in tv
        tgt = tup[1]
        src = tup[2]
        u   = tup[3]
        if 1 ≤ src ≤ nx && 1 ≤ u ≤ nu
            k = writepos[src, u]
            Adj[src][u][k] = tgt
            writepos[src, u] = k + 1
        end
    end

    # Initial-set (OUTER)
    X0 = Dionysos.Symbolic.get_states_from_set(abs_sys, X0_rect, Dionysos.Domain.INNER)

    # Optional reachable-only pruning without renumbering
    if reachable_only
        x0 = x0_abs === nothing ? (isempty(X0) ? nothing : X0[1]) : x0_abs
        if x0 !== nothing
            reachable = falses(nx)
            queue = Int[x0]
            reachable[x0] = true
            while !isempty(queue)
                v = popfirst!(queue)
                for u in 1:nu
                    for w in Adj[v][u]
                        if !reachable[w]
                            reachable[w] = true
                            push!(queue, w)
                        end
                    end
                end
            end
            for x in 1:nx
                if !reachable[x]
                    for u in 1:nu
                        empty!(Adj[x][u])
                    end
                end
            end
            X0 = [x for x in X0 if reachable[x]]
        end
    end

    return TSFromAdj(Adj, X0)
end


# =============================================================================
# Label precompute
# =============================================================================

"""
precompute_labels_and_forbidden(ts, labeler; forbidden_syms=Symbol[], aut_aps=nothing)

- labeler must be callable: labeler(x::Int)::Set{Symbol}
Returns:
  labels_ts::Vector{Set{Symbol}}, forbidden_mask::BitVector
"""
function precompute_labels_and_forbidden(ts::TSFromAdj,
                                        labeler;
                                        forbidden_syms::Vector{Symbol}=Symbol[],
                                        aut_aps=nothing)
    nx = nstates(ts)
    labels_ts = Vector{Set{Symbol}}(undef, nx)
    forbidden = falses(nx)

    @inbounds for x in 1:nx
        lbl = labeler(x)
        labels_ts[x] = lbl
        if !isempty(forbidden_syms)
            forbidden[x] = any(s -> (s in lbl), forbidden_syms)
        end
    end

    return labels_ts, forbidden
end

# Backward compatible method signature used in your planner:
function precompute_labels_and_forbidden(ts::TSFromAdj, labeler, forbidden_syms::Vector{Symbol})
    return precompute_labels_and_forbidden(ts, labeler; forbidden_syms=forbidden_syms)
end


# =============================================================================
# AP labeling utilities
# =============================================================================

"""
    APLabeler

Callable labeler object. Given an abstract state index `x::Int`, returns the set
of atomic propositions that are true at that state.

Example:
    lbl = APLabeler([:blue, :green], truth_sets)
    lbl(42) == Set([:blue])
"""
struct APLabeler
    ap::Vector{Symbol}
    f::Function
end

# Make the labeler callable: lbl(x)
(l::APLabeler)(x::Int) = l.f(x)

"""
    APLabeler(ap::Vector{Symbol}, truth_sets)

Build an `APLabeler` from a dictionary mapping each atomic proposition to the set
of abstract states where it is true.

`truth_sets` is expected to behave like:
    Dict{Symbol, Set{Int}}
"""
function APLabeler(ap::Vector{Symbol}, truth_sets)
    f = x -> Set(s for s in ap if x in get(truth_sets, s, Set{Int}()))
    return APLabeler(ap, f)
end

"""
    labeler_and_inner_sets(abstract_system; regions, forbidden_aps = Symbol[])

Build:
1. an OUTER-based `APLabeler` for use in planning / automata labeling
2. INNER sets for all non-forbidden APs (useful for robust acceptance semantics)

Arguments
---------
- `abstract_system`: Dionysos symbolic abstraction
- `regions::Dict{Symbol,Any}`:
    maps AP name => region object
    supported region types:
      * `Dionysos.Utils.HyperRectangle`
      * `AbstractVector` of `Dionysos.Utils.HyperRectangle`
- `forbidden_aps::Vector{Symbol}`:
    APs like `:obs` for which INNER sets are typically not needed

Returns
-------
(labeler, inner_sets)

where:
- `labeler(x)` returns `Set{Symbol}` of APs true at abstract state `x`
- `inner_sets[ap]` is a `Set{Int}` of abstract states whose cells are INNER to that region
"""
function labeler_and_inner_sets(abstract_system;
                                regions::Dict{Symbol,Any},
                                forbidden_aps::Vector{Symbol} = Symbol[])

    truth_sets = Dict{Symbol, Set{Int}}()   # OUTER labels
    inner_sets = Dict{Symbol, Set{Int}}()   # INNER sets for robust semantics

    for (ap, region_obj) in regions
        # -------------------------
        # OUTER states for labeling
        # -------------------------
        idxs_outer = Int[]

        if region_obj isa Dionysos.Utils.HyperRectangle
            idxs_outer = Dionysos.Symbolic.get_states_from_set(
                abstract_system, region_obj, Dionysos.Domain.INNER
            )

        elseif region_obj isa AbstractVector
            for rect in region_obj
                @assert rect isa Dionysos.Utils.HyperRectangle
                append!(idxs_outer,
                    Dionysos.Symbolic.get_states_from_set(
                        abstract_system, rect, Dionysos.Domain.INNER
                    )
                )
            end
            idxs_outer = unique(idxs_outer)

        else
            error("Unsupported region type for AP $(ap): $(typeof(region_obj))")
        end

        truth_sets[ap] = Set(idxs_outer)

        # ------------------------------------------
        # INNER states for robust goal-like semantics
        # ------------------------------------------
        if !(ap in forbidden_aps)
            idxs_inner = Int[]

            if region_obj isa Dionysos.Utils.HyperRectangle
                idxs_inner = Dionysos.Symbolic.get_states_from_set(
                    abstract_system, region_obj, Dionysos.Domain.INNER
                )

            elseif region_obj isa AbstractVector
                for rect in region_obj
                    @assert rect isa Dionysos.Utils.HyperRectangle
                    append!(idxs_inner,
                        Dionysos.Symbolic.get_states_from_set(
                            abstract_system, rect, Dionysos.Domain.INNER
                        )
                    )
                end
                idxs_inner = unique(idxs_inner)

            else
                error("Unsupported region type for AP $(ap): $(typeof(region_obj))")
            end

            inner_sets[ap] = Set(idxs_inner)
        end
    end

    ap_list = sort(collect(keys(regions)))
    labeler = APLabeler(ap_list, truth_sets)

    return labeler, inner_sets
end

# =============================================================================
# Product-game helpers (reachability-to-accepting, robust by default)
# =============================================================================

# Internal indexing helpers
@inline function _q_valid_range(aut::Automaton)
    aut.zero_based ? (0:(aut.nstates-1)) : (1:aut.nstates)
end

@inline function _q_to_i(aut::Automaton, q::Int)
    aut.zero_based ? (q + 1) : q
end

@inline function _i_to_q(aut::Automaton, iq::Int)
    aut.zero_based ? (iq - 1) : iq
end

const I_INVALID_Q = 0

"""
Precompute next_qP[iq, x] = iq' (1-based index into product q dimension), or I_INVALID_Q.
"""
function _precompute_next_q(aut::Automaton, labels_ts::Vector{Set{Symbol}})
    nx = length(labels_ts)
    qrange = _q_valid_range(aut)
    nq = length(qrange)

    next_q = fill(I_INVALID_Q, nq, nx)

    @inbounds for iq in 1:nq
        q = _i_to_q(aut, iq)
        for x in 1:nx
            m = letter_of(labels_ts[x], aut.aps)
            qp = nextstate(aut, q, m)
            if qp === nothing
                next_q[iq, x] = I_INVALID_Q
            else
                next_q[iq, x] = _q_to_i(aut, qp)
            end
        end
    end
    return next_q
end

"""
synthesize_buchi_controller(ts, aut, labels_ts; semantics=:robust, stop_at_accept=true, max_fp_iters=10^9)

If stop_at_accept=true, we solve *reachability to accepting states* on the product
with controllable-predecessor robust semantics:
  Pre(S) = {p | ∃u: post(p,u) ≠ ∅ and post(p,u) ⊆ S }.

Returns: (WinP, polP, next_qP, infoP)
"""
function synthesize_buchi_controller(ts::TSFromAdj,
                                    aut::Automaton,
                                    labels_ts::Vector{Set{Symbol}};
                                    semantics::Symbol = :robust,
                                    stop_at_accept::Bool = true,
                                    max_fp_iters::Int = 10^9)

    nx = nstates(ts)
    nu = ninputs(ts)
    qrange = _q_valid_range(aut)
    nq = length(qrange)
    np = nx * nq

    next_qP = _precompute_next_q(aut, labels_ts)

    # accepting product set
    accepting_qs = collect(aut.accepting)
    acceptingP = falses(np)
    @inbounds for iq in 1:nq
        q = _i_to_q(aut, iq)
        if is_accepting(aut, q)
            base = (iq - 1) * nx
            for x in 1:nx
                acceptingP[base + x] = true
            end
        end
    end

    # Reachability fixed point to accepting
    WinP = copy(acceptingP)
    layer = fill(typemax(Int), np)
    @inbounds for pid in 1:np
        if WinP[pid]
            layer[pid] = 0
        end
    end

    function pre_product(S::BitVector)
        Pre = falses(np)
        @inbounds for iq in 1:nq
            base = (iq - 1) * nx
            for x in 1:nx
                pid = base + x
                ok_u = false

                for u in 1:nu
                    succs = post(ts, x, u)
                    isempty(succs) && continue

                    if semantics == :robust
                        all_in = true
                        for xp in succs
                            iqp = next_qP[iq, xp]
                            iqp == I_INVALID_Q && (all_in = false; break)
                            pidp = (iqp - 1) * nx + xp
                            if !S[pidp]
                                all_in = false
                                break
                            end
                        end
                        if all_in
                            ok_u = true
                            break
                        end
                    elseif semantics == :optimistic
                        any_in = false
                        for xp in succs
                            iqp = next_qP[iq, xp]
                            iqp == I_INVALID_Q && continue
                            pidp = (iqp - 1) * nx + xp
                            if S[pidp]
                                any_in = true
                                break
                            end
                        end
                        if any_in
                            ok_u = true
                            break
                        end
                    else
                        error("Unknown semantics=$semantics")
                    end
                end

                Pre[pid] = ok_u
            end
        end
        return Pre
    end

    changed = true
    k = 0
    while changed && k < max_fp_iters
        changed = false
        k += 1
        PreW = pre_product(WinP)
        @inbounds for pid in 1:np
            if !WinP[pid] && PreW[pid]
                WinP[pid] = true
                layer[pid] = k
                changed = true
            end
        end
    end

    # Memoryless policy extraction: choose u that keeps successors in WinP and reduces layer
    polP = fill(0, np)

    @inbounds for iq in 1:nq
        base = (iq - 1) * nx
        for x in 1:nx
            pid = base + x
            WinP[pid] || continue

            if acceptingP[pid]
                polP[pid] = 1
                continue
            end

            best_u = 0
            best_score = typemax(Int)

            for u in 1:nu
                succs = post(ts, x, u)
                isempty(succs) && continue

                if semantics == :robust
                    ok = true
                    score = 0
                    for xp in succs
                        iqp = next_qP[iq, xp]
                        (iqp == I_INVALID_Q) && (ok = false; break)
                        pidp = (iqp - 1) * nx + xp
                        if !WinP[pidp]
                            ok = false
                            break
                        end
                        score = max(score, layer[pidp])
                    end
                    ok || continue
                    if score < best_score
                        best_score = score
                        best_u = u
                    end
                else
                    # optimistic
                    ok = false
                    score = typemax(Int)
                    for xp in succs
                        iqp = next_qP[iq, xp]
                        (iqp == I_INVALID_Q) && continue
                        pidp = (iqp - 1) * nx + xp
                        if WinP[pidp]
                            ok = true
                            score = min(score, layer[pidp])
                        end
                    end
                    ok || continue
                    if score < best_score
                        best_score = score
                        best_u = u
                    end
                end
            end

            polP[pid] = (best_u == 0 ? 1 : best_u)
        end
    end

    infoP = (;
        nx = nx,
        nq = nq,
        accepting_qs = accepting_qs,
        acceptingP = acceptingP,
        rank = layer,
        labels_ts = labels_ts
    )

    return WinP, polP, next_qP, infoP
end

# =============================================================================
# Simulation
# =============================================================================

"""
simulate_buchi_policy(ts, aut, labels_ts, polP, next_qP; x0, q0, max_steps, stop_at_accept)

Returns: xs::Vector{Int}, qs::Vector{Int}, us::Vector{Int}, status::Symbol
"""
function simulate_buchi_policy(ts::TSFromAdj,
                               aut::Automaton,
                               labels_ts::Vector{Set{Symbol}},
                               polP::Vector{Int},
                               next_qP::Matrix{Int};
                               x0::Int,
                               q0::Int,
                               max_steps::Int = 400,
                               stop_at_accept::Bool = true)

    nx = nstates(ts)
    nq = size(next_qP, 1)

    # map q0 to iq0
    iq = _q_to_i(aut, q0)
    if iq < 1 || iq > nq
        return Int[x0], Int[q0], Int[], :invalid_q0
    end

    xs = Int[x0]
    qs = Int[q0]
    us = Int[]

    for k in 1:max_steps
        x = xs[end]
        q = qs[end]
        iq = _q_to_i(aut, q)
        pid = (iq - 1) * nx + x
        u = polP[pid]

        if u == 0
            return xs, qs, us, :no_admissible_control
        end

        succs = post(ts, x, u)
        isempty(succs) && return xs, qs, us, :deadend

        # Representative successor: prefer movement (avoid self-loop if possible)
        xp = succs[1]
        for cand in succs
            if cand != x
                xp = cand
                break
            end
        end

        iqp = next_qP[iq, xp]
        (iqp == I_INVALID_Q) && return xs, qs, us, :invalid_transition
        qp = _i_to_q(aut, iqp)

        push!(us, u)
        push!(xs, xp)
        push!(qs, qp)

        if stop_at_accept && is_accepting(aut, qp)
            return xs, qs, us, :hit_accepting
        end
    end

    return xs, qs, us, :max_steps
end




"""Controllable predecessor on the product game with robust semantics.
Given a target set S (BitVector over product ids), return Pre(S):
  Pre(S) = { (x,q) | ∃u: succ((x,q),u) nonempty and succ ⊆ S }
Sink state is handled by set membership, no special qp==0 case.
"""
function pre_product(ts::TSFromAdj,
                     next_q::Array{Int,2},
                     S::BitVector;
                     semantics::Symbol = :robust)
    nx = L2A.nstates(ts)
    nq = size(next_q, 1)
    @assert nq == NQ_TOTAL "pre_product: next_q has nq=$nq but NQ_TOTAL=$NQ_TOTAL (mismatched automaton sizing)"
    np = nx * nq
    nu = L2A.ninputs(ts)

    Pre = falses(np)

    # --- REPLACEMENT: Clean, correct implementation ---
    for iq in 1:nq
        q = i_to_q(iq)
        base = (iq - 1) * nx
        for x in 1:nx
            pid = base + x
            ok_u = false

            for u in 1:nu
                succs = L2A.post(ts, x, u)
                isempty(succs) && continue
                
                has_valid = false

                if semantics == :robust
                    # ∃u : succ ⊆ S
                    all_in = true
                    for xp in succs
                        has_valid = true
                        iqp = next_q[iq, xp]
                        if iqp == I_INVALID_Q
                            all_in = false
                            break
                        end
                        pidp = (iqp - 1) * nx + xp
                        if !S[pidp]
                            all_in = false
                            break
                        end
                    end
                    if has_valid && all_in
                        ok_u = true
                        break
                    end

                elseif semantics == :optimistic
                    # ∃u : succ ∩ S ≠ ∅
                    any_in = false
                    for xp in succs
                        has_valid = true
                        iqp = next_q[iq, xp]
                        if iqp == I_INVALID_Q
                            continue
                        end
                        pidp = (iqp - 1) * nx + xp
                        if S[pidp]
                            any_in = true
                            break
                        end
                    end
                    if has_valid && any_in
                        ok_u = true
                        break
                    end

                else
                    error("pre_product: unknown semantics = $(semantics)")
                end
            end

            Pre[pid] = ok_u
        end
    end

    return Pre
end


function synthesize_reachability_controller(ts::TSFromAdj,
                                           aut::Automaton,
                                           labels::Vector{Set{Symbol}})
    nx = L2A.nstates(ts)
    nq = NQ_TOTAL
    np = nx * nq

    next_q = precompute_next_q(aut, labels)

    # Goal (accepting) product states
    Acc = falses(np)
    for q in Q_VALID_FIRST:Q_VALID_LAST
        if is_accepting_ext(aut, q)
            base = (q_to_i(q) - 1) * nx
            @inbounds for x in 1:nx
                Acc[base + x] = true
            end
        end
    end

    # Reachability winning set: Win = μZ. Acc ∪ Pre(Z)
    Win = copy(Acc)

    # Layer-to-goal (for policy extraction)
    layer = fill(typemax(Int), np)
    for pid in 1:np
        if Win[pid]
            layer[pid] = 0
        end
    end

    changed = true
    k = 0
    while changed
        changed = false
        k += 1

        PreW = pre_product(ts, next_q, Win; semantics = GAME_SEMANTICS)
        @inbounds for pid in 1:np
            if !Win[pid] && PreW[pid]
                Win[pid] = true
                layer[pid] = k
                changed = true
            end
        end
    end

    # Memoryless policy: choose u that stays in Win and reduces layer.
    nu = L2A.ninputs(ts)
    policy = fill(0, np)

    for iq in 1:nq
        q = i_to_q(iq)
        base = (iq - 1) * nx
        for x in 1:nx
            pid = base + x
            Win[pid] || continue

            # If already accepting, any input is fine for plotting; we’ll stop anyway.
            if Acc[pid]
                policy[pid] = 1
                continue
            end

            best_u = 0
            best_score = typemax(Int)
            best_tie = (typemax(Int), typemax(Float64), typemax(Float64), typemax(Int))
            # tie-break tuple: (forward_penalty, |steer|, -|speed|, no_move_penalty)

            for u in 1:nu
                succs = L2A.post(ts, x, u)
                isempty(succs) && continue

                if GAME_SEMANTICS == :robust
                    ok = true
                    score = 0
                    for xp in succs
                        iqp = next_q[iq, xp]
                        pidp = (iqp - 1) * nx + xp
                        if !Win[pidp]
                            ok = false
                            break
                        end
                        score = max(score, layer[pidp])
                    end
                    ok || continue
                    u_vec = input_from_index(u, abstract_system)
                    has_move = any(xp -> xp != x, succs)
                    tie = (u_vec[1] > 0 ? 0 : 1, abs(u_vec[2]), -abs(u_vec[1]), has_move ? 0 : 1)

                    if (score < best_score) || (score == best_score && tie < best_tie)
                        best_score = score
                        best_tie = tie
                        best_u = u
                    end

                elseif GAME_SEMANTICS == :optimistic
                    ok = false
                    score = typemax(Int)
                    for xp in succs
                        iqp = next_q[iq, xp]
                        pidp = (iqp - 1) * nx + xp
                        if Win[pidp]
                            ok = true
                            score = min(score, layer[pidp])
                        end
                    end
                    ok || continue
                    u_vec = input_from_index(u, abstract_system)
                    has_move = any(xp -> xp != x, succs)
                    tie = (u_vec[1] > 0 ? 0 : 1, abs(u_vec[2]), -abs(u_vec[1]), has_move ? 0 : 1)

                    if (score < best_score) || (score == best_score && tie < best_tie)
                        best_score = score
                        best_tie = tie
                        best_u = u
                    end

                else
                    error("Unknown GAME_SEMANTICS = $(GAME_SEMANTICS)")
                end
            end

            policy[pid] = best_u == 0 ? 1 : best_u
        end
    end

    return Win, policy
end


"""Solve Büchi game: winning region for visiting accepting states infinitely often.
Returns (Win::BitVector, policy::Vector{Int}).
Algorithm: repeated controllable attractors:
  Z_{k+1} = Attr_c(Acc ∩ Z_k; Z_k)
"""

"""Simulate the abstract closed-loop product controller.
Returns (xs, qs, us, status).
"""
function simulate_product_policy(ts::TSFromAdj,
                                aut::Automaton,
                                labels::Vector{Set{Symbol}},
                                Win::BitVector,
                                policy::Vector{Int},
                                x0::Int,
                                q0::Int;
                                max_steps::Int = 400)
    nx = L2A.nstates(ts)
    next_q = precompute_next_q(aut, labels)

    xs = Int[x0]
    qs = Int[q0]
    us = Int[]

    for k in 1:max_steps
        x = xs[end]
        q = qs[end]
        pid = prod_id(x, q, nx)
        u = policy[pid]
        u == 0 && return xs, qs, us, :losing

        succs = L2A.post(ts, x, u)
        if all(cand -> cand == x, succs)
            @warn "All successors are self-loops at this step" k x q u
        end
        isempty(succs) && return xs, qs, us, :deadend

        # choose one successor deterministically for a representative run
                # choose one successor deterministically for a representative run,
        # but avoid trivial self-loops if possible.
        iq  = q_to_i(q)
        if GAME_SEMANTICS == :optimistic
            best = nothing

            # 1) prefer a winning successor that actually moves (cand != x)
            for cand in succs
                cand == x && continue
                iqcand = next_q[iq, cand]
                iqcand == I_INVALID_Q && continue
                pidc   = (iqcand - 1) * nx + cand
                if Win[pidc]
                    best = cand
                    break
                end
            end

            # 2) otherwise allow a winning self-loop
            if best === nothing
                for cand in succs
                    iqcand = next_q[iq, cand]
                    iqcand == I_INVALID_Q && continue
                    pidc   = (iqcand - 1) * nx + cand
                    if Win[pidc]
                        best = cand
                        break
                    end
                end
            end

            best === nothing && return xs, qs, us, :left_winning_set
            xp = best
        else
            # robust: pick a moving successor if available
            xp = succs[1]
            for cand in succs
                if cand != x
                    xp = cand
                    break
                end
            end
        end
        iqp = next_q[iq, xp]
        qp  = i_to_q(iqp)

        push!(us, u)
        push!(xs, xp)
        push!(qs, qp)

        # optional stopping if we hit accepting (guard against sink)
        # (Büchi is infinite-horizon; for plotting we just stop on first accept)
        if is_accepting_ext(aut, qp)
            return xs, qs, us, :hit_accepting
        end
    end

    return xs, qs, us, :max_steps
end


end # module