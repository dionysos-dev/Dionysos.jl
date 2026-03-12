# ============================================================
# Safety Control
# ============================================================

function compute_largest_invariant_set(autom::AbstractAutomatonList, safelist;)
    contr_tab = SymbolicControlTable(get_n_state(autom))
    nstates = get_n_state(autom)
    nsymbols = get_n_input(autom)
    pairstable = [false for i in 1:nstates, j in 1:nsymbols]

    _compute_pairstable(pairstable, autom)
    nsymbolslist = sum(pairstable; dims = 2)

    # Remove unsafe states
    safeset = Set(safelist)
    for source in safeset
        if nsymbolslist[source] == 0
            delete!(safeset, source)
        end
    end

    unsafeset = Set(1:nstates)
    setdiff!(unsafeset, safeset)

    for source in unsafeset
        for symbol in 1:nsymbols
            pairstable[source, symbol] = false
        end
    end
    nextunsafeset = Set{Int}()

    # Iterate until convergence
    while true
        for target in unsafeset
            for soursymb in pre(autom, target)
                if pairstable[soursymb[1], soursymb[2]]
                    pairstable[soursymb[1], soursymb[2]] = false
                    nsymbolslist[soursymb[1]] -= 1
                    if nsymbolslist[soursymb[1]] == 0
                        push!(nextunsafeset, soursymb[1])
                    end
                end
            end
        end

        if isempty(nextunsafeset)
            break
        end

        setdiff!(safeset, nextunsafeset)
        unsafeset, nextunsafeset = nextunsafeset, unsafeset
        empty!(nextunsafeset)
    end

    # Populate controller
    for source in safeset
        for symbol in 1:nsymbols
            if pairstable[source, symbol]
                set_control!(contr_tab, source, symbol)
            end
        end
    end
    unsafeset = setdiff(Set(safelist), safeset)
    abstract_controller = to_ms_controller(contr_tab)
    return abstract_controller, safeset, unsafeset
end

function _compute_pairstable(pairstable, autom)
    for target in enum_states(autom)
        for soursymb in pre(autom, target)
            pairstable[soursymb[1], soursymb[2]] = true
        end
    end
end