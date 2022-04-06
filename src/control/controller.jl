abstract type Controller end

NewControllerList() = UT.SortedTupleSet{2,Int}()

function _compute_num_targets_unreachable(num_targets_unreachable, autom)
    for target in 1:autom.nstates
        for soursymb in SY.pre(autom, target)
            num_targets_unreachable[soursymb[1], soursymb[2]] += 1
        end
    end
end

# Assumes contr is "empty"
function _compute_controller_reach!(contr, autom, init_set, target_set, num_targets_unreachable, current_targets, next_targets)
    num_init_unreachable = length(init_set)
    while !isempty(current_targets) && !iszero(num_init_unreachable)
        empty!(next_targets)
        for target in current_targets
            for (source, symbol) in SY.pre(autom, target)
                if !(source in target_set) && iszero(num_targets_unreachable[source, symbol] -= 1)
                    push!(target_set, source)
                    push!(next_targets, source)
                    UT.push_new!(contr, (source, symbol))
                    if source in init_set
                        num_init_unreachable -= 1
                    end
                end
            end
        end
        current_targets, next_targets = next_targets, current_targets
    end
    return iszero(num_init_unreachable)
end
function _data(contr, autom, initlist, targetlist)
    num_targets_unreachable = zeros(Int, autom.nstates, autom.nsymbols)
    _compute_num_targets_unreachable(num_targets_unreachable, autom)
    initset = BitSet(initlist)
    targetset = BitSet(targetlist)
    current_targets = copy(targetlist)
    next_targets = Int[]
    return initset, targetset, num_targets_unreachable, current_targets, next_targets
end
function compute_controller_reach!(contr, autom, initlist, targetlist::Vector{Int})
    println("compute_controller_reach! started")
    # TODO: try to infer whether num_targets_unreachable is sparse or not,
    # and if sparse, use a dictionary instead
    if !_compute_controller_reach!(contr, autom, _data(contr, autom, initlist, targetlist)...)
        println("\ncompute_controller_reach! terminated without covering init set")
        # ProgressMeter.finish!(prog)
        return
    end
    # ProgressMeter.finish!(prog)
    println("\ncompute_controller_reach! terminated with success")
end

function _compute_pairstable(pairstable, autom)
    for target in 1:autom.nstates
        for soursymb in SY.pre(autom, target)
            pairstable[soursymb[1], soursymb[2]] = true
        end
    end
end

function compute_controller_safe!(contr, autom, initlist, safelist)
    println("compute_controller_safe! started")
    nstates = autom.nstates
    nsymbols = autom.nsymbols
    pairstable = [false for i in 1:nstates, j in 1:nsymbols]
    _compute_pairstable(pairstable, autom)
    nsymbolslist = sum(pairstable, dims = 2)
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

    # prog = ProgressUnknown("# iterations computing controller:")
    while true
        # ProgressMeter.next!(prog)
        for target in unsafeset
            for soursymb in SY.pre(autom, target)
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
    # ProgressMeter.finish!(prog)

    for source in safeset
        for symbol in 1:nsymbols
            if pairstable[source, symbol]
                UT.push_new!(contr, (source, symbol))
            end
        end
    end

    if ⊆(initlist, safeset)
        println("\ncompute_controller_safe! terminated with success")
    else
        println("\ncompute_controller_safe! terminated without covering init set")
    end
end
