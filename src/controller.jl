abstract type Controller end

mutable struct ControllerList
    pairs::Vector{Tuple{Int,Int}}
    issorted::Bool
end

function NewControllerList()
    return ControllerList(Tuple{Int, Int}[], true)
end

function ensure_sorted!(contr::ControllerList)
    if !contr.issorted
        sort!(contr.pairs)
        contr.issorted = true
    end
end

# Assumes not add twice same pair...
function add_pair!(contr::ControllerList, source, symbol)
    push!(contr.pairs, (source, symbol))
    contr.issorted = false
end

function Base.empty!(contr::ControllerList)
    empty!(contr.pairs)
    contr.issorted = true
end

function get_npairs(contr::ControllerList)
    return length(contr.pairs)
end

function compute_enabled_symbols!(symbollist, contr::ControllerList, source)
    ensure_sorted!(contr)
    idxlist = searchsorted(contr.pairs, source, by = x -> x[1])
    for idx in idxlist
        push!(symbollist, contr.pairs[idx][2])
    end
end

function _compute_num_targets_unreachable(num_targets_unreachable, autom)
    soursymblist = Tuple{Int,Int}[]
    for target in 1:autom.nstates
        empty!(soursymblist)
        compute_pre!(soursymblist, autom, target)
        for soursymb in soursymblist
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
            for (source, symbol) in pre(autom, target)
                if !(source in target_set) && iszero(num_targets_unreachable[source, symbol] -= 1)
                    push!(target_set, source)
                    push!(next_targets, source)
                    add_pair!(contr, source, symbol)
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
    soursymblist = Tuple{Int, Int}[]
    for target in 1:autom.nstates
        empty!(soursymblist)
        compute_pre!(soursymblist, autom, target)
        for soursymb in soursymblist
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
    soursymblist = Tuple{Int, Int}[]

    # prog = ProgressUnknown("# iterations computing controller:")
    while true
        # ProgressMeter.next!(prog)
        for target in unsafeset
            empty!(soursymblist)
            compute_pre!(soursymblist, autom, target)
            for soursymb in soursymblist
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
                add_pair!(contr, source, symbol)
            end
        end
    end

    if ⊆(initlist, safeset)
        println("\ncompute_controller_safe! terminated with success")
    else
        println("\ncompute_controller_safe! terminated without covering init set")
    end
end

# # "Set" assumes symmodel_contr is empty initially... May fail if not respected
# function set_controller_reach!(symmodel_contr, symmodel_sys, X_init, X_target)
#   @assert X_init.grid_space === X_target.grid_space ===
#       symmodel_contr.Xgrid === symmodel_contr.Ygrid
#       symmodel_sys.Xgrid === symmodel_sys.Ygrid
#   println("set_controller_reach! started")
#   # Commented sizehints because seems not to improve performances
#   # sizehint_symmodel!(symmodel_contr, get_symmodel_size(symmodel_sys))
#   Xgrid = symmodel_sys.Xgrid
#   Ugrid = symmodel_sys.Ugrid
#   X_remain = NewSubSet(Xgrid)
#   for x_ref in enumerategridspace_ref(Xgrid)
#       if is_xref_controllable(symmodel_sys, x_ref)
#           add_set_by_new_ref!(X_remain, x_ref)
#       end
#   end
#   # size_remain = get_ncells(X_remain)
#   X_init2 = NewSubSet(Xgrid)
#   unionsubsets!(X_init2, X_init)
#   X_target2 = NewSubSet(Xgrid)
#   # sizehintsubset!(X_target2, size_remain)
#   unionsubsets!(X_target2, X_target)
#   setdiffsubsets!(X_remain, X_target)
#   setdiffsubsets!(X_init2, X_target)
#   xref_new_coll = getgridspace_reftype(symmodel_sys.Xgrid)[]
#   uref_enabled_coll = getgridspace_reftype(symmodel_sys.Ugrid)[]
#
#   prog = ProgressUnknown("# iterations computing controller:")
#   while !issubset_empty(X_init2)
#       ProgressMeter.next!(prog)
#       empty!(xref_new_coll)
#       for x_ref in enum_pos(X_remain)
#           empty!(uref_enabled_coll)
#           add_inputs_by_xref_ysub!(uref_enabled_coll, symmodel_sys, x_ref, X_target2)
#           if isempty(uref_enabled_coll)
#               continue
#           end
#           push!(xref_new_coll, x_ref)
#           add_to_symmodel_by_new_refs_coll!(symmodel_contr,
#               (x_ref, u_ref, x_ref) for u_ref in uref_enabled_coll)
#       end
#       if isempty(xref_new_coll)
#           println("\nset_controller_reach! terminated without covering init set")
#           return
#       end
#       add_set_by_new_ref_coll!(X_target2, xref_new_coll)
#       remove_fromsubset_by_ref_coll!(X_remain, xref_new_coll)
#       remove_fromsubset_by_ref_coll!(X_init2, xref_new_coll)
#   end
#   ProgressMeter.finish!(prog)
#   println("\nset_controller_reach! terminated with success")
# end

# # "Set" assumes symmodel_contr is empty initially... May fail if not respected
# function set_controller_safe!(symmodel_contr, symmodel_sys, X_init, X_safe)
#   @assert X_init.grid_space === X_safe.grid_space ===
#       symmodel_contr.Xgrid === symmodel_contr.Ygrid
#       symmodel_sys.Xgrid === symmodel_sys.Ygrid
#   println("set_controller_reach! started")
#   # Commented sizehints because seems not to improve performances
#   # sizehint_symmodel!(symmodel_contr, get_symmodel_size(symmodel_sys))
#   Xgrid = symmodel_sys.Xgrid
#   Ugrid = symmodel_sys.Ugrid
#   X_safe2 = NewSubSet(Xgrid)
#   for x_ref in enum_pos(X_safe)
#       if is_xref_controllable(symmodel_sys, x_ref)
#           add_set_by_new_ref!(X_safe2, x_ref)
#       end
#   end
#   xref_remove_coll = getgridspace_reftype(symmodel_sys.Xgrid)[]
#   uref_enabled_coll = getgridspace_reftype(symmodel_sys.Ugrid)[]
#
#   while true
#       print(".")
#       display(get_ncells(X_safe2))
#       empty!(xref_remove_coll)
#       for x_ref in enum_pos(X_safe2)
#           empty!(uref_enabled_coll)
#           add_inputs_by_xref_ysub!(uref_enabled_coll, symmodel_sys, x_ref, X_safe2)
#           if isempty(uref_enabled_coll)
#               push!(xref_remove_coll, x_ref)
#           end
#       end
#       if isempty(xref_remove_coll)
#           break
#       end
#       remove_fromsubset_by_ref_coll!(X_safe2, xref_remove_coll)
#   end
#
#   for x_ref in enum_pos(X_safe2)
#       empty!(uref_enabled_coll)
#       add_inputs_by_xref_ysub!(uref_enabled_coll, symmodel_sys, x_ref, X_safe2)
#       add_to_symmodel_by_new_refs_coll!(symmodel_contr,
#           (x_ref, u_ref, x_ref) for u_ref in uref_enabled_coll)
#   end
#
#   if issubset1_insubset2(X_init, X_safe2)
#       println("\nset_controller_safe! terminated with success")
#   else
#       println("\nset_controller_safe! terminated without covering init set")
#   end
# end
