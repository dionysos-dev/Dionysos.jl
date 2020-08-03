abstract type Controller end

mutable struct ControllerList
    pairs::Vector{Tuple{Int, Int}}
    issorted::Bool
end

function NewControllerList()
    return ControllerList(Tuple{Int, Int}[], true)
end

function ensure_sorted!(contr::ControllerList)
    if !contr.issorted
        # display("contr not sorted")
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

function _compute_npoststable(npoststable, autom)
	soursymblist = Tuple{Int, Int}[]
    for target in 1:autom.nstates
		empty!(soursymblist)
		compute_pre!(soursymblist, autom, target)
		for soursymb in soursymblist
			npoststable[soursymb[1], soursymb[2]] += 1
		end
	end
end

# Assumes contr is "empty"
function compute_controller_reach!(contr, autom, initlist, targetlist)
	println("compute_controller_reach! started")
	nstates = autom.nstates
	nsymbols = autom.nsymbols
	npoststable = [0 for i = 1:nstates, j = 1:nsymbols]
	_compute_npoststable(npoststable, autom)
	initset = Set(initlist)
	targetlist = Set(targetlist)
	nexttargetlist = Set{Int}()
	soursymblist = Tuple{Int, Int}[]

	prog = ProgressUnknown("# iterations computing controller:")
	while !isempty(initset)
		for source in targetlist
			for symbol = 1:nsymbols
				npoststable[source, symbol] = -1
			end
		end
		ProgressMeter.next!(prog)
		for target in targetlist
			empty!(soursymblist)
			compute_pre!(soursymblist, autom, target)
			for soursymb in soursymblist
				npoststable[soursymb[1], soursymb[2]] -= 1
				if npoststable[soursymb[1], soursymb[2]] == 0
					push!(nexttargetlist, soursymb[1])
					add_pair!(contr, soursymb[1], soursymb[2])
				end
			end
		end
		if isempty(nexttargetlist)
			println("\ncompute_controller_reach! terminated without covering init set")
			return
		end
		setdiff!(initset, nexttargetlist)
		temp = targetlist
		targetlist = nexttargetlist
		nexttargetlist = temp
		empty!(nexttargetlist)
	end

	ProgressMeter.finish!(prog)
	println("\ncompute_controller_reach! terminated with success")
end
