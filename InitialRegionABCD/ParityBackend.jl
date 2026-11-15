module ParityBackend
using JSON3, SHA

export BuchiUnavailable, ParityController, translate_deterministic_buchi,
    export_parity_plant, synthesize_parity, verify_parity_controller, simulate_parity_controller

struct BuchiUnavailable <: Exception
    message::String
end
Base.showerror(io::IO,e::BuchiUnavailable)=print(io,e.message)

"""Try a deterministic, state-based Büchi translation; reject only nondeterminism.

The -D preference can obtain a DBA even when a default Büchi translation is
nondeterministic. Spot does not promise every LTL language has a DBA.
"""
function translate_deterministic_buchi(formula::AbstractString,L2A::Module)
    return mktempdir() do directory
        path=joinpath(directory,"buchi.hoa")
        open(path,"w") do io
            run(pipeline(`ltl2tgba -D -B -S --hoa -f $formula`,stdout=io))
        end
        check=run(pipeline(ignorestatus(`autfilt --is-deterministic --count $path`),stdout=devnull))
        check.exitcode in (0,1) || error("Spot determinism check failed")
        check.exitcode==1 && throw(BuchiUnavailable("Spot did not produce a deterministic Büchi automaton"))
        hoa=read(path,String)
        occursin(r"Acceptance:\s*1\s+Inf\(0\)",hoa) || error("Unexpected acceptance in state-based Büchi output")
        return L2A.hoa_to_automaton(hoa)
    end
end

function write_u32(io,x)
    write(io,htol(UInt32(x)))
end
function read_u32(io)
    Int(ltoh(read(io,UInt32)))
end
function label_mask(label,aps)
    mask=UInt64(0)
    for ap in label
        i=findfirst(==(String(ap)),aps)
        isnothing(i) && error("Plant label absent from frozen alphabet: $ap")
        mask |= UInt64(1)<<(i-1)
    end
    return mask
end

"""Export an exact action-labelled nondeterministic plant to the native bridge."""
function export_parity_plant(ts,labels,path,L2A::Module;
    allowed_aps=["blue","brown","green","purple","yellow","obs"])
    aps=sort!(unique(vcat(String.(allowed_aps),[String(a) for l in labels for a in l])))
    length(aps)<=63 || error("The plant binary supports at most 63 labelled propositions")
    length(labels)==L2A.nstates(ts) || error("Plant labels/state count mismatch")
    mkpath(dirname(abspath(path)))
    open(path*".tmp","w") do io
        write(io,"ABCDTS01")
        for n in (L2A.nstates(ts),L2A.ninputs(ts),length(L2A.initset(ts)),length(aps))
            write_u32(io,n)
        end
        for ap in aps
            write_u32(io,ncodeunits(ap));write(io,ap)
        end
        for x in L2A.initset(ts);write_u32(io,x-1);end
        for label in labels;write(io,htol(label_mask(label,aps)));end
        for x in 1:L2A.nstates(ts),u in 1:L2A.ninputs(ts)
            successors=L2A.post(ts,x,u)
            write_u32(io,length(successors))
            for xp in successors;write_u32(io,xp-1);end
        end
    end
    mv(path*".tmp",path;force=true)
    return aps
end

struct ParityController
    policy::Dict{Tuple{Int,Int},Int} # (one-based plant state, zero-based DPA memory) -> input
    delta::Dict{Tuple{Int,UInt64},Tuple{Int,Int}} # (memory,label) -> (next memory,max-odd priority)
    starts::Vector{Tuple{Int,Int}} # memory has already consumed initial label
    aps::Vector{String}
    automaton_states::Int
    priorities::Int
end

function read_policy(path,aps)
    open(path,"r") do io
        String(read(io,8))=="ABCDPG01" || error("Invalid parity policy binary")
        nq,colors,nm,nc=read_u32(io),read_u32(io),read_u32(io),read_u32(io)
        policy=Dict{Tuple{Int,Int},Int}()
        for _ in 1:nc
            x,q,u=read_u32(io)+1,read_u32(io),read_u32(io)
            u>0 && (policy[(x,q)]=u)
        end
        masks=[ltoh(read(io,UInt64)) for _ in 1:nm]
        delta=Dict{Tuple{Int,UInt64},Tuple{Int,Int}}()
        for q in 0:(nq-1),mask in masks
            delta[(q,mask)]=(read_u32(io),read_u32(io))
        end
        starts=[(read_u32(io)+1,read_u32(io)) for _ in 1:read_u32(io)]
        eof(io) || error("Trailing parity policy data")
        return ParityController(policy,delta,starts,aps,nq,colors)
    end
end

"""Synthesize over the full plant product, not free controllable AP valuations."""
function synthesize_parity(formula,ts,labels,L2A::Module;
    bridge,plant_binary,aps,directory,max_steps=1200,max_nodes=2_000_000,max_edges=40_000_000)
    max_steps>=1 || throw(ArgumentError("max_steps must be positive"))
    mkpath(directory)
    hoa=joinpath(directory,"specification.hoa")
    open(hoa,"w") do io
        run(pipeline(`ltl2tgba -D -C --parity=max\ odd -S --hoa -f $formula`,stdout=io))
    end
    result=joinpath(directory,"spot-result.json")
    policy=joinpath(directory,"spot-policy.bin")
    errors=IOBuffer()
    process=run(pipeline(ignorestatus(`$bridge $plant_binary $hoa $result $policy $max_steps $max_nodes $max_edges`),stderr=errors))
    process.exitcode==0 || error("Spot parity bridge failed: "*strip(String(take!(errors))))
    info=JSON3.read(read(result,String))
    controller=info.winning ? read_policy(policy,aps) : nothing
    return (;controller,info,hoa)
end

"""Independently check every adversarial run of a fixed max-odd parity policy.

For each even priority c, reject any reachable cycle whose maximum is c.
This checks subcycles rather than just a whole SCC's maximum priority.
"""
function verify_parity_controller(ts,labels,controller::ParityController,L2A::Module)
    keys=copy(controller.starts)
    isempty(keys) && return Dict("satisfied"=>false,"error"=>"empty initial set")
    index=Dict(k=>i for (i,k) in enumerate(keys))
    graph=Vector{Vector{Tuple{Int,Int}}}()
    masks=[label_mask(label,controller.aps) for label in labels]
    head=1
    while head<=length(keys)
        x,q=keys[head]
        u=get(controller.policy,(x,q),0)
        1<=u<=L2A.ninputs(ts) || return Dict("satisfied"=>false,"error"=>"undefined input")
        successors=L2A.post(ts,x,u)
        isempty(successors) && return Dict("satisfied"=>false,"error"=>"blocking policy")
        edges=Tuple{Int,Int}[]
        for xp in successors
            transition=get(controller.delta,(q,masks[xp]),nothing)
            isnothing(transition) && return Dict("satisfied"=>false,"error"=>"undefined memory update")
            qp,color=transition
            0<=qp<controller.automaton_states && 0<=color<controller.priorities || error("Invalid parity memory or priority")
            key=(xp,qp)
            if !haskey(index,key);push!(keys,key);index[key]=length(keys);end
            push!(edges,(index[key],color))
        end
        push!(graph,edges);head+=1
    end
    n=length(keys)
    for c in 0:2:(controller.priorities-1)
        visited=falses(n);order=Int[]
        for seed in 1:n
            visited[seed] && continue
            stack=Tuple{Int,Int}[(seed,1)];visited[seed]=true
            while !isempty(stack)
                v,i=last(stack)
                while i<=length(graph[v]) && graph[v][i][2]>c;i+=1;end
                if i>length(graph[v])
                    push!(order,v);pop!(stack)
                else
                    stack[end]=(v,i+1);w=graph[v][i][1]
                    if !visited[w];visited[w]=true;push!(stack,(w,1));end
                end
            end
        end
        reverse=[Int[] for _ in 1:n]
        for v in 1:n,(w,color) in graph[v];color<=c && push!(reverse[w],v);end
        components=zeros(Int,n);component=0
        for seed in Iterators.reverse(order)
            components[seed]!=0 && continue
            component+=1;stack=[seed];components[seed]=component
            while !isempty(stack)
                v=pop!(stack)
                for w in reverse[v]
                    if components[w]==0;components[w]=component;push!(stack,w);end
                end
            end
        end
        for v in 1:n,(w,color) in graph[v]
            if color==c && components[v]==components[w]
                return Dict("satisfied"=>false,"rejecting_cycle"=>true,"rejecting_max_priority"=>c,"reachable_product_states"=>n)
            end
        end
    end
    return Dict("satisfied"=>true,"rejecting_cycle"=>false,"reachable_product_states"=>n,
        "acceptance"=>"parity max odd","verification"=>"universal infinite horizon")
end

function simulate_parity_controller(ts,labels,controller::ParityController,L2A::Module;max_steps=1200)
    max_steps>=1 || throw(ArgumentError("max_steps must be positive"))
    x,q=first(controller.starts)
    xs=[x];qs=[q];us=Int[]
    for _ in 1:max_steps
        u=get(controller.policy,(x,q),0)
        u>0 || error("Undefined parity policy during bounded simulation")
        successors=L2A.post(ts,x,u)
        isempty(successors) && error("Blocking parity policy during bounded simulation")
        i=findfirst(!=(x),successors);xp=isnothing(i) ? first(successors) : successors[i]
        qp,_=controller.delta[(q,label_mask(labels[xp],controller.aps))]
        push!(xs,xp);push!(qs,qp);push!(us,u);x,q=xp,qp
    end
    return (;states=xs,memory=qs,inputs=us,status="horizon_reached",steps=max_steps)
end
end # module
