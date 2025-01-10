function analyze_non_determinism(a, sys)
    count = Dict{Tuple{Int, Int}, Int}()
    self_loops = Dict{Tuple{Int, Int}, Int}()

    write_histo = false
    write_transi = false
    path = "C:/Users/adrie/OneDrive - UCL/Master 2/mémoire visus/data/"

    for (t, s, sym) ∈ a.transitions.data
        # s(ource) -> t(arget) with input-symbol sym 
        # parce que non déterministe si plusieurs edge partant du même state avec même input
        if (s, sym) ∈ keys(count)
            count[s, sym] += 1
        else
            count[s, sym] = 1
        end
        if t == s
            self_loops[t, sym] = 1
        end
    end

    if write_transi
        # Write transitions to a text file
        open(joinpath(path, "transitions2.txt"), "w") do file
            for (s, sym) ∈ keys(count)
                pos = sys.xint2pos[s] # SY.get_state_by_xpos(abstract_system, pos)
                sl = 0
                if (s, sym) ∈ keys(self_loops)
                    sl = 1
                end
                println(file, "$s $pos $sym $(count[s, sym]) $sl")
            end
        end
    end

    for (t, s, sym) ∈ a.transitions.data
        if (s,sym) ∈ keys(self_loops)
            self_loops[s, sym] = count[s, sym]
        end
    end

    all_counts = collect(values(count))
    all_self_loops = collect(values(self_loops))

    if write_histo
        # Write all_counts and self_loop_values to a text file
        open(joinpath(path,"all_counts.txt"), "w") do file
            # Write all_counts
            for value in all_counts
                println(file, value)
            end
        end
        open(joinpath(path, "all_self_loops.txt"), "w") do file
            # Write self_loop_values
            for value in all_self_loops
                println(file, value)
            end
        end
    end

    b = range(0, stop=maximum(all_counts), length=maximum(all_counts)+1)
    # Plot the data
    p = histogram(all_counts, color = :blue, bins = b, label="total")
    histogram!(all_self_loops, bins=b, label="with self loops")
    display(p)
    return
end

function analyze_self_loops(a)
    count = 0
    for (t, s, _) ∈ a.transitions.data
        if t == s
            count += 1
        end
    end
    return count
end