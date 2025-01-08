function analyze_non_determinism(a)
    count = Dict{Tuple{Int, Int}, Int}()
    self_loops = Dict{Tuple{Int, Int}, Int}()
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

    for (t, s, sym) ∈ a.transitions.data
        if (s,sym) ∈ keys(self_loops)
            self_loops[s, sym] = count[s, sym]
        end
    end

    all_counts = collect(values(count))
    all_self_loops = collect(values(self_loops))

    path = "C:/Users/adrie/OneDrive - UCL/Master 2"

    # Write all_counts and self_loop_values to a text file
    open(joinpath(path,"all_counts2.txt"), "w") do file
        # Write all_counts
        for value in all_counts
            println(file, value)
        end
    end
    open(joinpath(path, "all_self_loops2.txt"), "w") do file
        # Write self_loop_values
        for value in all_self_loops
            println(file, value)
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