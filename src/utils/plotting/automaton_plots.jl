function analyze_non_determinism(a)
    count = Dict{Tuple{Int, Int}, Int}()
    for (_, s, sym) ∈ a.transitions.data
        if (s, sym) ∈ keys(count)
            count[s, sym] += 1
        else
            count[s, sym] = 1
        end
    end
    p = histogram(collect(values(count)), legend = false)
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