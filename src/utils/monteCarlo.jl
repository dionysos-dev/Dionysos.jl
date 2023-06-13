import Distributions

function sample_from_rec(rec, N)
    L = []
    n = length(rec.lb)
    for i in 1:n
        unifdist = Distributions.Uniform(rec.lb[i], rec.ub[i])
        push!(L, Distributions.rand(unifdist, N))
    end
    points = []
    for i in 1:N
        vec = Float64[L[j][i] for j in 1:n]
        push!(points, SVector{n, Float64}(vec))
    end
    return points
end

function count_occurences(tab)
    sort!(tab)
    n = length(tab)
    occurences = Tuple{Int, Int}[]
    val = tab[1]
    count = 1
    for i in 2:n
        if tab[i] == val
            count += 1
        else
            push!(occurences, (val, count))
            val = tab[i]
            count = 1
        end
    end
    push!(occurences, (val, count))
    return occurences
end
