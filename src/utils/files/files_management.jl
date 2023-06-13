# returns the index of the k-th occurrence of a substring substring in string p:
function findk_idx(substring::AbstractString, string::AbstractString, k::Int)
    indices = Int[]
    idx = 0
    while length(indices) < k && idx < length(string)
        idx = findnext(substring, string, idx + 1)
        if idx === nothing
            return 0
        end
        idx = first(idx)
        push!(indices, idx)
    end
    return indices[end]
end

function replace_occurrence!(
    file_path::AbstractString,
    old::AbstractString,
    new::AbstractString,
    k,
)
    file_string = read(file_path, String)
    start_idx = findk_idx(old, file_string, k)
    if start_idx === nothing
        return
    end
    start_idx = first(start_idx)
    old_len = length(old)
    new_string =
        string(file_string[1:(start_idx - 1)], new, file_string[(start_idx + old_len):end])
    open(file_path, "w") do file
        return write(file, new_string)
    end
end

# for occurence k
function write_after_string!(
    file::AbstractString,
    str_before::AbstractString,
    str_after::AbstractString,
    k,
)
    return replace_occurrence!(file, str_before, str_before * str_after, k)
end
