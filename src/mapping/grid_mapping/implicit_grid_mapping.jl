"""
ImplicitGridMapping:
- universe is a rectangular index box (min_pos..max_pos) possibly periodic
- id computed by row-major linearization (no dict)
- pos recovered by inverse mapping
"""
struct ImplicitGridMapping{N, T, G} <: GridMapping{N, T}
    grid::G
    min_pos::NTuple{N, Int}
    max_pos::NTuple{N, Int}
    n_per_dim::NTuple{N, Int}
    strides::NTuple{N, Int}
end

function ImplicitGridMapping(
    grid::G,
    min_pos::NTuple{N, Int},
    max_pos::NTuple{N, Int},
) where {N, G}
    n_per_dim = ntuple(d -> max_pos[d] - min_pos[d] + 1, N)
    @assert all(n_per_dim .>= 1) "Invalid implicit box: some max_pos < min_pos"

    # row-major strides: stride[1]=1, stride[2]=n1, stride[3]=n1*n2, ...
    strides = ntuple(d -> begin
        s = 1
        for k in 1:(d - 1)
            s *= n_per_dim[k]
        end
        s
    end, N)

    T = eltype(get_origin(grid))
    return ImplicitGridMapping{N, T, G}(grid, min_pos, max_pos, n_per_dim, strides)
end

function ImplicitGridMapping(
    grid::G,
    rect::UT.HyperRectangle{N, T};
    incl_mode = OUTER,
) where {N, T, G}
    rectI = get_pos_lims(grid, rect, incl_mode)
    min_pos = rectI.lb isa NTuple{N, Int} ? rectI.lb : NTuple{N, Int}(rectI.lb)
    max_pos = rectI.ub isa NTuple{N, Int} ? rectI.ub : NTuple{N, Int}(rectI.ub)
    return ImplicitGridMapping(grid, min_pos, max_pos)
end

get_grid(m::ImplicitGridMapping) = m.grid
get_n_state(m::ImplicitGridMapping) = prod(m.n_per_dim)
enum_states(m::ImplicitGridMapping) = 1:get_n_state(m)

is_valid_pos(m::ImplicitGridMapping{N}, pos::NTuple{N, Int}) where {N} =
    all(d -> (m.min_pos[d] <= pos[d] <= m.max_pos[d]), 1:N)

function get_state_by_pos(m::ImplicitGridMapping{N}, pos::NTuple{N, Int}) where {N}
    is_valid_pos(m, pos) || throw(DomainError(pos, "pos out of implicit mapping"))
    q = 1
    @inbounds for d in 1:N
        q += (pos[d] - m.min_pos[d]) * m.strides[d]
    end
    return q
end

function get_pos_by_state(m::ImplicitGridMapping{N}, q::Int) where {N}
    (1 <= q <= get_n_state(m)) || throw(DomainError(q, "state out of implicit universe"))
    z = q - 1
    return ntuple(d -> begin
        span = m.n_per_dim[d]
        digit = (z ÷ m.strides[d]) % span
        m.min_pos[d] + digit
    end, N)
end
