# To have type-stable Iterators.product
function _make_iterator_from_lims(lb::NTuple{1, Int}, ub::NTuple{1, Int})
    return Iterators.product(UnitRange(lb[1], ub[1]))
end

function _make_iterator_from_lims(lb::NTuple{2, Int}, ub::NTuple{2, Int})
    return Iterators.product(UnitRange(lb[1], ub[1]),
        UnitRange(lb[2], ub[2]))
end

function _make_iterator_from_lims(lb::NTuple{3, Int}, ub::NTuple{3, Int})
    return Iterators.product(UnitRange(lb[1], ub[1]),
        UnitRange(lb[2], ub[2]),
        UnitRange(lb[3], ub[3]))
end

function _make_iterator_from_lims(lb::NTuple{4, Int}, ub::NTuple{4, Int})
    return Iterators.product(UnitRange(lb[1], ub[1]),
        UnitRange(lb[2], ub[2]),
        UnitRange(lb[3], ub[3]),
        UnitRange(lb[4], ub[4]))
end

function _make_iterator_from_lims(lb::NTuple{5, Int}, ub::NTuple{5, Int})
    return Iterators.product(UnitRange(lb[1], ub[1]),
        UnitRange(lb[2], ub[2]),
        UnitRange(lb[3], ub[3]),
        UnitRange(lb[4], ub[4]),
        UnitRange(lb[5], ub[5]))
end

function _make_iterator_from_lims(lb::NTuple{6, Int}, ub::NTuple{6, Int})
    return Iterators.product(UnitRange(lb[1], ub[1]),
        UnitRange(lb[2], ub[2]),
        UnitRange(lb[3], ub[3]),
        UnitRange(lb[4], ub[4]),
        UnitRange(lb[5], ub[5]),
        UnitRange(lb[6], ub[6]))
end
