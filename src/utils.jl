# To have type-stable Iterators.product
function _make_iterator_from_lims(rect::HyperRectangle{NTuple{1,Int}})
    return Iterators.product(UnitRange(rect.lb[1], rect.ub[1]))
end

function _make_iterator_from_lims(rect::HyperRectangle{NTuple{2,Int}})
    return Iterators.product(UnitRange(rect.lb[1], rect.ub[1]),
        UnitRange(rect.lb[2], rect.ub[2]))
end

function _make_iterator_from_lims(rect::HyperRectangle{NTuple{3,Int}})
    return Iterators.product(UnitRange(rect.lb[1], rect.ub[1]),
        UnitRange(rect.lb[2], rect.ub[2]),
        UnitRange(rect.lb[3], rect.ub[3]))
end

function _make_iterator_from_lims(rect::HyperRectangle{NTuple{4,Int}})
    return Iterators.product(UnitRange(rect.lb[1], rect.ub[1]),
        UnitRange(rect.lb[2], rect.ub[2]),
        UnitRange(rect.lb[3], rect.ub[3]),
        UnitRange(rect.lb[4], rect.ub[4]))
end

function _make_iterator_from_lims(rect::HyperRectangle{NTuple{5,Int}})
    return Iterators.product(UnitRange(rect.lb[1], rect.ub[1]),
        UnitRange(rect.lb[2], rect.ub[2]),
        UnitRange(rect.lb[3], rect.ub[3]),
        UnitRange(rect.lb[4], rect.ub[4]),
        UnitRange(rect.lb[5], rect.ub[5]))
end

function _make_iterator_from_lims(rect::HyperRectangle{NTuple{6,Int}})
    return Iterators.product(UnitRange(rect.lb[1], rect.ub[1]),
        UnitRange(rect.lb[2], rect.ub[2]),
        UnitRange(rect.lb[3], rect.ub[3]),
        UnitRange(rect.lb[4], rect.ub[4]),
        UnitRange(rect.lb[5], rect.ub[5]),
        UnitRange(rect.lb[6], rect.ub[6]))
end
