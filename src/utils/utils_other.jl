module Utils
"""
    List of functions that could be added elsewhere
"""

using ..Abstraction
const AB = Abstraction

using ..DomainList
const D = DomainList
using Plots


## add to rectangle.jl file
# function center(rect::AB.HyperRectangle)
#     return (rect.lb+rect.ub)/2
# end
# function is_intersection(a::AB.HyperRectangle, b::AB.HyperRectangle)
#     c = AB.intersect(a, b)
#     return !AB.isempty(c)
# end
# function dims(rect::AB.HyperRectangle)
#     return length(rect.lb)
# end
# function scale(rect::AB.HyperRectangle,α)
#     return AB.HyperRectangle(rect.lb*α,rect.ub*α)
# end
## Plotting
# function rectangle(lb,ub)
#     Shape([(lb[1],lb[2]),(lb[1],ub[2]),(ub[1],ub[2]),(ub[1],lb[2])])
# end
# function rectangle2(c,r)
#     Shape(c[1].-r[1] .+ [0,2*r[1],2*r[1],0], c[2].-r[2] .+ [0,0,2*r[2],2*r[2]])
# end
#
# function plot_domain!(Xdom;dims=[1,2],opacity=0.2,color=:red)
#     grid = Xdom.grid
#     r = grid.h/2
#     dict = Dict{NTuple{2,Int}, Any}()
#     for pos in AB.enum_pos(Xdom)
#         center = AB.get_coord_by_pos(grid, pos)
#         if !haskey(dict,pos[dims])
#             dict[pos[dims]] = true
#             plot!(rectangle2(center[dims],r[dims]), opacity=opacity,color=color)
#         end
#     end
# end

# ## add to symbolicmodel.jl
# function get_symbol(symmodel,subset,incl_mode::AB.INCL_MODE)
#     Xdom = symmodel.Xdom
#     grid = Xdom.grid
#     posL = AB.get_subset_pos(Xdom,subset,incl_mode)
#     symbolsList = [AB.get_state_by_xpos(symmodel, pos) for pos in posL]
#     return symbolsList
# end
#
# function get_symbols(symmodel,subsetList,incl_mode::AB.INCL_MODE)
#     symbols = Int[]
#     for subset in subsetList
#         append!(symbols,get_symbol(symmodel,subset,incl_mode))
#     end
#     return symbols
# end

## add to symbolicmodel.jl

mutable struct SymbolicModelList2# <: AB.SymbolicModel{N,M}
    Xdom
    Udom
    autom
    xint2elem
    elem2int
    uint2elem
end
function SymbolicModelList2(Xdom,Udom)
    xint2elem = [elem for elem in AB.enum_pos(Xdom)]
    elem2int = Dict((elem, i) for (i, elem) in enumerate(AB.enum_pos(Xdom)))
    uint2elem = [elem for elem in AB.enum_pos(Udom)]
    return SymbolicModelList2(Xdom, Udom, nothing, xint2elem, elem2int, uint2elem)
end
## new type of domain struct add to domain.jl
# struct CustomList{N,T} <: AB.Domain{N,T}
#     elems::Vector{T}
# end
#
# function CustomList(elems::Vector{T}) where T
#     return CustomList{length(elems),T}(elems)
# end
#
# function AB.enum_pos(Xdom::CustomList)
#     return Xdom.elems
# end
# function AB.get_ncells(domain::CustomList)
#     return length(domain.elems)
# end
# function plot_domain!(Xdom::CustomList)
#     for rec in Xdom.elems
#         plot!(rectangle(rec.lb,rec.ub), opacity=.2,color=:red)
#     end
# end
end # end module
