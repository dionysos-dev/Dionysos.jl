module Abstraction

# For info, https://github.com/JuliaLang/julia/issues/23269

using StaticArrays
using Base.Cartesian
using CUDD

const INT_MAX = typemax(Int)
const INT_MIN = typemin(Int)
@enum INCL_MODE INNER OUTER

abstract type Manager end
abstract type GridManager <: Manager end

abstract type XType end
abstract type UType end

xcelltype(mng::Manager) = _celltype(typeof(mng.XDisc))
ucelltype(mng::Manager) = _celltype(typeof(mng.UDisc))
get_disc(mng::Manager, ::Type{XType}) = mng.XDisc
get_disc(mng::Manager, ::Type{UType}) = mng.UDisc
get_grid(mng::GridManager, XUType) = get_disc(mng, XUType).grid
xspacetype(mng::Manager) = _spacetype(typeof(mng.XDisc))
uspacetype(mng::Manager) = _spacetype(typeof(mng.UDisc))
get_map(mng::Manager, ::Type{XType}) = mng.XMap
get_map(mng::Manager, ::Type{UType}) = mng.UMap
statetype(mng::Manager) = _symboltype(typeof(mng.XMap))
labeltype(mng::Manager) = _symboltype(typeof(mng.UMap))
get_nstates(mng::Manager) = _nsymbols(mng.XMap)
get_nlabels(mng::Manager) = _nsymbols(mng.UMap)
enum_states(mng::Manager) = enum_symbols(mng, mng.XMap)
enum_labels(mng::Manager) = enum_symbols(mng, mng.UMap)
transitiontype(mng::Manager) = _transitiontype(typeof(mng), typeof(mng.XMap), typeof(mng.UMap))
controltype(mng::Manager) = _controltype(typeof(mng), typeof(mng.XMap), typeof(mng.UMap))

abstract type Discretization end

abstract type Domain end
abstract type XDomain <: Domain end
abstract type UDomain <: Domain end

_XUType(::XDomain) = XType
_XUType(::UDomain) = UType

abstract type Mapping end

abstract type SymbolSet end
abstract type StateSet <: SymbolSet end
abstract type LabelSet <: SymbolSet end

_XUType(::StateSet) = XType
_XUType(::LabelSet) = UType

abstract type Automaton end

abstract type Controller end

abstract type Specification end

mutable struct ListGridManager{XD,UD,XM,UM} <: GridManager
    XDisc::XD
    UDisc::UD
    XMap::XM
    UMap::UM
end

function ListGridManager(
        xorig::XV, xh::XV,
        uorig::UV, uh::UV
        ) where {XV<:SVector,UV<:SVector}
    XDisc = ListGridDisc(xorig, xh)
    UDisc = ListGridDisc(uorig, uh)
    XCT = _celltype(typeof(XDisc))
    UCT = _celltype(typeof(UDisc))
    XMap = ListGridMap(XCT)
    UMap = ListGridMap(UCT)
    return ListGridManager(
        XDisc, UDisc,
        XMap, UMap)
end

mutable struct BDDGridManager{XD,UD,XM,UM,BM} <: GridManager
    XDisc::XD
    UDisc::UD
    XMap::XM
    UMap::UM
    BDD::BM
end

function BDDGridManager(
        xorig::XV, xh::XV, xmin::XV,
        uorig::UV, uh::UV, umin::UV
        ) where {XV<:SVector,UV<:SVector}
    N = length(XV)
    M = length(UV)
    XGrid = TheGrid(xorig, xh)
    UGrid = TheGrid(uorig, uh)
    xposmin = _grid_coord2pos(XGrid, xmin)
    uposmin = _grid_coord2pos(UGrid, umin)
    XDisc = BDDGridDisc(XGrid, xposmin)
    UDisc = BDDGridDisc(UGrid, uposmin)
    XCT = _celltype(typeof(XDisc))
    UCT = _celltype(typeof(UDisc))
    XMap = 0 # ListGridMap(XCT)
    UMap = 0 # ListGridMap(UCT)
    return ListGridManager(
        XDisc, UDisc,
        XMap, UMap)
end

include("BDD.jl")
include("rectangle.jl")
include("indexedcollection.jl")
include("controlsystems.jl")
include("discretization.jl")
include("domain.jl")
include("mapping.jl")
include("symbolic.jl")
include("automaton.jl")
include("controller.jl")

end  # module Abstraction
