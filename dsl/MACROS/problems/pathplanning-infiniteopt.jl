# See https://infiniteopt.github.io/InfiniteOpt.jl/stable/examples/Optimal%20Control/hovercraft/
using InfiniteOpt
model = InfiniteModel()
@infinite_parameter(model, t in [0, 60])
@variable(model, x[1:3], Infinite(t))
@variable(model, u[1:2], Infinite(t))
@objective(model, Min, integral((x[1] - 3.3)^2 + (x[2] - 0.5)^2, t))
@expression(model, α, atan(tan(u[2]) / 2))
@constraint(model, ∂(x[1], t) == u[1] * cos(α + x[3]) * sec(α))
@constraint(model, ∂(x[2], t) == u[1] * sin(α + x[3]) * sec(α))
@constraint(model, ∂(x[3], t) == u[1] * tan(u[2]))

@variable(model, mode[1:2], Infinite(t), Int)

@enum(VariableType, INPUT, STATE, MODE)

import OrderedCollections
variables = OrderedCollections.OrderedDict{GeneralVariableRef,VariableType}()

for var in all_variables(model)
    variables[var] = INPUT
end

function parse_constraint!(_, c, _, _)
    error("Constraint $c not supported")
end

function parse_constraint!(variables, _, v::GeneralVariableRef, ::MOI.Integer)
    if variables[v] != INPUT
        error("Cannot specify $v to be integer")
    end
    variables[v] = MODE
end

function _underiv(d::GeneralVariableRef)
    return GeneralVariableRef(d.model, d.raw_index, InfiniteOpt.InfiniteVariableIndex, d.param_index)
end

function parse_constraint!(variables, _, f::NLPExpr, ::MOI.EqualTo)
    if f.tree_root.data == NodeData(:-) &&
        f.tree_root.child.data.value.index_type == InfiniteOpt.DerivativeIndex
        v = _underiv(f.tree_root.child.data.value)
        variables[v] = STATE
    else
        error("Constraint $c not understood")
    end
end

for con_ref in all_constraints(model)
    con_obj = constraint_object(con_ref)
    parse_constraint!(variables, con_ref, con_obj.func, con_obj.set)
end
