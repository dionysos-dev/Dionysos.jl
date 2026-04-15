module DionysosMathOptSymbolicADExt

import Dionysos
import Dionysos: Optimizer

import Symbolics
import MathOptSymbolicAD

const DI = Dionysos
const UT = DI.Utils
const ST = DI.System
const PR = DI.Problem
const OP = DI.Optim

include("../src/MOI_wrapper.jl")

Optimizer() = SymbolicsOptimizer()

end
