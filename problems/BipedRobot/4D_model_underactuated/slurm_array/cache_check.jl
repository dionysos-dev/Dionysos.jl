# Shared by precompile.jl and compute_chunk.jl (a plain `include`): which packages of the experiment have no
# valid precompilation cache in the Julia process that runs it. `Base.isprecompiled` applies the test that
# `using` applies, including the one for the CPU the cache was built for, without loading anything.

# What the experiment loads (examples/BipedRobot/biped_4d_normal_form.jl and the files it includes), and the
# package extensions of Dionysos that this loads.
const CACHE_PACKAGES = [
    "StaticArrays",
    "MathematicalSystems",
    "LazySets",
    "MathOptInterface",
    "Dionysos",
    "Symbolics",
    "IntervalArithmetic",
    "JLD2",
    "Plots",
]
const CACHE_EXTENSIONS = ["DionysosSymbolicsExt", "DionysosPlotsExt"]

function stale_packages()
    stale = String[]
    for name in CACHE_PACKAGES
        id = Base.identify_package(name)
        (id === nothing || !Base.isprecompiled(id)) && push!(stale, name)
    end
    dionysos = Base.identify_package("Dionysos")
    for name in CACHE_EXTENSIONS
        Base.isprecompiled(Base.PkgId(Base.uuid5(dionysos.uuid, name), name)) ||
            push!(stale, name)
    end
    return stale
end
