import Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()

using PackageCompiler

# Adjust this list to what you really want baked into the sysimage.
# Only include packages that are stable and frequently used.
packages = [
    "Dionysos",
    "MathematicalSystems",
    "StaticArrays",
    "LinearAlgebra",
    "JuMP",
    "MathOptInterface",
]

# Optional: compile representative calls too
precompile_file = joinpath(@__DIR__, "precompile_execution.jl")

create_sysimage(
    packages;
    sysimage_path = joinpath(@__DIR__, "dionysos_sysimage.dll"),   # Windows
    project = @__DIR__,
    # precompile_execution_file = isfile(precompile_file) ? precompile_file : nothing,
    incremental = true,
)