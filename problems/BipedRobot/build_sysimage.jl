import Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()

using PackageCompiler

packages = ["Dionysos", "RigidBodyDynamics", "MathematicalSystems", "StaticArrays"]

precompile_file = joinpath(@__DIR__, "precompile_execution.jl")

create_sysimage(
    packages;
    sysimage_path = joinpath(@__DIR__, "dionysos_robot_sysimage.dll"),
    project = @__DIR__,
    precompile_execution_file = isfile(precompile_file) ? precompile_file : nothing,
    incremental = true,
)
