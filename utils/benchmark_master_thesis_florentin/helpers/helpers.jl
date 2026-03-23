import IntervalArithmetic as IA

if !isdefined(@__MODULE__, :_BENCHMARK_MASTER_THESIS_HELPERS_LOADED)
    const _BENCHMARK_MASTER_THESIS_HELPERS_LOADED = true
    include(joinpath(@__DIR__, "periodic.jl"))
    include(joinpath(@__DIR__, "pipeline.jl"))
    include(joinpath(@__DIR__, "plot.jl"))
end
