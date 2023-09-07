include(joinpath(@__DIR__, "benchmarking_definitions.jl"))

function benchmark(optimizer::MOI.AbstractOptimizer)
    @time MOI.optimize!(optimizer)
end

for (i, ((method, problem), optimizer)) in enumerate(bench)
    println("$i) ($method, $problem)")
    benchmark(optimizer)
    println()
end
