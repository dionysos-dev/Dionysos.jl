using Dionysos
using Documenter, Literate

const EXAMPLES_DIR = joinpath(@__DIR__, "src", "examples")
const OUTPUT_DIR   = joinpath(@__DIR__, "src/generated")

const EXAMPLES = readdir(EXAMPLES_DIR)

for example in EXAMPLES
    example_filepath = joinpath(EXAMPLES_DIR, example)
    Literate.markdown(example_filepath, OUTPUT_DIR)
    Literate.notebook(example_filepath, OUTPUT_DIR)
    Literate.script(example_filepath, OUTPUT_DIR)
end

makedocs(
    sitename = "Dionysos",
    # See https://github.com/JuliaDocs/Documenter.jl/issues/868
    format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true"),
    # See https://github.com/jump-dev/JuMP.jl/issues/1576
    strict = true,
    pages = [
        "Index" => "index.md",
        "Examples" => map(EXAMPLES) do jl_file
            # Need `string` as Documenter fails if `name` is a `SubString{String}`.
            name = string(split(jl_file, ".")[1])
            return name => "generated/$name.md"
        end
    ],
)

# deploydocs(
#     repo   = "github.com/dionysos-dev/Dionysos.jl.git",
#     push_preview = true,
# )
