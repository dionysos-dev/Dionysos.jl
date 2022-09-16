using Dionysos
using Documenter, Literate

const EXAMPLES_DIR = joinpath(@__DIR__, "src", "examples")
const REFERENCE_DIR = joinpath(@__DIR__, "src", "reference")
const OUTPUT_DIR   = joinpath(@__DIR__, "src/generated")

const EXAMPLES = readdir(EXAMPLES_DIR)
const REFERENCE = readdir(REFERENCE_DIR)

for example in EXAMPLES
    example_filepath = joinpath(EXAMPLES_DIR, example)
    Literate.markdown(example_filepath, OUTPUT_DIR)
    Literate.notebook(example_filepath, OUTPUT_DIR)
    Literate.script(example_filepath, OUTPUT_DIR)
end

const _PAGES = [
    "Index" => "index.md",
    "Examples" => map(EXAMPLES) do jl_file
        # Need `string` as Documenter fails if `name` is a `SubString{String}`.
        name = string(split(jl_file, ".")[1])
        return name => "generated/$name.md"
    end,
    "API Reference" => map(REFERENCE) do jl_file
        # Need `string` as Documenter fails if `name` is a `SubString{String}`.
        name = string(split(jl_file, ".")[1])
        return name => "reference/$name.md"
    end,
    "Developer Docs" => [
        "Set up" => "developers/setup.md",
        "Git" => "developers/git.md",
    ],
]

makedocs(
    sitename = "Dionysos",
    # See https://github.com/JuliaDocs/Documenter.jl/issues/868
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
        assets = ["assets/extra_styles.css"]
    ),
    # See https://github.com/jump-dev/JuMP.jl/issues/1576
    strict = true,
    pages = _PAGES,
    # The following ensures that we only include the docstrings from
    # this module for functions defined in Base that we overwrite.
    # It also errors in case we don't include a docstring in the docs
    modules = [Dionysos],
)

deploydocs(
    repo   = "github.com/dionysos-dev/Dionysos.jl.git",
    push_preview = true,
)