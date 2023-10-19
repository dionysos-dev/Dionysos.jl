using Dionysos
using Documenter, Literate

const EXAMPLES_DIR = joinpath(@__DIR__, "src", "examples")
const EXAMPLES_SOLVERS_DIR = joinpath(@__DIR__, "src", "examples", "solvers")
const EXAMPLES_UTILS_DIR = joinpath(@__DIR__, "src", "examples", "utils")

const REFERENCE_DIR = joinpath(@__DIR__, "src", "reference")
const OUTPUT_DIR = joinpath(@__DIR__, "src", "generated")

const EXAMPLES = readdir(EXAMPLES_DIR)
const EXAMPLES_SOLVERS = readdir(EXAMPLES_SOLVERS_DIR)
const EXAMPLES_UTILS = readdir(EXAMPLES_UTILS_DIR)
const REFERENCE = readdir(REFERENCE_DIR)

for example in EXAMPLES_SOLVERS
    example_filepath = joinpath(EXAMPLES_SOLVERS_DIR, example)
    Literate.markdown(example_filepath, OUTPUT_DIR)
    Literate.notebook(example_filepath, OUTPUT_DIR)
    Literate.script(example_filepath, OUTPUT_DIR)
end
for example in EXAMPLES_UTILS
    example_filepath = joinpath(EXAMPLES_UTILS_DIR, example)
    Literate.markdown(example_filepath, OUTPUT_DIR)
    Literate.notebook(example_filepath, OUTPUT_DIR)
    Literate.script(example_filepath, OUTPUT_DIR)
end

const _PAGES = [
    "Index" => "index.md",
    "Manual" => ["manual/abstraction-based-control.md", "manual/manual.md"],
    "Getting Started" => "generated/Getting Started.md",
    "Solvers" => map(EXAMPLES_SOLVERS) do jl_file
        # Need `string` as Documenter fails if `name` is a `SubString{String}`.
        name = string(split(jl_file, ".")[1])
        return name => "generated/$name.md"
    end,
    "Utils" => map(EXAMPLES_UTILS) do jl_file
        # Need `string` as Documenter fails if `name` is a `SubString{String}`.
        name = string(split(jl_file, ".")[1])
        return name => "generated/$name.md"
    end,
    "API Reference" => map(REFERENCE) do jl_file
        # Need `string` as Documenter fails if `name` is a `SubString{String}`.
        name = string(split(jl_file, ".")[1])
        return name => "reference/$name.md"
    end,
    "Developer Docs" =>
        ["Set up" => "developers/setup.md", "Git" => "developers/git.md"],
]

makedocs(;
    sitename = "Dionysos",
    # See https://github.com/JuliaDocs/Documenter.jl/issues/868
    format = Documenter.HTML(;
        prettyurls = get(ENV, "CI", nothing) == "true",
        assets = ["assets/extra_styles.css"],
    ),
    pages = _PAGES,
    # The following ensures that we only include the docstrings from
    # this module for functions defined in Base that we overwrite.
    # It also errors in case we don't include a docstring in the docs
    modules = [Dionysos],
)

deploydocs(; repo = "github.com/dionysos-dev/Dionysos.jl.git", push_preview = true)
