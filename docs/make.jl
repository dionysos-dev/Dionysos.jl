using Dionysos
using Documenter, Literate
import DocumenterCitations

const SRC_DIR = joinpath(@__DIR__, "src")
const EXAMPLES_DIR = joinpath(SRC_DIR, "examples")
const OUTPUT_DIR = joinpath(SRC_DIR, "generated")
const REFERENCE_DIR = joinpath(SRC_DIR, "reference")

# Example tiers, in nav order. `jump/` shows how a problem is *stated* through the JuMP
# front-end; `solvers/` shows how many algorithms accept it through the MOI contract.
# Adding a tier is one entry here; adding an example is just a file in the folder.
const TIERS = ["jump" => "Examples", "solvers" => "Solver families"]

# A tier may open with a hand-written landing page, listed before its generated pages.
const TIER_INDEX = Dict("jump" => "examples.md")

# Ordering within a tier: the basenames listed here come first, in this order; anything
# else is appended alphabetically with a notice. This keeps a new example working with no
# edit to this file, while still letting the reading order be curated.
const ORDER = Dict(
    "jump" => [
        "path_planning",
        "unicycle_robot",
        "dcdc_converter",
        "thermostat",
        "integrator_ltl",
    ],
    "solvers" => [
        "gol_lazar_belta",
        "lazy_ellipsoids_abstraction",
        "trajectory_certification",
        "biped_velocity_control",
        "state_feedback_pwa",
    ],
)

# Executing every Literate example dominates the build time; skip it for fast local
# Documenter iterations with DIONYSOS_SKIP_LITERATE=true. Pages whose markdown is not
# already in `src/generated` are then dropped from the nav rather than breaking the build.
const SKIP_LITERATE = get(ENV, "DIONYSOS_SKIP_LITERATE", "false") == "true"

# Public-API predicate for the `@autodocs` reference blocks: keep everything except
# underscore-prefixed internals. Referenced by name (`Filter = _is_public`) from the
# reference pages, so the API reference stays complete and never breaks when a new
# public docstring is added (`checkdocs = :all` enforces coverage).
_is_public(x) =
    let n = try
            string(nameof(x))
        catch
            return true
        end
        !startswith(n, "_")
    end

# The nav entry is the example's own `# # Title: what it achieves` line, cut at the colon.
# That decouples the sidebar text ("DC-DC converter") from the URL (`dcdc_converter.html`),
# so filenames can stay lowercase and space-free without the nav going with them.
function example_title(path)
    for line in eachline(path)
        m = match(r"^#\s+#\s+(.*\S)\s*$", line)
        if m !== nothing
            return String(first(split(String(m.captures[1]), ':')))
        end
    end
    return titlecase(replace(basename(path)[1:(end - 3)], '_' => ' '))
end

function tier_stems(tier)
    dir = joinpath(EXAMPLES_DIR, tier)
    isdir(dir) || return String[]
    stems = sort([f[1:(end - 3)] for f in readdir(dir) if endswith(f, ".jl")])
    listed = get(ORDER, tier, String[])
    known = filter(in(stems), listed)
    rest = sort(setdiff(stems, known))
    isempty(rest) ||
        @info "docs: examples missing from ORDER[\"$tier\"] — appended alphabetically" rest
    stale = setdiff(listed, stems)
    isempty(stale) || @warn "docs: ORDER[\"$tier\"] lists files that do not exist" stale
    return vcat(known, rest)
end

literate_actions(file) =
    (Literate.markdown(file, OUTPUT_DIR); Literate.script(file, OUTPUT_DIR))

const GETTING_STARTED = joinpath(EXAMPLES_DIR, "getting_started.jl")

if !SKIP_LITERATE
    literate_actions(GETTING_STARTED)
    for (tier, _) in TIERS, stem in tier_stems(tier)
        literate_actions(joinpath(EXAMPLES_DIR, tier, stem * ".jl"))
    end
end

# In skip mode only pages already present in `src/generated` can be listed.
has_page(stem) = !SKIP_LITERATE || isfile(joinpath(OUTPUT_DIR, stem * ".md"))

function tier_pages(tier)
    pages = Pair{String, String}[]
    index = get(TIER_INDEX, tier, nothing)
    index === nothing || push!(pages, "Overview" => index)
    for stem in tier_stems(tier)
        has_page(stem) || continue
        title = example_title(joinpath(EXAMPLES_DIR, tier, stem * ".jl"))
        push!(pages, title => "generated/$stem.md")
    end
    return pages
end

const _PAGES = Any[
    "Home" => "index.md",
    "Manual" => ["manual/abstraction-based-control.md", "manual/overview.md"],
]

has_page("getting_started") &&
    push!(_PAGES, "Getting Started" => "generated/getting_started.md")

for (tier, label) in TIERS
    pages = tier_pages(tier)
    isempty(pages) || push!(_PAGES, label => pages)
end

push!(_PAGES, "API Reference" => map(sort(readdir(REFERENCE_DIR))) do file
    # `string` because Documenter fails on a `SubString{String}` name.
    name = string(split(file, ".")[1])
    return name => "reference/$name.md"
end)
push!(
    _PAGES,
    "Developer Docs" => [
        "Set up" => "developers/setup.md",
        "Conventions" => "developers/conventions.md",
        "Adding an example" => "developers/examples.md",
        "Git" => "developers/git.md",
    ],
)
push!(_PAGES, "Bibliography" => "bibliography.md")

makedocs(;
    sitename = "Dionysos",
    # See https://github.com/JuliaDocs/Documenter.jl/issues/868
    format = Documenter.HTML(;
        prettyurls = get(ENV, "CI", nothing) == "true",
        assets = ["assets/extra_styles.css", "assets/citations.css"],
    ),
    pages = _PAGES,
    # The following ensures that we only include the docstrings from
    # this module for functions defined in Base that we overwrite.
    # It also errors in case we don't include a docstring in the docs
    modules = [Dionysos],
    plugins = [
        DocumenterCitations.CitationBibliography(
            joinpath(SRC_DIR, "references.bib");
            style = :authoryear,
        ),
    ],
)

deploydocs(; repo = "github.com/dionysos-dev/Dionysos.jl.git", push_preview = true)
