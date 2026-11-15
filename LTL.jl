"""
Utilities for parsing and analyzing LTL formulas.

The main entry point is [`analyze_ltl`](@ref).  No Julia packages are required,
but Spot's `ltl2tgba` executable must be installed to compute automaton sizes.
"""

abstract type LTLNode end

struct LTLAtom <: LTLNode
    name::String
    is_constant::Bool
end

struct LTLUnary <: LTLNode
    operator::String
    child::LTLNode
end

struct LTLBinary <: LTLNode
    operator::String
    left::LTLNode
    right::LTLNode
end

struct LTLToken
    kind::Symbol
    text::String
    position::Int
end

mutable struct LTLParser
    tokens::Vector{LTLToken}
    current::Int
end

const UNARY_OPERATORS = Set(["!", "X", "F", "G"])
const TEMPORAL_OPERATORS = Set(["X", "F", "G", "U", "R", "W", "M"])

# Larger numbers bind more tightly.  These precedences follow Spot's syntax.
const BINARY_OPERATOR_INFO = Dict(
    "<->" => (precedence = 10, associativity = :left),
    "xor" => (precedence = 15, associativity = :left),
    "->"  => (precedence = 20, associativity = :right),
    "|"   => (precedence = 30, associativity = :left),
    "&"   => (precedence = 40, associativity = :left),
    "U"   => (precedence = 50, associativity = :right),
    "R"   => (precedence = 50, associativity = :right),
    "W"   => (precedence = 50, associativity = :right),
    "M"   => (precedence = 50, associativity = :right),
)

const OPERATOR_ALIASES = Dict(
    "~" => "!",
    "&&" => "&",
    "||" => "|",
    "=>" => "->",
    "<=>" => "<->",
    "[]" => "G",
    "<>" => "F",
    "^" => "xor",
    "V" => "R",
)

"""Return the character at `position`, or `nothing` after the end of a string."""
function character_at(input::AbstractString, position::Int)
    position > lastindex(input) && return nothing
    return input[position]
end

"""Test whether `input` starts with `text` at the supplied string index."""
function starts_at(input::AbstractString, text::AbstractString, position::Int)
    position > lastindex(input) && return false
    return startswith(SubString(input, position), text)
end

"""
Tokenize the subset of Spot's LTL syntax used by this project.

Uppercase `X`, `F`, and `G` at the start of a token are unary operators, so
compact formulas such as `GFprop_1` are read as `G(F(prop_1))`.  Atomic
propositions beginning with those letters should be double-quoted, as in
`"G_enabled"`.
"""
function tokenize_ltl(formula::AbstractString)
    input = String(formula)
    tokens = LTLToken[]
    position = firstindex(input)

    while position <= lastindex(input)
        character = input[position]

        if isspace(character)
            position = nextind(input, position)
            continue
        end

        if character == '('
            push!(tokens, LTLToken(:left_parenthesis, "(", position))
            position = nextind(input, position)
            continue
        elseif character == ')'
            push!(tokens, LTLToken(:right_parenthesis, ")", position))
            position = nextind(input, position)
            continue
        end

        matched_operator = nothing
        for candidate in ("<=>", "<->", "&&", "||", "=>", "->", "[]", "<>")
            if starts_at(input, candidate, position)
                matched_operator = candidate
                break
            end
        end
        if matched_operator !== nothing
            normalized = get(OPERATOR_ALIASES, matched_operator, matched_operator)
            kind = normalized in UNARY_OPERATORS ? :unary_operator : :binary_operator
            push!(tokens, LTLToken(kind, normalized, position))
            for _ in matched_operator
                position = nextind(input, position)
            end
            continue
        end

        if character in ('!', '~')
            push!(tokens, LTLToken(:unary_operator, get(OPERATOR_ALIASES, string(character), string(character)), position))
            position = nextind(input, position)
            continue
        elseif character in ('&', '|', '^')
            operator = get(OPERATOR_ALIASES, string(character), string(character))
            push!(tokens, LTLToken(:binary_operator, operator, position))
            position = nextind(input, position)
            continue
        elseif character in ('X', 'F', 'G')
            push!(tokens, LTLToken(:unary_operator, string(character), position))
            position = nextind(input, position)
            continue
        elseif character in ('U', 'R', 'V', 'W', 'M')
            operator = get(OPERATOR_ALIASES, string(character), string(character))
            push!(tokens, LTLToken(:binary_operator, operator, position))
            position = nextind(input, position)
            continue
        end

        if character == '"'
            start = position
            position = nextind(input, position)
            escaped = false
            closed = false
            while position <= lastindex(input)
                current = input[position]
                if escaped
                    escaped = false
                elseif current == '\\'
                    escaped = true
                elseif current == '"'
                    position = nextind(input, position)
                    closed = true
                    break
                end
                position = nextind(input, position)
            end
            closed || throw(ArgumentError("Unterminated quoted proposition at character $start"))
            push!(tokens, LTLToken(:atom, input[start:prevind(input, position)], start))
            continue
        end

        start = position
        while position <= lastindex(input)
            current = input[position]
            if isspace(current) || current in ('(', ')', '!', '~', '&', '|', '^')
                break
            end
            if any(candidate -> starts_at(input, candidate, position),
                   ("<=>", "<->", "&&", "||", "=>", "->", "[]", "<>"))
                break
            end
            position = nextind(input, position)
        end
        position == start && throw(ArgumentError("Unexpected character '$(input[position])' at character $position"))

        text = input[start:prevind(input, position)]
        lowercase_text = lowercase(text)
        if lowercase_text == "xor"
            push!(tokens, LTLToken(:binary_operator, "xor", start))
        else
            is_constant = lowercase_text in ("true", "false", "tt", "ff") || text in ("0", "1")
            push!(tokens, LTLToken(is_constant ? :constant : :atom, text, start))
        end
    end

    push!(tokens, LTLToken(:end_of_input, "", lastindex(input) + 1))
    return tokens
end

peek(parser::LTLParser) = parser.tokens[parser.current]

function consume!(parser::LTLParser)
    token = peek(parser)
    parser.current += 1
    return token
end

function parse_primary!(parser::LTLParser)
    token = peek(parser)

    if token.kind == :unary_operator
        consume!(parser)
        return LTLUnary(token.text, parse_primary!(parser))
    elseif token.kind == :atom || token.kind == :constant
        consume!(parser)
        return LTLAtom(token.text, token.kind == :constant)
    elseif token.kind == :left_parenthesis
        consume!(parser)
        expression = parse_expression!(parser, 0)
        closing = peek(parser)
        closing.kind == :right_parenthesis || throw(ArgumentError(
            "Expected ')' at character $(closing.position), found '$(closing.text)'"
        ))
        consume!(parser)
        return expression
    end

    displayed = token.kind == :end_of_input ? "end of formula" : "'$(token.text)'"
    throw(ArgumentError("Expected an atomic proposition, unary operator, or '(' at character $(token.position); found $displayed"))
end

function parse_expression!(parser::LTLParser, minimum_precedence::Int)
    left = parse_primary!(parser)

    while true
        token = peek(parser)
        token.kind == :binary_operator || break
        information = BINARY_OPERATOR_INFO[token.text]
        information.precedence < minimum_precedence && break

        consume!(parser)
        next_minimum = information.associativity == :left ? information.precedence + 1 : information.precedence
        right = parse_expression!(parser, next_minimum)
        left = LTLBinary(token.text, left, right)
    end

    return left
end

"""Parse an LTL formula and return its AST."""
function parse_ltl(formula::AbstractString)
    isempty(strip(formula)) && throw(ArgumentError("The LTL formula cannot be empty"))
    parser = LTLParser(tokenize_ltl(formula), 1)
    root = parse_expression!(parser, 0)
    trailing = peek(parser)
    trailing.kind == :end_of_input || throw(ArgumentError(
        "Unexpected token '$(trailing.text)' at character $(trailing.position)"
    ))
    return root
end

function collect_metrics(node::LTLAtom)
    return (
        ast_size = 1,
        temporal_depth = 0,
    )
end

function collect_metrics(node::LTLUnary)
    child = collect_metrics(node.child)
    return (
        ast_size = child.ast_size + 1,
        temporal_depth = child.temporal_depth + (node.operator in TEMPORAL_OPERATORS ? 1 : 0),
    )
end

function collect_metrics(node::LTLBinary)
    left = collect_metrics(node.left)
    right = collect_metrics(node.right)
    return (
        ast_size = left.ast_size + right.ast_size + 1,
        temporal_depth = max(left.temporal_depth, right.temporal_depth) +
                         (node.operator in TEMPORAL_OPERATORS ? 1 : 0),
    )
end

"""Compute structural metrics from a parsed LTL AST."""
ltl_ast_metrics(root::LTLNode) = collect_metrics(root)

"""
    automaton_size(formula; buchi = false, ltl2tgba = "ltl2tgba")

Build a Büchi automaton using Spot and return one size value, defined as the
number of states plus the number of edges.  With `buchi=false`, this runs:

    ltl2tgba -f FORMULA --stats=%s,%e

With `buchi=true`, `-B` is added before `-f`.

Arguments are passed directly to the process, without a shell.
"""
function automaton_size(
    formula::AbstractString;
    buchi::Bool = false,
    ltl2tgba::AbstractString = "ltl2tgba",
)
    executable = Sys.which(ltl2tgba)
    executable === nothing && throw(ArgumentError(
        "Could not find '$ltl2tgba'. Install Spot or pass the executable path using the ltl2tgba keyword."
    ))

    arguments = buchi ?
        [executable, "-B", "-f", String(formula), "--stats=%s,%e"] :
        [executable, "-f", String(formula), "--stats=%s,%e"]
    command = Cmd(arguments)
    output = try
        strip(read(command, String))
    catch error
        throw(ErrorException("ltl2tgba failed for formula '$formula': $(sprint(showerror, error))"))
    end

    values = split(output, ',')
    length(values) == 2 || throw(ErrorException(
        "Unexpected ltl2tgba statistics output '$output'; expected 'states,edges'"
    ))

    states = tryparse(Int, strip(values[1]))
    edges = tryparse(Int, strip(values[2]))
    (states === nothing || edges === nothing) && throw(ErrorException(
        "Could not parse ltl2tgba statistics output '$output'"
    ))

    return states + edges
end

"""
    analyze_ltl(formula; ltl2tgba = "ltl2tgba")

Return structural metrics and both Spot automaton sizes for an LTL formula.
`automaton_size` uses `ltl2tgba -f`, while `buchi_automaton_size` uses
`ltl2tgba -B -f`.  Each size is the number of states plus edges.

# Example

```julia
include("LTL.jl")
info = analyze_ltl("G(F(prop_1))")
println(info.ast_size)             # 3
println(info.temporal_depth)       # 2
println(info.automaton_size)       # states + edges
println(info.buchi_automaton_size) # Büchi states + edges
```
"""
function analyze_ltl(formula::AbstractString; ltl2tgba::AbstractString = "ltl2tgba")
    root = parse_ltl(formula)
    structural = ltl_ast_metrics(root)
    standard_size = automaton_size(formula; ltl2tgba = ltl2tgba)
    buchi_size = automaton_size(formula; buchi = true, ltl2tgba = ltl2tgba)
    return merge(structural, (
        automaton_size = standard_size,
        buchi_automaton_size = buchi_size,
    ))
end

"""
    are_ltl_equivalent(first_formula, second_formula; ltlfilt = "ltlfilt")

Return `true` when Spot determines that two LTL formulas are equivalent and
`false` otherwise.  This runs the equivalent of:

    ltlfilt --equivalent-to=FIRST_FORMULA -f SECOND_FORMULA --count

Arguments are passed directly to the process, without a shell.
"""
function are_ltl_equivalent(
    first_formula::AbstractString,
    second_formula::AbstractString;
    ltlfilt::AbstractString = "ltlfilt",
)
    executable = Sys.which(ltlfilt)
    executable === nothing && throw(ArgumentError(
        "Could not find '$ltlfilt'. Install Spot or pass the executable path using the ltlfilt keyword."
    ))

    command = Cmd([
        executable,
        "--equivalent-to=$(String(first_formula))",
        "-f",
        String(second_formula),
        "--count",
    ])
    standard_output = IOBuffer()
    standard_error = IOBuffer()
    process = run(pipeline(ignorestatus(command); stdout = standard_output, stderr = standard_error))

    process.exitcode == 0 && return true
    process.exitcode == 1 && return false

    error_message = strip(String(take!(standard_error)))
    throw(ErrorException(
        isempty(error_message) ?
        "ltlfilt failed with exit code $(process.exitcode)" :
        "ltlfilt failed: $error_message",
    ))
end
