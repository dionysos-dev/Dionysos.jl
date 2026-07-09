function llm_to_spot_ltl(s::AbstractString)::String
   # 1) normalize whitespace
   t = strip(replace(String(s), r"\s+" => " "))


   # 2) tokenize words and punctuation we care about
   # tokens: identifiers (prop_1), operators, parentheses
   toks = String[]
   i = 1
   while i <= lastindex(t)
       c = t[i]
       if c in ('(', ')')
           push!(toks, string(c))
           i = nextind(t, i)
       elseif c in ('&','|','!')
           push!(toks, string(c))
           i = nextind(t, i)
       elseif c == '-'
           # possibly "->"
           if i < lastindex(t) && t[nextind(t,i)] == '>'
               push!(toks, "->")
               i = nextind(t, nextind(t,i))
           else
               push!(toks, "-")
               i = nextind(t,i)
           end
       elseif isspace(c)
           i = nextind(t, i)
       else
           # read a word/identifier
           j = i
           while j <= lastindex(t)
               cj = t[j]
               if isspace(cj) || cj in ('(', ')', '&', '|', '!', '-')
                   break
               end
               j = nextind(t, j)
           end
           push!(toks, t[i:prevind(t,j)])
           i = j
       end
   end


   # 3) map English-like keywords to Spot tokens
   function map_token(tok::String)::String
       w = lowercase(tok)
       if w == "globally" || w == "always"
           return "G"
       elseif w == "finally" || w == "eventually"
           return "F"
       elseif w == "next"
           return "X"
       elseif w == "until"
           return "U"
       elseif w == "release" || w == "releases"
           return "R"
       elseif w == "and"
           return "&"
       elseif w == "or"
           return "|"
       elseif w == "negation" || w == "not" || w == "negate"
           return "!"
       elseif w == "imply" || w == "implies" || w == "implie" || w == "implicate" || w == "implies," || w == "imply,"
           return "->"
       elseif w == "true" || w == "tt"
           return "true"
       elseif w == "false" || w == "ff"
           return "false"
       else
           return tok
       end
   end


   toks = [map_token(tok) for tok in toks if !isempty(tok)]


   # 4) rebuild string with sensible spacing:
   #    - no spaces after unary temporal ops before '(' when present: G( ... )
   #    - otherwise keep single spaces between tokens
   out = IOBuffer()
   prev = ""
   for tok in toks
       if tok == ")"
           print(out, tok)
       elseif tok == "("
           # attach "G(" and "F(" etc.
           if prev in ("G","F","X")
               # remove trailing space if any
               str = String(take!(out))
               str = replace(str, r"\s+$" => "")
               out = IOBuffer(); print(out, str)
               print(out, tok)
           else
               if position(out) > 0 && !endswith(String(take!(IOBuffer(String(take!(out))))), " ")
                   # no-op; we handle spacing below
               end
               # restore out
               # (simpler: ensure a space if needed)
               #
               # We'll just add directly with a preceding space if last char is alnum or ')'
               str = String(take!(out))
               out = IOBuffer(); print(out, str)
               if !isempty(str)
                   lastc = str[end]
                   if isletter(lastc) || isnumeric(lastc) || lastc == ')'
                       print(out, ' ')
                   end
               end
               print(out, tok)
           end
       else
           # normal token
           str = String(take!(out))
           out = IOBuffer(); print(out, str)
           if !isempty(str)
               lastc = str[end]
               if lastc != '(' && lastc != ' '
                   print(out, ' ')
               end
           end
           print(out, tok)
       end
       prev = tok
   end


   res = strip(String(take!(out)))


   # 5) fix common LLM artifact: extra closing parentheses
   n_open  = count(==( '(' ), res)
   n_close = count(==( ')' ), res)
   if n_close > n_open
       # remove extra ')' from the end only
       extra = n_close - n_open
       while extra > 0 && endswith(res, ")")
           res = chop(res)
           extra -= 1
           res = rstrip(res)
       end
   end


   # 6) final whitespace cleanup
   res = replace(res, r"\s+" => " ")
   return res
end


"""Replace prop_i tokens in text with their original AP names using invmap."""
function props_to_colors(txt::AbstractString, invmap::Dict{String,String})
   out = String(txt)
   # Replace higher indices first to avoid prop_1 touching prop_10-like tokens (if any)
   props = sort(collect(keys(invmap)), by=length, rev=true)
   for p in props
       w = invmap[p]
       rx = Regex("\\b" * p * "\\b")
       out = replace(out, rx => w)
   end
   return out
end


#############################################
# NL <-> prop_i mapping (colors/obstacles)
#############################################


"""Replace task-relevant region names in NL with prop_i tokens (for LLM input).


Returns:
 nl_mapped::String, invmap::Dict{String,String}


Where invmap maps "prop_1" => "blue" etc. We map obstacles to "obs" (to match your AP name).


Notes:
 - Matching is case-insensitive.
 - Replacements are done on word boundaries to avoid partial matches.
"""
function nl_colors_to_props(nl::AbstractString)
   # Canonical AP surface names used in your code/labels
   vocab = [
       "blue",
       "brown",
       "purple",
       "green",
       "yellow",
       # obstacle synonyms -> obs
       "obstacle",
       "obstacles",
       "gray",
       "grey",
   ]


   # Find which tokens appear (case-insensitive)
   present = String[]
   low = lowercase(String(nl))
   for w in vocab
       # crude presence test, refined replacement uses regex word boundaries
       occursin(w, low) && push!(present, w)
   end
   present = unique(present)


   # Deterministic ordering so prop indices are stable
   order = Dict(
       "blue"=>1, "brown"=>2, "purple"=>3, "green"=>4, "yellow"=>5,
       "obstacle"=>6, "obstacles"=>6, "gray"=>6, "grey"=>6,
   )
   sort!(present, by = w -> get(order, w, 10_000))


   # Build mapping word -> prop_k
   word_to_prop = Dict{String,String}()
   invmap = Dict{String,String}()
   k = 1
   for w in present
       # merge all obstacle-like words into the same logical AP surface name "obs"
       canon = (w in ("obstacle","obstacles","gray","grey")) ? "obs" : w
       # If canon already assigned, reuse its prop id
       if haskey(word_to_prop, canon)
           word_to_prop[w] = word_to_prop[canon]
           continue
       end
       p = "prop_$(k)"
       word_to_prop[w] = p
       word_to_prop[canon] = p
       invmap[p] = canon
       k += 1
   end


   # Apply replacements on word boundaries, longest-first to be safe
   # (e.g., "obstacles" before "obstacle")
   keys_sorted = sort(collect(keys(word_to_prop)), by=length, rev=true)
   out = String(nl)
   for w in keys_sorted
       p = word_to_prop[w]
       # word boundary regex (case-insensitive)
       rx = Regex("\\b" * w * "\\b", "i")
       out = replace(out, rx => p)
   end


   return out, invmap
end



function fallback_nl_to_spot_formula(nl_sentence::AbstractString)
    s = lowercase(strip(String(nl_sentence)))

    # normalize punctuation/whitespace a bit
    s = replace(s, r"[\.;,]+" => " ")
    s = replace(s, r"\s+" => " ")

    # canonicalize obstacle names first
    s = replace(s, "obstacles" => "obs")
    s = replace(s, "obstacle" => "obs")

    aps = ("yellow", "blue", "brown", "purple", "green", "obs")

    # collect APs in order of appearance in the sentence
    found = Tuple{Int,String}[]
    for ap in aps
        idx = findfirst(ap, s)
        if idx !== nothing
            start_idx = idx isa UnitRange ? first(idx) : idx
            push!(found, (start_idx, ap))
        end
    end
    sort!(found, by = first)
    ordered_aps = [ap for (_, ap) in found]

    isempty(ordered_aps) && return ""

    first_ap = ordered_aps[1]
    second_ap = length(ordered_aps) >= 2 ? ordered_aps[2] : nothing

    # --- Binary temporal patterns first ---
    # Examples:
    #   "brown until blue" -> "brown U blue"
    #   "negate brown until blue" -> "! brown U blue"
    #   "not brown until blue" -> "! brown U blue"
    if occursin("until", s) && second_ap !== nothing
        left_ap = first_ap
        right_ap = second_ap

        neg_left = occursin("negate " * left_ap, s) ||
                   occursin("not " * left_ap, s) ||
                   occursin("avoid " * left_ap, s) ||
                   occursin("never " * left_ap, s)

        left = neg_left ? "! " * left_ap : left_ap
        return left * " U " * right_ap
    end

    # --- Unary wrappers ---
    if occursin("eventually", s) || occursin("finally", s) || occursin("reach", s)
        if occursin("negate " * first_ap, s) || occursin("not " * first_ap, s)
            return "F ! " * first_ap
        end
        return "F " * first_ap
    elseif occursin("always", s) || occursin("globally", s)
        if occursin("negate " * first_ap, s) || occursin("not " * first_ap, s)
            return "G ! " * first_ap
        end
        return "G " * first_ap
    elseif occursin("not ", s) || occursin("avoid", s) || occursin("never", s) || occursin("negate", s)
        return "G ! " * first_ap
    else
        return first_ap
    end
end

function is_obviously_incomplete_spot_formula(s::AbstractString)
    t = strip(String(s))
    isempty(t) && return true

    # Bare operators with no operand(s)
    t in ("X", "F", "G", "!", "U", "R", "&", "|", "->", "<->") && return true

    # Ends with an operator / incomplete connective
    endswith(t, "!") && return true
    endswith(t, "U") && return true
    endswith(t, "R") && return true
    endswith(t, "&") && return true
    endswith(t, "|") && return true
    endswith(t, "->") && return true
    endswith(t, "<->") && return true

    # Unbalanced parentheses
    nopen = count(==( '(' ), t)
    nclose = count(==( ')' ), t)
    nopen != nclose && return true

    return false
end

"""Normalize obvious punctuation artifacts from LLM output before Spot parsing."""
function sanitize_spot_formula_candidate(s::AbstractString)
    t = String(s)
    t = strip(t)
    # Remove common leading/trailing punctuation artifacts from generation.
    t = replace(t, r"^[,;:\.\s]+" => "")
    t = replace(t, r"[,;:\.\s]+$" => "")
    t = replace(t, r"\s+" => " ")
    return strip(t)
end

"""Heuristic guard for obviously invalid Spot formulas.
This catches punctuation-only strings, commas, and characters outside a small
safe LTL/AP alphabet used in this project.
"""
function is_obviously_invalid_spot_formula(s::AbstractString)
    t = sanitize_spot_formula_candidate(s)
    isempty(t) && return true

    # Only reject punctuation-only garbage here.
    # Do not over-filter valid formulas; Spot is the real parser.
    has_alnum = any(c -> (isletter(c) || isdigit(c)), t)
    !has_alnum && return true

    return false
end