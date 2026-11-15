Task-specific prompt preview (version 2)

These are the exact prompts saved by the temperature runner. The 40 selected record/field pairs match the repaired run. No API calls have been made for this prepared run.

Example without legacy labels:

```text
You translate one natural-language requirement into a Linear Temporal Logic (LTL) formula.

Task atomic propositions:
AP = {blue, brown, green, purple, yellow}
Each color proposition is true exactly when the robot is in that color's region at the current time.
Use only these five color names as atomic propositions in the output.

Semantics:

- Use standard LTL over infinite traces, evaluated from the current time.
- `U` is strong Until: `(a U b)` requires `b` to eventually hold, while `a` holds at every preceding time.
- “Whenever,” “every time,” and “globally” introduce a recurring global requirement.
- “From then on” means that the consequent must hold globally from the triggering time.
- Do not add a global operator unless the requirement expresses a global or recurring condition.

Allowed syntax:

- Atomic propositions: `blue`, `brown`, `green`, `purple`, `yellow`.
- Negation: `!`
- Conjunction: `&`
- Disjunction: `|`
- Implication: `->`
- Equivalence: `<->`
- Next: `X`
- Eventually: `F`
- Globally: `G`
- Until: `U`

Translation rules:

- Preserve the complete logical and temporal meaning.
- Use the exact lowercase color names listed in AP.
- Do not introduce new propositions.
- Do not introduce domain assumptions, such as mutual exclusion, reachability, or ordering, unless explicitly stated.
- Do not use every available proposition unless the requirement needs it.
- Use parentheses to make operator scope unambiguous.
- Before answering, silently check temporal scope, ordering, eventuality, and whether each Until obligation is strong.

Output rules:

- Output exactly one LTL formula.
- Do not output explanations, labels, comments, or code fences.

Examples:

Natural-language requirement:
Whenever the robot is in the blue region, it must remain in the yellow region from then on.

LTL:
G((blue) -> G(yellow))

Natural-language requirement:
Whenever the robot is in the blue region, it must stay in the green region until it reaches the yellow region.

LTL:
G((blue) -> ((green) U (yellow)))

Natural-language requirement:
Eventually reach the blue region, and from that time eventually reach the yellow region.

LTL:
F((blue) & F(yellow))

Natural-language requirement:
The robot must never enter the purple region.

LTL:
G(!(purple))

Natural-language requirement:
Stay in the blue region until reaching the brown region.

LTL:
(blue) U (brown)

Natural-language requirement:
At the next step, the robot must be in the green region.

LTL:
X(green)

Target input:

Available atomic proposition set:
AP = {blue, brown, green, purple, yellow}

Natural-language requirement:
<requirement>
Stay in yellow until either (you are in both purple and green, or in neither, and remain in blue forever from then on) or (you are in exactly one of purple and green and eventually leave blue).
</requirement>

LTL:
```

Alias sections for requirements containing legacy labels:

Record 615, field natural_paraphrase_gpt54mini_rewritten:

```text
These labels are aliases for the color propositions; output their color names.
- `prop_2` means `yellow` (the robot is in the yellow region).
```

Record 913, field natural_paraphrase_gpt54mini_rewritten:

```text
These labels are aliases for the color propositions; output their color names.
- `prop_2` means `brown` (the robot is in the brown region).
```

Record 620, field paraphrase_deepseek_gpt54mini_rewritten:

```text
These labels are aliases for the color propositions; output their color names.
- `prop_1` means `yellow` (the robot is in the yellow region).
```

Record 1064, field paraphrase_gemini-2.5-flash_gpt54mini_rewritten:

```text
These labels are aliases for the color propositions; output their color names.
- `prop_2` means `green` (the robot is in the green region).
```

Record 1002, field natural_paraphrase_gpt54mini_rewritten:

```text
These labels are aliases for the color propositions; output their color names.
- `prop_2` means `blue` (the robot is in the blue region).
```

