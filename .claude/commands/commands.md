---
description: List available project commands with concise usage and examples
allowed-tools: Glob, Read
---

Produce a concise reference of this project's custom Claude commands.

Build it **dynamically** — never hardcode command names, descriptions, or examples. Always
read the current contents of `.claude/commands/` so the help stays in sync as commands are
added or changed.

## Steps

1. `Glob` `.claude/commands/*.md` to find the command files. Each file `foo.md` maps to the
   command `/foo`.
2. For each file, read it (frontmatter + body) and extract:
   - **description** — the frontmatter `description:` field (fall back to the first prose
     line if absent).
   - **usage** — the frontmatter `argument-hint:` field, if any.
   - **options** — any enumerated choices in the body used to build realistic examples: an
     actions/scopes table, `Options:`/`Scopes:` lists, or bullet lists of flags.
3. Render the output below with the commands sorted **alphabetically**. Keep the summary
   table to one line per command; put the examples in the per-command sections.
4. For each command, write **3–5 realistic examples** derived from its own options and
   usage — the no-argument form plus the most common argument/flag variants. Do not invent
   options a command does not define.

## Output format

```
## Dionysos — Claude commands

| Command | Description |
|---------|-------------|
| `/name` | one-line description |
| …       | …            |

### /name
> one-line description

**Usage:** `/name <argument-hint>`

**Examples:**
- `/name` — what it does with no arguments
- `/name <option>` — what it does with that option
- `/name <arg> --flag` — what the flag changes

### Docs & guidelines
- Operating guide & commit-message format — `CLAUDE.md`
- Coding conventions — `docs/src/developers/conventions.md`
- Git workflow — `docs/src/developers/git.md`
```
