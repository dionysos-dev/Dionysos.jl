---
description: Commit changes following the project's "[ACTION] module: description" convention
argument-hint: [optional hint about what to commit or how to split]
allowed-tools: Bash(git status:*), Bash(git diff:*), Bash(git add:*), Bash(git commit:*), Bash(git restore:*), Bash(git rev-parse:*)
---

Commit the current changes following this project's commit convention.

## Message format

`[ACTION] module: description`

- **description** — concise, lowercase, no trailing period, ≤ 60 characters.
- Do **not** add a `Co-Authored-By` line.

### Actions

| Action | Meaning |
| :----- | :------ |
| `ADD`  | new feature |
| `IMP`  | improvement to an existing feature |
| `FIX`  | bug fix |
| `REF`  | refactor (no behaviour change) |
| `REM`  | removal |
| `MOV`  | move / rename |
| `REV`  | revert |

### Module

The subsystem the change belongs to, inferred from the touched paths:

- `src/utils/` → `utils`, `src/system/` → `system`, `src/problem/` → `problem`,
  `src/mapping/` → `mapping`, `src/symbolic/` → `symbolic`, `src/optim/` → `optim`
  (use the solver family, e.g. `optim/uniform-grid`, when it sharpens the message).
- `ext/` → the extension (e.g. `ext/csv`), `test/` → `test`, `docs/` → `docs`,
  `problems/<name>/` → `problems/<name>`, `scripts/<name>/` → `scripts/<name>`.
- Repo config / CI / tooling → `meta`.

Example: `FIX ext/csv: export controller via the controller protocol`

## Steps

1. Re-read the **Git workflow** section of the project-root `CLAUDE.md` and honour it
   (never commit to `master` — branch first; format before committing).
2. `git status` — see the full working tree (staged, unstaged, untracked).
3. `git diff --cached` and `git diff` — review staged and unstaged changes.
4. Decide **one commit or several**. Split when changes are logically independent
   (a bug fix plus an unrelated improvement, or edits across different modules with
   different purposes); keep as one when everything serves a single purpose.
   - Changes across unrelated modules → separate commits.
   - If one commit would stage **more than 15 files**, warn and suggest splitting
     unless every file clearly serves the same purpose.
5. For each commit:
   a. `git add` only the files of this logical unit — include the relevant unstaged and
      untracked files. Never stage unrelated files or secrets (`.env`, credentials, key
      material).
   b. Pick the `ACTION` and `module` from the nature and location of the changes.
   c. Write the description, then commit with the formatted message.
6. Repeat for each logical commit when splitting.
7. `git status` again to verify nothing is left behind. If leftovers belong to the work
   just committed, amend or add a follow-up commit; if they are unrelated, state
   explicitly what was left out and why.

User hint: $ARGUMENTS
