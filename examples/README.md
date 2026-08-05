# Examples

Runnable driver scripts showing how to use Dionysos on concrete control problems — one folder per
problem, mirroring [`../problems/`](../problems/) (which each script `include`s for the system +
specification).

These are **reference drivers**, not library code (that lives in [`../src/`](../src/)). For the
**polished, documented** tutorials that render on the website, see
[`../docs/src/examples/`](../docs/src/examples/); the scripts here are the fuller set of runnable
variants (extra specs, backends, and experiments per problem).

Each script loads its problem via `joinpath(dirname(dirname(pathof(Dionysos))), "problems", …)`, so
it runs from any working directory. See the first lines of each file for what it does.
