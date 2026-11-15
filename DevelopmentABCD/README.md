# Frozen development controller-synthesis experiment

Primary: DeepSeek V4 Pro, T=0.7, 40 development prompts, 10 existing outputs
per prompt, k=1..5. No model API calls. Mistral T=0.7, k_max=4 is reserved.

Candidates were parsed, checked for allowed APs, checked individually for SAT,
and clustered by equivalence in the frozen upstream report. This runner verifies
its provenance against all 40 original response checkpoints. Clusters rank by
frequency with earliest-sample ties. Prefixes saturate at min(k,m).

Each of the 200 prompt/budget cases receives exactly one outcome:

1. `controller_issued`: configured initial product states are winning, universal
   infinite-horizon verification passes, and the policy has been exported.
2. `logically_inconsistent`: Spot proves the conjunction UNSAT.
3. `plant_unrealizable`: the initial label has no automaton transition, or the
   configured initial product state is outside the converged robust winning set.
4. `computational_failure`: timeout, unsupported automaton, worker/resource/tool
   failure, verification failure, or controller-export failure.
5. `no_valid_candidate`: no cluster survives the upstream filter.

Only outcome 1 counts as issuance. Outcomes 2–5 abstain. Tool failures never
become UNSAT or plant-unrealizability. All rates use all 40 prompts per k.
Pending cases are separate; an incomplete run is never labelled complete.

There is one shared vehicle plant for the 40 instances because the benchmark
supplies no per-prompt plants. A **fresh** abstraction matches the configuration
in `InteractiveWithFeedBack.jl`: grid origin (0,0,pi/2), state spacing
(0.2,0.2,0.3), input spacing (0.3,0.3), DT=0.3, obstacle margin 0.12, initial
state (4.8,5.0,pi/2), and the existing region labels. The legacy abstraction.jld2
has a different grid origin and is not substituted. The fresh plant is cached
and fingerprinted. Label semantics are exactly those implemented in the local
`LTL2Automata.labeler_and_inner_sets`.

"Plant-unrealizable" is relative to this labelled finite abstraction and initial
state, not a proof of physical impossibility for the continuous system. Controller
guarantees likewise refer to the labelled abstraction and its abstraction assumptions.

The current solver requires deterministic state-based Büchi automata. A
nondeterministic result from its translator is recorded as unsupported under
computational failure. No claim of unrealizability follows from that limitation.

The per-distinct-formula wall-clock limit is frozen at 120 seconds, including SAT,
automaton translation, game solving, verification and export. Backend startup
and a small JIT warm-up are excluded, with a separate 300-second startup limit.
Duplicate saturated prefixes share the same job and artifact but remain separate
prompt/budget observations. A timed-out Julia process and its children are killed;
a fresh worker loads the fingerprinted plant for the next job.

From the Dionysos checkout:

```sh
python3 DevelopmentABCD/run.py --dionysos-root . \
  --source Results/task-specific-temperature-prompt-v2 \
  --output Results/deepseek-development-abcd-v1
```

Repeat to resume. A changed source/code/configuration requires a new output
directory. The completed experiment is `Results/deepseek-development-abcd-v1/`.
The unsuffixed directory retains a failed initial backend startup, with no
development cases evaluated; it is not the completed experiment.

- `frozen-development.json`: immutable settings, prompt cases, source/code hashes.
- `plant.jls`, `plant-identity.json`: plant, labels, physical input vectors and
  representative state centers, plus identity.
- `jobs/`, `job-results/`: unique formula jobs and checkpointed outcomes.
- `controllers/*.jls`: product policy, winning set, automaton/update, configured
  initial product states and universal verification certificate, linked to plant.
- `summary.json`, `summary.csv`, `comparison.md`: all 200 cases and outcome table.
- `worker.log`: backend and failure diagnostics.

Python 3 standard library, the existing Dionysos Julia environment, and Spot CLI
are required. The GUI and interactive LLM pipeline are not started.

```sh
python3 DevelopmentABCD/test_runner.py
julia --project=. test/ltl/test_buchi_synthesis.jl
julia --project=. DevelopmentABCD/audit_controllers.jl . Results/deepseek-development-abcd-v1
```
