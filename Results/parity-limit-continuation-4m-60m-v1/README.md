# Remaining parity-limit cases, with cooldown breaks

159 remaining ground-truth parity cases, plus six resolved pilot decisions
reused with provenance. Same production Spot bridge and independent verifier.
Four-million-node and sixty-million-edge caps; 600 seconds per new case,
including independent verification and export. Cases run one at a time.
Julia verification uses one thread. No LLM/model requests.

Monitored process-group RSS guard: 3 GiB. System-wide free-memory floor: 15%.
These are memory guards, not temperature sensors or OS allocation quotas.
Every new completed case is followed by 60 seconds with no synthesis process.
After roughly 30 minutes of case-processing time, use a 5-minute break instead.
The final case needs no cooldown because all processing has finished.
Cool-downs reduce sustained load; CPU temperature is not measured.

Pause by sending SIGINT to this runner. Its active child process group is
stopped in finally; interrupted cases are not saved as failed decisions.
Resume with the same command; completed case/result.json files are reused.
progress.json distinguishes active synthesis, verification, cooldown and
waiting for memory. summary.json complete=true marks final report completion.
Original audits and pilot files remain unchanged. Controller artifacts retain
the frozen plant-cache path and hash.
