# Transition-level certificate audit

Monte Carlo stress tests are empirical diagnostics, not formal proofs. A transition is treated as locally formal only when the LMI solved, the source ellipsoid is inside the local linearization box, the affine feedback maps the full source ellipsoid inside U, and the direct transition-level replay has no violations at the selected tolerance.

## Indexing diagnostics
- k=1: ||E_k.c - x_k||=0.0, ||E_{k+1}.c - x_{k+1}||=0.0, ||kappa(E_k.c)-u_k||=1.0413170815559502e-7
- k=2: ||E_k.c - x_k||=0.0, ||E_{k+1}.c - x_{k+1}||=0.0, ||kappa(E_k.c)-u_k||=2.0731526175935266e-5
- k=3: ||E_k.c - x_k||=0.0, ||E_{k+1}.c - x_{k+1}||=0.0, ||kappa(E_k.c)-u_k||=2.877374635420756e-7
- k=4: ||E_k.c - x_k||=0.0, ||E_{k+1}.c - x_{k+1}||=0.0, ||kappa(E_k.c)-u_k||=1.322869940700193e-6
- k=5: ||E_k.c - x_k||=0.0, ||E_{k+1}.c - x_{k+1}||=0.0, ||kappa(E_k.c)-u_k||=1.87591424809747e-8

## Verdicts
- none__box_1.0: NOT FORMAL, global_success=1.0, box_all=false, input_all=true, empirical_pass=true, sim_ok=true
- none__box_2.0: FORMAL, global_success=1.0, box_all=true, input_all=true, empirical_pass=true, sim_ok=true
- none__box_3.0: FORMAL, global_success=1.0, box_all=true, input_all=true, empirical_pass=true, sim_ok=true
- trajectory_std__box_1.0: NOT FORMAL, global_success=1.0, box_all=false, input_all=true, empirical_pass=true, sim_ok=true
- trajectory_std__box_2.0: FORMAL, global_success=1.0, box_all=true, input_all=true, empirical_pass=true, sim_ok=true
- trajectory_std__box_3.0: FORMAL, global_success=1.0, box_all=true, input_all=true, empirical_pass=true, sim_ok=true
- trajectory_range__box_1.0: NOT FORMAL, global_success=1.0, box_all=false, input_all=true, empirical_pass=true, sim_ok=true
- trajectory_range__box_2.0: NOT FORMAL, global_success=1.0, box_all=false, input_all=true, empirical_pass=true, sim_ok=true
- trajectory_range__box_3.0: FORMAL, global_success=1.0, box_all=true, input_all=true, empirical_pass=true, sim_ok=true
