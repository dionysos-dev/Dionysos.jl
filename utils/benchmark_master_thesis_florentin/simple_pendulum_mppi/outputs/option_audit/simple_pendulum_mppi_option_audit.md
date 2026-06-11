# Simple pendulum MPPI option audit

The MPPI/certification trajectory was generated once from `run_simple_pendulum_mppi_audit_copy.jl`; each row below recertifies that same candidate with one option change.

| option | cert | failed k | ok steps | V initial | V min | V median | r0 min | r0 max | X margin | U margin | statuses |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|
| lambda_more_volume | true |  | 66 | 0.00246208 | 5.40999e-06 | 0.000132074 | 0.00726311 | 0.107902 | 1.00341 | 1.00043 | base_fallback:55;max_volume_selected:11 |
| current_defaults | true |  | 66 | 0.00246073 | 5.42258e-06 | 0.000132178 | 0.00730069 | 0.107288 | 1.00336 | 1.00045 | base_fallback:55;max_volume_selected:11 |
| terminal_radius_fallback_0_50 | true |  | 66 | 0.00246073 | 5.42878e-06 | 0.000132183 | 0.00730068 | 0.107288 | 1.00353 | 1.00039 | base_fallback:54;max_volume_selected:12 |
| terminal_shrink_small | true |  | 66 | 0.00246073 | 5.36724e-06 | 0.000132163 | 0.00730068 | 0.107288 | 1.00348 | 1.00043 | base_fallback:56;max_volume_selected:10 |
| terminal_shrink_full_box | true |  | 66 | 0.00246073 | 5.45484e-06 | 0.000132188 | 0.00730068 | 0.107288 | 1.00345 | 1.00044 | base_fallback:55;max_volume_selected:11 |
| wider_adaptive_box_search | true |  | 66 | 0.00245904 | 5.42878e-06 | 0.000132128 | 0.00729697 | 0.107269 | 1.00352 | 1.00046 | base_fallback:58;max_volume_selected:8 |
| first_consistent_adaptive_box | true |  | 66 | 0.002459 | 5.37044e-06 | 0.000132004 | 0.0072969 | 0.107268 | 1.0036 | 1.00049 | accepted:66 |
| lambda_more_cost | true |  | 66 | 0.00238246 | 5.37039e-06 | 0.000133215 | 0.0076314 | 0.0993735 | 1.00354 | 1.00011 | base_fallback:55;max_volume_selected:11 |
| fixed_linearization_boxes | true |  | 66 | 0.00189902 | 5.22436e-06 | 0.000129821 | 0.00653197 | 0.0925414 | 0.0441787 | 0.0286467 | fixed:66 |
| trajectory_std_state_scaling | false | 57 | 9 | 0.000883136 | 0.000883136 | 0.0692312 | 0.00152083 | 0.18484 | 1.02892 | 1.04723 | base_fallback:7;lmi_infeasible_at_max_box:1;max_volume_selected:2 |
| no_state_scaling | false | 52 | 14 | 6.02075e-05 | 6.02075e-05 | 0.00290765 | 0.000556892 | 0.0344136 | 1.00487 | 1.00071 | base_fallback:10;lmi_infeasible_at_max_box:1;max_volume_selected:4 |
| double_default_state_scaling | false | 52 | 14 | 3.42048e-05 | 3.42048e-05 | 0.00138874 | 0.000543771 | 0.0200226 | 1.00125 | 1.00866 | base_fallback:12;lmi_infeasible_at_max_box:1;max_volume_selected:2 |
