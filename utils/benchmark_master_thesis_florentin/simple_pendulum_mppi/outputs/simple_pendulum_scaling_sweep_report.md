# Simple pendulum scaling sweep

Adaptive linearization boxes were enabled with objective `:max_volume`.
Scaling alone is not accepted unless post-hoc containment checks pass.

Best valid configuration: trajectory_std

| candidate | cert | boxes | all target | V0 | Vmin | Vmed | X margin | U margin | samples |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| trajectory_std | true | true | true | 0.001605 | 0.0002126 | 0.0007666 | 1.003 | 1.082 | 1000 |
| double_current | true | true | true | 0.001043 | 0.000155 | 0.000412 | 1.001 | 1.085 | 1000 |
| no_scaling | true | true | true | 0.0008068 | 7.051e-05 | 0.000345 | 1.005 | 1.07 | 1000 |
| current | true | true | true | 0.0005119 | 7.104e-05 | 0.0001799 | 1.008 | 1.005 | 1000 |
| half_current | true | true | true | 0.0001924 | 2.227e-05 | 6.844e-05 | 1.019 | 1.063 | 1000 |
| trajectory_range | false | false | missing | missing | 0.2466 | 0.2507 | 1.149 | 2.645 | 1000 |
