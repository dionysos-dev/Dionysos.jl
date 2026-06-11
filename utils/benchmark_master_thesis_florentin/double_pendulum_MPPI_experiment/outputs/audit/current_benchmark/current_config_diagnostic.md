# Double Pendulum Current Benchmark Diagnostic

## Resultat global

- Configuration auditee: defaults actuels de `DoublePendulumMPPIConfig`, trajectoire sauvegardee, `state_scaling = diag([2, 2, 4, 4])`, boites fixes (`adaptive_linearization_boxes = false`).
- Horizon candidat: 66 inputs. La certification echoue a `k=34`.
- Transitions OK: 32. Tentative suivante: 34. Si on compte le terminal comme ellipsoide, cela fait 33 ellipsoides stockes.
- Plus petit ellipsoide certifie: `k=35`, profondeur 32, log-volume `-68.713`, rayons axes `[2.7304333e-08, 7.9797822e-08, 1.9554165e-07, 3.2734429e-07]`.
- Fichiers: `utils/benchmark_master_thesis_florentin/double_pendulum_MPPI_experiment/outputs/audit/current_benchmark/current_config_step_diagnostics.csv`, `utils/benchmark_master_thesis_florentin/double_pendulum_MPPI_experiment/outputs/audit/current_benchmark/current_config_phase_summary.csv`.

## Phases observees

| phase | k range | depth | count | logvol min | logvol max | max Lx | max reqX/box | max reqU/box |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| terminal-large | 66-64 | 1-3 | 3 | -3.5025 | -0.65187 | 7.1857 | 14.776 | 0.28063 |
| box-inconsistent | 63-59 | 4-8 | 5 | -11.135 | -4.7964 | 7.998 | 4.689 | 0.17952 |
| large-remainder | 58-52 | 9-15 | 7 | -21.686 | -12.584 | 12.777 | 0.69116 | 0.47403 |
| box-inconsistent | 51-51 | 16-16 | 1 | -21.684 | -21.684 | 13.599 | 0.099837 | 1.3496 |
| large-remainder | 50-38 | 17-29 | 13 | -24.862 | -22.165 | 18.927 | 0.083953 | 0.54727 |
| box-inconsistent | 37-36 | 30-31 | 2 | -30.852 | -26.489 | 6.1339 | 0.03986 | 1.959 |
| large-remainder | 35-35 | 32-32 | 1 | -68.713 | -68.713 | 4.9685 | 3.2734e-06 | 0.53606 |
| failed | 34-34 | 33-33 | 1 | nan | nan | 6.3415 | -inf | -inf |

## Points dominants

### Plus fortes bornes de linearisation `L`

| k | depth | phase | max Lx | logvol | reqU/box |
|---:|---:|---|---:|---:|---:|
| 46 | 21 | large-remainder | 18.927 | -23.227 | 0.14144 |
| 45 | 22 | large-remainder | 18.122 | -23.54 | 0.41873 |
| 47 | 20 | large-remainder | 17.408 | -24.322 | 0.51994 |
| 44 | 23 | large-remainder | 15.78 | -24.003 | 0.16083 |
| 48 | 19 | large-remainder | 13.991 | -23.198 | 0.54371 |
| 51 | 16 | box-inconsistent | 13.599 | -21.684 | 1.3496 |

### Plus gros depassements de boite controle

| k | depth | phase | reqU/box | reqX/box | logvol | rayons axes |
|---:|---:|---|---:|---:|---:|---|
| 37 | 30 | box-inconsistent | 1.959 | 0.03986 | -26.489 | `[0.00056611816, 0.0014056317, 0.0032254872, 0.0039860128]` |
| 36 | 31 | box-inconsistent | 1.8462 | 0.014185 | -30.852 | `[0.00012918631, 0.00043211425, 0.0011638165, 0.0014185393]` |
| 51 | 16 | box-inconsistent | 1.3496 | 0.099837 | -21.684 | `[0.0034245502, 0.0041205616, 0.0099837348, 0.0093891147]` |
| 50 | 17 | large-remainder | 0.54727 | 0.083953 | -22.165 | `[0.0032580953, 0.0040067325, 0.0083953425, 0.0061837338]` |
| 48 | 19 | large-remainder | 0.54371 | 0.060634 | -23.198 | `[0.0032852869, 0.0038557803, 0.0060634448, 0.0028386937]` |
| 38 | 29 | large-remainder | 0.54123 | 0.047343 | -24.862 | `[0.00085712248, 0.0017646379, 0.0047343425, 0.0044092267]` |

## Zone de rupture

| k | status | phase | logvol | min axis | max axis | max Lx | reqU/box | cost |
|---:|---|---|---:|---:|---:|---:|---:|---:|
| 42 | ok | large-remainder | -24.083 | 0.0018278 | 0.0041044 | 10.483 | 0.22467 | 401.18 |
| 41 | ok | large-remainder | -24.07 | 0.0016137 | 0.0044301 | 7.7601 | 0.22751 | 271.22 |
| 40 | ok | large-remainder | -24.071 | 0.0013774 | 0.0048886 | 6.687 | 0.42183 | 207.72 |
| 39 | ok | large-remainder | -24.381 | 0.0010936 | 0.0047289 | 6.9861 | 0.2294 | 287.49 |
| 38 | ok | large-remainder | -24.862 | 0.00085712 | 0.0047343 | 6.7579 | 0.54123 | 109.59 |
| 37 | ok | box-inconsistent | -26.489 | 0.00056612 | 0.003986 | 6.1339 | 1.959 | 312.97 |
| 36 | ok | box-inconsistent | -30.852 | 0.00012919 | 0.0014185 | 5.295 | 1.8462 | 302.71 |
| 35 | ok | large-remainder | -68.713 | 2.7304e-08 | 3.2734e-07 | 4.9685 | 0.53606 | 137.12 |
| 34 | infeasible | failed | nan | nan | nan | 6.3415 | nan | nan |

## Lecture technique

- Les premiers steps partent d’un terminal large: a `k=66`, les rayons requis en etat montent a environ `14.8 x` la boite fixe `DeltaX = 0.1`. En mode fixe, ce depassement est seulement diagnostique; il ne rend pas le step infeasible dans le code actuel.
- Les phases `large-remainder` sont dominees par les composantes vitesse de `L` (`max Lx` jusqu’a `18.93` a `k=46`). C’est le signal principal de non-linearite/couplage du double pendule le long de cette trajectoire.
- Les depassements de boite controle sont localises: `reqU/box = 1.35` a `k=51`, puis `1.96` et `1.85` a `k=37..36`. Avec `DeltaU = 0.01`, une boite controle autour de `0.02` couvrirait ces steps OK, mais elle recalculerait aussi des bornes de Taylor plus grandes.
- La vraie rupture est la contraction brutale a `k=35`: log-volume `-68.7`, rayons axes `2.7e-8..3.3e-7`. Le step suivant `k=34` devient infeasible non pas parce que `L` y est maximal, mais parce que l’ellipsoide cible issu de `k=35` est quasiment degenere.
- Donc le benchmark est complexe par combinaison de terminal large, dynamique fortement non lineaire en vitesse, trajectoire avec commandes saturees puis changement brutal, et absence de verification de coherence des boites en mode fixe.

## Conseils immediats

- Ne conclus pas simplement “boites plus grandes”. Pour une certification sound, active plutot `adaptive_linearization_boxes=true` et laisse le diagnostic rejeter les steps incoherents.
- Si tu restes en boites fixes, `DeltaU = 0.02` est le premier test raisonnable; `DeltaX` plus grand que `0.1` risque surtout d’augmenter `L`, sauf pour rendre les premiers steps terminaux coherents, ou il faudrait en fait des boites enormes.
- Teste un terminal plus petit (`terminal_john_shrink < 1`, par exemple `0.5` ou `0.7`) ou une trajectoire terminale moins rapide; le terminal actuel autorise des vitesses jusqu’a `+-5`, ce qui donne des premiers ellipsoides tres larges.
- Le scaling `diag([2,2,4,4])` n’est pas manifestement absurde, mais les bornes `L` restent dominees par les vitesses. Un scaling domaine `diag([pi, pi, 6, 6])` vaut un sweep cible, pas une hypothese a adopter sans comparaison de volume et `failed_k`.
