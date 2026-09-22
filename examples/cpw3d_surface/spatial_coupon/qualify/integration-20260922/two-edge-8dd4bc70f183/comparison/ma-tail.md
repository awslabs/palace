# Sharp-edge MA extrapolation of two-edge-8dd4bc70f183

MA_sharp = MA_raw + MA_tail per source; MA_raw = the sum of the per-ring MA shells; MA_tail = the `Consistent` estimator of radial_ma_profile.py: the top (sharp 90-degree) edge's -2/3 law anchored on ring 2 (Theory@2) plus the bottom (metal / trench) edge's own law fitted over rings 2-4 (Fit2-4), each remainder = the law's ring-1 energy minus the resolved ring 1; the ring-1 factor is measured on this run at each order

## p4 (two-edge-8dd4bc70f183-p4; 78 sources)

- alpha (top edge, rings 2-4): median -0.648 (quartiles -0.657 / -0.613), se median 0.006; strongest-0 median n/a; theory -0.667
- ring-1 factor (resolved / -2/3 anchored on ring 2, top edge): median 0.665, range 0.648-0.668

| deficit estimator | median (quartiles) | strongest median | at 53 | at 58 |
|---|---|---:|---:|---:|
| Fit2-4 | +2.21% (+1.34% / +3.20%) | n/a | +5.40% | +5.40% |
| Theory@2 | +4.88% (+3.93% / +5.23%) | n/a | +5.20% | +5.20% |
| Consistent (tail) | +3.07% (+2.54% / +3.56%) | n/a | +5.11% | +5.11% |

