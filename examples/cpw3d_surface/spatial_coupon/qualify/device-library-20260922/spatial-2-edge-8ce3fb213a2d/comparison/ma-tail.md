# Sharp-edge MA extrapolation of spatial-2-edge-8ce3fb213a2d

MA_sharp = MA_raw + MA_tail per source; MA_raw = the sum of the per-ring MA shells; MA_tail = the `Consistent` estimator of radial_ma_profile.py: the top (sharp 90-degree) edge's -2/3 law anchored on ring 2 (Theory@2) plus the bottom (metal / trench) edge's own law fitted over rings 2-4 (Fit2-4), each remainder = the law's ring-1 energy minus the resolved ring 1; the ring-1 factor is measured on this run at each order

## p4 (spatial-2-edge-8ce3fb213a2d-p4; 134 sources)

- alpha (top edge, rings 2-4): median -0.663 (quartiles -0.665 / -0.661), se median 0.001; strongest-0 median n/a; theory -0.667
- ring-1 factor (resolved / -2/3 anchored on ring 2, top edge): median 0.667, range 0.650-0.668

| deficit estimator | median (quartiles) | strongest median | at 53 | at 58 |
|---|---|---:|---:|---:|
| Fit2-4 | +3.18% (+2.70% / +3.47%) | n/a | +3.01% | +3.18% |
| Theory@2 | +4.69% (+4.16% / +4.91%) | n/a | +4.90% | +4.76% |
| Consistent (tail) | +3.33% (+2.89% / +3.55%) | n/a | +3.18% | +3.32% |

