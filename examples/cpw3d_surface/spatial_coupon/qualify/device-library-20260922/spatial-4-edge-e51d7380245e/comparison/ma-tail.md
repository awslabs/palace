# Sharp-edge MA extrapolation of spatial-4-edge-e51d7380245e

MA_sharp = MA_raw + MA_tail per source; MA_raw = the sum of the per-ring MA shells; MA_tail = the `Consistent` estimator of radial_ma_profile.py: the top (sharp 90-degree) edge's -2/3 law anchored on ring 2 (Theory@2) plus the bottom (metal / trench) edge's own law fitted over rings 2-4 (Fit2-4), each remainder = the law's ring-1 energy minus the resolved ring 1; the ring-1 factor is measured on this run at each order

## p4 (spatial-4-edge-e51d7380245e-p4; 100 sources)

- alpha (top edge, rings 2-4): median -0.659 (quartiles -0.663 / -0.645), se median 0.002; strongest-0 median n/a; theory -0.667
- ring-1 factor (resolved / -2/3 anchored on ring 2, top edge): median 0.667, range 0.644-0.668

| deficit estimator | median (quartiles) | strongest median | at 53 | at 58 |
|---|---|---:|---:|---:|
| Fit2-4 | +2.66% (+1.42% / +3.23%) | n/a | +3.32% | +24.10% |
| Theory@2 | +4.68% (+4.38% / +4.95%) | n/a | +4.72% | +25.03% |
| Consistent (tail) | +3.15% (+2.83% / +3.39%) | n/a | +3.45% | +24.10% |

