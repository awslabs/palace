# Sharp-edge MA extrapolation of spatial-3-edge-9cd906b81cbe

MA_sharp = MA_raw + MA_tail per source; MA_raw = the sum of the per-ring MA shells; MA_tail = the `Consistent` estimator of radial_ma_profile.py: the top (sharp 90-degree) edge's -2/3 law anchored on ring 2 (Theory@2) plus the bottom (metal / trench) edge's own law fitted over rings 2-4 (Fit2-4), each remainder = the law's ring-1 energy minus the resolved ring 1; the ring-1 factor is measured on this run at each order

## p4 (spatial-3-edge-9cd906b81cbe-p4; 131 sources)

- alpha (top edge, rings 2-4): median -0.659 (quartiles -0.663 / -0.653), se median 0.002; strongest-0 median n/a; theory -0.667
- ring-1 factor (resolved / -2/3 anchored on ring 2, top edge): median 0.667, range 0.648-0.668

| deficit estimator | median (quartiles) | strongest median | at 53 | at 58 |
|---|---|---:|---:|---:|
| Fit2-4 | +2.38% (+1.57% / +3.21%) | n/a | +4.46% | +2.65% |
| Theory@2 | +4.84% (+3.83% / +5.64%) | n/a | +10.08% | +4.98% |
| Consistent (tail) | +2.63% (+2.28% / +3.35%) | n/a | +4.56% | +2.88% |

