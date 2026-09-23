# Sharp-edge MA extrapolation of spatial-10-edge-65450ff47b9b

MA_sharp = MA_raw + MA_tail per source; MA_raw = the sum of the per-ring MA shells; MA_tail = the `Consistent` estimator of radial_ma_profile.py: the top (sharp 90-degree) edge's -2/3 law anchored on ring 2 (Theory@2) plus the bottom (metal / trench) edge's own law fitted over rings 2-4 (Fit2-4), each remainder = the law's ring-1 energy minus the resolved ring 1; the ring-1 factor is measured on this run at each order

## p4 (spatial-10-edge-65450ff47b9b-p4; 191 sources)

- alpha (top edge, rings 2-4): median -0.657 (quartiles -0.663 / -0.656), se median 0.002; strongest-0 median n/a; theory -0.667
- ring-1 factor (resolved / -2/3 anchored on ring 2, top edge): median 0.666, range 0.653-0.668

| deficit estimator | median (quartiles) | strongest median | at 53 | at 58 |
|---|---|---:|---:|---:|
| Fit2-4 | +2.17% (+0.22% / +2.24%) | n/a | +2.47% | +2.21% |
| Theory@2 | +5.62% (+0.24% / +5.75%) | n/a | +5.29% | +5.69% |
| Consistent (tail) | +2.39% (+0.24% / +2.46%) | n/a | +2.67% | +2.42% |

