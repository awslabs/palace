# Sharp-edge MA extrapolation of spatial-3-edge-5d3b5e644745

MA_sharp = MA_raw + MA_tail per source; MA_raw = the sum of the per-ring MA shells; MA_tail = the `Consistent` estimator of radial_ma_profile.py: the top (sharp 90-degree) edge's -2/3 law anchored on ring 2 (Theory@2) plus the bottom (metal / trench) edge's own law fitted over rings 2-4 (Fit2-4), each remainder = the law's ring-1 energy minus the resolved ring 1; the ring-1 factor is measured on this run at each order

## p4 (spatial-3-edge-5d3b5e644745-p4; 159 sources)

- alpha (top edge, rings 2-4): median -0.660 (quartiles -0.663 / -0.656), se median 0.001; strongest-0 median n/a; theory -0.667
- ring-1 factor (resolved / -2/3 anchored on ring 2, top edge): median 0.667, range 0.652-0.668

| deficit estimator | median (quartiles) | strongest median | at 53 | at 58 |
|---|---|---:|---:|---:|
| Fit2-4 | +2.28% (+1.59% / +3.07%) | n/a | +5.87% | +3.19% |
| Theory@2 | +4.90% (+2.30% / +5.67%) | n/a | +10.27% | +4.76% |
| Consistent (tail) | +2.51% (+2.11% / +3.23%) | n/a | +6.00% | +3.33% |

