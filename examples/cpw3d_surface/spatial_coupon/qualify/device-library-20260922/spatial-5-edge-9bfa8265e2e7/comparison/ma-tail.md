# Sharp-edge MA extrapolation of spatial-5-edge-9bfa8265e2e7

MA_sharp = MA_raw + MA_tail per source; MA_raw = the sum of the per-ring MA shells; MA_tail = the `Consistent` estimator of radial_ma_profile.py: the top (sharp 90-degree) edge's -2/3 law anchored on ring 2 (Theory@2) plus the bottom (metal / trench) edge's own law fitted over rings 2-4 (Fit2-4), each remainder = the law's ring-1 energy minus the resolved ring 1; the ring-1 factor is measured on this run at each order

## p4 (spatial-5-edge-9bfa8265e2e7-p4; 119 sources)

- alpha (top edge, rings 2-4): median -0.656 (quartiles -0.659 / -0.645), se median 0.004; strongest-0 median n/a; theory -0.667
- ring-1 factor (resolved / -2/3 anchored on ring 2, top edge): median 0.666, range 0.642-0.668

| deficit estimator | median (quartiles) | strongest median | at 53 | at 58 |
|---|---|---:|---:|---:|
| Fit2-4 | +0.01% (+0.00% / +0.32%) | n/a | +0.00% | -0.00% |
| Theory@2 | +0.01% (+0.00% / +0.61%) | n/a | +0.00% | +0.00% |
| Consistent (tail) | +0.01% (+0.00% / +0.40%) | n/a | +0.00% | +0.00% |

