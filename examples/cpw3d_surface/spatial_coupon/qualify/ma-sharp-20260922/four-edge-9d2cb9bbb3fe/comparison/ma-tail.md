# Sharp-edge MA extrapolation of four-edge-9d2cb9bbb3fe

MA_sharp = MA_raw + MA_tail per source; MA_raw = the sum of the per-ring MA shells; MA_tail = the `Consistent` estimator of radial_ma_profile.py: the top (sharp 90-degree) edge's -2/3 law anchored on ring 2 (Theory@2) plus the bottom (metal / trench) edge's own law fitted over rings 2-4 (Fit2-4), each remainder = the law's ring-1 energy minus the resolved ring 1; the ring-1 factor is measured on this run at each order

## p4 (four-edge-9d2cb9bbb3fe-p4; 60 sources)

- alpha (top edge, rings 2-4): median -0.655 (quartiles -0.661 / -0.639), se median 0.003; strongest-20 median -0.660; theory -0.667
- ring-1 factor (resolved / -2/3 anchored on ring 2, top edge): median 0.666, range 0.644-0.668

| deficit estimator | median (quartiles) | strongest median | at 53 | at 58 |
|---|---|---:|---:|---:|
| Fit2-4 | +2.27% (+1.33% / +2.69%) | +1.99% | +0.00% | +1.48% |
| Theory@2 | +4.68% (+3.21% / +5.45%) | +2.44% | +0.01% | +4.39% |
| Consistent (tail) | +2.86% (+2.40% / +3.21%) | +2.14% | +0.01% | +4.39% |

## Reference

- reference unextrapolated (no ring / edge sizing recorded): its sharp-edge deficit is modelled per source as deficit_run(Consistent) x (eps_ref / eps_run)^(1/3) from its edge size eps_ref and the run's innermost ring radius eps_run (the eps^(1/3) law of the edge MA inside a cutoff, decision 55)
- reference unextrapolated; sharp-edge deficit modelled at 3.64% (median) from its 0.5 nm edge size via the eps^(1/3) law

