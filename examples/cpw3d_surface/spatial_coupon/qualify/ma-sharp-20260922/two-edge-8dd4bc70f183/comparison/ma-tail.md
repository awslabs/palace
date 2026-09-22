# Sharp-edge MA extrapolation of two-edge-8dd4bc70f183

MA_sharp = MA_raw + MA_tail per source; MA_raw = the sum of the per-ring MA shells; MA_tail = the `Consistent` estimator of radial_ma_profile.py: the top (sharp 90-degree) edge's -2/3 law anchored on ring 2 (Theory@2) plus the bottom (metal / trench) edge's own law fitted over rings 2-4 (Fit2-4), each remainder = the law's ring-1 energy minus the resolved ring 1; the ring-1 factor is measured on this run at each order

## p4 (two-edge-8dd4bc70f183-p4; 78 sources)

- alpha (top edge, rings 2-4): median -0.648 (quartiles -0.657 / -0.613), se median 0.006; strongest-20 median -0.646; theory -0.667
- ring-1 factor (resolved / -2/3 anchored on ring 2, top edge): median 0.665, range 0.648-0.668

| deficit estimator | median (quartiles) | strongest median | at 53 | at 58 |
|---|---|---:|---:|---:|
| Fit2-4 | +2.21% (+1.34% / +3.20%) | +2.18% | +5.40% | +5.40% |
| Theory@2 | +4.88% (+3.93% / +5.23%) | +4.82% | +5.20% | +5.20% |
| Consistent (tail) | +3.07% (+2.54% / +3.56%) | +3.02% | +5.11% | +5.11% |

## p5 (two-edge-8dd4bc70f183-p5; 78 sources)

- alpha (top edge, rings 2-4): median -0.648 (quartiles -0.658 / -0.609), se median 0.006; strongest-20 median -0.645; theory -0.667
- ring-1 factor (resolved / -2/3 anchored on ring 2, top edge): median 0.706, range 0.687-0.710

| deficit estimator | median (quartiles) | strongest median | at 53 | at 58 |
|---|---|---:|---:|---:|
| Fit2-4 | +1.89% (+0.96% / +2.85%) | +1.76% | +4.75% | +4.75% |
| Theory@2 | +4.44% (+3.42% / +4.90%) | +4.44% | +4.52% | +4.52% |
| Consistent (tail) | +2.68% (+2.22% / +3.21%) | +2.65% | +4.44% | +4.44% |

## Reference

- reference unextrapolated (no ring / edge sizing recorded): its sharp-edge deficit is modelled per source as deficit_run(Consistent) x (eps_ref / eps_run)^(1/3) from its edge size eps_ref and the run's innermost ring radius eps_run (the eps^(1/3) law of the edge MA inside a cutoff, decision 55)
- reference unextrapolated; sharp-edge deficit modelled at 3.38% (median) from its 0.5 nm edge size via the eps^(1/3) law

