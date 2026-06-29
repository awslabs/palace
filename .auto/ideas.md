# Ideas backlog — volume/domain postprocessing

- Establish `PALACE_VOLUME_PROFILE=1` baseline first: compare `volume_paraview_seconds`, `boundary_paraview_seconds`, `final_volume_paraview_seconds`, and full `paraview_seconds`.
- Fusing domain derived point fields is the main structural candidate: `U_e`, `U_m`, and `S` currently have separate `DomainPointFieldEvaluator` operators and writer callbacks. A multi-field evaluator could reuse E/B interpolation, geometry, and material lookup and fill multiple VTK arrays in one pass.
- If fusion is too large, try a narrower `U_m + S` fusion because both use B/H and inverse permeability; `S` also needs E.
- If the volume split shows base GridFunction domain output dominates, investigate MFEM `RegisterField` VTU output for base E/B rather than libCEED derived fields.
- If `operator_construction_seconds` becomes a target, revisit buffer-only domain evaluator assembly. Prior run #15 improved construction/HWM but did not improve old full `paraview_seconds`.
- Avoid repeating failed micro-optimizations from the previous segment: refined-count metadata caching, trace zero-fill/IR caching, thread-local payload buffers, and larger ofstream buffers.
