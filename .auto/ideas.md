# Ideas backlog — volume/domain postprocessing

- Established `PALACE_VOLUME_PROFILE=1` baseline: volume save was 15.84-17.07s before base-field changes. Temporary detail profiling showed base domain E/B GridFunction VTU output dominated; derived `U_e/U_m/S` payloads were much smaller.
- Fusing domain derived point fields is the main structural candidate: `U_e`, `U_m`, and `S` currently have separate `DomainPointFieldEvaluator` operators and writer callbacks. A multi-field evaluator could reuse E/B interpolation, geometry, and material lookup and fill multiple VTK arrays in one pass.
- If fusion is too large, try a narrower `U_m + S` fusion because both use B/H and inverse permeability; `S` also needs E.
- Base GridFunction domain output dominated and was improved structurally: current best uses libCEED `FIELD_E/FIELD_B` probe evaluators for E/B plus point-major domain point buffers, reducing `volume_paraview_seconds` to ~1.92s. Future work should validate field ordering/values against MFEM on a small case.
- If `operator_construction_seconds` becomes a target, revisit buffer-only domain evaluator assembly. Prior run #15 and the volume-segment retry improved construction/HWM but worsened primary `volume_paraview_seconds`.
- Avoid repeating failed micro-optimizations from the previous segment and volume segment: refined-count metadata caching, trace zero-fill/IR caching, thread-local payload buffers, larger ofstream buffers, domain IR pointer caching, chunked point-major float conversion, and buffer-only evaluator assembly for this primary metric.
