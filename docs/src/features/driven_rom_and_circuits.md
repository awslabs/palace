```@raw html
<!---
Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
SPDX-License-Identifier: Apache-2.0
--->
```

# Feature: Reduced-Order Modelling and Circuit Extraction

The *Palace* driven solver has the option of building a reduced-order model (ROM) that can
accelerate the computation of a frequency sweep. Using this machinery, *Palace* can also perform a
certain type of circuit synthesis, returning $L^{-1}$, $R^{-1}$, $C$-like matrices to the user.

These features are described in two tutorials:

  - [Adaptive Driven Solver for Fast-Frequency Sweeps](../examples/tutorial_driven_uniform_v_adaptive.md)
  - [Circuit Synthesis from AC Simulations](../examples/tutorial_circuit_extraction.md)

!!! warning

    Please note that both of these features are considered "advanced". They require some care in their configuration and interpretation of the results. They are also more sensitive to algorithmic choices and implementation details that are considered "internal" to *Palace*. These may change if better numerical algorithms become available. Please proceed with caution.

## Configuration options

For reference, options related to the above two features are:

  - In [`config["Solver"]["Driven"]`](../config/solver.md#solver%5B%22Driven%22%5D):
      + `"AdaptiveTol"`
      + `"AdaptiveMaxSamples"`
      + `"AdaptiveConvergenceMemory"`
      + `"AdaptiveGSOrthogonalization"`
      + `"AdaptiveCircuitSynthesis"`
      + `"AdaptiveCircuitSynthesisDomainOrthogonalization"`
  - In [`config["Solver"]["Driven"]["Samples"]`](../config/solver.md#solver%5B%22Driven%22%5D%5B%22Samples%22%5D):
      + `"AddToPROM"`
