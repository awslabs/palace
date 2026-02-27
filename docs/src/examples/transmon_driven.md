```@raw html
<!---
Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
SPDX-License-Identifier: Apache-2.0
--->
```

# Driven Solver: Uniform vs Adaptive

In this tutorial, we will discusses the driven solver (LINK) in more detail. This computes the
frequency-domain response (steady-state) of a system driven by external excitations. We will
especially focus on properties of the “adaptive” driven solves, which uses reduced-order modelling
(ROM) techniques. These ROM techniques can dramatically increase the speed for a frequency-domain,
but require more careful usage.

The uniform and adaptive driven solvers were also discussed in the example on [cross-talk between
coplanar waveguides](cpw.md).

In this tutorial, we will go into more of the algorithmic detail of the adaptive solvers. These are
helpful not only to use this solver mode more effectively, but also to understand other features
based on it like the [circuit extraction](transmon_circuit.md). We will consider the [cpw
example](cpw.md) in more detail as well as apply the driven solver to the [single transmon
model](transmon.md). The adaptive solver works particularly well for the transmon case.

!!! warning

    This tutorial is advanced. The algorithmic choices described here are considered "internal" to
    *Palace* and may change at any time if the developers decide that alternative approaches are better.
    Please proceed with caution.

!!! note

    The simulations below were performed with git has ... on ...

    TODO: Where to get files and folders.


## Reminder of the Transmon Model

We consider the same model of a transmon discussed in the tutorial here LINK. We encourage the reader to work through that example first, before continuing below.

The system consists of a transmon qubit, a readout-resonator, and a feed-line with two $50~ \Omega$ ports at the end. The qubit is modelled as a lumped with a capacitive and inductance surface impedance. The model has two eigenmodes in the frequency region XX - YY, corresponding to the mode on the transmon and one on the resonator.

Picture of GDS

## Uniform Driven Solver

We use exactly the same mesh and basic set up as in the tutorial. But we modify our configuration files as appropriate for a driven solver.

Specifically 


When we run this configuration, *Palace* will iterate over each requested frequency point and solve the full linear equation $$[]$$, see the reference section. 


Driven Solver Config file: Uniform solve (what does it do)
- Single-Excitation & Multi-excitation

We can look at — for example — the scattering parameters:



We see resonance dips corresponding to the eigenmodes as well as the background signal.

These scattering parameters are currently only calculated on purely dissipative ports with the assumption that the reference impedance $Z_R$ and surface impedance values $R$ specified the configuration file.

Of course, we can analyze this data

Fitting the uniform solve using Vector Fit and AAA with determinant surrugate, plot. 



Can get response and eigenmodes, but requires more detailed checks and tuning.
Validate against eigenmodes. We will return to this point in our tutorial on 
synthisizing circuits (LINK).

## Adaptive Driven Solver and Reduced Order Model

- The uniform solver is great baseline, but quite slow since it solves indepedantly for every frequenccy.
- If we make a finer mesh, this becomes really bad
- As we can see from our vector fitting examples, we can extract much of the information from a much smaller sample of data. This leads us to the idea of reduced-order modeling.
- The idea is to solve the full system only at a few key frequency points and the recustruct the rest of the response at other frequncies bases on that. This is implemented in the "adaptive" driven solver of *Palace*.
- We will first show and example of running with this solver and then return to discuss the
  algorithmic details.


### A first run of the adaptive tol


### Background reduced order model based on solutions

- Adaptive Solve based on Reduced order modeling
	- Give good references at end
	- How does it work? Mode shapes
- Only manages to fit what it can see — think of it as an rational matrix interpolation on the real axis

### More details reduced order modeling solver

- Tolerances and what they mean
- How frequency points are chose
- How does it pick frequencies: manual "AddToPROM" and "AdaptiveTol"
- AdaptiveConvergenceMemory

### Tuning tolerances and convergence

What is the error tolerance.


## Summary and take-aways

We have now reached a natural break if this tutorial 

- Works really well since there are few poles in or close to the region of interest.
- Takeaway: rom_tol >> solver_tol but there are limits.
- It is a greedy algorithm — default AdaptiveConvergenceMemory is low at two — if it is a difficult case might need to increase for challenging problem.
- When not to use this — (small number of samples compared to number of poles needed). High precision in difficult region.

A more difficult example is here.

---

## Obtaining a circuit out of the Adaptive Driven Solver

- How to interpret the circuits
- Is this a normal circuit?
- How to post-process these circuits

## CPW Line with and LC Port

!! 

- Convergence: "AdaptiveTol"


- Getting the circuit parameters
- Changing circuit parameters — warning

- CPW Line with port — a more difficult example
	- Why?


## Literature & References


Complex conventions of stability