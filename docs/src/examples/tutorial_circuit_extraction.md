```@raw html
<!---
Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
SPDX-License-Identifier: Apache-2.0
--->
```

# Circuit Synthesis from AC Simulations

In this tutorial we discuss *Palace*'s circuit synthesis feature, which extracts a lumped-circuit
model of a device directly from the full-wave FEM. The synthesis feature is based on the rational
interpolation of the adaptive driven solver, which is discussed in [Driven Solver: Uniform vs
Adaptive](tutorial_driven_uniform_v_adaptive.md). We assume familiarity with the details therein.
We will continue to use the transmon model as our reference example, which features in the
[eigenmode](transmon.md) and [driven](tutorial_driven_uniform_v_adaptive.md) tutorials.

!!! warning "Warning: More Algorithmic Details Ahead!"

    The synthesis feature prints out circuit matrices, whose detailed form depends on "internal" algorithmic choices of the [adaptive driven solver](tutorial_driven_uniform_v_adaptive.md) and on *Palace*'s conventions.

    The internal algorithmic choices may change if better numerical algorithms become available. Furthermore, the interpretation of the output results requires some care to get right. Please proceed with caution.

!!! note

    The data can be generated with the script `examples/transmon/transmon_tutorial_circuit.jl` and the result plots generated with `examples/transmon/transmon_tutorial_circuit_plots.py`.

!!! note "Naming Convention"

    There are a few cases where variables can cause naming confusion. *Palace* uses the notation ``\bm{A}(\omega) = \bm{K} + i\omega \bm{C} - \omega^2 \bm{M}`` for the matrix at a given ``\omega``. Here ``\bm{C}`` is the loss matrix in the finite element space, which is not the same as the capacitance matrix. The ``R`` matrix of ``QR`` orthogonalization is not the circuit resistance matrix.

    For this tutorial we will consistently label circuit matrices with a hat: the inverse inductance ``\widehat{\bm{L}}^{-1}``, inverse resistance ``\widehat{\bm{R}}^{-1}``, capacitance ``\widehat{\bm{C}}``, etc. Similarly, we will label the voltage vector ``\widehat{\bm{V}}`` and current ``\widehat{\bm{I}}``. When referring to to scalar quantities, like voltage $V$, we will drop the hat when there is no risk of confusion.

## Circuit Synthesis Quick-Start

The circuit extraction requires using the adaptive driven solver (`"AdaptiveTol" > 0`) and setting
the flag [`"AdaptiveCircuitSynthesis": true`](../config/solver.md#solver%5B%22Driven%22%5D). A
`"Solver"/"Driven"` configuration might look like:

```json
"Driven": {
  "Samples": [ {"Type": "Linear", "MinFreq": 3.0, "MaxFreq": 7.0, "FreqStep": 0.01} ],
  "AdaptiveTol": 1e-5, // Use adaptive solver
  "AdaptiveCircuitSynthesis": true // Enable synthesis
}
```

Synthesis is opt-in since it alters the internal construction of the reduced-order model (ROM) of
the adaptive solver. We will discuss this further below. *Palace* will proceed to perform the
adaptive driven solver with the new ROM. At the end, it will print the standard measurement output
of the driven solver as well as the following additional files:

  - `rom-Linv-re.csv`, `rom-Rinv-re.csv`, `rom-C-re.csv`. These are the real parts of the
    synthesized inverse inductance ``\mathrm{Re}~\widehat{\bm{L}}^{-1}``, inverse resistance
    ``\mathrm{Re}~\widehat{\bm{R}}^{-1}``, and capacitance ``\mathrm{Re}~\widehat{\bm{C}}`` matrices. Since driven
    simulations require a resistive port, all three matrices will be present. The csv header
    describes the type of the node (port or synthesized). The matrices are in SI units.
  - (Optional): `rom-Linv-im.csv`, `rom-Rinv-im.csv`, `rom-C-im.csv`. The imaginary component of the
    synthesized matrices ``\mathrm{Im}~\widehat{\bm{L}}^{-1}``, ``\mathrm{Im}~\widehat{\bm{R}}^{-1}``,
    ``\mathrm{Im}~\widehat{\bm{C}}``. Each matrix is only printed when the *Palace* simulation contains a
    non-zero contribution to that matrix. For example, a [material loss
    tangent](../config/domains.md#domains%5B%22Materials%22%5D) results in a contribution to
    ``\mathrm{Im}~\widehat{\bm{C}}``. These terms may seem unfamiliar, since in "textbook" circuits ``\widehat{\bm{L}}``,
    ``\widehat{\bm{R}}``, ``\widehat{\bm{C}}`` are real.
  - `rom-orthogonalization-matrix-R.csv`: the Gram–Schmidt ``R`` factor of the synthesized circuit
    modes. This is very useful in advanced circuit postprocessing, but can be ignored by most users.

There are several constraints and considerations for using this feature:

  - Currently, the circuit synthesis only supports models that have a pure quadratic frequency
    dependence. This means `WavePort`, `WavePortPEC`, `Conductivity`, and second-order `Absorbing`
    boundary conditions are not supported.
  - All `LumpedPort` attributes must be orthogonal to each other, since these are separated out as
    individual rows and columns in the circuit matrix. For *Palace*'s Nédélec meshes, this means
    that that lumped ports cannot share parallel edges, since the degree of freedom on the edge
    contributes to both ports.
  - The
    [`"AdaptiveCircuitSynthesisDomainOrthogonalization"`](../config/solver.md#solver%5B%22Driven%22%5D)
    option controls how the non-port basis vectors are orthogonalized and therefore determines the
    voltage normalization of synthesized nodes. In general, a user should not need to switch from
    the default value of `"Energy"`.
  - `"AdaptiveCircuitSynthesis": true` requires `"AdaptiveTol" > 0`; *Palace* reports an error
    otherwise.
  - All the guidance and caveats on using the adaptive solver discussed in [Driven Solver: Uniform
    vs Adaptive](tutorial_driven_uniform_v_adaptive.md) still apply.

## Circuit Theory and Conventions

The general question of "What exactly is a circuit?" and how to connect Maxwell's equation to
circuit representation is covered in many electrical engineering textbooks. In particular, there are
subtleties and choices of convention in defining circuits when going from DC (electro- and
magnetostatic) to AC simulations. We will discuss our approach below, but refer to references [1-5]
for in-depth discussions.

The circuit synthesis of *Palace* is currently based on the adaptive driven solver — [see the driven
solver tutorial](tutorial_driven_uniform_v_adaptive.md) — so is AC *only*. The result will be an
effective (synthesized) circuit that will only accurately reproduce the response in the domain it
was trained on.

### Projective Construction

Let us recall the basics of the [ROM construction](tutorial_driven_uniform_v_adaptive.md). The
linear equation that Palace solves when evaluating an driven simulation is ``\bm{A}(\omega) \bm{x} = i \omega \bm{b}``, where

```math
\bm{A}(\omega) = \left[\bm{K} + i\omega \bm{C} - \omega^2 \bm{M}\right].
```

Here, we are neglecting the possibility of terms that non-quadratic in ``\omega`` (
``\bm{A}_2(\omega)`` and ``\bm{b}_2(\omega)``) since *Palace* cannot currently synthesize them. The
adaptive solver creates a ROM by projecting this large system onto a set of
orthogonal vectors ``\bm{Q}``. This forms projected matrices ``\bm{K}_r = \bm{Q}^T\bm{K}\bm{Q}``,
``\bm{C}_r = \bm{Q}^T\bm{C}\bm{Q}``, ``\bm{M}_r = \bm{Q}^T\bm{M}\bm{Q}``, and projected vector
``\bm{b}_r = \bm{Q}^T\bm{b}``. The projection onto a smaller basis accurately reproduces the
response of the system in the user specified frequency interval ``[f_{\mathrm{min}}, f_{\mathrm{max}}]``.

The linear equation ``[\bm{A}_r(\omega) / i \omega] \bm{x}_r = \bm{b}_r``, looks like a circuit
admittance equation ``\widehat{\bm{Y}}(\omega)\widehat{\bm{V}} = \widehat{\bm{I}}``, suggesting that
we can just identify $\bm{K}_r$ with the inverse inductance $\widehat{\bm{L}}^{-1}$, $\bm{C}_r$ with
the inverse resistance $\widehat{\bm{R}}^{-1}$, and $\bm{M}_r$ with the capacitance
$\widehat{\bm{C}}$. We could have said the exact same thing about the non-projected system.

This identification is “morally speaking“ correct [1-5]. The matrix ``\bm{K}`` arises from
discretizing the magnetic energy term ``\tfrac{1}{2}\int dV\, (\nabla\times\bm{E}^*)\mu^{-1}(\nabla\times\bm{E})/{\omega^2}`` which we identify as
``\tfrac{1}{2}\widehat{\bm{I}}^{*} \widehat{\bm{L}}^{-1} \widehat{\bm{I}}`` in a circuit. The matrix
``\bm{M}`` come from discretizing the electric energy ``\tfrac{1}{2}\int dV \, \bm{E}^* \varepsilon \bm{E}`` which we identify as ``\tfrac{1}{2}\widehat{\bm{V}}^{*} \widehat{\bm{C}} \widehat{\bm{V}}``
. And the matrix ``\bm{C}`` from discretizing the Rayleigh dissipation term.

However, just printing ``\bm{K}_r``, ``\bm{C}_r``, and ``\bm{M}_r`` matrices is neither a precise
nor a useful circuit identification, for two related reasons :

 1. The basis $\bm{Q}$ used to construct these matrices has lost the notion of system ports. There
    is no way to calculate the scattering matrix of this circuit or to cascade it with other
    circuits.
 2. The ``\bm{K}_r``, ``\bm{C}_r``, ``\bm{M}_r`` matrices have no consistent definition of voltage
    ``V`` and current ``I``. This is certainly true because these matrices are interpolation
    coefficients of a finite element basis that discretize Maxwell's equations. For example, if you
    change the polynomial `Order` of the FEM basis, the basis changes so do the matrix entries. Even
    without the finite element interpolation, however, ``V`` and ``I`` at AC are only generically
    defined by convention at ports [1,3-5]. There is no a priori notion of ``V`` and ``I`` for a 3D
    electric field solution $\bm{x}$ that make up $\bm{Q}$.

In order to make the connection to circuits, we will change both the content of the basis $\bm{Q}$
as well as its orthogonalisation rule. This means that running *Palace* with
`"AdaptiveCircuitSynthesis"` with `true` and `false` use different ROMs.

### Adding Ports to the ROM Basis

The circuit-synthesis ROM adds the port modes as the first ``N_p`` columns of the basis matrix
``\bm{Q}``. Before orthogonalization, the vectors added to the ROM basis look like:

```math
\bm{W} = \big[
\underbrace{\bm{e}_1 \;\; \bm{e}_2 \;\; \cdots \;\; \bm{e}_{N_p}}_{\text{port nodes}} \;\;
\underbrace{\mathrm{Re}\,\bm{x}(\omega^*_1) \;\; \mathrm{Im}\,\bm{x}(\omega^*_1) \;\; \cdots}_{\text{synthesised interior nodes}}\big].
```

Here ``\bm{e}_j`` are the port mode fields and ``\bm{x}(\omega^*_k)`` are high-dimensional model
(HDM) solutions at the sample frequencies. Because we add lumped ports first and demand that the
lumped ports do not overlap, the orthogonalization does does not alter the port structure. We can
sensibly interpret the top left ``N_p \times N_p`` block of the synthesized matrices as the physical
port block.  When driving the circuit with an external excitation at a port corresponds to exciting
the appropriate row ``j``. The adaptive solver appends HDM samples one frequency at a time until the
ROM meets the prescribed `"AdaptiveTol"`. The orthogonalization routine will now change these
vectors and, most notably, always remove the port mode contribution since these appear earlier in
the basis.

!!! note "Ports in the output CSV files"

    The header of every `rom-*.csv` file lists the node names in the order they appear in
    ``\bm{Q}``. Nodes `port_1_re`, `port_2_re`, …, `port_Np_re` appear first. Here the number
    is the `"Index"` values in the `"LumpedPort"` configuration. Lumped ports are currently always
    real fields, so they occupy a single column and have the label `_re` at the end.

### Orthogonalization

The vectors in ``\bm{W}`` above are coefficients in the finite element basis that discretize
Maxwell's equations. To interpret them as sensible microwave circuits, where the ``N_p \times N_p``
port block recovers our expected circuit formulae, we have to proceed a little carefully, in
defining how to orthogonalize these fields.

#### Lumped port boundary

*Palace*'s implementation of lumped ports, voltage convention choices, and definition of the
scattering matrix are all defined in the [Reference Documentation](../reference.md). We assume
familiarity with that section. Here will also assume that ports are made of a single element ``e``,
although all formulae generalize.

On the lumped port ``\Gamma_j``, the correct orthogonalization of an electric field ``E_1``  against
the port ``\bm{e}_j`` is via the inner product defined by the power flow:

```math
\langle \bm{E}_1, \bm{e}_j \rangle_{\Gamma_j} = \int_{\Gamma_j} dS \, \bm{n} \cdot (\bm{E}_1 \times \frac{[\bm{e}_j]^*}{Z_s})
```

This matches the requirements from the Lorentz reciprocity theorem. Intuitively, it means that if
every boundary mode ``\bm{e}`` carries energy independently, which is essential to recover the
correct excitation behaviour. When ``E = \bm{e}_j``, this is the defines the normalization of the
port mode. [TODO: WHAT IS THE NORAMLIZATION].

The [reference](../reference.md) discusses the relationship between physical surface impedance on a
port ``Z_s`` and the circuit characteristic impedance ``Z``. For ports that are not purely resistive
(such as junctions or lossy waveguide ports), it is generally more helpful to normalize port field
with respect to a different reference impedance ``Z_R`` which is pure resistive. For the port mode
added to the basis ``\bm{W}``, this is essential to have a frequency independent voltage
normalization and we pick the impedance of free space ``Z_R = Z_0 \approx 376.73~\mathrm{\Omega}``.
However, this will also mean that we will have to post-process the *Palace* synthesized circuit
matrices matrices carefully to recover the values (voltage, current, scattering parameters) with
respect to the conventional circuit characteristic impedance ``Z``.

#### Bulk Volume

After we have enforced port orthogonalize rule above, we need a rule for bulk degrees of freedom.
Since these will only contribute to effective synthesized nodes, we do not necessarily need to pick
a physical normalization choice — i.e. there is no "natural" voltage choice for so generic AC mode
shape in 3D. Different orthogonalization choices correspond to rotations or scaling in the
synthesized space, which cancel out once we calculated physical quantities like energy, port
scattering parameters, or eigenmode frequencies.

By default, *Palace* normalizes the field according to the system mass matrix ``\bm{M}``. For two
bulk modes this corresponds to a orthogonalization according to the inner product

```math
\langle \bm{E}_1, \bm{E}_2 \rangle_\Omega = \int_\Omega dV \bm{E}_1 \varepsilon \bm{E}_2,
```

which is the electric domain energy, up to a factor 2. An appealing aspect of this inner product is
not only that it corresponds to a sensible physical quantity and so converges to stable values with
changing finite element order or mesh refinement. We refer to the orthogonalization as `"Energy"`.

*Palace* offers the user to change the orthogonalization rule with the flag
[`"AdaptiveCircuitSynthesisDomainOrthogonalization"`](../config/solver.md#solver%5B%22Driven%22%5D),
although most users should not need this.

## Running the Transmon Model with Circuit Synthesis

We continue with the transmon model from the [eigenmode](transmon.md) and [driven
solver](tutorial_driven_uniform_v_adaptive.md) tutorials. We enable circuit synthesis with the
`"AdaptiveCircuitSynthesis": true` flag:

```json
{
  "Problem"   : {
    "Verbose": 2,
    "Output" : "postpro/transmon_tutorial_driven_rom/driven_adaptive_1e-1",
    "Type"   : "Driven"
  },
  "Boundaries": {
    "LumpedPort": [
      { "Attributes": [6], "Index": 1, "Direction": "+X", "Excitation": 1, "R": 50 },
      { "Attributes": [7], "Index": 2, "Direction": "-X", "Excitation": 2, "R": 50 },
      { "Attributes": [4], "Index": 3, "Direction": "+Y", "C": 5.5e-15, "L": 1.486e-8 }
    ],
    "PEC"       : { "Attributes": [5] },
    "Absorbing" : { "Order": 1, "Attributes": [3] }
  },
  "Model"     : { "Refinement": {"MaxIts": 0}, "Mesh": "mesh/transmon.msh2", "L0": 1.0e-6 },
  "Domains"   : {
    "Postprocessing": { "Energy": [ { "Attributes": [1], "Index": 1 } ] },
    "Materials"     : [
      { "Permittivity": 1.0, "Attributes": [2], "Permeability": 1.0 },
      {
        "LossTan"     : [ 3.0e-5,          3.0e-5,           8.6e-5          ],
        "Permittivity": [ 9.3,             9.3,              11.5            ],
        "Attributes"  : [ 1                                                  ],
        "Permeability": [ 0.99999975,      0.99999975,       0.99999979      ],
        "MaterialAxes": [ [0.8, 0.6, 0.0], [-0.6, 0.8, 0.0], [0.0, 0.0, 1.0] ]
      }
    ]
  },
  "Solver"    : {
    "Driven": {
      "Samples"    : [ {"Type": "Linear", "MinFreq": 3.5, "MaxFreq": 6.5, "FreqStep": 0.025} ],
      "AdaptiveTol": 1e-1,
      "AdaptiveCircuitSynthesis": true
    },
    "Order" : 2,
    "Linear": {"Type": "Default", "Tol": 1.0e-12, "MaxIts": 1000}
  }
}
```

In addition to the standard output (`port-S.csv`, `domain-E.csv`, etc.), *Palace* writes the circuit
files listed in the quick-start section, above.

Let us look at the difference in the standard output `domain-E.csv` between the simulations with
`"AdaptiveCircuitSynthesis": true` and `"AdaptiveCircuitSynthesis": false`.

```@raw html
<br/><p align="center">
  <img src="../../assets/examples/circuit_transmon_domain_e_comparison.svg" width="100%" />
</p><br/>
```

The left column above shows the domain electric energy ``E_\mathrm{elec}`` from the uniform driven
solver, which we use as our reference baseline. This is exactly the same data shown previously in
the [driven solver tutorial](tutorial_driven_uniform_v_adaptive.md#driving-the-transmon-model). The
middle column shows the RMS normalized absolute error between the electric energy calculated via the
adaptive solver with circuit synthesis turned on and the uniform reference. The meaning of this plot
was described [previously](tutorial_driven_uniform_v_adaptive.md#driving-the-transmon-model). The
dashed lines indicate the adaptive tolerance `"AdaptiveTol"`; however, these are *not* bounds on
the error of the domain energy. As discussed in the driven solver tutorial, the adaptive tolerance
controls the error in the electric field coefficients and not the error in a derived quantity such
as energy, so the dashed line is only a rough order-of-magnitude heuristic.

The right column shows the RMS scaled absolute error between ``E_\mathrm{elec}`` calculated through
the adaptive driven solver with circuit synthesis turned on ``E_\mathrm{circuit}`` and turned off
``E_\mathrm{adaptive}``. It is substantially smaller than the error relative to the uniform
reference data in the middle column, but still larger than solver precision. The choice of HDM
frequency samples is the same for both — e.g. at `AdaptiveTol` of `1e-1` the combined pivot
frequencies from both excitations are ``\{3.5, 5.02, 5.899, 5.9, 6.5\}~\mathrm{GHz}``. How can we
understand the discrepancy?

As discussed in the theory section above, during circuit synthesis *Palace* builds a different ROM
by (a) adding extra port modes to the basis $\bm{Q}$, so it spans a larger space and (b) using a
different orthogonalization weight, so the basis of this larger space is different. Adding port
modes will change the condition number ``\kappa[A_r(\omega)]`` of the ROM system matrix. We already
know that linear solves will have a worse conditioning close to the poles and errors will spike. We
can estimate this error with matrix stability analysis [6], but that is beyond the scope of this
tutorial. Furthermore, the different orthogonalization routine used with synthesis adds many more
floating-point operations when calculating the basis ``\bm{Q}``. However, we expect the additional
error from these operations to be mostly frequency-independent and contribute to the
noise floor on the plots above.

For completeness, we show the analogous plot for the scattering parameters from `port-S.csv`. This
plot has a similar interpretation to the `domain-E.csv` plot above.

```@raw html
<br/><p align="center">
  <img src="../../assets/examples/circuit_transmon_port_s_comparison.svg" width="100%" />
</p><br/>
```

Let us now look at the `rom-*.csv` output. For this model, there are `5` output matrices of size
``19 \times 19``: `rom-Linv-re.csv`, `rom-Rinv-re.csv`, `rom-C-re.csv`, `rom-C-im.csv`, and
`rom-orthogonalization-matrix-R.csv`. The header of the `csv` files describes the origin of the
nodes as described above (shortened for brevity):

```csv
port_1_re,
port_2_re,
port_3_re,
sample_e1_s0_re,
sample_e1_s0_im,
sample_e1_s1_re,
sample_e1_s1_im,
# ...
sample_e2_s0_re,
sample_e2_s0_im,
# ...
sample_e2_s3_re,
sample_e2_s3_im
```

The first three columns correspond to the port modes, which are all real for lumped ports. The
remaining columns correspond to synthesized nodes coming from samples during the adaptive solve. The
label `e1` / `e2` corresponds to the excitation and `s0`, `s1`, ... are the sample indices within
each excitation. As we tighten the `"AdaptiveTol"` of the driven solver that gave rise to this
model, the number of samples increases and so does the size of the synthesized matrices.

Let us look at the first 5 rows and columns of `rom-Rinv-re.csv`:

```csv
 port_1_re, port_2_re,  port_3_re,  sample_e1_s0_re,  sample_e1_s0_im,  ...
 2.000e-02, 0.000e+00,  0.000e+00,       -1.799e-23,       -2.878e-22,  
 0.000e+00, 2.000e-02,  0.000e+00,        1.799e-23,        2.878e-22,  
 0.000e+00, 0.000e+00,  0.000e+00,        0.000e+00,        0.000e+00,  
-1.799e-23, 1.799e-23,  0.000e+00,        6.091e-10,        1.320e-08,  
-2.878e-22, 2.878e-22,  0.000e+00,        1.320e-08,        3.293e-07,  
...
```

The diagonal elements of `port_1` and `port_2` are the expected ``1 / (50~\Omega)`` port values. The
full row and column of `port_3` are zero, since it is pure $L$, $C$. The other rows and columns
correspond to the orthogonalized basis of samples. For the first sample shown above, we can think of
this as the driven solve with the lumped port modes removed. Why does each HDM sample have a small
diagonal dissipation? First, the model above has `Absorbing` boundary conditions, which means there
is dissipation beyond what happens at the resistive ports. Second, the electric field of a HDM
sample at a given frequency $\omega$ may have a different shape at a lumped port than the port mode
(i.e., mode shape is not enforced at port). This means that the HDM could pick up more or less than
the full $1 / (50~\Omega)$ dissipative term.

Now consider `rom-C-re.csv`:

```csv
 port_1_re  port_2_re  port_3_re  sample_e1_s0_re  sample_e1_s0_im, ...
 4.689e-16  0.000e+00  0.000e+00       -5.844e-17       -9.542e-17,
 0.000e+00  4.282e-16  0.000e+00       -4.632e-17        8.473e-17,
 0.000e+00  0.000e+00  5.529e-15       -6.781e-22       -2.415e-20,
-5.844e-17 -4.632e-17 -6.781e-22        3.542e-14       -2.962e-20,
-9.542e-17  8.473e-17 -2.415e-20       -2.962e-20        3.542e-14,
...
```

The diagonal value of port 3 is `5.529e-15`, which is larger than the `5.5e-15` term we added to `C`
in the `LumpedPort` configuration. Additionally, ports 1 and 2 have small diagonal contributions. In
all these cases, this is because the boundary port term will still pick up some of the capacitive
contribution from the volume term of the mesh elements neighbouring the port surface. If we were to
refine the mesh elements close to the port, the size of this term would decrease.

Finally, *Palace* outputs a `rom-C-im.csv` for this model. This is because of the presence of the
`"LossTan"` in our model configuration (see [theory
reference](../reference.md#mathematical-background)).

## Literature & References

[1] R. B. Adler, L. J. Chu, and R. M. Fano, Electromagnetic energy transmission and radiation, 5. print. in The MIT Press classics. Cambridge, Mass.: M.I.T. Press, 1969.

[2] J.-M. Jin, The finite element method in electromagnetics, 3. ed. Piscataway, NJ: IEEE Press, 2014.

[3] R. B. Marks and D. F. Williams, "A general waveguide circuit theory," J. Res. Natl. Inst. Stan., vol. 97, no. 5, p. 533, Sep. 1992, doi: 10.6028/jres.097.024.

[4] D. M. Pozar, Microwave engineering, Fourth edition. Hoboken, NJ: John Wiley & Sons, Inc, 2012.

[5] G. Wendt et al., Elektrische Felder und Wellen / Electric Fields and Waves, vol. 4 / 16. in Handbuch der Physik / Encyclopedia of Physics, vol. 4 / 16. Berlin, Heidelberg: Springer Berlin Heidelberg, 1958. doi: 10.1007/978-3-642-45895-8.

[6] N. J. Higham, Accuracy and stability of numerical algorithms, 2nd ed. in Other titles in applied mathematics, no. 80. Philadelphia, Pa: Society for Industrial and Applied Mathematics (SIAM, 3600 Market Street, Floor 6, Philadelphia, PA 19104), 2002. doi: 10.1137/1.9780898718027.

[7] Z. K. Minev et al., "Energy-participation quantization of Josephson circuits," *npj Quantum Information*, vol. 7, p. 131, 2021. doi: 10.1038/s41534-021-00461-8.
