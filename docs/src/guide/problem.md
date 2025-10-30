```@raw html
<!--- Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved. --->
<!--- SPDX-License-Identifier: Apache-2.0 --->
```

# Problem Types

## Eigenmode problems

For eigenmode simulations,
[`config["Problem"]["Type"]: "Eigenmode"`](../config/problem.md#config%5B%22Problem%22%5D),
the user should specify a nonzero (but arbitrarily small) frequency above which to search
for eigenmodes. The computed eigenvalues are written to an ASCII file named `eig.csv`, in
the directory specified by
[`config["Problem"]["Output"]`](../config/problem.md#config%5B%22Problem%22%5D). Also in
this file are the mode quality factors and errors (absolute and backward) computed for each
eigenpair.

Calculations related to
[energy-participation ratio (EPR) quantization](https://www.nature.com/articles/s41534-021-00461-8)
can be performed with *Palace* when the user specifies lumped ports corresponding to the
linearized lumped circuit elements in the model. In this case, the participation matrix for
inductive elements is automatically generated for the specified number of modes and number
of inductive lumped ports. The participation matrix is output in an ASCII file named
`port-EPR.csv`.

The EPR framework can be used to characterize the dissipative elements in the model as well.
In particular, lumped ports with nonzero resistance in the model will trigger coupling rate
and quality factor calculations based on input-output (I-O) line coupling loss: By
specifying resistive lumped ports in the model, the mode coupling quality factors will be
computed as ``Q_{ml} = \omega_m/\kappa_{ml}``. The output file `port-Q.csv` will be created
in the output directory containing these mode quality factor contributions. For bulk and
interface dielectric loss calculations, which are not unique to the eigenmode simulation
type, see the sections [Domain postprocessing](postprocessing.md#Domain-postprocessing) and
[Boundary postprocessing](postprocessing.md#Boundary-postprocessing) of this guide.

## Driven problems in the frequency domain

For frequency domain driven simulations,
[`config["Problem"]["Type"]: "Driven"`](../config/problem.md#config%5B%22Problem%22%5D), the
model is excited by a time harmonic incident field (port boundary) or surface current.
The user can specify a port excitation using
[lumped ports or numeric wave ports](boundaries.md#Lumped-and-wave-port-excitation).

The default frequency sweep behavior for frequency domain driven simulations is to perform a
uniform sampling from the minimum to the maximum specified frequency of interest, using the
user specified step size. An adaptive fast frequency sweep strategy can also be used,
activated by specifying a nonzero value for `"AdaptiveTol"` under the
[`config["Solver"]["Driven"]`](../config/solver.md#solver%5B%22Driven%22%5D) object. In this
case, using the high-dimensional model solution computed at a few automatically selected
frequency samples, a low-cost model is constructed and used to compute the frequency
response over the entire frequency range of interest. The specified error tolerance ensures
that the approximate low-cost model is reliably accurate relative to the high-dimensional
model within the frequency band of interest. This is particularly useful for
fine-resolution sweeps containing many sample points, where it can yield a significant
speedup over the default strategy.

Port scattering parameters, or S-parameters, are postprocessed for the column of the
scattering matrix corresponding to the driven port index automatically for this simulation
type and stored in an ASCII file named `port-S.csv`, in the directory specified by
[`config["Problem"]["Output"]`](../config/problem.md#config%5B%22Problem%22%5D). Both the
``\text{dB}`` magnitude (``20\log_{10}(|S_{ij}|)``) and the phase ``\angle(S_{ij})``
(in degrees) are written to the file. In the case that more than a single lumped or wave
port is excited or surface current excitations are used, scattering parameter output will
be disabled for the simulation (though other quantities of interest are still
postprocessed). When lumped ports are present, the peak complex lumped port voltages and
currents computed for each excitation frequency are written to ASCII files named
`port-V.csv` and `port-I.csv`, respectively, Additionally, the surface current excitations
are written to `surface-I.csv`.

It is often the case that a user wants to compute the entire scattering matrix rather than
just a single column. In this case, each column can be computed in parallel by running
*Palace* multiple times. For example, consider the following short Python code which
modifies a base configuration file `config.json` to generate a complete 4x4 scattering
matrix by running 4 *Palace* simulations, each with 2 MPI processes:

```python
import json
import os
import subprocess

# Base configuration file
config_path = "config.json"

for i in range(4):
    # Prepare configuration file for simulation
    with open(config_path, "r") as f:
        config_json = json.loads(f.read())
    for port in config_json["Boundaries"]["LumpedPort"]:
        port["Excitation"] = (1+i == port["Index"])

    # Write new config file
    config_path_i = os.path.splitext(config_path)[0] + f"-{1+i}.json"
    with open(config_path_i, "w") as f:
        f.write(json.dumps(config_json))

    # Run Palace simulation (alternatively, use Popen and wait)
    subprocess.run(["palace", "-np", 2, config_path_i])
```

## Driven problems in the time domain

The previous simulation types describe simulations based on frequency domain formulations of
Maxwell's equations. Time domain simulations are also possible through the transient
simulation type:
[`config["Problem"]["Type"]: "Transient"`](../config/problem.md#config%5B%22Problem%22%5D).

Similar to the driven simulation type in the frequency domain, transient simulations involve
simulating the response of the system to a time-dependent excitation field specified at
lumped ports or surface current excitations in the model. The system is always started from
rest with zero initial conditions and time-integrated for a user specified duration, given
in nanoseconds. There are several available excitation types which define the time
dependence of the pulse or excitation waveform. These are specified under the
[`config["Solver"]["Transient"]`](../config/solver.md#solver%5B%22Transient%22%5D) object
using the `"Excitation"` keyword.

The time histories of the lumped port voltages and currents are postprocessed and
automatically written to ASCII files named `port-V.csv` and `port-I.csv`, respectively, in
the directory specified by
[`config["Problem"]["Output"]`](../config/problem.md#config%5B%22Problem%22%5D).
Additionally, surface current excitation time histories are written to `surface-I.csv`.

## Electrostatic problems

For electrostatic simulations,
([`config["Problem"]["Type"]: "Electrostatic"`](../config/problem.md#config%5B%22Problem%22%5D),
the user should specify a number of terminal boundaries
([`config["Boundaries"]["Terminal"]`](../config/boundaries.md#boundaries%5B%22Terminal%22%5D))
as well as boundaries which are grounded
([`config["Boundaries"]["Ground"]`](../config/boundaries.md#boundaries%5B%22Ground%22%5D)).
For each terminal, an electrostatic field is computed by assigning the terminal of interest
a positive unit voltage and all other terminals and grounded boundaries a zero voltage. The
resulting fields are then used to compute the Maxwell capacitance matrix and its inverse,
which are written to an ASCII file named `terminal-C.csv` and `terminal-Cinv.csv`,
respectively, in the directory specified by
[`config["Problem"]["Output"]`](../config/problem.md#config%5B%22Problem%22%5D). The mutual
capacitance matrix is also computed and written to `terminal-Cm.csv` in the same directory.

## Magnetostatic problems

For magnetostatic simulations,
([`config["Problem"]["Type"]: "Magnetostatic"`](../config/problem.md#config%5B%22Problem%22%5D),
the user should specify a number of source current boundaries. For each current source, a
magnetostatic field is computed by applying a unit current to the source index of interest,
leaving all other sources open with no excitation. Surfaces which are expected to carry
current should be labeled as perfectly conducting, which prescribes a zero magnetic flux, or
[magnetic insulation](https://doc.comsol.com/5.5/doc/com.comsol.help.comsol/comsol_ref_acdc.17.74.html),
boundary condition. The resulting fields are used to compute the inductance matrix and its
inverse, which are written to an ASCII file named `terminal-M.csv` and `terminal-Minv.csv`,
respectively, in the directory specified by
[`config["Problem"]["Output"]`](../config/problem.md#config%5B%22Problem%22%5D). A "mutual"
inductance matrix which has the same form as the mutual capacitance matrix (its entries are
based on current differences between ports rather than absolute currents) is computed and
written to `terminal-Mm.csv` in the same directory.
