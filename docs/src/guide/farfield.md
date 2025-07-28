```@raw html
<!--- Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved. --->
<!--- SPDX-License-Identifier: Apache-2.0 --->
```

# Extracting Fields in the Radiation Zone

When run in the [`Driven` mode](problem.md), *Palace* can extrapolate fields
from the near-field region (being simulated) to the far-field zone. This
capability can be used to study the radiative properties of a system.

The mathematical details on how this is accomplished are available in the
[reference](../reference.md#Far-field-extraction). The key points are as follows:

  - Computing fields requires evaluation of surface integrals, so users must
    specify a surface of integration (more on this later).
  - For the result to be physically accurate, it is important to properly model
    the propagation of waves to infinity. This can be accomplished by enclosing
    the system inside a sphere or box and applying [`config["Boundaries"]["Absorbing"]` boundary
    conditions](../config/boundaries.md#boundaries%5B%22Absorbing%22%5D).
  - The result is provided as complex vectors ``r \mathbf{E}(\theta, \phi)``, where
    ``(\theta, \phi)`` identify a point on sphere at infinite distance and the
    result is defined up to a global phase.

*Palace* outputs ``r\mathbf{E}`` because this is a well-defined, finite quantity
(``\mathbf{E}`` itself goes to zero at infinity as 1/r). You can also think of
this as the electric field measured at a distance of one meter. With this
output, you can immediately compute various quantities of interest. For
instance, ``|r\mathbf{E}|^2`` gives the relative radiative power.

!!! warning "Limitations"
    
    *Palace* only supports propagation of fields to infinity when the integration
    surface
    
     1. is an external boundary
     2. does not cross anisotropic materials

## Setup

A typical setup consists of starting from the system under consideration and
enclosing the system inside an outer boundary (typically a sphere or a box).
Then, we use this boundary for the Absorbing boundary conditions and as the
surface for the integration.

!!! tip
    
    For best accuracy, it is a good idea to make sure that this outer boundary is
    meshed finely enough to resolve the expected wavelength.

Next, we need to activate this feature. To do so, we look at the `"FarField"`
section under `"Postprocessing"` in `"Boundaries"`. Here, we need to specify the
identifier of the integration surface in `"Attributes"` and specify how we want
to sample the sphere at infinity.

The simplest way to do this is by setting `"NSample"`, so that the far-field
sphere is uniformly sampled with `NSample` points. (Note that uniform on a
sphere means that the polar angle ``\theta`` of the sampled points is not
uniformly distributed, to avoid bunching of points on the poles.)

Often, however, uniform sampling on the far-field sphere is not the best idea to
accurately capture the radiative properties because the radiation properties are
highly directional. In this case, you can specify at what angles ``\theta, \phi``
you want to evaluate the fields. This is done by passing 2-vectors to the
`"ThetaPhis"` array with the angular coordinates of your choosing (in degrees).

This can be combined with scripts to target specific regions. For instance, the
following Python code produces the required `"ThetaPhis"` section to sample over
the xy plane:

```python
import json
dphi = 1  # degree
print(json.dumps({"ThetaPhis": [[90.0, p] for p in range(0, 361, dphi)]}))
```

Both `"NSample"` and `"ThetaPhis"` can be specified simultaneously and the results
will be combined and duplicates removed.

## Output

Once a simulation is run, *Palace* generates a CSV file named
`farfield-E-<frequency>.csv` for each frequency point, where `<frequency>` is
encoded in scientific notation in units of GHz (e.g., `farfield-E-1.000e+00.csv`
for 1 GHz).

!!! tip "Extracting the frequency from the filename"
    
    The following regex pattern matches the filenames and can be used to extract the frequency:
    
    ```
    farfield-E-([+-]?\d+\.\d+e[+-]?\d+)\.csv
    ```

The CSV contains a header line that describes the columns, and one row for each
``(\theta, \phi)`` pair. The columns are:

  - `theta (deg.)` : Polar angle in degrees
  - `phi (deg.)` : Azimuthal angle in degrees
  - `r*Re{E_x} (V)`, `r*Im{E_x} (V)` : Real and imaginary parts of ``r E_x`` component
  - `r*Re{E_y} (V)`, `r*Im{E_y} (V)` : Real and imaginary parts of ``r E_y`` component
  - `r*Re{E_z} (V)`, `r*Im{E_z} (V)` : Real and imaginary parts of ``r E_z`` component

All field values are in SI units.

For example, this Python code prints the maxima for all frequencies in the
`postpro` directory:

```python
import re, pathlib, pandas

postpro_dir = pathlib.Path("postpro")
pattern = r'farfield-E-([+-]?\d+\.\d+e[+-]?\d+)\.csv'

for csv_file in postpro_dir.glob("farfield-E-*.csv"):
    match_ = re.search(pattern, csv_file.name)
    if match_:
        frequency = float(match_.group(1))
        df = pandas.read_csv(csv_file)
        df.columns = df.columns.str.strip()

        # Calculate field magnitude |E|
        E_mag = (df['r*Re{E_x} (V)']**2 + df['r*Im{E_x} (V)']**2 + 
                 df['r*Re{E_y} (V)']**2 + df['r*Im{E_y} (V)']**2 + 
                 df['r*Re{E_z} (V)']**2 + df['r*Im{E_z} (V)']**2)**0.5 
        
        max_field = E_mag.max()
        
        print(f"Frequency: {frequency:.3e} Hz")
        print(f"Max |rE|: {max_field:.3e} V")
```

See also the [antenna example](../examples/antenna.md) for more complex
postprocessing and for a script that directly plots antenna patterns. In
particular, the example includes a Python script that converts CSV files to PVD
files that can be immediately visualized with
[ParaView](https://www.paraview.org/).
