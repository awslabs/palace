#!/usr/bin/env pvpython

# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""Extract farfield data from Palace simulations.

This script interpolates electric field components (E_real and E_imag) onto a
spherical surface and exports the results to CSV format.

Usage:
    pvpython extract_farfield.py --radius 0.75

Requires ParaView Python bindings (pvpython).

"""

import argparse
import csv
import math
import os
import sys

try:
    import paraview.simple as pv
except ImportError as e:
    print(f"ERROR: ParaView Python module not found ({e}).")
    print("Run with 'pvpython extract_farfield.py' instead of 'python'.")
    sys.exit(1)


def interpolate_on_sphere(
    data, radius: float, resolution: int = 360
) -> list[dict[str, float]]:
    """Interpolate E-field components on a sphere of given radius with the given azimuthal resolution.

    The meriodional resolution is assumed to be half of the azimuthal one.

    Returns list of dictionaries with theta, phi, and field components.
    """

    # Create a ParaView sphere
    sphere = pv.Sphere()
    sphere.Radius = radius
    sphere.Center = [0.0, 0.0, 0.0]
    sphere.PhiResolution = resolution  # Around equator
    sphere.ThetaResolution = resolution // 2  # From pole to pole
    sphere.UpdatePipeline()

    data.UpdatePipeline()

    # Interpolate data onto the sphere
    interpolator = pv.ResampleWithDataset(SourceDataArrays=data, DestinationMesh=sphere)
    interpolator.UpdatePipeline()

    # Get interpolated data
    data = pv.servermanager.Fetch(interpolator)
    points = data.GetPoints()

    # Extract field components and convert to spherical coordinates
    results = []
    for i in range(points.GetNumberOfPoints()):
        x, y, z = points.GetPoint(i)

        # Convert Cartesian to spherical coordinates
        r = math.sqrt(x * x + y * y + z * z)
        theta = math.degrees(math.acos(z / r))
        phi = math.degrees(math.atan2(y, x))
        if phi < 0:
            phi += 360

        # Get field components from vector fields
        e_real_x = get_field_value(data, "E_real", 0, i)
        e_real_y = get_field_value(data, "E_real", 1, i)
        e_real_z = get_field_value(data, "E_real", 2, i)
        e_imag_x = get_field_value(data, "E_imag", 0, i)
        e_imag_y = get_field_value(data, "E_imag", 1, i)
        e_imag_z = get_field_value(data, "E_imag", 2, i)

        results.append(
            {
                "theta": theta,
                "phi": phi,
                "E_real_x": e_real_x,
                "E_real_y": e_real_y,
                "E_real_z": e_real_z,
                "E_imag_x": e_imag_x,
                "E_imag_y": e_imag_y,
                "E_imag_z": e_imag_z,
            }
        )

    return results


def get_field_value(data, field_name: str, component: int, point_index: int) -> float:
    """Extract the vector field component values"""
    array = data.GetPointData().GetArray(field_name)
    if array:
        return array.GetTuple(point_index)[component]

    # array is None
    print(f"Field '{field_name}' not found!")
    sys.exit(1)


def save_to_csv(data: list[dict[str, float]], filename: str) -> None:
    """Save interpolated data to CSV file following Palace's output format"""

    with open(filename, "w", newline="") as f:
        writer = csv.writer(f)
        headers = [
            "theta (deg.)",
            "phi (deg.)",
            "Re{E_x} (V/m)",
            "Re{E_y} (V/m)",
            "Re{E_z} (V/m)",
            "Im{E_x} (V/m)",
            "Im{E_y} (V/m)",
            "Im{E_z} (V/m)",
        ]
        writer.writerow([f"{header:>20}" for header in headers])
        for row in data:
            writer.writerow(
                [
                    f"{row['theta']:20.8e}",
                    f"{row['phi']:20.8e}",
                    f"{row['E_real_x']:20.8e}",
                    f"{row['E_real_y']:20.8e}",
                    f"{row['E_real_z']:20.8e}",
                    f"{row['E_imag_x']:20.8e}",
                    f"{row['E_imag_y']:20.8e}",
                    f"{row['E_imag_z']:20.8e}",
                ]
            )


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Extract data defined on a sphere of a specified radius."
    )
    parser.add_argument(
        "--postpro",
        default="postpro",
        help="Path to postprocessing folder (default: postpro)",
    )
    parser.add_argument(
        "--radius",
        type=float,
        default=0.75,
        help="Radius fraction of maximum domain extent (default: 0.75)",
    )
    parser.add_argument(
        "--resolution",
        type=int,
        default=360,
        help="Interpolation resolution (default: 360)",
    )
    parser.add_argument(
        "--output",
        help="Output filename (default: farfield_sphere_r{radius_fraction}.csv)",
    )

    args = parser.parse_args()

    if not (0 < args.radius <= 1):
        print(f"Radius has to be between 0 and 1 (passed {args.radius})")
        sys.exit(1)

    # Set default filename if not provided
    args.output = args.output or f"farfield_r{args.radius:.2f}.csv"

    pvd_path = os.path.join(args.postpro, "paraview", "driven", "driven.pvd")

    if not os.path.exists(pvd_path):
        print(f"Error: PVD file not found at {pvd_path}")
        sys.exit(1)

    # Load the Palace data
    reader = pv.PVDReader(FileName=pvd_path)
    reader.UpdatePipeline()

    # Select the first frame (which contains the actual data)
    data = pv.ExtractTimeSteps()
    data.Input = reader
    data.TimeStepIndices = [0]
    data.UpdatePipeline()

    # Calculate domain bounds and actual radius
    bounds = data.GetDataInformation().GetBounds()
    max_extent = max(
        bounds[1] - bounds[0], bounds[3] - bounds[2], bounds[5] - bounds[4]
    )
    actual_radius = args.radius * max_extent * 0.5

    print(f"Loading data from: {pvd_path}")
    print(
        f"Domain bounds: x=[{bounds[0]:.2f}, {bounds[1]:.2f}], y=[{bounds[2]:.2f}, {bounds[3]:.2f}], z=[{bounds[4]:.2f}, {bounds[5]:.2f}]"
    )
    print(f"Maximum domain extent: {max_extent:.2f}")
    print(f"Radius fraction: {args.radius}")
    print(f"Actual interpolation radius: {actual_radius:.2f}")
    print(f"Interpolation resolution: {args.resolution}")

    # Perform interpolation
    print("Interpolating... ", end="")
    sphere_data = interpolate_on_sphere(data, actual_radius, args.resolution)
    save_to_csv(sphere_data, args.output)

    print("Completed!")
    print(f"Saved file {args.output} ({len(sphere_data)} points)")


if __name__ == "__main__":
    main()
