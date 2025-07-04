#!/usr/bin/env python3
"""
Basic usage examples for Palace Python interface.

This script demonstrates the fundamental ways to interact with Palace
through the Python API.
"""

import os
import sys

import numpy as np

# Add the parent directory to the path to import palace
try:
    # Try relative import first
    sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))
    from palace import PalaceSolver, run_palace
    from palace.utils import create_basic_config, read_csv_results, save_config
except ImportError:
    # Fallback for installed package
    from palace import PalaceSolver, run_palace
    from palace.utils import create_basic_config, read_csv_results, save_config

# Try to import matplotlib, handle gracefully if not available
try:
    import matplotlib.pyplot as plt

    HAS_MATPLOTLIB = True
except ImportError:
    print("Warning: matplotlib not available. Plotting will be disabled.")
    HAS_MATPLOTLIB = False

    # Create a dummy plt module for compatibility
    class DummyPlt:
        @staticmethod
        def figure(*args, **kwargs):
            pass

        @staticmethod
        def plot(*args, **kwargs):
            pass

        @staticmethod
        def xlabel(*args, **kwargs):
            pass

        @staticmethod
        def ylabel(*args, **kwargs):
            pass

        @staticmethod
        def title(*args, **kwargs):
            pass

        @staticmethod
        def grid(*args, **kwargs):
            pass

        @staticmethod
        def legend(*args, **kwargs):
            pass

        @staticmethod
        def show(*args, **kwargs):
            print("(Plot would be displayed here if matplotlib was available)")

        @staticmethod
        def savefig(*args, **kwargs):
            print(
                f"(Plot would be saved as {args[0] if args else 'plot.png'} if matplotlib was available)"
            )

        @staticmethod
        def close(*args, **kwargs):
            pass

        @staticmethod
        def subplot(*args, **kwargs):
            pass

        @staticmethod
        def tight_layout(*args, **kwargs):
            pass

    plt = DummyPlt()


def example_1_simple_run():
    """Example 1: Simple simulation run using existing config."""
    print("=" * 60)
    print("Example 1: Running a simple Palace simulation")
    print("=" * 60)

    # Check if we have example configs available
    config_file = "../../examples/cylinder/cavity_pec.json"
    if not os.path.exists(config_file):
        print(f"Config file not found: {config_file}")
        print("This example requires the Palace repository examples.")
        return

    try:
        # Method 1: Using the convenience function
        print("Running simulation using convenience function...")
        result = run_palace(config_file, num_procs=2)

        if result.returncode == 0:
            print("✓ Simulation completed successfully!")
        else:
            print(f"✗ Simulation failed with exit code: {result.returncode}")
            print("STDERR:", result.stderr)

    except Exception as e:
        print(f"Error running simulation: {e}")


def example_2_solver_class():
    """Example 2: Using the PalaceSolver class for more control."""
    print("\n" + "=" * 60)
    print("Example 2: Using PalaceSolver class")
    print("=" * 60)

    try:
        # Initialize solver
        solver = PalaceSolver()
        print(f"Palace executable found at: {solver.executable}")

        # Check if we have an example config
        config_file = "../../examples/rings/rings.json"
        if not os.path.exists(config_file):
            print(f"Config file not found: {config_file}")
            return

        # Validate configuration first
        if solver.validate_config(config_file):
            print("✓ Configuration file is valid")
        else:
            print("✗ Configuration file has issues")
            return

        # Run simulation with custom settings
        print("Running simulation with custom output directory...")
        result = solver.run(config_file, output_dir="./output_rings", num_procs=1)

        if result.returncode == 0:
            print("✓ Simulation completed successfully!")
            print("Results saved to: ./output_rings")
        else:
            print(f"✗ Simulation failed: {result.stderr}")

    except Exception as e:
        print(f"Error: {e}")


def example_3_create_config():
    """Example 3: Creating a configuration programmatically."""
    print("\n" + "=" * 60)
    print("Example 3: Creating configuration programmatically")
    print("=" * 60)

    # Look for an example mesh file
    mesh_files = [
        "../../examples/coaxial/mesh/coaxial.msh",
        "../../examples/cpw/mesh/cpw_coax.msh",
        "../../test/unit/mesh/star-tri.mesh",
    ]

    mesh_file = None
    for mf in mesh_files:
        if os.path.exists(mf):
            mesh_file = mf
            break

    if not mesh_file:
        print("No example mesh files found. Creating a placeholder config.")
        mesh_file = "path/to/your/mesh.msh"

    # Create different types of configurations
    configs = {
        "eigenmode": create_basic_config(mesh_file, "Eigenmode"),
        "driven": create_basic_config(mesh_file, "Driven"),
    }

    # Customize the eigenmode config
    configs["eigenmode"]["Solver"]["Eigenmode"]["Target"] = 2.5e9  # 2.5 GHz
    configs["eigenmode"]["Solver"]["Eigenmode"]["MaxSize"] = 20

    # Add material properties
    configs["eigenmode"]["Model"]["Domain"] = [
        {
            "Index": 1,
            "Material": {
                "Permeability": 1.0,
                "Permittivity": 4.0,  # Relative permittivity
                "LossTan": 1e-4,  # Loss tangent
            },
        }
    ]

    # Save configurations
    for name, config in configs.items():
        output_file = f"example_{name}_config.json"
        save_config(config, output_file)
        print(f"✓ Created {name} configuration: {output_file}")

        # Show key settings
        if "Eigenmode" in config.get("Solver", {}):
            target_freq = config["Solver"]["Eigenmode"]["Target"]
            print(f"  Target frequency: {target_freq / 1e9:.1f} GHz")
        elif "Driven" in config.get("Solver", {}):
            min_freq = config["Solver"]["Driven"]["MinFreq"]
            max_freq = config["Solver"]["Driven"]["MaxFreq"]
            print(f"  Frequency range: {min_freq / 1e9:.1f} - {max_freq / 1e9:.1f} GHz")


def example_4_postprocess_results():
    """Example 4: Post-processing simulation results."""
    print("\n" + "=" * 60)
    print("Example 4: Post-processing results")
    print("=" * 60)

    # Look for example result files
    result_dirs = [
        "../../test/examples/ref/coaxial/matched",
        "../../test/examples/ref/cpw/lumped_uniform",
        "./output_rings",
    ]

    found_results = False

    for result_dir in result_dirs:
        if os.path.exists(result_dir):
            print(f"Processing results from: {result_dir}")

            # Look for CSV files
            csv_files = [f for f in os.listdir(result_dir) if f.endswith(".csv")]

            for csv_file in csv_files[:2]:  # Process first 2 CSV files
                csv_path = os.path.join(result_dir, csv_file)
                try:
                    freq, data = read_csv_results(csv_path)
                    print(f"  ✓ {csv_file}: {len(freq)} data points")

                    if len(data.shape) > 1 and data.shape[1] > 0:
                        print(
                            f"    Frequency range: {freq[0] / 1e9:.2f} - {freq[-1] / 1e9:.2f} GHz"
                        )
                        print(f"    Data columns: {data.shape[1]}")

                        # Simple plot if matplotlib is available
                        try:
                            plt.figure(figsize=(10, 6))
                            for i in range(
                                min(3, data.shape[1])
                            ):  # Plot first 3 columns
                                plt.plot(
                                    freq / 1e9,
                                    np.abs(data[:, i]),
                                    label=f"Column {i + 1}",
                                )
                            plt.xlabel("Frequency (GHz)")
                            plt.ylabel("Magnitude")
                            plt.title(f"Results from {csv_file}")
                            plt.legend()
                            plt.grid(True)

                            output_plot = f"plot_{csv_file.replace('.csv', '.png')}"
                            plt.savefig(output_plot, dpi=150, bbox_inches="tight")
                            plt.close()
                            print(f"    ✓ Plot saved: {output_plot}")

                        except ImportError:
                            print("    (matplotlib not available for plotting)")

                except Exception as e:
                    print(f"  ✗ Error reading {csv_file}: {e}")

            found_results = True
            break

    if not found_results:
        print("No result files found. Run a simulation first to generate results.")

        # Create some dummy data for demonstration
        print("\nCreating example data for demonstration...")
        freq = np.linspace(1e9, 10e9, 100)  # 1-10 GHz
        s11_mag = np.abs(1 / (1 + 1j * freq / 5e9))  # Simple resonance
        s11_phase = np.angle(1 / (1 + 1j * freq / 5e9))

        # Save example data
        data = np.column_stack([freq, s11_mag, s11_phase])
        np.savetxt(
            "example_s11.csv",
            data,
            delimiter=",",
            header="Frequency(Hz),S11_Magnitude,S11_Phase(rad)",
        )

        # Plot example data
        try:
            plt.figure(figsize=(12, 5))

            plt.subplot(1, 2, 1)
            plt.plot(freq / 1e9, 20 * np.log10(s11_mag))
            plt.xlabel("Frequency (GHz)")
            plt.ylabel("|S11| (dB)")
            plt.title("S11 Magnitude")
            plt.grid(True)

            plt.subplot(1, 2, 2)
            plt.plot(freq / 1e9, s11_phase * 180 / np.pi)
            plt.xlabel("Frequency (GHz)")
            plt.ylabel("S11 Phase (degrees)")
            plt.title("S11 Phase")
            plt.grid(True)

            plt.tight_layout()
            plt.savefig("example_s11_plot.png", dpi=150, bbox_inches="tight")
            plt.close()
            print("✓ Example plot saved: example_s11_plot.png")

        except ImportError:
            print("matplotlib not available for plotting")


def main():
    """Run all examples."""
    print("Palace Python Interface Examples")
    print("=" * 60)

    # Check if Palace is available
    try:
        solver = PalaceSolver()
        print(f"Palace executable found: {solver.executable}")
    except Exception as e:
        print(f"Warning: Palace not found in PATH: {e}")
        print("Some examples may not work without Palace installed.")

    print()

    # Run examples
    example_1_simple_run()
    example_2_solver_class()
    example_3_create_config()
    example_4_postprocess_results()

    print("\n" + "=" * 60)
    print("Examples completed!")
    print("Check the generated files:")
    print("  - example_*_config.json: Configuration files")
    print("  - example_*.csv: Sample data files")
    print("  - *.png: Result plots (if matplotlib available)")
    print("=" * 60)


if __name__ == "__main__":
    main()
