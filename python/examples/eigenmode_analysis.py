#!/usr/bin/env python3
"""
Eigenmode Analysis Example for Palace Python Interface.

This script demonstrates how to set up and analyze eigenmode simulations
for cavity resonators and other electromagnetic structures.
"""

import json
import os
import sys

import matplotlib.pyplot as plt
import numpy as np

# Add parent directory to path
sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))

from palace.utils import create_basic_config, save_config


def create_cavity_eigenmode_config(
    mesh_file="cavity.msh", target_freq=5e9, num_modes=10
):
    """
    Create eigenmode configuration for cavity resonator analysis.

    Args:
        mesh_file: Path to mesh file
        target_freq: Target frequency for eigenmode search (Hz)
        num_modes: Number of eigenmodes to find

    Returns:
        Configuration dictionary
    """
    config = create_basic_config(mesh_file, "Eigenmode")

    # Configure eigenmode solver
    config["Solver"]["Eigenmode"] = {
        "Target": target_freq,
        "Tol": 1e-9,
        "MaxIts": 1000,
        "MaxSize": max(50, num_modes * 3),
        "Save": num_modes,
    }

    # Set up cavity boundary conditions (PEC walls)
    config["Model"]["Boundary"] = [
        {
            "Index": [1, 2, 3, 4, 5, 6],  # All boundaries
            "PEC": {},
        }
    ]

    # Material properties (vacuum/air)
    config["Model"]["Domain"] = [
        {"Index": 1, "Material": {"Permeability": 1.0, "Permittivity": 1.0}}
    ]

    # Post-processing options
    config["Model"]["PostProcessing"] = {
        "Energy": [
            {
                "Index": 1  # Compute energy in domain 1
            }
        ]
    }

    return config


def create_dielectric_resonator_config(
    mesh_file="dielectric_resonator.msh", epsilon_r=10.0, loss_tan=1e-4
):
    """
    Create configuration for dielectric resonator eigenmode analysis.

    Args:
        mesh_file: Path to mesh file
        epsilon_r: Relative permittivity of dielectric
        loss_tan: Loss tangent

    Returns:
        Configuration dictionary
    """
    config = create_basic_config(mesh_file, "Eigenmode")

    # Configure for dielectric resonator
    config["Solver"]["Eigenmode"].update(
        {
            "Target": 10e9,  # Higher frequency for dielectric resonators
            "Tol": 1e-9,
            "MaxIts": 1000,
            "MaxSize": 50,
            "Save": 10,
        }
    )

    # Material regions
    config["Model"]["Domain"] = [
        {
            "Index": 1,  # Dielectric region
            "Material": {
                "Permeability": 1.0,
                "Permittivity": epsilon_r,
                "LossTan": loss_tan,
            },
        },
        {
            "Index": 2,  # Air region
            "Material": {"Permeability": 1.0, "Permittivity": 1.0},
        },
    ]

    # Boundary conditions
    config["Model"]["Boundary"] = [
        {
            "Index": [1, 2, 3, 4],  # Side walls - PEC
            "PEC": {},
        },
        {
            "Index": [5, 6],  # Top/bottom - open (radiation)
            "Impedance": {
                "Rs": 377.0,  # Free space impedance
                "Xs": 0.0,
            },
        },
    ]

    return config


def analyze_eigenmode_results(result_file="domain-E.csv"):
    """
    Analyze eigenmode simulation results.

    Args:
        result_file: Path to results CSV file

    Returns:
        Analysis results dictionary
    """
    if not os.path.exists(result_file):
        # Create example data for demonstration
        print(f"Results file {result_file} not found. Creating example data...")

        # Simulate cavity resonator modes (TE/TM modes)
        modes = np.arange(1, 11)
        # Example: rectangular cavity a=2cm, b=1cm, c=1cm
        a, b, c = 0.02, 0.01, 0.01  # meters
        c0 = 3e8  # speed of light

        frequencies = []
        mode_types = []

        for m in range(3):
            for n in range(3):
                for p in range(2):
                    if m + n + p > 0:  # Exclude TEM mode
                        f = (c0 / (2 * np.pi)) * np.sqrt(
                            (m * np.pi / a) ** 2
                            + (n * np.pi / b) ** 2
                            + (p * np.pi / c) ** 2
                        )
                        frequencies.append(f)
                        if p == 0:
                            mode_types.append(f"TE_{m}{n}{p}")
                        else:
                            mode_types.append(f"TM_{m}{n}{p}")

        # Sort by frequency and take first 10
        sorted_indices = np.argsort(frequencies)[:10]
        frequencies = np.array(frequencies)[sorted_indices]
        mode_types = np.array(mode_types)[sorted_indices]

        # Add some realistic Q factors
        q_factors = 1000 + 500 * np.random.rand(len(frequencies))

        results = {
            "frequencies": frequencies,
            "q_factors": q_factors,
            "mode_types": mode_types,
            "num_modes": len(frequencies),
        }

    else:
        # Load actual results from file
        data = np.loadtxt(result_file, delimiter=",", skiprows=1)
        results = {
            "frequencies": data[:, 1],  # Assuming freq in column 1
            "q_factors": data[:, 2] if data.shape[1] > 2 else None,
            "mode_types": [f"Mode_{i + 1}" for i in range(len(data))],
            "num_modes": len(data),
        }

    return results


def plot_eigenmode_analysis(results, title="Eigenmode Analysis"):
    """
    Create comprehensive plots for eigenmode analysis.

    Args:
        results: Results dictionary from analyze_eigenmode_results
        title: Plot title
    """
    fig, axes = plt.subplots(2, 2, figsize=(15, 12))
    fig.suptitle(title, fontsize=16)

    frequencies = results["frequencies"]
    q_factors = results["q_factors"]
    mode_types = results["mode_types"]

    # Plot 1: Frequency spectrum
    axes[0, 0].bar(range(1, len(frequencies) + 1), frequencies / 1e9)
    axes[0, 0].set_xlabel("Mode Number")
    axes[0, 0].set_ylabel("Frequency (GHz)")
    axes[0, 0].set_title("Eigenmode Frequencies")
    axes[0, 0].grid(True, alpha=0.3)

    # Add mode labels
    for i, (f, mt) in enumerate(zip(frequencies, mode_types)):
        axes[0, 0].text(i + 1, f / 1e9 + 0.1, mt, rotation=45, ha="left", fontsize=8)

    # Plot 2: Q factors (if available)
    if q_factors is not None:
        axes[0, 1].bar(range(1, len(q_factors) + 1), q_factors)
        axes[0, 1].set_xlabel("Mode Number")
        axes[0, 1].set_ylabel("Q Factor")
        axes[0, 1].set_title("Quality Factors")
        axes[0, 1].grid(True, alpha=0.3)

        # Add Q values as text
        for i, q in enumerate(q_factors):
            axes[0, 1].text(
                i + 1, q + max(q_factors) * 0.01, f"{q:.0f}", ha="center", fontsize=8
            )
    else:
        axes[0, 1].text(
            0.5,
            0.5,
            "Q factors not available",
            ha="center",
            va="center",
            transform=axes[0, 1].transAxes,
        )
        axes[0, 1].set_title("Quality Factors")

    # Plot 3: Mode density
    axes[1, 0].hist(
        frequencies / 1e9,
        bins=min(10, len(frequencies) // 2 + 1),
        alpha=0.7,
        edgecolor="black",
    )
    axes[1, 0].set_xlabel("Frequency (GHz)")
    axes[1, 0].set_ylabel("Number of Modes")
    axes[1, 0].set_title("Mode Density")
    axes[1, 0].grid(True, alpha=0.3)

    # Plot 4: Q vs Frequency (if Q factors available)
    if q_factors is not None:
        scatter = axes[1, 1].scatter(
            frequencies / 1e9, q_factors, s=100, alpha=0.7, c=range(len(frequencies))
        )
        axes[1, 1].set_xlabel("Frequency (GHz)")
        axes[1, 1].set_ylabel("Q Factor")
        axes[1, 1].set_title("Q Factor vs Frequency")
        axes[1, 1].grid(True, alpha=0.3)

        # Add colorbar
        cbar = plt.colorbar(scatter, ax=axes[1, 1])
        cbar.set_label("Mode Number")

        # Add trend line
        if len(frequencies) > 2:
            z = np.polyfit(frequencies / 1e9, q_factors, 1)
            p = np.poly1d(z)
            axes[1, 1].plot(
                frequencies / 1e9,
                p(frequencies / 1e9),
                "r--",
                alpha=0.8,
                label=f"Trend: Q = {z[0]:.0f}f + {z[1]:.0f}",
            )
            axes[1, 1].legend()
    else:
        axes[1, 1].text(
            0.5,
            0.5,
            "Q factors not available",
            ha="center",
            va="center",
            transform=axes[1, 1].transAxes,
        )
        axes[1, 1].set_title("Q Factor vs Frequency")

    plt.tight_layout()
    plt.show()

    # Print summary
    print("\nEigenmode Analysis Summary:")
    print(f"{'=' * 50}")
    print(f"Number of modes found: {len(frequencies)}")
    print(
        f"Frequency range: {frequencies[0] / 1e9:.3f} - {frequencies[-1] / 1e9:.3f} GHz"
    )
    print(f"Mode spacing (avg): {np.mean(np.diff(frequencies)) / 1e6:.1f} MHz")

    if q_factors is not None:
        print(f"Q factor range: {np.min(q_factors):.0f} - {np.max(q_factors):.0f}")
        print(f"Average Q factor: {np.mean(q_factors):.0f}")

    print("\nMode Details:")
    print(f"{'Mode':<6} {'Type':<10} {'Freq (GHz)':<12} {'Q Factor':<10}")
    print(f"{'-' * 48}")
    for i, (f, mt) in enumerate(zip(frequencies, mode_types)):
        q_str = f"{q_factors[i]:.0f}" if q_factors is not None else "N/A"
        print(f"{i + 1:<6} {mt:<10} {f / 1e9:<12.3f} {q_str:<10}")


def main():
    """Main function demonstrating eigenmode analysis."""
    print("Palace Eigenmode Analysis Example")
    print("=" * 50)

    # Example 1: Cavity resonator
    print("\n1. Creating cavity resonator configuration...")
    cavity_config = create_cavity_eigenmode_config(
        mesh_file="rectangular_cavity.msh", target_freq=5e9, num_modes=10
    )

    save_config(cavity_config, "cavity_eigenmode.json")
    print("✓ Saved cavity_eigenmode.json")

    # Example 2: Dielectric resonator
    print("\n2. Creating dielectric resonator configuration...")
    dr_config = create_dielectric_resonator_config(
        mesh_file="dielectric_resonator.msh", epsilon_r=10.0, loss_tan=1e-4
    )

    save_config(dr_config, "dielectric_resonator_eigenmode.json")
    print("✓ Saved dielectric_resonator_eigenmode.json")

    # Example 3: Analyze results
    print("\n3. Analyzing eigenmode results...")
    results = analyze_eigenmode_results("domain-E.csv")

    # Create plots
    plot_eigenmode_analysis(results, "Rectangular Cavity Eigenmode Analysis")

    # Save results
    results_summary = {
        "frequencies_GHz": (results["frequencies"] / 1e9).tolist(),
        "q_factors": results["q_factors"].tolist()
        if results["q_factors"] is not None
        else None,
        "mode_types": results["mode_types"].tolist(),
        "analysis_notes": {
            "num_modes": results["num_modes"],
            "freq_range_GHz": [
                float(results["frequencies"][0] / 1e9),
                float(results["frequencies"][-1] / 1e9),
            ],
            "avg_mode_spacing_MHz": float(
                np.mean(np.diff(results["frequencies"])) / 1e6
            ),
        },
    }

    with open("eigenmode_analysis_results.json", "w") as f:
        json.dump(results_summary, f, indent=2)
    print("✓ Saved eigenmode_analysis_results.json")

    # Example 4: Parameter study
    print("\n4. Parameter study: Effect of dielectric constant...")
    epsilon_values = [1.0, 2.2, 4.0, 9.8, 12.0]
    material_names = ["Air", "PTFE", "Quartz", "Silicon", "GaAs"]

    param_results = []
    for eps_r, mat_name in zip(epsilon_values, material_names):
        config = create_dielectric_resonator_config(
            mesh_file="dielectric_resonator.msh", epsilon_r=eps_r, loss_tan=1e-4
        )

        config_file = f"dr_eps_{eps_r:.1f}.json"
        save_config(config, config_file)

        # Simulate frequency scaling: f ∝ 1/√εᵣ
        f0 = 10e9  # Reference frequency at εᵣ = 1
        scaled_freq = f0 / np.sqrt(eps_r)
        param_results.append((mat_name, eps_r, scaled_freq))

        print(
            f"  {mat_name} (εᵣ={eps_r}): f₀ ≈ {scaled_freq / 1e9:.2f} GHz → {config_file}"
        )

    # Plot parameter study results
    materials, eps_vals, freqs = zip(*param_results)

    plt.figure(figsize=(12, 5))

    plt.subplot(1, 2, 1)
    plt.bar(materials, [f / 1e9 for f in freqs])
    plt.ylabel("Fundamental Frequency (GHz)")
    plt.title("Frequency vs Material")
    plt.xticks(rotation=45)
    plt.grid(True, alpha=0.3)

    plt.subplot(1, 2, 2)
    plt.loglog(eps_vals, [f / 1e9 for f in freqs], "bo-", markersize=8)
    plt.xlabel("Relative Permittivity")
    plt.ylabel("Fundamental Frequency (GHz)")
    plt.title("Frequency Scaling: f ∝ 1/√εᵣ")
    plt.grid(True)

    # Add theoretical line
    eps_theory = np.logspace(0, 1.2, 50)
    f_theory = (f0 / np.sqrt(eps_theory)) / 1e9
    plt.loglog(eps_theory, f_theory, "r--", alpha=0.7, label="f ∝ 1/√εᵣ")
    plt.legend()

    plt.tight_layout()
    plt.savefig("dielectric_parameter_study.png", dpi=150, bbox_inches="tight")
    plt.show()

    print("\n✓ Parameter study completed")
    print("✓ Plot saved as dielectric_parameter_study.png")

    print(f"\n{'=' * 50}")
    print("Eigenmode analysis example completed!")
    print("Generated files:")
    print("  - cavity_eigenmode.json")
    print("  - dielectric_resonator_eigenmode.json")
    print("  - dr_eps_*.json (parameter sweep configs)")
    print("  - eigenmode_analysis_results.json")
    print("  - dielectric_parameter_study.png")


if __name__ == "__main__":
    main()
