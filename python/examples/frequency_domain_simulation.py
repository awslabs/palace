#!/usr/bin/env python3
"""
Frequency Domain Simulation Example for Palace Python Interface.

This script demonstrates how to set up and analyze frequency domain driven
simulations for S-parameter extraction, filter analysis, and antenna modeling.
"""

import json
import os
import sys

import matplotlib.pyplot as plt
import numpy as np

# Add parent directory to path
sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))

from palace.utils import create_basic_config, read_csv_results, save_config


def create_two_port_network_config(
    mesh_file="two_port.msh",
    freq_min=1e9,
    freq_max=10e9,
    freq_step=0.1e9,
    port_impedance=50.0,
):
    """
    Create configuration for two-port network S-parameter extraction.

    Args:
        mesh_file: Path to mesh file
        freq_min: Minimum frequency (Hz)
        freq_max: Maximum frequency (Hz)
        freq_step: Frequency step (Hz)
        port_impedance: Port characteristic impedance (Ohms)

    Returns:
        Configuration dictionary
    """
    config = create_basic_config(mesh_file, "Driven")

    # Configure frequency sweep
    config["Solver"]["Driven"] = {
        "MinFreq": freq_min,
        "MaxFreq": freq_max,
        "FreqStep": freq_step,
        "SaveStep": max(1, int(freq_step / (0.1e9))),  # Save every ~100 MHz
        "Restart": 0,
    }

    # Define ports
    config["Model"]["Boundary"] = [
        {
            "Index": 1,  # Port 1 (excitation)
            "LumpedPort": {"R": port_impedance, "L": 0.0, "C": 0.0, "Excitation": True},
        },
        {
            "Index": 2,  # Port 2 (matched load)
            "LumpedPort": {"R": port_impedance, "L": 0.0, "C": 0.0},
        },
    ]

    # Material properties
    config["Model"]["Domain"] = [
        {
            "Index": 1,  # Dielectric substrate
            "Material": {
                "Permeability": 1.0,
                "Permittivity": 4.3,  # FR-4 substrate
                "LossTan": 0.02,
            },
        },
        {
            "Index": 2,  # Air regions
            "Material": {"Permeability": 1.0, "Permittivity": 1.0},
        },
    ]

    # Add surface impedance for conductor losses
    config["Model"]["Boundary"].append(
        {
            "Index": [10, 11, 12],  # Conductor surfaces
            "Impedance": {
                "Rs": 0.01,  # Surface resistance (Ohms)
                "Ls": 0.0,  # Surface inductance
                "Cs": 0.0,  # Surface capacitance
            },
        }
    )

    # Post-processing options
    config["Model"]["PostProcessing"] = {
        "Probe": [
            {
                "Index": 1,
                "Center": [0.0, 0.0, 0.001],  # E-field probe
            }
        ]
    }

    return config


def create_filter_analysis_config(
    mesh_file="bandpass_filter.msh", center_freq=5e9, bandwidth_factor=0.2
):
    """
    Create configuration for microwave filter analysis.

    Args:
        mesh_file: Path to mesh file
        center_freq: Center frequency (Hz)
        bandwidth_factor: Fractional bandwidth for frequency sweep

    Returns:
        Configuration dictionary
    """
    # Calculate frequency range based on center frequency and bandwidth
    freq_span = center_freq * bandwidth_factor
    freq_min = center_freq - freq_span
    freq_max = center_freq + freq_span
    freq_step = freq_span / 100  # 100 points across bandwidth

    config = create_two_port_network_config(
        mesh_file=mesh_file,
        freq_min=freq_min,
        freq_max=freq_max,
        freq_step=freq_step,
        port_impedance=50.0,
    )

    # Add specific post-processing for filter analysis
    config["Model"]["PostProcessing"]["Surface"] = [
        {
            "Index": [20, 21, 22],  # Resonator surfaces
            "Type": "Capacitive",
        }
    ]

    return config


def create_antenna_simulation_config(
    mesh_file="antenna.msh", freq_center=2.4e9, freq_bandwidth=0.5e9
):
    """
    Create configuration for antenna simulation with far-field analysis.

    Args:
        mesh_file: Path to mesh file
        freq_center: Center frequency (Hz)
        freq_bandwidth: Bandwidth around center frequency (Hz)

    Returns:
        Configuration dictionary
    """
    freq_min = freq_center - freq_bandwidth / 2
    freq_max = freq_center + freq_bandwidth / 2
    freq_step = freq_bandwidth / 50

    config = create_basic_config(mesh_file, "Driven")

    # Configure frequency sweep
    config["Solver"]["Driven"] = {
        "MinFreq": freq_min,
        "MaxFreq": freq_max,
        "FreqStep": freq_step,
        "SaveStep": 1,
    }

    # Antenna feed
    config["Model"]["Boundary"] = [
        {
            "Index": 1,  # Feed point
            "LumpedPort": {"R": 50.0, "Excitation": True},
        }
    ]

    # Far-field boundary (absorbing boundary condition)
    config["Model"]["Boundary"].append(
        {
            "Index": [10, 11, 12, 13, 14, 15],  # Far-field boundaries
            "Absorbing": {
                "Order": 2  # Second-order ABC
            },
        }
    )

    # Add far-field post-processing
    config["Model"]["PostProcessing"]["FarField"] = {
        "Theta": {"Min": 0, "Max": 180, "Step": 5},
        "Phi": {"Min": 0, "Max": 360, "Step": 10},
    }

    return config


def analyze_s_parameters(s_param_file="port-S.csv"):
    """
    Analyze S-parameter results from Palace simulation.

    Args:
        s_param_file: Path to S-parameter CSV file

    Returns:
        Analysis results dictionary
    """
    if not os.path.exists(s_param_file):
        print(f"S-parameter file {s_param_file} not found. Creating example data...")

        # Create example S-parameter data for a bandpass filter
        freq = np.linspace(4e9, 6e9, 201)  # 4-6 GHz
        f0 = 5e9  # Center frequency
        Q = 50  # Quality factor

        # Simple bandpass filter model
        s21_mag = 1 / np.sqrt(1 + (2 * Q * (freq - f0) / f0) ** 2)
        s21_phase = -np.arctan(2 * Q * (freq - f0) / f0)

        s11_mag = np.sqrt(1 - s21_mag**2)
        s11_phase = s21_phase + np.pi

        # Add some ripple for realism
        ripple = 0.1 * np.sin(20 * np.pi * (freq - freq[0]) / (freq[-1] - freq[0]))
        s21_mag = s21_mag * (1 + ripple)
        s11_mag = s11_mag * (1 - ripple * 0.5)

        # Convert to complex
        s11 = s11_mag * np.exp(1j * s11_phase)
        s21 = s21_mag * np.exp(1j * s21_phase)
        s12 = s21  # Reciprocal network
        s22 = s11  # Symmetric network

        # Save example data
        s_data = np.column_stack(
            [
                freq,
                np.real(s11),
                np.imag(s11),
                np.real(s12),
                np.imag(s12),
                np.real(s21),
                np.imag(s21),
                np.real(s22),
                np.imag(s22),
            ]
        )

        np.savetxt(
            "example_s_parameters.csv",
            s_data,
            delimiter=",",
            header="Freq(Hz),Re(S11),Im(S11),Re(S12),Im(S12),Re(S21),Im(S21),Re(S22),Im(S22)",
        )

        s_param_file = "example_s_parameters.csv"

    # Load S-parameter data
    try:
        freq, data = read_csv_results(s_param_file)

        # Parse S-parameters (assuming Re/Im format)
        if data.shape[1] >= 8:
            s11 = data[:, 0] + 1j * data[:, 1]
            s12 = data[:, 2] + 1j * data[:, 3]
            s21 = data[:, 4] + 1j * data[:, 5]
            s22 = data[:, 6] + 1j * data[:, 7]
        else:
            raise ValueError("Insufficient columns in S-parameter file")

    except Exception as e:
        print(f"Error reading S-parameters: {e}")
        return None

    # Analysis calculations
    results = {
        "frequency": freq,
        "s11": s11,
        "s12": s12,
        "s21": s21,
        "s22": s22,
        "return_loss_db": 20 * np.log10(np.abs(s11)),
        "insertion_loss_db": 20 * np.log10(np.abs(s21)),
        "vswr": (1 + np.abs(s11)) / (1 - np.abs(s11)),
    }

    # Find key metrics
    min_s11_idx = np.argmin(np.abs(s11))
    max_s21_idx = np.argmax(np.abs(s21))

    results["best_match_freq"] = freq[min_s11_idx]
    results["best_match_s11_db"] = results["return_loss_db"][min_s11_idx]
    results["max_transmission_freq"] = freq[max_s21_idx]
    results["max_transmission_s21_db"] = results["insertion_loss_db"][max_s21_idx]

    # Bandwidth calculation (3dB points)
    s21_db = results["insertion_loss_db"]
    max_s21_db = np.max(s21_db)
    bw_3db_level = max_s21_db - 3

    bw_indices = np.where(s21_db >= bw_3db_level)[0]
    if len(bw_indices) > 0:
        bw_freq_low = freq[bw_indices[0]]
        bw_freq_high = freq[bw_indices[-1]]
        results["bandwidth_3db"] = bw_freq_high - bw_freq_low
        results["fractional_bandwidth"] = (
            results["bandwidth_3db"] / results["max_transmission_freq"]
        )
    else:
        results["bandwidth_3db"] = None
        results["fractional_bandwidth"] = None

    return results


def plot_s_parameter_analysis(results, title="S-Parameter Analysis"):
    """
    Create comprehensive S-parameter plots.

    Args:
        results: Results dictionary from analyze_s_parameters
        title: Plot title
    """
    if results is None:
        print("No results to plot")
        return

    freq = results["frequency"]
    s11 = results["s11"]
    s21 = results["s21"]

    fig, axes = plt.subplots(2, 3, figsize=(18, 12))
    fig.suptitle(title, fontsize=16)

    # S11 magnitude
    axes[0, 0].plot(freq / 1e9, results["return_loss_db"])
    axes[0, 0].set_xlabel("Frequency (GHz)")
    axes[0, 0].set_ylabel("|S11| (dB)")
    axes[0, 0].set_title("Return Loss")
    axes[0, 0].grid(True)
    axes[0, 0].axhline(-10, color="r", linestyle="--", alpha=0.7, label="-10 dB")
    axes[0, 0].legend()

    # S21 magnitude
    axes[0, 1].plot(freq / 1e9, results["insertion_loss_db"])
    axes[0, 1].set_xlabel("Frequency (GHz)")
    axes[0, 1].set_ylabel("|S21| (dB)")
    axes[0, 1].set_title("Insertion Loss")
    axes[0, 1].grid(True)
    axes[0, 1].axhline(-3, color="r", linestyle="--", alpha=0.7, label="-3 dB")
    axes[0, 1].legend()

    # VSWR
    axes[0, 2].plot(freq / 1e9, results["vswr"])
    axes[0, 2].set_xlabel("Frequency (GHz)")
    axes[0, 2].set_ylabel("VSWR")
    axes[0, 2].set_title("Voltage Standing Wave Ratio")
    axes[0, 2].grid(True)
    axes[0, 2].axhline(2, color="r", linestyle="--", alpha=0.7, label="VSWR = 2")
    axes[0, 2].set_ylim([1, min(10, np.max(results["vswr"]))])
    axes[0, 2].legend()

    # Smith chart (S11)
    theta = np.linspace(0, 2 * np.pi, 100)
    axes[1, 0].plot(np.cos(theta), np.sin(theta), "k-", alpha=0.3)  # Unit circle
    axes[1, 0].plot(np.real(s11), np.imag(s11), "b-", linewidth=2)
    axes[1, 0].plot(np.real(s11[0]), np.imag(s11[0]), "go", markersize=8, label="Start")
    axes[1, 0].plot(np.real(s11[-1]), np.imag(s11[-1]), "ro", markersize=8, label="End")
    axes[1, 0].set_xlabel("Real(S11)")
    axes[1, 0].set_ylabel("Imag(S11)")
    axes[1, 0].set_title("S11 Smith Chart")
    axes[1, 0].axis("equal")
    axes[1, 0].grid(True)
    axes[1, 0].legend()

    # Phase response
    axes[1, 1].plot(freq / 1e9, np.angle(s11) * 180 / np.pi, label="S11")
    axes[1, 1].plot(freq / 1e9, np.angle(s21) * 180 / np.pi, label="S21")
    axes[1, 1].set_xlabel("Frequency (GHz)")
    axes[1, 1].set_ylabel("Phase (degrees)")
    axes[1, 1].set_title("Phase Response")
    axes[1, 1].grid(True)
    axes[1, 1].legend()

    # Group delay
    phase_s21 = np.unwrap(np.angle(s21))
    group_delay = -np.gradient(phase_s21) / np.gradient(2 * np.pi * freq)
    axes[1, 2].plot(freq / 1e9, group_delay * 1e9)  # Convert to ns
    axes[1, 2].set_xlabel("Frequency (GHz)")
    axes[1, 2].set_ylabel("Group Delay (ns)")
    axes[1, 2].set_title("Group Delay")
    axes[1, 2].grid(True)

    plt.tight_layout()
    plt.show()

    # Print summary
    print("\nS-Parameter Analysis Summary:")
    print(f"{'=' * 50}")
    print(f"Frequency range: {freq[0] / 1e9:.2f} - {freq[-1] / 1e9:.2f} GHz")
    print(f"Best match frequency: {results['best_match_freq'] / 1e9:.3f} GHz")
    print(f"Best return loss: {results['best_match_s11_db']:.1f} dB")
    print(
        f"Max transmission frequency: {results['max_transmission_freq'] / 1e9:.3f} GHz"
    )
    print(f"Max transmission: {results['max_transmission_s21_db']:.1f} dB")

    if results["bandwidth_3db"] is not None:
        print(f"3-dB bandwidth: {results['bandwidth_3db'] / 1e6:.1f} MHz")
        print(f"Fractional bandwidth: {results['fractional_bandwidth'] * 100:.1f}%")

    # Quality metrics
    avg_return_loss = np.mean(results["return_loss_db"])
    avg_insertion_loss = np.mean(results["insertion_loss_db"])
    print(f"Average return loss: {avg_return_loss:.1f} dB")
    print(f"Average insertion loss: {avg_insertion_loss:.1f} dB")


def main():
    """Main function demonstrating frequency domain simulations."""
    print("Palace Frequency Domain Simulation Examples")
    print("=" * 60)

    # Example 1: Two-port network (general S-parameter extraction)
    print("\n1. Creating two-port network configuration...")
    two_port_config = create_two_port_network_config(
        mesh_file="microstrip_line.msh",
        freq_min=1e9,
        freq_max=10e9,
        freq_step=0.1e9,
        port_impedance=50.0,
    )

    save_config(two_port_config, "two_port_network.json")
    print("✓ Saved two_port_network.json")

    # Example 2: Microwave filter analysis
    print("\n2. Creating bandpass filter configuration...")
    filter_config = create_filter_analysis_config(
        mesh_file="bandpass_filter.msh", center_freq=5e9, bandwidth_factor=0.4
    )

    save_config(filter_config, "bandpass_filter.json")
    print("✓ Saved bandpass_filter.json")

    # Example 3: Antenna simulation
    print("\n3. Creating antenna simulation configuration...")
    antenna_config = create_antenna_simulation_config(
        mesh_file="patch_antenna.msh", freq_center=2.4e9, freq_bandwidth=0.5e9
    )

    save_config(antenna_config, "antenna_simulation.json")
    print("✓ Saved antenna_simulation.json")

    # Example 4: S-parameter analysis
    print("\n4. Analyzing S-parameter results...")
    s_param_results = analyze_s_parameters("port-S.csv")

    if s_param_results:
        plot_s_parameter_analysis(
            s_param_results, "Bandpass Filter S-Parameter Analysis"
        )

        # Save detailed results
        analysis_summary = {
            "best_match_freq_GHz": float(s_param_results["best_match_freq"] / 1e9),
            "best_return_loss_dB": float(s_param_results["best_match_s11_db"]),
            "max_transmission_freq_GHz": float(
                s_param_results["max_transmission_freq"] / 1e9
            ),
            "max_transmission_dB": float(s_param_results["max_transmission_s21_db"]),
            "bandwidth_3dB_MHz": float(s_param_results["bandwidth_3db"] / 1e6)
            if s_param_results["bandwidth_3db"]
            else None,
            "fractional_bandwidth_percent": float(
                s_param_results["fractional_bandwidth"] * 100
            )
            if s_param_results["fractional_bandwidth"]
            else None,
        }

        with open("s_parameter_analysis.json", "w") as f:
            json.dump(analysis_summary, f, indent=2)
        print("✓ Saved s_parameter_analysis.json")

    # Example 5: Frequency sweep optimization
    print("\n5. Frequency sweep optimization example...")

    # Demonstrate adaptive frequency sampling
    def adaptive_frequency_points(f_min, f_max, initial_points=51):
        """Generate adaptive frequency points with higher density near resonances."""
        # Start with uniform sampling
        f_uniform = np.linspace(f_min, f_max, initial_points)

        # Add extra points around suspected resonance (example)
        f_resonance = 5e9  # Expected resonance
        f_span = 0.5e9  # Span around resonance
        f_dense = np.linspace(f_resonance - f_span / 2, f_resonance + f_span / 2, 21)

        # Combine and sort
        f_combined = np.unique(np.concatenate([f_uniform, f_dense]))
        return f_combined

    adaptive_freqs = adaptive_frequency_points(1e9, 10e9)

    print("Uniform sampling: 51 points")
    print(f"Adaptive sampling: {len(adaptive_freqs)} points")
    print(
        f"Extra resolution around 5 GHz: {np.sum((adaptive_freqs >= 4.5e9) & (adaptive_freqs <= 5.5e9))} points"
    )

    # Create adaptive config
    adaptive_config = create_two_port_network_config(
        mesh_file="resonant_structure.msh",
        freq_min=1e9,
        freq_max=10e9,
        freq_step=0.1e9,  # Will be overridden by adaptive sampling
    )

    # Add custom frequency points (this would require Palace support)
    adaptive_config["Solver"]["Driven"]["FrequencyList"] = adaptive_freqs.tolist()

    save_config(adaptive_config, "adaptive_frequency_sweep.json")
    print("✓ Saved adaptive_frequency_sweep.json")

    # Plot frequency sampling comparison
    plt.figure(figsize=(12, 6))

    plt.subplot(1, 2, 1)
    f_uniform = np.linspace(1e9, 10e9, 51)
    plt.plot(
        f_uniform / 1e9, np.ones_like(f_uniform), "bo", markersize=4, label="Uniform"
    )
    plt.xlabel("Frequency (GHz)")
    plt.ylabel("Sampling Density")
    plt.title("Uniform Frequency Sampling")
    plt.grid(True, alpha=0.3)
    plt.ylim([0.5, 1.5])

    plt.subplot(1, 2, 2)
    plt.plot(
        adaptive_freqs / 1e9,
        np.ones_like(adaptive_freqs),
        "ro",
        markersize=4,
        label="Adaptive",
    )
    plt.xlabel("Frequency (GHz)")
    plt.ylabel("Sampling Density")
    plt.title("Adaptive Frequency Sampling")
    plt.grid(True, alpha=0.3)
    plt.ylim([0.5, 1.5])

    # Highlight dense region
    dense_region = (adaptive_freqs >= 4.5e9) & (adaptive_freqs <= 5.5e9)
    plt.plot(
        adaptive_freqs[dense_region] / 1e9,
        np.ones_like(adaptive_freqs[dense_region]),
        "go",
        markersize=6,
        label="High resolution",
    )
    plt.legend()

    plt.tight_layout()
    plt.savefig("frequency_sampling_comparison.png", dpi=150, bbox_inches="tight")
    plt.show()

    print("✓ Frequency sampling comparison saved as frequency_sampling_comparison.png")

    print(f"\n{'=' * 60}")
    print("Frequency domain simulation examples completed!")
    print("Generated files:")
    print("  - two_port_network.json")
    print("  - bandpass_filter.json")
    print("  - antenna_simulation.json")
    print("  - adaptive_frequency_sweep.json")
    print("  - s_parameter_analysis.json")
    print("  - example_s_parameters.csv")
    print("  - frequency_sampling_comparison.png")
    print("\nTo run simulations:")
    print("  palace two_port_network.json -o results_two_port/")
    print("  palace bandpass_filter.json -o results_filter/")
    print("  palace antenna_simulation.json -o results_antenna/")


if __name__ == "__main__":
    main()
