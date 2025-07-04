#!/usr/bin/env python3
"""
Time Domain Simulation Example for Palace Python Interface.

This script demonstrates how to set up and analyze time domain transient
simulations for pulse propagation, time-domain reflectometry, and broadband analysis.
"""

import json
import os
import sys

import matplotlib.pyplot as plt
import numpy as np
from scipy import signal
from scipy.fft import fft, fftfreq

# Add parent directory to path
sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))

from palace.utils import create_basic_config, read_csv_results, save_config


def create_transient_config(
    mesh_file="transmission_line.msh",
    time_final=5e-9,
    time_step=1e-12,
    excitation_type="gaussian",
):
    """
    Create configuration for time domain transient simulation.

    Args:
        mesh_file: Path to mesh file
        time_final: Final simulation time (s)
        time_step: Time step size (s)
        excitation_type: Type of excitation pulse

    Returns:
        Configuration dictionary
    """
    config = create_basic_config(mesh_file, "Transient")

    # Configure time domain solver
    config["Solver"]["Transient"] = {
        "MaxTime": time_final,
        "TimeStep": time_step,
        "Order": 2,  # Second-order time integration
        "SaveStep": max(1, int(1e-11 / time_step)),  # Save every 10 ps
        "Type": "GeneralizedAlpha",  # Time integration scheme
    }

    # Define excitation
    if excitation_type == "gaussian":
        # Gaussian pulse parameters
        pulse_width = 100e-12  # 100 ps width
        pulse_center = 200e-12  # Center at 200 ps

        config["Model"]["Excitation"] = {
            "Type": "Gaussian",
            "Amplitude": 1.0,
            "Width": pulse_width,
            "Center": pulse_center,
        }

    elif excitation_type == "modulated_gaussian":
        # Modulated Gaussian (for specific frequency content)
        carrier_freq = 5e9  # 5 GHz carrier
        pulse_width = 200e-12
        pulse_center = 300e-12

        config["Model"]["Excitation"] = {
            "Type": "ModulatedGaussian",
            "Amplitude": 1.0,
            "Width": pulse_width,
            "Center": pulse_center,
            "Frequency": carrier_freq,
        }

    elif excitation_type == "step":
        # Step function
        config["Model"]["Excitation"] = {
            "Type": "Step",
            "Amplitude": 1.0,
            "RiseTime": 10e-12,  # 10 ps rise time
        }

    # Port configuration for transient
    config["Model"]["Boundary"] = [
        {
            "Index": 1,  # Excitation port
            "LumpedPort": {"R": 50.0, "Excitation": True},
        },
        {
            "Index": 2,  # Observation port
            "LumpedPort": {"R": 50.0},
        },
    ]

    # Material properties
    config["Model"]["Domain"] = [
        {
            "Index": 1,  # Dielectric
            "Material": {
                "Permeability": 1.0,
                "Permittivity": 2.2,  # Low-loss dielectric
                "LossTan": 0.001,
            },
        },
        {
            "Index": 2,  # Air
            "Material": {"Permeability": 1.0, "Permittivity": 1.0},
        },
    ]

    # Time domain post-processing
    config["Model"]["PostProcessing"] = {
        "Probe": [
            {
                "Index": 1,
                "Center": [0.001, 0.0, 0.0005],  # Near input
                "Type": "Electric",
            },
            {
                "Index": 2,
                "Center": [0.005, 0.0, 0.0005],  # Middle
                "Type": "Electric",
            },
            {
                "Index": 3,
                "Center": [0.009, 0.0, 0.0005],  # Near output
                "Type": "Electric",
            },
        ]
    }

    return config


def create_tdr_config(
    mesh_file="device_under_test.msh", rise_time=20e-12, simulation_time=10e-9
):
    """
    Create configuration for Time Domain Reflectometry (TDR) simulation.

    Args:
        mesh_file: Path to mesh file
        rise_time: Step function rise time (s)
        simulation_time: Total simulation time (s)

    Returns:
        Configuration dictionary
    """
    # Use small time step for accurate TDR
    time_step = rise_time / 20  # 20 points per rise time

    config = create_transient_config(
        mesh_file=mesh_file,
        time_final=simulation_time,
        time_step=time_step,
        excitation_type="step",
    )

    # Modify for TDR-specific settings
    config["Model"]["Excitation"]["RiseTime"] = rise_time
    config["Model"]["Excitation"]["Amplitude"] = 0.5  # TDR typically uses 0.5V step

    # Single port for TDR (incident and reflected from same port)
    config["Model"]["Boundary"] = [
        {
            "Index": 1,  # TDR port
            "LumpedPort": {"R": 50.0, "Excitation": True},
        }
    ]

    # Additional probes along transmission line
    config["Model"]["PostProcessing"]["Probe"] = [
        {
            "Index": i,
            "Center": [i * 0.001, 0.0, 0.0005],  # Probes every 1mm
            "Type": "Electric",
        }
        for i in range(1, 11)  # 10 probes along the line
    ]

    return config


def analyze_transient_results(time_file="probe-E.csv", port_file="port-V.csv"):
    """
    Analyze time domain simulation results.

    Args:
        time_file: Path to time domain probe data
        port_file: Path to port voltage data

    Returns:
        Analysis results dictionary
    """
    results = {}

    # Create example data if files don't exist
    if not os.path.exists(time_file) or not os.path.exists(port_file):
        print("Result files not found. Creating example data...")

        # Generate example time domain data
        t_max = 5e-9  # 5 ns
        dt = 1e-12  # 1 ps
        t = np.arange(0, t_max, dt)

        # Gaussian pulse excitation
        pulse_width = 100e-12
        pulse_center = 200e-12
        excitation = np.exp(-(((t - pulse_center) / pulse_width) ** 2))

        # Simulate transmission line response with reflections
        c = 3e8 / np.sqrt(2.2)  # Wave speed in dielectric
        line_length = 0.01  # 1 cm line
        delay = line_length / c

        # Direct transmission
        transmitted = np.zeros_like(t)
        delay_samples = int(delay / dt)
        if delay_samples < len(t):
            transmitted[delay_samples:] = 0.8 * excitation[: len(t) - delay_samples]

        # Reflection from load
        reflection_delay = 2 * delay
        reflected = np.zeros_like(t)
        valid_idx = t >= reflection_delay
        if np.sum(valid_idx) > 0:
            refl_start = int(reflection_delay / dt)
            refl_signal = 0.2 * excitation[: len(t) - refl_start]
            reflected[refl_start : refl_start + len(refl_signal)] = refl_signal

        # Input voltage (incident + reflected)
        v_input = excitation + reflected

        # Output voltage (transmitted)
        v_output = transmitted

        # Add some noise for realism
        noise_level = 0.01
        v_input += noise_level * np.random.randn(len(t))
        v_output += noise_level * np.random.randn(len(t))

        # Save example data
        probe_data = np.column_stack([t, v_input, v_output])
        np.savetxt(
            "example_probe_data.csv",
            probe_data,
            delimiter=",",
            header="Time(s),E_field_input,E_field_output",
        )

        port_data = np.column_stack([t, v_input, v_output])
        np.savetxt(
            "example_port_data.csv",
            port_data,
            delimiter=",",
            header="Time(s),V_port1,V_port2",
        )

        time_file = "example_probe_data.csv"
        port_file = "example_port_data.csv"

    # Load data
    try:
        time_data, probe_values = read_csv_results(time_file)
        port_time, port_values = read_csv_results(port_file)

        results["time"] = time_data
        results["probe_data"] = probe_values
        results["port_data"] = port_values

    except Exception as e:
        print(f"Error loading results: {e}")
        return None

    # Perform frequency domain analysis via FFT
    dt = time_data[1] - time_data[0]
    n_samples = len(time_data)

    # FFT of input and output signals
    if probe_values.shape[1] >= 2:
        input_signal = probe_values[:, 0]
        output_signal = probe_values[:, 1]

        # Window the data to reduce spectral leakage
        window = signal.windows.hann(n_samples)
        input_windowed = input_signal * window
        output_windowed = output_signal * window

        # Compute FFT
        input_fft = fft(input_windowed)
        output_fft = fft(output_windowed)
        freqs = fftfreq(n_samples, dt)

        # Only keep positive frequencies
        pos_freq_idx = freqs > 0
        freqs_pos = freqs[pos_freq_idx]
        input_fft_pos = input_fft[pos_freq_idx]
        output_fft_pos = output_fft[pos_freq_idx]

        # Calculate transfer function
        transfer_function = output_fft_pos / (
            input_fft_pos + 1e-12
        )  # Avoid division by zero

        results["frequencies"] = freqs_pos
        results["input_spectrum"] = input_fft_pos
        results["output_spectrum"] = output_fft_pos
        results["transfer_function"] = transfer_function
        results["input_signal"] = input_signal
        results["output_signal"] = output_signal

    # Time domain metrics
    if "input_signal" in results:
        # Find pulse arrival times
        threshold = 0.1 * np.max(np.abs(results["input_signal"]))
        input_arrival = np.where(np.abs(results["input_signal"]) > threshold)[0]
        output_arrival = np.where(np.abs(results["output_signal"]) > threshold)[0]

        if len(input_arrival) > 0 and len(output_arrival) > 0:
            propagation_delay = (output_arrival[0] - input_arrival[0]) * dt
            results["propagation_delay"] = propagation_delay

        # Peak values and times
        input_max_idx = np.argmax(np.abs(results["input_signal"]))
        output_max_idx = np.argmax(np.abs(results["output_signal"]))

        results["input_peak_time"] = time_data[input_max_idx]
        results["input_peak_value"] = results["input_signal"][input_max_idx]
        results["output_peak_time"] = time_data[output_max_idx]
        results["output_peak_value"] = results["output_signal"][output_max_idx]

        # Transmission coefficient
        results["transmission_coefficient"] = (
            results["output_peak_value"] / results["input_peak_value"]
        )

    return results


def plot_transient_analysis(results, title="Time Domain Analysis"):
    """
    Create comprehensive time domain analysis plots.

    Args:
        results: Results dictionary from analyze_transient_results
        title: Plot title
    """
    if results is None:
        print("No results to plot")
        return

    fig, axes = plt.subplots(2, 3, figsize=(18, 12))
    fig.suptitle(title, fontsize=16)

    time = results["time"]
    input_signal = results["input_signal"]
    output_signal = results["output_signal"]

    # Time domain signals
    axes[0, 0].plot(time * 1e9, input_signal, label="Input", linewidth=2)
    axes[0, 0].plot(time * 1e9, output_signal, label="Output", linewidth=2)
    axes[0, 0].set_xlabel("Time (ns)")
    axes[0, 0].set_ylabel("Amplitude (V)")
    axes[0, 0].set_title("Time Domain Signals")
    axes[0, 0].grid(True)
    axes[0, 0].legend()

    # Zoomed view of pulse arrival
    if "propagation_delay" in results:
        delay_ns = results["propagation_delay"] * 1e9
        center_time = results["input_peak_time"] * 1e9
        time_window = 0.5  # ±0.5 ns window

        time_mask = (time * 1e9 >= center_time - time_window) & (
            time * 1e9 <= center_time + time_window
        )
        axes[0, 1].plot(
            time[time_mask] * 1e9, input_signal[time_mask], label="Input", linewidth=2
        )
        axes[0, 1].plot(
            time[time_mask] * 1e9, output_signal[time_mask], label="Output", linewidth=2
        )
        axes[0, 1].set_xlabel("Time (ns)")
        axes[0, 1].set_ylabel("Amplitude (V)")
        axes[0, 1].set_title(f"Pulse Detail (Delay: {delay_ns:.3f} ns)")
        axes[0, 1].grid(True)
        axes[0, 1].legend()
    else:
        axes[0, 1].text(
            0.5,
            0.5,
            "Propagation delay\nnot calculated",
            ha="center",
            va="center",
            transform=axes[0, 1].transAxes,
        )
        axes[0, 1].set_title("Pulse Detail")

    # Frequency domain - Input spectrum
    if "frequencies" in results:
        freqs = results["frequencies"]
        input_spectrum = results["input_spectrum"]

        axes[0, 2].semilogx(freqs / 1e9, 20 * np.log10(np.abs(input_spectrum)))
        axes[0, 2].set_xlabel("Frequency (GHz)")
        axes[0, 2].set_ylabel("Magnitude (dB)")
        axes[0, 2].set_title("Input Signal Spectrum")
        axes[0, 2].grid(True)
        axes[0, 2].set_xlim([0.1, 50])

    # Transfer function magnitude
    if "transfer_function" in results:
        transfer_mag = np.abs(results["transfer_function"])
        transfer_phase = np.angle(results["transfer_function"])

        axes[1, 0].semilogx(freqs / 1e9, 20 * np.log10(transfer_mag))
        axes[1, 0].set_xlabel("Frequency (GHz)")
        axes[1, 0].set_ylabel("|H(f)| (dB)")
        axes[1, 0].set_title("Transfer Function Magnitude")
        axes[1, 0].grid(True)
        axes[1, 0].set_xlim([0.1, 50])

        # Transfer function phase
        axes[1, 1].semilogx(freqs / 1e9, transfer_phase * 180 / np.pi)
        axes[1, 1].set_xlabel("Frequency (GHz)")
        axes[1, 1].set_ylabel("∠H(f) (degrees)")
        axes[1, 1].set_title("Transfer Function Phase")
        axes[1, 1].grid(True)
        axes[1, 1].set_xlim([0.1, 50])

    # Group delay
    if "transfer_function" in results and len(results["transfer_function"]) > 10:
        phase_unwrapped = np.unwrap(np.angle(results["transfer_function"]))
        group_delay = -np.gradient(phase_unwrapped) / (2 * np.pi * np.gradient(freqs))

        # Smooth the group delay to remove noise
        from scipy.ndimage import gaussian_filter1d

        group_delay_smooth = gaussian_filter1d(group_delay, sigma=2)

        axes[1, 2].semilogx(
            freqs[10:-10] / 1e9, group_delay_smooth[10:-10] * 1e9
        )  # Convert to ns
        axes[1, 2].set_xlabel("Frequency (GHz)")
        axes[1, 2].set_ylabel("Group Delay (ns)")
        axes[1, 2].set_title("Group Delay")
        axes[1, 2].grid(True)
        axes[1, 2].set_xlim([0.1, 50])
    else:
        axes[1, 2].text(
            0.5,
            0.5,
            "Group delay calculation\nrequires more data points",
            ha="center",
            va="center",
            transform=axes[1, 2].transAxes,
        )
        axes[1, 2].set_title("Group Delay")

    plt.tight_layout()
    plt.show()

    # Print summary
    print("\nTime Domain Analysis Summary:")
    print(f"{'=' * 50}")
    print(f"Simulation time: {time[-1] * 1e9:.1f} ns")
    print(f"Time resolution: {(time[1] - time[0]) * 1e12:.1f} ps")

    if "propagation_delay" in results:
        print(f"Propagation delay: {results['propagation_delay'] * 1e12:.1f} ps")

    if "transmission_coefficient" in results:
        print(f"Transmission coefficient: {results['transmission_coefficient']:.3f}")
        print(
            f"Transmission (dB): {20 * np.log10(abs(results['transmission_coefficient'])):.1f} dB"
        )

    if "input_peak_value" in results:
        print(
            f"Peak input: {results['input_peak_value']:.3f} V at {results['input_peak_time'] * 1e9:.3f} ns"
        )
        print(
            f"Peak output: {results['output_peak_value']:.3f} V at {results['output_peak_time'] * 1e9:.3f} ns"
        )

    if "frequencies" in results:
        # Find -3dB bandwidth
        transfer_mag_db = 20 * np.log10(np.abs(results["transfer_function"]))
        max_transfer_db = np.max(transfer_mag_db)
        bw_3db = freqs[transfer_mag_db >= max_transfer_db - 3]
        if len(bw_3db) > 0:
            bandwidth_3db = np.max(bw_3db) - np.min(bw_3db)
            print(f"3-dB bandwidth: {bandwidth_3db / 1e9:.2f} GHz")


def main():
    """Main function demonstrating time domain simulations."""
    print("Palace Time Domain Simulation Examples")
    print("=" * 60)

    # Example 1: Basic transient simulation
    print("\n1. Creating basic transient simulation configuration...")
    transient_config = create_transient_config(
        mesh_file="microstrip_line.msh",
        time_final=5e-9,
        time_step=1e-12,
        excitation_type="gaussian",
    )

    save_config(transient_config, "transient_basic.json")
    print("✓ Saved transient_basic.json")

    # Example 2: Modulated Gaussian pulse
    print("\n2. Creating modulated Gaussian pulse configuration...")
    modulated_config = create_transient_config(
        mesh_file="bandpass_structure.msh",
        time_final=3e-9,
        time_step=5e-13,
        excitation_type="modulated_gaussian",
    )

    save_config(modulated_config, "transient_modulated.json")
    print("✓ Saved transient_modulated.json")

    # Example 3: Time Domain Reflectometry (TDR)
    print("\n3. Creating TDR simulation configuration...")
    tdr_config = create_tdr_config(
        mesh_file="connector_dut.msh", rise_time=20e-12, simulation_time=10e-9
    )

    save_config(tdr_config, "tdr_simulation.json")
    print("✓ Saved tdr_simulation.json")

    # Example 4: Analyze transient results
    print("\n4. Analyzing time domain results...")
    transient_results = analyze_transient_results("probe-E.csv", "port-V.csv")

    if transient_results:
        plot_transient_analysis(
            transient_results, "Transmission Line Transient Analysis"
        )

        # Save analysis results
        analysis_summary = {
            "simulation_time_ns": float(transient_results["time"][-1] * 1e9),
            "time_resolution_ps": float(
                (transient_results["time"][1] - transient_results["time"][0]) * 1e12
            ),
            "propagation_delay_ps": float(transient_results["propagation_delay"] * 1e12)
            if "propagation_delay" in transient_results
            else None,
            "transmission_coefficient": float(
                transient_results["transmission_coefficient"]
            )
            if "transmission_coefficient" in transient_results
            else None,
            "peak_input_V": float(transient_results["input_peak_value"])
            if "input_peak_value" in transient_results
            else None,
            "peak_output_V": float(transient_results["output_peak_value"])
            if "output_peak_value" in transient_results
            else None,
        }

        with open("transient_analysis.json", "w") as f:
            json.dump(analysis_summary, f, indent=2)
        print("✓ Saved transient_analysis.json")

    # Example 5: Eye diagram analysis (for digital signals)
    print("\n5. Creating eye diagram analysis example...")

    def create_eye_diagram(bit_stream, bit_rate=1e9, samples_per_bit=100):
        """Create an eye diagram from a digital bit stream."""
        # Generate time axis
        bit_duration = 1 / bit_rate
        dt = bit_duration / samples_per_bit
        t_bit = np.arange(0, bit_duration, dt)

        # Create full signal
        signal_full = []
        for bit in bit_stream:
            if bit == 1:
                # Simple NRZ encoding with some rise/fall time
                rise_samples = int(0.1 * samples_per_bit)
                bit_signal = np.concatenate(
                    [
                        np.linspace(0, 1, rise_samples),
                        np.ones(samples_per_bit - 2 * rise_samples),
                        np.linspace(1, 0, rise_samples),
                    ]
                )
            else:
                bit_signal = np.zeros(samples_per_bit)

            signal_full.extend(bit_signal)

        # Add noise
        noise_level = 0.05
        signal_full = np.array(signal_full) + noise_level * np.random.randn(
            len(signal_full)
        )

        # Create eye diagram by overlapping bit periods
        eye_data = []
        for i in range(len(bit_stream) - 1):
            start_idx = i * samples_per_bit
            end_idx = (i + 2) * samples_per_bit  # Two bit periods for eye
            if end_idx <= len(signal_full):
                eye_data.append(signal_full[start_idx:end_idx])

        return np.array(eye_data), np.arange(0, 2 * bit_duration, dt)

    # Generate example digital signal
    np.random.seed(42)  # For reproducible results
    bit_sequence = np.random.randint(0, 2, 20)  # 20 random bits
    eye_traces, eye_time = create_eye_diagram(
        bit_sequence, bit_rate=1e9, samples_per_bit=100
    )

    # Plot eye diagram
    plt.figure(figsize=(12, 8))

    plt.subplot(2, 2, 1)
    # Plot several traces
    for i in range(min(10, len(eye_traces))):
        plt.plot(eye_time * 1e9, eye_traces[i], "b-", alpha=0.3)
    plt.xlabel("Time (ns)")
    plt.ylabel("Amplitude (V)")
    plt.title("Eye Diagram")
    plt.grid(True)
    plt.axvline(1.0, color="r", linestyle="--", alpha=0.7, label="Bit boundary")
    plt.legend()

    # Original bit stream
    plt.subplot(2, 2, 2)
    t_stream = np.arange(len(bit_sequence))
    plt.step(t_stream, bit_sequence, where="post", linewidth=2)
    plt.xlabel("Bit Number")
    plt.ylabel("Bit Value")
    plt.title("Original Bit Stream")
    plt.grid(True)
    plt.ylim([-0.5, 1.5])

    # Eye diagram statistics
    plt.subplot(2, 2, 3)
    # Calculate eye opening
    mid_point = len(eye_time) // 2  # Middle of eye
    eye_opening = []
    for trace in eye_traces:
        eye_opening.append(
            np.max(trace[mid_point - 10 : mid_point + 10])
            - np.min(trace[mid_point - 10 : mid_point + 10])
        )

    plt.hist(eye_opening, bins=15, alpha=0.7, edgecolor="black")
    plt.xlabel("Eye Opening (V)")
    plt.ylabel("Count")
    plt.title("Eye Opening Distribution")
    plt.grid(True)

    # Jitter analysis
    plt.subplot(2, 2, 4)
    # Find zero crossings for jitter measurement
    zero_crossings = []
    for trace in eye_traces:
        # Find where signal crosses 0.5V (midpoint)
        crossings = np.where(np.diff(np.sign(trace - 0.5)))[0]
        if len(crossings) > 0:
            zero_crossings.extend(eye_time[crossings])

    if zero_crossings:
        jitter_ps = (
            np.array(zero_crossings) - 1e-9
        ) * 1e12  # Convert to ps relative to ideal
        plt.hist(jitter_ps, bins=15, alpha=0.7, edgecolor="black")
        plt.xlabel("Timing Jitter (ps)")
        plt.ylabel("Count")
        plt.title("Timing Jitter Distribution")
        plt.grid(True)
    else:
        plt.text(
            0.5,
            0.5,
            "No zero crossings found",
            ha="center",
            va="center",
            transform=plt.gca().transAxes,
        )
        plt.title("Timing Jitter Distribution")

    plt.tight_layout()
    plt.savefig("eye_diagram_analysis.png", dpi=150, bbox_inches="tight")
    plt.show()

    print("✓ Eye diagram analysis saved as eye_diagram_analysis.png")

    # Example 6: Pulse compression analysis
    print("\n6. Pulse compression analysis...")

    # Simulate chirped pulse compression
    def chirped_pulse(t, duration, bandwidth):
        """Generate a linear chirped pulse."""
        # Linear chirp
        chirp_rate = bandwidth / duration
        envelope = np.exp(-((t / (duration / 4)) ** 2))  # Gaussian envelope
        phase = np.pi * chirp_rate * t**2
        return envelope * np.cos(phase)

    t_pulse = np.linspace(-1e-9, 1e-9, 2000)  # 2 ns window, 1 ps resolution
    pulse_duration = 500e-12  # 500 ps
    pulse_bandwidth = 2e9  # 2 GHz bandwidth

    original_pulse = chirped_pulse(t_pulse, pulse_duration, pulse_bandwidth)

    # Simulate dispersion compensation (pulse compression)
    # This would normally be done by the Palace simulation
    compressed_pulse = chirped_pulse(
        t_pulse, pulse_duration / 10, pulse_bandwidth
    )  # 10x compression

    plt.figure(figsize=(15, 10))

    # Time domain comparison
    plt.subplot(2, 3, 1)
    plt.plot(t_pulse * 1e12, original_pulse, label="Original chirped pulse")
    plt.xlabel("Time (ps)")
    plt.ylabel("Amplitude")
    plt.title("Original Chirped Pulse")
    plt.grid(True)
    plt.legend()

    plt.subplot(2, 3, 2)
    plt.plot(t_pulse * 1e12, compressed_pulse, label="Compressed pulse", color="red")
    plt.xlabel("Time (ps)")
    plt.ylabel("Amplitude")
    plt.title("Compressed Pulse")
    plt.grid(True)
    plt.legend()

    # Frequency domain
    original_fft = fft(original_pulse)
    compressed_fft = fft(compressed_pulse)
    freqs_pulse = fftfreq(len(t_pulse), t_pulse[1] - t_pulse[0])

    plt.subplot(2, 3, 3)
    pos_freq = freqs_pulse > 0
    plt.plot(
        freqs_pulse[pos_freq] / 1e9, np.abs(original_fft[pos_freq]), label="Original"
    )
    plt.plot(
        freqs_pulse[pos_freq] / 1e9,
        np.abs(compressed_fft[pos_freq]),
        label="Compressed",
    )
    plt.xlabel("Frequency (GHz)")
    plt.ylabel("Magnitude")
    plt.title("Frequency Spectra")
    plt.grid(True)
    plt.legend()
    plt.xlim([0, 5])

    # Pulse width analysis
    def pulse_width_fwhm(signal, time):
        """Calculate Full Width at Half Maximum."""
        max_val = np.max(np.abs(signal))
        half_max_indices = np.where(np.abs(signal) >= max_val / 2)[0]
        if len(half_max_indices) > 0:
            return time[half_max_indices[-1]] - time[half_max_indices[0]]
        return 0

    original_width = pulse_width_fwhm(original_pulse, t_pulse)
    compressed_width = pulse_width_fwhm(compressed_pulse, t_pulse)

    plt.subplot(2, 3, 4)
    widths = [original_width * 1e12, compressed_width * 1e12]
    labels = ["Original", "Compressed"]
    colors = ["blue", "red"]
    plt.bar(labels, widths, color=colors, alpha=0.7)
    plt.ylabel("Pulse Width (ps)")
    plt.title("Pulse Width Comparison")
    plt.grid(True, alpha=0.3)

    # Compression ratio
    compression_ratio = original_width / compressed_width
    plt.text(
        0.5,
        max(widths) * 0.8,
        f"Compression Ratio: {compression_ratio:.1f}:1",
        ha="center",
        fontweight="bold",
    )

    # Time-bandwidth product
    plt.subplot(2, 3, 5)
    tb_original = original_width * pulse_bandwidth
    tb_compressed = compressed_width * pulse_bandwidth

    tb_products = [tb_original, tb_compressed]
    plt.bar(labels, tb_products, color=colors, alpha=0.7)
    plt.ylabel("Time-Bandwidth Product")
    plt.title("Time-Bandwidth Product")
    plt.grid(True, alpha=0.3)
    plt.axhline(
        0.441,
        color="black",
        linestyle="--",
        alpha=0.7,
        label="Transform limit (Gaussian)",
    )
    plt.legend()

    # Peak power enhancement
    plt.subplot(2, 3, 6)
    peak_original = np.max(np.abs(original_pulse))
    peak_compressed = np.max(np.abs(compressed_pulse))
    peak_enhancement = peak_compressed / peak_original

    peaks = [peak_original, peak_compressed]
    plt.bar(labels, peaks, color=colors, alpha=0.7)
    plt.ylabel("Peak Amplitude")
    plt.title("Peak Power Enhancement")
    plt.grid(True, alpha=0.3)
    plt.text(
        0.5,
        max(peaks) * 0.8,
        f"Enhancement: {peak_enhancement:.1f}×",
        ha="center",
        fontweight="bold",
    )

    plt.tight_layout()
    plt.savefig("pulse_compression_analysis.png", dpi=150, bbox_inches="tight")
    plt.show()

    print("✓ Pulse compression analysis saved as pulse_compression_analysis.png")
    print(f"  Original pulse width: {original_width * 1e12:.1f} ps")
    print(f"  Compressed pulse width: {compressed_width * 1e12:.1f} ps")
    print(f"  Compression ratio: {compression_ratio:.1f}:1")
    print(f"  Peak enhancement: {peak_enhancement:.1f}×")

    print(f"\n{'=' * 60}")
    print("Time domain simulation examples completed!")
    print("Generated files:")
    print("  - transient_basic.json")
    print("  - transient_modulated.json")
    print("  - tdr_simulation.json")
    print("  - transient_analysis.json")
    print("  - example_probe_data.csv")
    print("  - example_port_data.csv")
    print("  - eye_diagram_analysis.png")
    print("  - pulse_compression_analysis.png")
    print("\nTo run simulations:")
    print("  palace transient_basic.json -o results_transient/")
    print("  palace tdr_simulation.json -o results_tdr/")


if __name__ == "__main__":
    main()
