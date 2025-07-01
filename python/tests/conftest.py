"""
Pytest configuration and shared fixtures for Palace Python tests.
"""

import os
import shutil
import tempfile

import numpy as np
import pytest


@pytest.fixture
def temp_dir():
    """Create a temporary directory for test files."""
    temp_dir = tempfile.mkdtemp()
    yield temp_dir
    shutil.rmtree(temp_dir)


@pytest.fixture
def sample_config():
    """Create a sample Palace configuration for testing."""
    return {
        "Problem": {"Type": "Eigenmode", "Verbose": 1},
        "Model": {"Mesh": "test_mesh.msh", "L0": 1e-3},
        "Solver": {
            "Eigenmode": {"Target": 5.0e9, "Tol": 1e-9, "MaxIts": 1000, "MaxSize": 50}
        },
    }


@pytest.fixture
def sample_mesh_file(temp_dir):
    """Create a minimal mesh file for testing."""
    mesh_content = """$MeshFormat
2.2 0 8
$EndMeshFormat
$Nodes
8
1 0 0 0
2 1 0 0
3 1 1 0
4 0 1 0
5 0 0 1
6 1 0 1
7 1 1 1
8 0 1 1
$EndNodes
$Elements
1
1 15 2 0 1 1
$EndElements
"""
    mesh_file = os.path.join(temp_dir, "test_mesh.msh")
    with open(mesh_file, "w") as f:
        f.write(mesh_content)
    return mesh_file


@pytest.fixture
def sample_csv_data(temp_dir):
    """Create sample CSV data files for testing."""
    # S-parameter data
    freq = np.linspace(1e9, 10e9, 101)
    s11_real = 0.1 * np.cos(2 * np.pi * freq / 5e9)
    s11_imag = 0.1 * np.sin(2 * np.pi * freq / 5e9)
    s21_real = 0.9 * np.cos(2 * np.pi * freq / 5e9 + np.pi / 4)
    s21_imag = 0.1 * np.sin(2 * np.pi * freq / 5e9 + np.pi / 4)

    s_param_data = np.column_stack([freq, s11_real, s11_imag, s21_real, s21_imag])
    s_param_file = os.path.join(temp_dir, "test_s_params.csv")
    np.savetxt(
        s_param_file,
        s_param_data,
        delimiter=",",
        header="Frequency(Hz),Re(S11),Im(S11),Re(S21),Im(S21)",
    )

    # Eigenmode data
    eigenfreqs = np.array([4.87e9, 5.23e9, 6.14e9, 7.89e9, 8.76e9])
    q_factors = np.array([15000, 12000, 18000, 8500, 11000])
    eigen_data = np.column_stack([eigenfreqs, q_factors])
    eigen_file = os.path.join(temp_dir, "test_eigenmodes.csv")
    np.savetxt(eigen_file, eigen_data, delimiter=",", header="Frequency(Hz),Q_Factor")

    # Time domain data
    t = np.linspace(0, 5e-9, 5000)
    pulse = np.exp(-(((t - 1e-9) / 200e-12) ** 2))  # Gaussian pulse
    transmitted = np.zeros_like(t)
    delay_samples = int(0.5e-9 / (t[1] - t[0]))  # 0.5 ns delay
    if delay_samples < len(transmitted):
        transmitted[delay_samples:] = 0.8 * pulse[: len(transmitted) - delay_samples]

    time_data = np.column_stack([t, pulse, transmitted])
    time_file = os.path.join(temp_dir, "test_time_domain.csv")
    np.savetxt(
        time_file, time_data, delimiter=",", header="Time(s),Input_Signal,Output_Signal"
    )

    return {
        "s_params": s_param_file,
        "eigenmodes": eigen_file,
        "time_domain": time_file,
    }


@pytest.fixture
def mock_palace_executable(temp_dir):
    """Create a mock Palace executable for testing."""
    mock_exe = os.path.join(temp_dir, "palace")

    # Create a simple script that mimics Palace behavior
    script_content = """#!/bin/bash
# Mock Palace executable for testing
echo "Palace Mock - Processing: $1"
echo "Simulation completed successfully"
exit 0
"""

    with open(mock_exe, "w") as f:
        f.write(script_content)

    os.chmod(mock_exe, 0o755)  # Make executable
    return mock_exe


@pytest.fixture
def palace_not_available():
    """Fixture to test behavior when Palace is not available."""
    # This fixture can be used to test graceful degradation
    return True


class MockSubprocessResult:
    """Mock subprocess result for testing."""

    def __init__(self, returncode=0, stdout="", stderr=""):
        self.returncode = returncode
        self.stdout = stdout
        self.stderr = stderr


@pytest.fixture
def mock_successful_run():
    """Mock a successful Palace run."""
    return MockSubprocessResult(
        returncode=0,
        stdout="Palace simulation completed successfully\nResults written to output directory",
        stderr="",
    )


@pytest.fixture
def mock_failed_run():
    """Mock a failed Palace run."""
    return MockSubprocessResult(
        returncode=1, stdout="", stderr="Error: Configuration file validation failed"
    )
