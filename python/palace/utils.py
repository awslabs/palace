"""
Utility functions for Palace simulations.

This module provides helper functions for mesh processing, configuration
management, and result analysis.
"""

import json
import os
from typing import Any, Dict, List, Tuple

import numpy as np


def load_config(config_file: str) -> Dict[str, Any]:
    """
    Load Palace configuration from JSON file.

    Args:
        config_file: Path to configuration file

    Returns:
        Configuration dictionary
    """
    with open(config_file) as f:
        return json.load(f)


def save_config(config: Dict[str, Any], config_file: str) -> None:
    """
    Save Palace configuration to JSON file.

    Args:
        config: Configuration dictionary
        config_file: Path to output file
    """
    with open(config_file, "w") as f:
        json.dump(config, f, indent=2)


def create_basic_config(
    mesh_file: str, problem_type: str = "Eigenmode"
) -> Dict[str, Any]:
    """
    Create a basic Palace configuration.

    Args:
        mesh_file: Path to mesh file
        problem_type: Type of simulation ("Eigenmode", "Driven", "Transient", etc.)

    Returns:
        Basic configuration dictionary
    """
    config = {
        "Problem": {"Type": problem_type, "Verbose": 1},
        "Model": {
            "Mesh": mesh_file,
            "L0": 1e-3,  # Length scale in meters
        },
    }

    if problem_type == "Eigenmode":
        config["Solver"] = {
            "Eigenmode": {
                "Target": 5.0e9,  # Target frequency in Hz
                "Tol": 1e-9,
                "MaxIts": 1000,
                "MaxSize": 50,
            }
        }
    elif problem_type == "Driven":
        config["Solver"] = {
            "Driven": {"MinFreq": 1.0e9, "MaxFreq": 10.0e9, "FreqStep": 1.0e8}
        }
    elif problem_type == "Transient":
        config["Solver"] = {
            "Transient": {
                "MaxTime": 5.0e-9,  # 5 ns simulation time
                "TimeStep": 1.0e-12,  # 1 ps time step
                "Order": 2,
            }
        }

    return config


def read_csv_results(csv_file: str) -> Tuple[np.ndarray, np.ndarray]:
    """
    Read Palace CSV output file.

    Args:
        csv_file: Path to CSV results file

    Returns:
        Tuple of (frequency/time array, data array)
    """
    try:
        data = np.loadtxt(csv_file, delimiter=",", skiprows=1)
        if data.ndim == 1:
            data = data.reshape(1, -1)
        return data[:, 0], data[:, 1:]
    except Exception as e:
        raise ValueError(f"Error reading CSV file {csv_file}: {e}")


def find_mesh_files(directory: str, extensions: List[str] = None) -> List[str]:
    """
    Find mesh files in a directory.

    Args:
        directory: Directory to search
        extensions: List of file extensions to look for

    Returns:
        List of mesh file paths
    """
    if extensions is None:
        extensions = [".msh", ".mesh", ".vtk", ".vtu", ".exo"]

    mesh_files = []
    for root, dirs, files in os.walk(directory):
        for file in files:
            if any(file.lower().endswith(ext) for ext in extensions):
                mesh_files.append(os.path.join(root, file))

    return sorted(mesh_files)


def get_example_configs() -> Dict[str, str]:
    """
    Get paths to example configuration files.

    Returns:
        Dictionary mapping example names to config file paths
    """
    package_dir = os.path.dirname(__file__)
    examples_dir = os.path.join(
        os.path.dirname(os.path.dirname(package_dir)), "examples"
    )

    examples = {}
    if os.path.exists(examples_dir):
        for item in os.listdir(examples_dir):
            item_path = os.path.join(examples_dir, item)
            if os.path.isdir(item_path):
                # Look for JSON config files
                for file in os.listdir(item_path):
                    if file.endswith(".json"):
                        examples[f"{item}_{file[:-5]}"] = os.path.join(item_path, file)

    return examples


def validate_mesh_file(mesh_file: str) -> bool:
    """
    Basic validation of mesh file.

    Args:
        mesh_file: Path to mesh file

    Returns:
        True if file appears to be a valid mesh
    """
    if not os.path.exists(mesh_file):
        return False

    # Basic file size check
    if os.path.getsize(mesh_file) < 100:  # Too small to be a real mesh
        return False

    # Check file extension
    valid_extensions = [".msh", ".mesh", ".vtk", ".vtu", ".exo", ".e"]
    return any(mesh_file.lower().endswith(ext) for ext in valid_extensions)


__all__ = [
    "load_config",
    "save_config",
    "create_basic_config",
    "read_csv_results",
    "find_mesh_files",
    "get_example_configs",
    "validate_mesh_file",
]
