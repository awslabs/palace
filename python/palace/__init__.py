"""
Palace: 3D Finite Element Solver for Computational Electromagnetics

This package provides Python bindings for Palace, a parallel finite element code
for full-wave 3D electromagnetic simulations in the frequency or time domain.

Key features:
- Eigenmode calculations with material/radiative loss
- Frequency domain driven simulations
- Time domain transient analysis
- Lumped parameter extraction
- Adaptive mesh refinement
- High-order finite elements
- GPU acceleration support
"""

__version__ = "0.1.0"
__author__ = "AWS Center for Quantum Computing"
__email__ = "palace-maint@amazon.com"

# Import main classes and functions
from .core import PalaceSolver, run_palace
from .utils import (
    create_basic_config,
    find_mesh_files,
    get_example_configs,
    load_config,
    read_csv_results,
    save_config,
    validate_mesh_file,
)

__all__ = [
    "__version__",
    "__author__",
    "__email__",
    "PalaceSolver",
    "run_palace",
    "load_config",
    "save_config",
    "create_basic_config",
    "read_csv_results",
    "find_mesh_files",
    "get_example_configs",
    "validate_mesh_file",
]
