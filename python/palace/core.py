"""
Core Palace functionality and bindings.

This module provides the main interface to Palace solver capabilities.
"""

import os
import subprocess
from typing import Optional


class PalaceSolver:
    """
    Main interface to Palace electromagnetic solver.

    This class provides a Python interface to run Palace simulations
    and manage configuration files.
    """

    def __init__(self, palace_executable: Optional[str] = None):
        """
        Initialize Palace solver.

        Args:
            palace_executable: Path to palace executable. If None, searches PATH.
        """
        self.executable = palace_executable or self._find_executable()
        if not self.executable:
            raise RuntimeError(
                "Palace executable not found. Please install Palace or specify path."
            )

    def _find_executable(self) -> Optional[str]:
        """Find palace executable in PATH."""
        try:
            result = subprocess.run(["which", "palace"], capture_output=True, text=True)
            if result.returncode == 0:
                return result.stdout.strip()
        except Exception:
            pass
        return None

    def run(
        self,
        config_file: str,
        output_dir: Optional[str] = None,
        num_procs: Optional[int] = None,
        **kwargs,
    ) -> subprocess.CompletedProcess:
        """
        Run Palace simulation.

        Args:
            config_file: Path to Palace configuration JSON file
            output_dir: Output directory for results
            num_procs: Number of MPI processes to use
            **kwargs: Additional command line arguments

        Returns:
            CompletedProcess object with simulation results
        """
        if not os.path.exists(config_file):
            raise FileNotFoundError(f"Configuration file not found: {config_file}")

        cmd = []
        if num_procs and num_procs > 1:
            cmd.extend(["mpirun", "-np", str(num_procs)])

        cmd.append(self.executable)
        cmd.append(config_file)

        if output_dir:
            cmd.extend(["-o", output_dir])

        # Add additional arguments
        for key, value in kwargs.items():
            if value is True:
                cmd.append(f"--{key}")
            elif value is not False and value is not None:
                cmd.extend([f"--{key}", str(value)])

        return subprocess.run(cmd, capture_output=True, text=True)

    def validate_config(self, config_file: str) -> bool:
        """
        Validate Palace configuration file.

        Args:
            config_file: Path to configuration file

        Returns:
            True if configuration is valid
        """
        # This would typically call Palace's built-in validation
        # For now, just check if file exists and is valid JSON
        try:
            import json

            with open(config_file) as f:
                json.load(f)
            return True
        except Exception:
            return False


def run_palace(config_file: str, **kwargs) -> subprocess.CompletedProcess:
    """
    Convenience function to run Palace simulation.

    Args:
        config_file: Path to Palace configuration file
        **kwargs: Additional arguments passed to PalaceSolver.run()

    Returns:
        CompletedProcess object with simulation results
    """
    solver = PalaceSolver()
    return solver.run(config_file, **kwargs)


__all__ = ["PalaceSolver", "run_palace"]
