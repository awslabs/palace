"""
Unit tests for Palace core functionality.
"""

import json
import os

# Import the modules to test
import sys
from unittest.mock import Mock, patch

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))

from palace.core import PalaceSolver, run_palace


class TestPalaceSolver:
    """Test the PalaceSolver class."""

    def test_init_with_executable_path(self, mock_palace_executable):
        """Test PalaceSolver initialization with explicit executable path."""
        solver = PalaceSolver(mock_palace_executable)
        assert solver.executable == mock_palace_executable

    def test_init_without_palace_in_path(self):
        """Test PalaceSolver initialization when Palace is not in PATH."""
        with patch("subprocess.run") as mock_run:
            mock_run.return_value = Mock(returncode=1, stdout="", stderr="")

            with pytest.raises(RuntimeError, match="Palace executable not found"):
                PalaceSolver()

    def test_init_with_palace_in_path(self, mock_palace_executable):
        """Test PalaceSolver initialization when Palace is found in PATH."""
        with patch("subprocess.run") as mock_run:
            mock_run.return_value = Mock(
                returncode=0, stdout=mock_palace_executable, stderr=""
            )

            solver = PalaceSolver()
            assert solver.executable == mock_palace_executable

    def test_validate_config_valid_json(self, temp_dir, sample_config):
        """Test configuration validation with valid JSON."""
        config_file = os.path.join(temp_dir, "valid_config.json")
        with open(config_file, "w") as f:
            json.dump(sample_config, f)

        solver = PalaceSolver("dummy_path")
        assert solver.validate_config(config_file) is True

    def test_validate_config_invalid_json(self, temp_dir):
        """Test configuration validation with invalid JSON."""
        config_file = os.path.join(temp_dir, "invalid_config.json")
        with open(config_file, "w") as f:
            f.write("{ invalid json content")

        solver = PalaceSolver("dummy_path")
        assert solver.validate_config(config_file) is False

    def test_validate_config_nonexistent_file(self):
        """Test configuration validation with nonexistent file."""
        solver = PalaceSolver("dummy_path")
        assert solver.validate_config("nonexistent.json") is False

    def test_run_basic(self, temp_dir, sample_config, mock_palace_executable):
        """Test basic simulation run."""
        config_file = os.path.join(temp_dir, "test_config.json")
        with open(config_file, "w") as f:
            json.dump(sample_config, f)

        solver = PalaceSolver(mock_palace_executable)

        with patch("subprocess.run") as mock_run:
            mock_run.return_value = Mock(returncode=0, stdout="Success", stderr="")

            result = solver.run(config_file)

            # Verify subprocess was called correctly
            mock_run.assert_called_once()
            call_args = mock_run.call_args[0][0]
            assert mock_palace_executable in call_args
            assert config_file in call_args

    def test_run_with_mpi(self, temp_dir, sample_config, mock_palace_executable):
        """Test simulation run with MPI."""
        config_file = os.path.join(temp_dir, "test_config.json")
        with open(config_file, "w") as f:
            json.dump(sample_config, f)

        solver = PalaceSolver(mock_palace_executable)

        with patch("subprocess.run") as mock_run:
            mock_run.return_value = Mock(returncode=0, stdout="Success", stderr="")

            result = solver.run(config_file, num_procs=4)

            # Verify MPI was included in command
            call_args = mock_run.call_args[0][0]
            assert "mpirun" in call_args
            assert "-np" in call_args
            assert "4" in call_args

    def test_run_with_output_dir(self, temp_dir, sample_config, mock_palace_executable):
        """Test simulation run with output directory."""
        config_file = os.path.join(temp_dir, "test_config.json")
        with open(config_file, "w") as f:
            json.dump(sample_config, f)

        output_dir = os.path.join(temp_dir, "results")

        solver = PalaceSolver(mock_palace_executable)

        with patch("subprocess.run") as mock_run:
            mock_run.return_value = Mock(returncode=0, stdout="Success", stderr="")

            result = solver.run(config_file, output_dir=output_dir)

            # Verify output directory was included
            call_args = mock_run.call_args[0][0]
            assert "-o" in call_args
            assert output_dir in call_args

    def test_run_with_kwargs(self, temp_dir, sample_config, mock_palace_executable):
        """Test simulation run with additional keyword arguments."""
        config_file = os.path.join(temp_dir, "test_config.json")
        with open(config_file, "w") as f:
            json.dump(sample_config, f)

        solver = PalaceSolver(mock_palace_executable)

        with patch("subprocess.run") as mock_run:
            mock_run.return_value = Mock(returncode=0, stdout="Success", stderr="")

            result = solver.run(config_file, verbose=True, threads=8)

            # Verify additional arguments were included
            call_args = mock_run.call_args[0][0]
            assert "--verbose" in call_args
            assert "--threads" in call_args
            assert "8" in call_args

    def test_run_nonexistent_config(self, mock_palace_executable):
        """Test run with nonexistent configuration file."""
        solver = PalaceSolver(mock_palace_executable)

        with pytest.raises(FileNotFoundError):
            solver.run("nonexistent_config.json")

    def test_run_simulation_failure(
        self, temp_dir, sample_config, mock_palace_executable
    ):
        """Test handling of simulation failure."""
        config_file = os.path.join(temp_dir, "test_config.json")
        with open(config_file, "w") as f:
            json.dump(sample_config, f)

        solver = PalaceSolver(mock_palace_executable)

        with patch("subprocess.run") as mock_run:
            mock_run.return_value = Mock(
                returncode=1, stdout="", stderr="Simulation failed"
            )

            result = solver.run(config_file)

            # Verify failure is properly returned
            assert result.returncode == 1
            assert "Simulation failed" in result.stderr


class TestRunPalaceFunction:
    """Test the run_palace convenience function."""

    def test_run_palace_basic(self, temp_dir, sample_config):
        """Test basic run_palace function call."""
        config_file = os.path.join(temp_dir, "test_config.json")
        with open(config_file, "w") as f:
            json.dump(sample_config, f)

        with patch("palace.core.PalaceSolver") as mock_solver_class:
            mock_solver = Mock()
            mock_solver.run.return_value = Mock(returncode=0)
            mock_solver_class.return_value = mock_solver

            result = run_palace(config_file)

            # Verify PalaceSolver was created and run was called
            mock_solver_class.assert_called_once()
            mock_solver.run.assert_called_once_with(config_file)

    def test_run_palace_with_kwargs(self, temp_dir, sample_config):
        """Test run_palace with keyword arguments."""
        config_file = os.path.join(temp_dir, "test_config.json")
        with open(config_file, "w") as f:
            json.dump(sample_config, f)

        with patch("palace.core.PalaceSolver") as mock_solver_class:
            mock_solver = Mock()
            mock_solver.run.return_value = Mock(returncode=0)
            mock_solver_class.return_value = mock_solver

            result = run_palace(config_file, num_procs=4, output_dir="results")

            # Verify arguments were passed through
            mock_solver.run.assert_called_once_with(
                config_file, num_procs=4, output_dir="results"
            )


class TestModuleImports:
    """Test that the module imports work correctly."""

    def test_import_palace_core(self):
        """Test that palace.core imports successfully."""
        import palace.core

        assert hasattr(palace.core, "PalaceSolver")
        assert hasattr(palace.core, "run_palace")

    def test_import_from_palace(self):
        """Test that imports from palace work."""
        from palace import PalaceSolver, run_palace

        assert PalaceSolver is not None
        assert run_palace is not None
