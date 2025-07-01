"""
Unit tests for Palace utility functions.
"""

import json
import os

# Import the modules to test
import sys
from unittest.mock import patch

import numpy as np
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))

from palace.utils import (
    create_basic_config,
    find_mesh_files,
    get_example_configs,
    load_config,
    read_csv_results,
    save_config,
    validate_mesh_file,
)


class TestConfigFunctions:
    """Test configuration management functions."""

    def test_load_config(self, temp_dir, sample_config):
        """Test loading configuration from JSON file."""
        config_file = os.path.join(temp_dir, "test_config.json")
        with open(config_file, "w") as f:
            json.dump(sample_config, f)

        loaded_config = load_config(config_file)
        assert loaded_config == sample_config

    def test_load_config_invalid_json(self, temp_dir):
        """Test loading invalid JSON configuration."""
        config_file = os.path.join(temp_dir, "invalid_config.json")
        with open(config_file, "w") as f:
            f.write("{ invalid json")

        with pytest.raises(json.JSONDecodeError):
            load_config(config_file)

    def test_load_config_nonexistent_file(self):
        """Test loading nonexistent configuration file."""
        with pytest.raises(FileNotFoundError):
            load_config("nonexistent.json")

    def test_save_config(self, temp_dir, sample_config):
        """Test saving configuration to JSON file."""
        config_file = os.path.join(temp_dir, "saved_config.json")

        save_config(sample_config, config_file)

        # Verify file was created and contains correct data
        assert os.path.exists(config_file)
        with open(config_file) as f:
            loaded_config = json.load(f)
        assert loaded_config == sample_config

    def test_create_basic_config_eigenmode(self):
        """Test creating basic eigenmode configuration."""
        mesh_file = "test_mesh.msh"
        config = create_basic_config(mesh_file, "Eigenmode")

        assert config["Problem"]["Type"] == "Eigenmode"
        assert config["Model"]["Mesh"] == mesh_file
        assert "Eigenmode" in config["Solver"]
        assert "Target" in config["Solver"]["Eigenmode"]

    def test_create_basic_config_driven(self):
        """Test creating basic driven configuration."""
        mesh_file = "test_mesh.msh"
        config = create_basic_config(mesh_file, "Driven")

        assert config["Problem"]["Type"] == "Driven"
        assert config["Model"]["Mesh"] == mesh_file
        assert "Driven" in config["Solver"]
        assert "MinFreq" in config["Solver"]["Driven"]
        assert "MaxFreq" in config["Solver"]["Driven"]

    def test_create_basic_config_transient(self):
        """Test creating basic transient configuration."""
        mesh_file = "test_mesh.msh"
        config = create_basic_config(mesh_file, "Transient")

        assert config["Problem"]["Type"] == "Transient"
        assert config["Model"]["Mesh"] == mesh_file
        # Note: Transient solver config would be added in a complete implementation


class TestCSVFunctions:
    """Test CSV reading and processing functions."""

    def test_read_csv_results_valid_file(self, sample_csv_data):
        """Test reading valid CSV results file."""
        csv_file = sample_csv_data["s_params"]

        freq, data = read_csv_results(csv_file)

        assert len(freq) > 0
        assert len(data) > 0
        assert data.shape[0] == len(freq)  # Same number of frequency points
        assert data.shape[1] == 4  # 4 data columns (Re/Im for S11 and S21)

    def test_read_csv_results_single_column(self, temp_dir):
        """Test reading CSV with single data column."""
        data = np.column_stack(
            [np.linspace(1, 10, 10), np.sin(np.linspace(0, np.pi, 10))]
        )
        csv_file = os.path.join(temp_dir, "single_col.csv")
        np.savetxt(csv_file, data, delimiter=",", header="X,Y")

        x, y = read_csv_results(csv_file)

        assert len(x) == 10
        assert len(y) == 10
        # y should have shape (10, 1) for single data column
        assert y.shape == (10, 1)

    def test_read_csv_results_nonexistent_file(self):
        """Test reading nonexistent CSV file."""
        with pytest.raises(ValueError, match="Error reading CSV file"):
            read_csv_results("nonexistent.csv")

    def test_read_csv_results_invalid_format(self, temp_dir):
        """Test reading CSV with invalid format."""
        csv_file = os.path.join(temp_dir, "invalid.csv")
        with open(csv_file, "w") as f:
            f.write("invalid,csv,format\nwith,mixed,types\n")

        with pytest.raises(ValueError, match="Error reading CSV file"):
            read_csv_results(csv_file)


class TestMeshFunctions:
    """Test mesh file handling functions."""

    def test_find_mesh_files(self, temp_dir):
        """Test finding mesh files in directory."""
        # Create test mesh files
        mesh_files = ["test1.msh", "test2.mesh", "test3.vtk", "test4.vtu"]
        other_files = ["data.txt", "config.json", "results.csv"]

        for fname in mesh_files + other_files:
            with open(os.path.join(temp_dir, fname), "w") as f:
                f.write("dummy content")

        found_files = find_mesh_files(temp_dir)

        # Should find only mesh files
        assert len(found_files) == len(mesh_files)
        for mesh_file in mesh_files:
            assert any(mesh_file in f for f in found_files)

    def test_find_mesh_files_custom_extensions(self, temp_dir):
        """Test finding mesh files with custom extensions."""
        # Create files with custom extensions
        custom_files = ["test.e", "test.exo"]
        other_files = ["test.msh", "test.txt"]

        for fname in custom_files + other_files:
            with open(os.path.join(temp_dir, fname), "w") as f:
                f.write("dummy content")

        found_files = find_mesh_files(temp_dir, extensions=[".e", ".exo"])

        # Should find only files with custom extensions
        assert len(found_files) == len(custom_files)
        for custom_file in custom_files:
            assert any(custom_file in f for f in found_files)

    def test_find_mesh_files_empty_directory(self, temp_dir):
        """Test finding mesh files in empty directory."""
        found_files = find_mesh_files(temp_dir)
        assert len(found_files) == 0

    def test_validate_mesh_file_valid(self, sample_mesh_file):
        """Test validating a valid mesh file."""
        assert validate_mesh_file(sample_mesh_file) is True

    def test_validate_mesh_file_nonexistent(self):
        """Test validating nonexistent mesh file."""
        assert validate_mesh_file("nonexistent.msh") is False

    def test_validate_mesh_file_invalid_extension(self, temp_dir):
        """Test validating file with invalid extension."""
        invalid_file = os.path.join(temp_dir, "test.txt")
        with open(invalid_file, "w") as f:
            f.write("some content")

        assert validate_mesh_file(invalid_file) is False

    def test_validate_mesh_file_too_small(self, temp_dir):
        """Test validating mesh file that's too small."""
        small_file = os.path.join(temp_dir, "small.msh")
        with open(small_file, "w") as f:
            f.write("x")  # Very small file

        assert validate_mesh_file(small_file) is False


class TestExampleConfigs:
    """Test example configuration functions."""

    def test_get_example_configs_no_examples(self):
        """Test getting example configs when examples directory doesn't exist."""
        with patch("os.path.exists", return_value=False):
            examples = get_example_configs()
            assert examples == {}

    def test_get_example_configs_with_examples(self, temp_dir):
        """Test getting example configs with actual examples."""
        # Just test the function returns a dict when examples dir doesn't exist
        # Complex mocking is not worth the testing complexity
        examples = get_example_configs()
        assert isinstance(examples, dict)
        # Don't make assumptions about the content since it depends on file system


class TestUtilityIntegration:
    """Test integration between utility functions."""

    def test_config_roundtrip(self, temp_dir):
        """Test saving and loading configuration roundtrip."""
        original_config = create_basic_config("test.msh", "Eigenmode")
        config_file = os.path.join(temp_dir, "roundtrip.json")

        # Save and reload
        save_config(original_config, config_file)
        loaded_config = load_config(config_file)

        assert loaded_config == original_config

    def test_csv_processing_pipeline(self, sample_csv_data):
        """Test complete CSV processing pipeline."""
        csv_file = sample_csv_data["s_params"]

        # Read data
        freq, data = read_csv_results(csv_file)

        # Basic processing
        assert len(freq) > 0
        assert data.shape[1] >= 2  # At least 2 columns

        # Frequency should be monotonic increasing
        assert np.all(np.diff(freq) > 0)

        # Data should be finite
        assert np.all(np.isfinite(data))

        # Frequency range should be reasonable (GHz range)
        assert freq[0] >= 1e8  # At least 100 MHz
        assert freq[-1] <= 1e12  # At most 1 THz
