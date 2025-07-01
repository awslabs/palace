# Palace Python Test Suite

This directory contains comprehensive tests for the Palace Python interface, ensuring that all functionality works correctly and examples are robust.

## Test Structure

```
tests/
├── conftest.py              # Pytest configuration and shared fixtures
├── pytest.ini               # Pytest settings
├── test_runner.py           # Standalone test runner script
├── unit/                    # Unit tests for core functionality
│   ├── test_core.py        # Tests for palace.core module
│   └── test_utils.py       # Tests for palace.utils module
├── integration/             # Integration tests for examples
│   └── test_examples.py    # Tests for example scripts
└── data/                   # Test data files (created during tests)
```

## Running Tests

### Method 1: Using the Test Runner (Recommended)

The test runner provides different levels of testing:

```bash
# Run quick example verification (fastest)
python tests/test_runner.py --examples

# Run unit tests only
python tests/test_runner.py --unit

# Run integration tests only
python tests/test_runner.py --integration

# Run all tests (default)
python tests/test_runner.py
```

### Method 2: Using pytest directly

If pytest is installed:

```bash
# Run all tests
pytest

# Run specific test categories
pytest tests/unit/           # Unit tests only
pytest tests/integration/    # Integration tests only

# Run with verbose output
pytest -v

# Run specific test file
pytest tests/unit/test_core.py
```

### Method 3: Using python -m pytest

```bash
# From the python/ directory
python -m pytest tests/

# With specific options
python -m pytest tests/ -v --tb=short
```

## Test Categories

### Unit Tests (`tests/unit/`)

Test individual components in isolation:

- **`test_core.py`**: Tests for `palace.core` module
  - `PalaceSolver` class functionality
  - `run_palace()` convenience function
  - Error handling and edge cases
  - Mock testing without requiring Palace installation

- **`test_utils.py`**: Tests for `palace.utils` module
  - Configuration file management
  - CSV result reading and processing
  - Mesh file utilities
  - Example configuration functions

### Integration Tests (`tests/integration/`)

Test complete workflows and example scripts:

- **`test_examples.py`**: Tests for all example scripts
  - Import verification for all examples
  - Configuration creation functionality
  - Analysis function testing
  - Main function execution
  - Error handling and graceful degradation

## Test Features

### Robust Testing Without Dependencies

The tests are designed to work even when:
- Palace is not installed
- Matplotlib is not available
- Example mesh files are missing

This is achieved through:
- Mock objects for Palace execution
- Dummy plotting classes when matplotlib is unavailable
- Generated example data when real results are missing

### Example Verification

The test suite includes comprehensive verification that:
- All example scripts can be imported successfully
- Configuration creation functions work correctly
- Analysis functions handle both real and simulated data
- Main functions execute without errors
- Error handling is robust

### Fixtures and Test Data

The test suite includes:
- **Sample configurations**: Valid Palace JSON configurations
- **Mock mesh files**: Minimal mesh files for testing
- **Sample CSV data**: Realistic simulation result data
- **Mock executables**: Fake Palace executables for testing

## Writing New Tests

### Adding Unit Tests

1. Create test functions starting with `test_`
2. Use fixtures from `conftest.py` for common test data
3. Mock external dependencies (Palace, matplotlib)
4. Test both success and failure cases

Example:
```python
def test_my_function(temp_dir, sample_config):
    """Test description."""
    # Setup
    config_file = os.path.join(temp_dir, "test.json")
    with open(config_file, 'w') as f:
        json.dump(sample_config, f)

    # Test
    result = my_function(config_file)

    # Assert
    assert result is not None
    assert result['status'] == 'success'
```

### Adding Integration Tests

1. Test complete workflows from import to execution
2. Verify examples work with and without dependencies
3. Check file creation and output generation
4. Test error handling and recovery

### Test Guidelines

- **Independence**: Tests should not depend on each other
- **Cleanup**: Use temporary directories for file operations
- **Mocking**: Mock external dependencies to ensure reliable testing
- **Error Testing**: Test both success and failure scenarios
- **Documentation**: Include clear docstrings for test functions

## Continuous Integration

The test suite is designed to work in CI environments:

- No GUI dependencies (matplotlib displays are mocked)
- No external software requirements (Palace is mocked)
- Fast execution (most tests complete in seconds)
- Clear failure reporting

## Test Data and Fixtures

### Generated Test Data

Tests automatically generate realistic data:
- S-parameter matrices with proper complex values
- Eigenmode frequencies and Q-factors
- Time-domain pulse signals
- Mesh file content

### Fixtures Available

From `conftest.py`:
- `temp_dir`: Temporary directory for test files
- `sample_config`: Valid Palace configuration
- `sample_mesh_file`: Minimal test mesh
- `sample_csv_data`: Realistic CSV result files
- `mock_palace_executable`: Fake Palace executable

## Debugging Test Failures

### Common Issues

1. **Import Errors**: Check that Python path includes the package
2. **File Not Found**: Verify test is using temporary directories
3. **Mock Failures**: Ensure external dependencies are properly mocked
4. **Platform Issues**: Some file operations may behave differently on Windows/Linux

### Debug Tips

```bash
# Run with maximum verbosity
python tests/test_runner.py --examples -v

# Run single test with debugging
pytest tests/unit/test_core.py::TestPalaceSolver::test_init_with_executable_path -v -s

# Check test discovery
pytest --collect-only
```

## Performance

The test suite is optimized for speed:
- **Example verification**: ~5 seconds
- **Unit tests**: ~10 seconds
- **Integration tests**: ~15 seconds
- **Full suite**: ~30 seconds

Most time is spent in example imports and matplotlib mocking rather than actual computation.

## Contributing

When adding new functionality:

1. Add corresponding unit tests
2. Update integration tests if examples are affected
3. Ensure tests pass without external dependencies
4. Update this README if new test categories are added

The goal is to maintain 100% test coverage of the Python interface while ensuring examples work reliably for users.
