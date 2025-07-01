# Palace Python Interface Examples

This directory contains comprehensive examples demonstrating how to use Palace through its Python interface for electromagnetic simulations.

## Overview

The Palace Python interface provides:
- **Easy configuration management**: Programmatically create and modify simulation setups
- **Simulation control**: Run Palace simulations from Python scripts
- **Result analysis**: Post-process and visualize simulation results
- **Integration**: Seamless integration with the Python scientific ecosystem

## Examples

### 1. Basic Usage (`basic_usage.py`)
**Beginner-friendly introduction to the Palace Python interface**

- Simple simulation execution using existing configurations
- Using the `PalaceSolver` class for more control
- Creating configurations programmatically
- Basic result post-processing and visualization

```bash
python basic_usage.py
```

**Key concepts:**
- `run_palace()` convenience function
- `PalaceSolver` class
- Configuration validation
- CSV result reading

### 2. Interactive Tutorial (`palace_tutorial.ipynb`)
**Comprehensive Jupyter notebook tutorial**

An interactive notebook covering:
- Getting started with Palace Python interface
- Configuration file creation and customization
- Running simulations programmatically
- Advanced post-processing techniques
- Parameter sweeps and automation
- Integration with matplotlib and scientific Python

```bash
jupyter notebook palace_tutorial.ipynb
```

### 3. Eigenmode Analysis (`eigenmode_analysis.py`)
**Cavity resonators and dielectric resonator analysis**

- Cavity resonator eigenmode simulations
- Dielectric resonator configurations
- Mode analysis and visualization
- Q-factor extraction
- Parameter studies with material properties

```bash
python eigenmode_analysis.py
```

**Features:**
- Automatic mode classification
- Energy participation analysis
- Frequency scaling laws
- Material property sweeps

### 4. Frequency Domain Simulations (`frequency_domain_simulation.py`)
**S-parameter extraction and filter analysis**

- Two-port network S-parameter extraction
- Microwave filter analysis
- Antenna simulation setup
- Advanced S-parameter post-processing
- Smith chart visualization
- Adaptive frequency sampling

```bash
python frequency_domain_simulation.py
```

**Applications:**
- Transmission line characterization
- Filter design verification
- Antenna impedance analysis
- Broadband characterization

### 5. Time Domain Simulations (`time_domain_simulation.py`)
**Transient analysis and pulse propagation**

- Time domain transient simulations
- Time Domain Reflectometry (TDR)
- Pulse propagation analysis
- Eye diagram generation
- Pulse compression studies

```bash
python time_domain_simulation.py
```

**Advanced features:**
- FFT-based frequency domain analysis
- Group delay calculation
- Digital signal integrity analysis
- Chirped pulse compression

## Quick Start

1. **Install Palace** following the [official documentation](https://awslabs.github.io/palace/dev/install/)

2. **Install the Python package**:
   ```bash
   cd /path/to/palace
   pip install -e .
   ```

3. **Run a basic example**:
   ```bash
   python basic_usage.py
   ```

4. **Try the interactive tutorial**:
   ```bash
   jupyter notebook palace_tutorial.ipynb
   ```

## Configuration Examples

### Eigenmode Simulation
```python
from palace.utils import create_basic_config, save_config

config = create_basic_config("cavity.msh", "Eigenmode")
config["Solver"]["Eigenmode"]["Target"] = 5e9  # 5 GHz target
config["Solver"]["Eigenmode"]["MaxSize"] = 20  # Find 20 modes

save_config(config, "cavity_eigenmode.json")
```

### Frequency Domain Simulation
```python
config = create_basic_config("filter.msh", "Driven")
config["Solver"]["Driven"] = {
    "MinFreq": 1e9,   # 1 GHz
    "MaxFreq": 10e9,  # 10 GHz
    "FreqStep": 0.1e9 # 100 MHz steps
}

# Add ports
config["Model"]["Boundary"] = [
    {"Index": 1, "LumpedPort": {"R": 50.0, "Excitation": True}},
    {"Index": 2, "LumpedPort": {"R": 50.0}}
]
```

### Time Domain Simulation
```python
config = create_basic_config("line.msh", "Transient")
config["Solver"]["Transient"] = {
    "MaxTime": 5e-9,   # 5 ns simulation
    "TimeStep": 1e-12, # 1 ps time step
    "Order": 2
}

# Gaussian pulse excitation
config["Model"]["Excitation"] = {
    "Type": "Gaussian",
    "Width": 100e-12,  # 100 ps width
    "Center": 200e-12  # Centered at 200 ps
}
```

## Running Simulations

### Method 1: Convenience Function
```python
from palace import run_palace

result = run_palace("config.json", num_procs=4, output_dir="results/")
if result.returncode == 0:
    print("Simulation completed successfully!")
```

### Method 2: Solver Class
```python
from palace import PalaceSolver

solver = PalaceSolver()
if solver.validate_config("config.json"):
    result = solver.run("config.json", output_dir="results/", num_procs=4)
```

## Post-Processing Results

### Reading CSV Results
```python
from palace.utils import read_csv_results
import matplotlib.pyplot as plt

# Read S-parameter results
freq, s_data = read_csv_results("results/port-S.csv")
s11_real, s11_imag = s_data[:, 0], s_data[:, 1]
s11_mag = np.sqrt(s11_real**2 + s11_imag**2)

# Plot return loss
plt.plot(freq/1e9, 20*np.log10(s11_mag))
plt.xlabel('Frequency (GHz)')
plt.ylabel('|S11| (dB)')
plt.show()
```

### Eigenmode Analysis
```python
# Read eigenfrequency results
freq, data = read_csv_results("results/domain-E.csv")
eigenfreqs = data[:, 0]  # Eigenfrequencies
q_factors = data[:, 1]   # Q factors

# Find highest Q mode
best_mode = np.argmax(q_factors)
print(f"Best mode: {eigenfreqs[best_mode]/1e9:.3f} GHz, Q = {q_factors[best_mode]:.0f}")
```

## Generated Files

Running the examples will create various output files:

**Configuration files:**
- `*_config.json` - Palace configuration files
- `*_eigenmode.json` - Eigenmode simulation setups
- `*_driven.json` - Frequency domain setups
- `*_transient.json` - Time domain setups

**Analysis results:**
- `*_analysis.json` - Analysis summaries and metrics
- `example_*.csv` - Sample data files for demonstration

**Visualizations:**
- `*.png` - Result plots and analysis figures
- Smith charts, S-parameter plots, eigenmode analysis

## Tips and Best Practices

### Configuration Management
- Use `create_basic_config()` as a starting point
- Validate configurations with `solver.validate_config()`
- Save configurations for reproducibility
- Use descriptive filenames for different scenarios

### Simulation Execution
- Always check return codes: `result.returncode == 0`
- Use appropriate number of MPI processes for your system
- Specify output directories to organize results
- Monitor simulation progress through log files

### Post-Processing
- Start with example data when developing analysis scripts
- Use appropriate frequency/time ranges for visualization
- Apply window functions for FFT analysis
- Save analysis results for documentation

### Performance Optimization
- Use adaptive frequency sampling for resonant structures
- Optimize time steps for transient simulations
- Consider mesh refinement for accuracy
- Parallel processing for parameter sweeps

## Advanced Topics

### Parameter Sweeps
```python
materials = ["silicon", "gaas", "quartz"]
permittivities = [11.7, 12.9, 3.8]

for mat, eps in zip(materials, permittivities):
    config = create_basic_config(f"{mat}_device.msh", "Eigenmode")
    config["Model"]["Domain"][0]["Material"]["Permittivity"] = eps
    save_config(config, f"sweep_{mat}.json")

    # Run simulation
    result = run_palace(f"sweep_{mat}.json", output_dir=f"results_{mat}/")
```

### Batch Processing
```python
import multiprocessing as mp

def run_single_simulation(config_file):
    return run_palace(config_file, output_dir=f"results_{config_file[:-5]}/")

config_files = ["config1.json", "config2.json", "config3.json"]

# Run simulations in parallel
with mp.Pool(processes=3) as pool:
    results = pool.map(run_single_simulation, config_files)
```

### Custom Analysis
```python
def analyze_filter_response(s_param_file):
    freq, data = read_csv_results(s_param_file)
    s21_mag = np.sqrt(data[:, 4]**2 + data[:, 5]**2)  # S21 magnitude

    # Find passband
    max_transmission = np.max(s21_mag)
    passband_3db = s21_mag >= max_transmission / np.sqrt(2)

    return {
        'center_freq': freq[np.argmax(s21_mag)],
        'bandwidth_3db': freq[passband_3db][-1] - freq[passband_3db][0],
        'insertion_loss': 20 * np.log10(max_transmission),
        'passband_ripple': 20 * np.log10(np.max(s21_mag[passband_3db]) / np.min(s21_mag[passband_3db]))
    }
```

## Troubleshooting

### Common Issues

1. **Palace not found**: Ensure Palace is installed and in your PATH
2. **Configuration errors**: Use `validate_config()` to check syntax
3. **Mesh file missing**: Verify mesh file paths in configurations
4. **Memory issues**: Reduce mesh size or use fewer MPI processes
5. **Convergence problems**: Adjust solver tolerances and iterations

### Getting Help

- Check the [Palace documentation](https://awslabs.github.io/palace/)
- Review configuration schema in `scripts/schema/config-schema.json`
- Examine existing examples in the main `examples/` directory
- Use `palace --help` for command-line options

## Contributing

To add new examples:

1. Follow the existing code structure and documentation style
2. Include comprehensive comments and docstrings
3. Provide both demonstration data and real simulation examples
4. Add appropriate visualization and analysis
5. Update this README with your new example

## References

- [Palace Documentation](https://awslabs.github.io/palace/)
- [Palace GitHub Repository](https://github.com/awslabs/palace)
- [MFEM Finite Element Library](https://mfem.org/)
- [Palace Configuration Reference](https://awslabs.github.io/palace/dev/config/)
