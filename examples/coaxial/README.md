# Coaxial example

This example simulates a coaxial cable with a 50 Ω port in one side and three cases on the other side: 1) open, 2) shorted and 3) matched.

## Running the example with Python

The example can be run with Python:

```bash
python coaxial.py
```

This should execute three simulations with *palace* and then produce plots for the voltage, the current and the energy stored in the electric and magnetic fields as a function of time comparing the three cases.

## Running the example with Julia

The example can be run with Julia by executing `coaxial.jl`.

## Running the example manually

You can simply run each of the simulations from the command line, for example for the open end case:

```bash
palace -np 1 coaxial_open.json
```

Then you will have to analyze the results yourself.
