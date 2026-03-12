## Molecular Dynamics and Thermostat Simulations

This directory contains three standalone C++ simulations covering diatomic molecules and thermostat-controlled particle dynamics, plus a notebook for plotting the generated results.

## Contents

- `molecules.cpp` - diatomic molecules in 3D with Lennard-Jones interactions
- `particles_w_thermostat.cpp` - particle simulation with a Nose-Hoover thermostat
- `particles_w_thermostat_extended.cpp` - longer thermostat run intended for equilibrium analysis
- `datafolder/` - generated CSV output from previous runs
- `plot_data.ipynb` - notebook for plots, diagnostics, and animations where supported
- `requirements.txt` - Python packages needed for the notebook
- `Molecular_&_Canonical_Dynamics.pdf` - accompanying report

## Building

Each simulation is compiled separately:

```bash
g++ -std=c++20 molecules.cpp -o molecules
g++ -std=c++20 particles_w_thermostat.cpp -o particles_w_thermostat
g++ -std=c++20 particles_w_thermostat_extended.cpp -o particles_w_thermostat_extended
```

## Running Simulations

### Molecules

`molecules` accepts an optional integer output index:

```bash
./molecules
./molecules 1
```

- Running without an index writes `molecule_data.csv`
- Running with an index writes `molecule_data_<index>.csv`

### Thermostat Particle Simulations

Both thermostat executables expect two arguments:

```bash
./particles_w_thermostat 1 0
./particles_w_thermostat_extended 2 1
```

Argument meanings:

- first argument: output file index, producing `particle_data_<index>.csv`
- second argument: exponent used for the thermostat mass parameter `Q = 10^exponent`

As with the particle project, most simulation parameters are declared near the top of each source file and are intended to be adjusted in code before recompiling.

## Plotting and Animation

Use `plot_data.ipynb` to inspect energies, temperatures, trajectories, and other observables from the CSV files in `datafolder/`. The notebook also supports animation workflows for the simulations that store full particle trajectories.

## Python Requirements

Install the notebook dependencies with:

```bash
pip install -r requirements.txt
```

The notebook currently depends on:

- `matplotlib`
- `numpy`
- `ipython`
- `scipy`
