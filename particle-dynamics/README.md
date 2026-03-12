## Particle Dynamics Simulation

This directory contains a standalone C++ particle simulation and a notebook for visualizing the saved trajectories and observables.

## Contents

- `particles.cpp` - simulation of particles in a 3D box with periodic boundary conditions and Lennard-Jones interactions
- `datafolder/` - generated CSV output from previous runs
- `plot_data.ipynb` - notebook for plotting observables and building animations from saved trajectory data
- `requirements.txt` - Python packages needed for the notebook
- `Particle_Dynamics.pdf` - accompanying report

## Building

The simulation is a single source file and can be compiled directly with `g++`:

```bash
g++ -std=c++20 particles.cpp -o particles
```

## Running Simulations

The executable accepts an optional integer output index:

```bash
./particles
./particles 1
```

- Running without an index writes `particle_data.csv`
- Running with an index writes `particle_data_<index>.csv`, for example `particle_data_1.csv`

Most physical and numerical parameters are defined as compile-time constants near the top of `particles.cpp`, so changing the simulation setup usually means editing the source and recompiling.

## Plotting and Animation

Use `plot_data.ipynb` to inspect trajectories, energies, and generated observables. The notebook is also the place where particle animations are created from the saved CSV output in `datafolder/`.

## Python Requirements

Install the notebook dependencies with:

```bash
pip install -r requirements.txt
```

The notebook currently depends on:

- `matplotlib`
- `numpy`
- `ipython`
