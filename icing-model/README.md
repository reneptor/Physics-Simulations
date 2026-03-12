## Ising Model Simulation

This directory contains a C++ implementation of the 2D Ising model, a Python plotting script, generated lattice data, and output figures used in the accompanying report.

## Contents

- `src/` - C++ source files for running simulations and exporting observables
- `include/` - model header files
- `scripts/` - Python plotting utilities
- `build/` - CMake build directory
- `lattice_data/` - generated CSV output from simulations
- `output/` - saved plots and lattice visualizations
- `Icing_Model.pdf` - accompanying report

## Building

Configure and build the executable from the project directory:

```bash
cmake -S . -B build
cmake --build build
```

This produces the executable `ICING` inside `build/`.

## Running Simulations

The executable supports two modes:

```bash
./build/ICING <folder_name> <lattice_size>
./build/ICING <folder_name> <lattice_size> <num_lattices>
```

Examples:

```bash
./build/ICING 10x10 10
./build/ICING multiple10 10 50
```

- the two-argument form runs a single lattice simulation and writes one output file
- the three-argument form runs multiple lattices and stores several outputs in the chosen folder

Some simulation settings, such as temperature ranges and sweep counts, are currently defined in the source code and should be adjusted there before rebuilding.

## Generating Plots

Use the Python plotting script from `scripts/` to visualize saved data:

```bash
python3 scripts/plotting.py 10x10 lattice spin_plot
```

The exact plotting arguments depend on the saved dataset and plot type you want to generate.
