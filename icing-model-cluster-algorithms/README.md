## Ising Model Simulation with Cluster Algorithms

This directory extends the Ising model project with multiple update strategies, including Metropolis, Wolff, and Swendsen-Wang style cluster methods. It contains the C++ simulation code, Python plotting utilities, generated datasets, and the accompanying report.

## Contents

- `src/` - C++ implementations of the simulation and update algorithms
- `include/` - shared model header files
- `scripts/` - Python plotting utilities
- `build/` - CMake build directory
- `lattice_data/` - generated CSV output from simulations
- `output/` - saved figures
- `Cluster_Algorithms_for_Ising_Model.pdf` - accompanying report

## Building

Configure and build the executable from the project directory:

```bash
cmake -S . -B build
cmake --build build
```

This produces the executable `ISING` inside `build/`.

## Running Simulations

Run the compiled executable with a task identifier and output folder name:

```bash
./build/ISING <task_id> <folder_name>
```

Example:

```bash
./build/ISING 1 observables
```

Simulation parameters such as lattice size, temperature range, and the chosen update method are currently set in `src/main.cpp` and should be edited there before rebuilding.

## Generating Plots

Use the plotting script to visualize saved observables:

```bash
python3 scripts/plotting.py <input_folder> <task_id> <output_folder>
```

Example:

```bash
python3 scripts/plotting.py observables 1 observables_plot
```
