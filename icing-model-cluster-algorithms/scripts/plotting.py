import numpy as np 
import matplotlib.pyplot as plt
from scipy.stats import norm
import sys
import os
import csv
from typing import List

T_C = 2/(np.log(1 + np.sqrt(2))) # Critical temperature.

class LatticeObservables:
    """Class to read and store lattice observables from simulations."""
    def __init__(self, name):
        self.name = name
        self.T = None
        self.E = None
        self.Cv = None
        self.M = None
        self.Chi = None
        self.lattice = None
        self.lattice_size = None


    def LoadFromFile(self, full_path):
        """Load lattice observables from file."""
        with open(full_path, "r") as file:
            csv_reader = csv.reader(file)

            # Skip header.
            for _ in range(9):
                next(csv_reader)
            
            self.T = np.array([float(c) for c in next(csv_reader)[:-1]])
            self.E = np.array([float(c) for c in next(csv_reader)[:-1]])
            self.Cv = np.array([float(c) for c in next(csv_reader)[:-1]])
            self.M = np.array([float(c) for c in next(csv_reader)[:-1]])
            self.Chi = np.array([float(c) for c in next(csv_reader)[:-1]])
            self.U_L = np.array([float(c) for c in next(csv_reader)[:-1]])

            lattice_row = [int(c) for c in next(csv_reader)[:-1]]
            self.lattice = [lattice_row]
            self.lattice_size = len(lattice_row)
            for _ in range(len(lattice_row) - 1):
                lattice_row = [int(c) for c in next(csv_reader)[:-1]]
                self.lattice.append(lattice_row)


class LatticeAutocorrelations:
    """Class to read and store lattice autocorrelation data from simulations."""
    def __init__(self, T, deltas, runtime, steps, lattice_size):
        self.T = T
        self.deltas = deltas
        self.runtime = runtime
        self.steps = steps
        self.lattice_size = lattice_size
        self.autocorr_time = np.amax(self.deltas)/self.deltas[0]
        self.mc_speed = steps / (self.autocorr_time * runtime)


def load_autocorrelations(full_path):
    print(full_path)
    """Load autocorrelation data for different temperatures/lattice sizes into array."""
    lattices = []
    with open(full_path, "r") as file:
        csv_reader = csv.reader(file)

        # Skip header.
        for _ in range(8):
            next(csv_reader)

        n_lattices = int(next(csv_reader)[0])
        
        for _ in range(n_lattices):
            T = float(next(csv_reader)[0])
            deltas = np.array([float(c) for c in next(csv_reader)[:-1]])
            runtime = float(next(csv_reader)[0])
            steps = int(next(csv_reader)[0])
            lattice_size = int(next(csv_reader)[0])

            lattice = LatticeAutocorrelations(T, deltas, runtime, steps, lattice_size)
            lattices.append(lattice)

    return lattices


def plot_autocorrelations_temps(lattices: List[LatticeAutocorrelations], method, output_path):
    temp_labels = ["$T_C-1.0$", "$T_C-0.5$", "$T_C$\t    ", "$T_C+0.5$", "$T_C+1.0$"]
    labels = [f"{temp_labels[i]}, $τ = $ {lattices[i].autocorr_time:.2f}" for i in range(len(lattices))]
    cmap = plt.get_cmap("viridis")

    fig, ax = plt.subplots(figsize=(8, 6))

    print(method)

    for i, lattice in enumerate(lattices):
        c = cmap((lattice.T - 1.27) / (2.27))  # Normalize color mapping
        l_range = np.arange(len(lattice.deltas))
        print(f"T = {lattice.T}, runtime = {lattice.runtime:.5f}, mc_speed = {lattice.mc_speed:.3f}")
        ax.plot(l_range, lattice.deltas, label=labels[i], color=c)
        ax.vlines(np.argmax(lattice.deltas), ymin=0, ymax=np.amax(lattice.deltas), linestyle="dashed", linewidth=1, color=c)
        ax.scatter(np.argmax(lattice.deltas), np.amax(lattice.deltas), color=c)


    ax.set_xlabel("level $L$")
    ax.set_ylabel("$\Delta_L$")
    ax.set_title("Autocorrelation via Binning Method for " + method)
    ax.legend()
    
    plt.savefig(output_path + method + "_autocorrelation.png", dpi=300)


def plot_autocorrelations_sizes(lattices: List[LatticeAutocorrelations], method, output_path):
    size_labels = ["10x10", "20x20", "30x30"]
    labels = [f"{size_labels[i]}, $τ = $ {lattices[i].autocorr_time:.2f}" for i in range(len(lattices))]
    cmap = plt.get_cmap("viridis")

    fig, ax = plt.subplots(figsize=(8, 6))

    print(method)

    for i, lattice in enumerate(lattices):
        c = cmap((lattice.lattice_size-10)/30)  # Normalize color mapping
        l_range = np.arange(len(lattice.deltas))
        print(f"L = {lattice.lattice_size}, runtime = {lattice.runtime:.5f}, mc_speed = {lattice.mc_speed:.3f}")
        ax.plot(l_range, lattice.deltas, label=labels[i], color=c)
        ax.vlines(np.argmax(lattice.deltas), ymin=0, ymax=np.amax(lattice.deltas), linestyle="dashed", linewidth=1, color=c)
        ax.scatter(np.argmax(lattice.deltas), np.amax(lattice.deltas), color=c)


    ax.set_xlabel("level $L$")
    ax.set_ylabel("$\Delta_L$")
    ax.set_title("Autocorrelation via Binning Method for " + method)
    ax.legend()
    
    plt.savefig(output_path + method + "_autocorrelation.png", dpi=300)


def plot_observables(lattices: List[LatticeObservables], output_path):
    """Plot observables against temperature in a combined 2x2 plot."""

    fig_M, ax_M = plt.subplots(figsize=(8, 6))
    for lattice in lattices:
        ax_M.plot(lattice.T, lattice.M, label=lattice.name)
    ax_M.set_xlabel("Temperature $T$")
    ax_M.set_ylabel("Absolute Magnetization $<|M|>$")
    ax_M.axvline(T_C, label="$T_C$", linestyle="dashed", color="red")
    ax_M.legend()
    ax_M.set_title(f"Magnetization Plot, {lattice.lattice_size}x{lattice.lattice_size} Lattices")
    fig_M.savefig(output_path + "magnetization.png", dpi=300)

    fig_Chi, ax_Chi = plt.subplots(figsize=(8, 6))
    for lattice in lattices:
        ax_Chi.plot(lattice.T, lattice.Chi, label=lattice.name)
    ax_Chi.set_xlabel("Temperature $T$")
    ax_Chi.set_ylabel("Suceptibility $χ$")
    ax_Chi.axvline(T_C, label="$T_C$", linestyle="dashed", color="red")
    ax_Chi.legend()
    ax_Chi.set_title(f"Susceptibility Plot, {lattice.lattice_size}x{lattice.lattice_size} Lattices")
    fig_Chi.savefig(output_path + "uceptibility.png", dpi=300)

    fig_U_L, ax_U_L = plt.subplots(figsize=(8, 6))
    for lattice in lattices:
        ax_U_L.plot(lattice.T, lattice.U_L, label=lattice.name)
    ax_U_L.set_xlabel("Temperature $T$")
    ax_U_L.set_ylabel("Binder cumulant $U_L$")
    ax_U_L.axvline(T_C, label="$T_C$", linestyle="dashed", color="red")
    ax_U_L.legend()
    ax_U_L.set_title(f"Binder cumulant Plot, {lattice.lattice_size}x{lattice.lattice_size} Lattices")
    fig_U_L.savefig(output_path + "binder_cumulant.png", dpi=300)



def main(input_path, task_id, output_path):

    os.system(f"rm -rf {output_path}")
    os.makedirs(output_path)

    if task_id == 1:
        metropolis_lattice = LatticeObservables("Metropolis")
        wolff_lattice = LatticeObservables("Wolff")
        swedsen_wang_lattice = LatticeObservables("Swedsen Wang")
        metropolis_lattice.LoadFromFile(input_path + "metropolis_sim_lattice.csv")
        wolff_lattice.LoadFromFile(input_path + "wolff_sim_lattice.csv")
        swedsen_wang_lattice.LoadFromFile(input_path + "swedsen_wang_sim_lattice.csv")
        plot_observables([metropolis_lattice, wolff_lattice, swedsen_wang_lattice], output_path)
    elif task_id == 2:
        metropolis_autocorrs = load_autocorrelations(input_path + "metropolis_autocorr.csv")    
        wolff_autocorrs = load_autocorrelations(input_path + "wolff_autocorr.csv")   
        swedsen_wang_autocorrs =  load_autocorrelations(input_path + "swedsen_wang_autocorr.csv")  
        plot_autocorrelations_temps(metropolis_autocorrs, "Metropolis", output_path)
        plot_autocorrelations_temps(wolff_autocorrs, "Wolff", output_path)
        plot_autocorrelations_temps(swedsen_wang_autocorrs, "Swedsen Wang", output_path)
    elif task_id == 3:
        pass
    elif task_id == 4:
        metropolis_autocorrs = load_autocorrelations(input_path + "metropolis_autocorr.csv")    
        wolff_autocorrs = load_autocorrelations(input_path + "wolff_autocorr.csv")   
        swedsen_wang_autocorrs =  load_autocorrelations(input_path + "swedsen_wang_autocorr.csv")  
        plot_autocorrelations_sizes(metropolis_autocorrs, "Metropolis", output_path)
        plot_autocorrelations_sizes(wolff_autocorrs, "Wolff", output_path)
        plot_autocorrelations_sizes(swedsen_wang_autocorrs, "Swedsen Wang", output_path)
    else:
        return # Unknown task_id.


if __name__ == "__main__":
    """Usage: python3 plotting.py (input_path) (task_id) (output_path)"""

    base_path = os.getcwd()
    folder_path = os.path.dirname(base_path)

    input_path = folder_path + "/lattice_data/" + sys.argv[1] + "/"
    task_id = int(sys.argv[2])
    output_path = folder_path + "/output/" + sys.argv[3] + "/"

    main(input_path, task_id, output_path)

