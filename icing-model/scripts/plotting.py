import numpy as np 
import matplotlib.pyplot as plt
from scipy.stats import norm
import sys
import os
import csv
from typing import List


class LatticeData:
    """Class to read and store lattice data from simulations."""
    def __init__(self):
        self.T = None
        self.E = None
        self.Cv = None
        self.M = None
        self.Chi = None
        self.lattice = None
        self.lattice_size = 0


    def load_from_file(self, full_path):
        """Load lattice observables from file."""
        with open(full_path, "r") as file:
            csv_reader = csv.reader(file)

            # Skip header.
            for _ in range(8):
                next(csv_reader)
            
            self.T = np.array([float(c) for c in next(csv_reader)[:-1]])
            self.E = np.array([float(c) for c in next(csv_reader)[:-1]])
            self.Cv = np.array([float(c) for c in next(csv_reader)[:-1]])
            self.M = np.array([float(c) for c in next(csv_reader)[:-1]])
            self.Chi = np.array([float(c) for c in next(csv_reader)[:-1]])

            lattice_row = [int(c) for c in next(csv_reader)[:-1]]
            self.lattice = [lattice_row]
            self.lattice_size = len(lattice_row)
            for _ in range(len(lattice_row) - 1):
                lattice_row = [int(c) for c in next(csv_reader)[:-1]]
                self.lattice.append(lattice_row)


def plot_lattice(lattice: LatticeData, output_folder_path):
    """Plots a 2D grid where 1 is red and -1 is blue."""

    lattice_size = len(lattice.lattice[0])
    colors = {1: 'red', -1: 'blue'}
    fig, ax = plt.subplots(figsize=(6, 6))

    for i in range(lattice_size):
        for j in range(lattice_size):
            color = colors[lattice.lattice[i][j]]
            ax.add_patch(plt.Rectangle((j, lattice_size-i-1), 1, 1, color=color))

    tick_positions = [i for i in range(0, lattice_size+1, lattice_size//10)]
    ax.set_xticks(tick_positions)
    ax.set_yticks(tick_positions)
    ax.set_xticklabels(tick_positions)
    ax.set_yticklabels(tick_positions)
    ax.set_title(f"{lattice_size}x{lattice_size} Lattice, T = {lattice.T[-1]}")
    ax.set_xlim(0, lattice_size)
    ax.set_ylim(0, lattice_size)
    
    plt.gca().set_aspect('equal')
    if not os.path.exists(output_folder_path):
        os.makedirs(output_folder_path)
    plt.savefig(output_folder_path + f"{lattice_size}x{lattice_size}_lattice_plot.png", dpi=300)


def plot_observables(lattice: LatticeData, output_folder_path):
    """Plot observables against temperature in a combined 2x2 plot."""

    fig, ax = plt.subplots(nrows=2, ncols=2, figsize=(10, 8))

    ax[0, 0].plot(lattice.T, lattice.E)
    ax[0, 0].set_title("Energy")
    ax[0, 1].plot(lattice.T, lattice.Cv)
    ax[0, 1].set_title("Heat Capacity")
    ax[1, 0].plot(lattice.T, lattice.M)
    ax[1, 0].set_title("Magnetization")
    ax[1, 1].plot(lattice.T, lattice.Chi)
    ax[1, 1].set_title("Susceptibility")

    plt.tight_layout()
    if not os.path.exists(output_folder_path):
        os.makedirs(output_folder_path)
    plt.savefig(output_folder_path + "observables_plot.png", dpi=300)


def plot_observables_varying_sizes(lattices: List[LatticeData], output_folder_path):
    """Plot observables for varying lattice sizes against temperature in a combined 2x2 plot."""

    fig, ax = plt.subplots(nrows=2, ncols=2, sharex=True, figsize=(10, 8))
    T_C = 2/np.log(1 + np.sqrt(2))

    for lattice in lattices:

        ax[0, 0].plot(lattice.T, lattice.E/lattice.lattice_size**2, label=f"{lattice.lattice_size}x{lattice.lattice_size}")
        ax[0, 0].set_title("Mean Energy")
        ax[0, 0].set_ylabel("<E>")
        ax[0, 1].plot(lattice.T, lattice.Cv/lattice.lattice_size**2, label=f"{lattice.lattice_size}x{lattice.lattice_size}")
        ax[0, 1].set_title("Heat Capacity")
        ax[0, 1].set_ylabel("Cv")
        ax[1, 0].plot(lattice.T, lattice.M, label=f"{lattice.lattice_size}x{lattice.lattice_size}")
        ax[1, 0].set_title("Mean Absolute Value Magnetization")
        ax[1, 0].set_ylabel("<|M|>")
        ax[1, 1].plot(lattice.T, lattice.Chi, label=f"{lattice.lattice_size}x{lattice.lattice_size}")
        ax[1, 1].set_title("Susceptibility")
        ax[1, 1].set_ylabel("χ")
        ax[1, 0].set_xlabel("Temperature T")
        ax[1, 1].set_xlabel("Temperature T")
    
    ax[0, 0].axvline(T_C, label="T_crit", color="purple")
    ax[0, 1].axvline(T_C, label="T_crit", color="purple")
    ax[1, 0].axvline(T_C, label="T_crit", color="purple")
    ax[1, 1].axvline(T_C, label="T_crit", color="purple")

    plt.legend()
    # plt.tight_layout()
    fig.suptitle("Icing Model Observables, Varying Sizes")
    if not os.path.exists(output_folder_path):
        os.makedirs(output_folder_path)
    plt.savefig(output_folder_path + "varying_sizes_plot.png", dpi=300)



def plot_avg_response_functions(lattices: List[LatticeData], output_folder_path):
    """Plot average graphs of multiple lattice simulations"""

    num_lattices = len(lattices)

    Temperatures = lattices[0].T

    Cv = np.zeros((num_lattices, lattices[0].Cv.shape[0]))
    Chi = np.zeros((num_lattices, lattices[0].Chi.shape[0]))

    for i in range(num_lattices):
        Cv[i, :] = lattices[i].Cv
        Chi[i, :] = lattices[i].Chi

    mean_Cv = np.mean(Cv, axis=0)
    std_Cv = np.std(Cv, axis=0)
    mean_Chi = np.mean(Chi, axis=0)
    std_Chi = np.std(Chi, axis=0)

    fig, ax = plt.subplots(ncols=2, figsize=(10, 6), sharex=True)

    # Mean
    ax[0].plot(Temperatures, mean_Cv, label="Mean Cv", color="blue")
    ax[1].plot(Temperatures, mean_Chi, label="Mean χ", color="blue")

    # Upper error bound
    ax[0].plot(Temperatures, mean_Cv+std_Cv, label="Error Margin Cv", color="orange")
    ax[1].plot(Temperatures, mean_Chi+std_Chi, label="Error Margin χ", color="orange")

    #Lower error bound
    ax[0].plot(Temperatures, mean_Cv-std_Cv, color="orange")
    ax[1].plot(Temperatures, mean_Chi-std_Chi, color="orange")

    ax[0].legend(loc='upper left')
    ax[1].legend(loc='upper left')

    ax[0].set_xlabel("Temperature T")
    ax[1].set_xlabel("Temperature T")
    ax[0].set_ylabel("Heat Capacity Cv")
    ax[1].set_ylabel("Susceptibility χ")

    fig.suptitle(f"Simulated Response Functions for {num_lattices} Lattices")

    if not os.path.exists(output_folder_path):
        os.makedirs(output_folder_path)
    plt.savefig(output_folder_path + "avg_plot.png", dpi=300)


def critical_temp_histogram(lattices: List[LatticeData], output_folder_path):
    """Plot to estimate critical temperature"""

    Temperatures = lattices[0].T

    temperature_indices_Cv = []
    temperature_indices_Chi = []

    for lattice in lattices:
        temperature_indices_Cv.append(np.argmax(lattice.Cv))
        temperature_indices_Chi.append(np.argmax(lattice.Chi))

    plt.clf()
    fig, ax = plt.subplots(ncols=2, figsize=(10, 6), sharex=True)

    Temperatures_Cv = Temperatures[temperature_indices_Cv]
    Temperatures_Chi = Temperatures[temperature_indices_Chi]

    mean_Cv_T, std_Cv_T = norm.fit(Temperatures_Cv)
    mean_Chi_T, std_Chi_T = norm.fit(Temperatures_Chi)

    gaussian_Cv = norm.pdf(Temperatures, mean_Cv_T, std_Cv_T)
    gaussian_Chi = norm.pdf(Temperatures, mean_Chi_T, std_Chi_T)

    ax[0].hist(Temperatures_Cv, label="Temperature Histogram for Cv")
    ax[1].hist(Temperatures_Chi, label="Temperature Histogram for χ")

    ax[0].plot(Temperatures, gaussian_Cv, label=f"Gaussian Fit, μ = {mean_Cv_T:.3f}, σ = {std_Cv_T:.3f}")
    ax[1].plot(Temperatures, gaussian_Chi, label=f"Gaussian Fit, μ = {mean_Chi_T:.3f}, σ = {std_Chi_T:.3f}")

    ax[0].set_xlabel("Temperature T")
    ax[1].set_xlabel("Temperature T")

    ax[0].set_ylabel("Frequency")
    ax[1].set_ylabel("Frequency")

    ax[0].legend()
    ax[1].legend()

    fig.suptitle("Temperature Histogram to estimate T_crit")

    if not os.path.exists(output_folder_path):
        os.makedirs(output_folder_path)
    plt.savefig(output_folder_path + "temperature_hist.png", dpi=300)


def main(lattices, plot_method, output_folder_path):

    os.system(f"rm -rf {output_folder_path}")

    if plot_method == "observables":
        plot_observables(lattices[0], output_folder_path)
    elif plot_method == "lattice":
        plot_lattice(lattices[0], output_folder_path)
    elif plot_method == "varying":
        plot_observables_varying_sizes(lattices, output_folder_path)
    elif plot_method == "multiple":
        plot_avg_response_functions(lattices, output_folder_path)
        critical_temp_histogram(lattices, output_folder_path)
    else:
        print("Unknown plot command.")


if __name__ == "__main__":
    """Usage: python3 plotting.py (input_path) (plot_method) (output_path)"""

    base_path = os.getcwd()
    folder_path = os.path.dirname(base_path)

    input_path = folder_path + "/lattice_data/" + sys.argv[1] + "/"
    plot_method = sys.argv[2]
    output_path = folder_path + "/output/" + sys.argv[3] + "/"

    print(input_path)
    print(output_path)

    lattices = []
    for filename in os.listdir(input_path):
        full_path = input_path + filename
        lattice = LatticeData()
        lattice.load_from_file(full_path)
        lattices.append(lattice)

    main(lattices, plot_method, output_path)

