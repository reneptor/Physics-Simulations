#include "ising_model.h"

std::atomic<int> task_counter(0);
int num_tasks;
unsigned int seed = 123456789;
std::vector<double> Temperatures;


namespace IsingModel{

    void SimulateAndWriteLattice(const std::string folder_name, int task_id, int lattice_size){
        // Function to run multiple lattice simulations in parallel.

        unsigned int this_seed = seed + task_id;
        IsingLattice Lattice = IsingModel::IsingLattice(lattice_size, this_seed);
        Lattice.Temperatures(Temperatures);
        Lattice.InitLatticeSpinUp();
        Lattice.Sweep(2000, 0); // Equilibrate to initial temperature.

        for (int T_index = 0; T_index < Temperatures.size(); T_index++){
            Lattice.Sweep(200, T_index);
            Lattice.MeasureForTemperature(100, T_index);
        }

        std::string filename = std::to_string(lattice_size) + "x" + std::to_string(lattice_size) 
        + "_lattice" + std::to_string(task_id) + ".csv";

        WriteDataToFile(folder_name, filename, Lattice, lattice_size);
    }

    void thread_task(const std::string folder_name, int thread_id, int lattice_size){
        while (true) {
            int task_id = task_counter.fetch_add(1); 
    
            // Stop when all tasks are assigned
            if (task_id >= num_tasks) {
                break;  
            }

            SimulateAndWriteLattice(folder_name, task_id, lattice_size);
            std::cout << "Simulated lattice " << task_id << " on thread " << thread_id << std::endl; 
        }
    }

    void SimulateMultipleLattices(const std::string folder_name, int lattice_size, int num_lattices){
    
        num_tasks = num_lattices;
        for (int i = 215; i <= 241; i++) {
            Temperatures.push_back(i * 0.01);
        }

        // Use maximum possible number of threads.
        const int num_threads = std::thread::hardware_concurrency();
        std::vector<std::thread> threads;

        for (int i = 0; i < num_threads; ++i) {
            threads.emplace_back(thread_task, folder_name, i, lattice_size);
        }

        for (auto& thread : threads) {
            thread.join();
        }
    
    }
    
}