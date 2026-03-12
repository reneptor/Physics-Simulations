
#ifndef ICING_H_
#define ICING_H_

#include <iostream>
#include <cstdlib>
#include <vector>
#include <map>
#include <random>
#include <algorithm>
#include <fstream>
#include <string>
#include <thread>
#include <atomic>
#include <filesystem>

namespace IcingModel {

    class IcingLattice {

        private:
            int lattice_size;  
            std::vector<double> T_list;
            std::vector<double> E_list;
            std::vector<double> M_list;
            int n_steps;
            int time_step;  
            int N;

            std::vector<std::vector<int>> lattice; 
            std::mt19937 rng;
            std::uniform_int_distribution<int> dist_init;  
            std::uniform_real_distribution<double> dist_flip; 
            std::uniform_int_distribution<int> dist_coord;   
            std::map<std::pair<int, int>, double> acceptance_probs; // Task 4
            
            // Functions.
            int SpinSum(int i, int j);
            std::pair<double, double> Observables();
            double AcceptanceProb(int i, int j, int T_index);
            
    
        public:
            double kB = 1.0;  
            double J = 1.0;  
            IcingLattice(int lattice_size, unsigned int rseed);
            void Temperatures(std::vector<double> Temps);
            std::vector<double> Temperatures() const;
            void InitLatticeRandomly();
            void InitLatticeSpinUp();
            void PrintLattice();
            void PrintLattice(std::ofstream& file) const;
            void Sweep(int iterations, int T_index);
            void MeasureForTemperature(int K, int T_index);

            std::vector<double> E_measured; // Energy
            std::vector<double> Cv_measured; // Heat capacity
            std::vector<double> M_measured; // Magnetization
            std::vector<double> Chi_measured; // Susceptibility
    };

    void WriteDataToFile(const std::string& folder_name, const std::string& filename, const IcingLattice& Lattice, int lattice_size);
    void SimulateAndWriteLattice(const std::string folder_name, int task_id, int lattice_size);
    void thread_task(const std::string folder_name, int thread_id, int lattice_size);
    void SimulateMultipleLattices(const std::string folder_name, int lattice_size, int num_lattices);

}  

#endif