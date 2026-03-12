
#ifndef Ising_H_
#define Ising_H_

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
#include <chrono>
#include <cmath>
#include <filesystem>



namespace IsingModel {

    class IsingLattice {

        private:
            int lattice_size;  
            double T; // Temperature
            int T_index; // For storing acceptance probabilities.
            double bond_prob; // For Wolff updates.
            int time_step;  
            int N;
            int M; // Live magnetization tracker, more efficient than Observables().
            int flip_method;

            std::vector<std::vector<int>> lattice;
            std::mt19937 rng;
            std::uniform_real_distribution<double> dist_prob; 
            std::uniform_int_distribution<int> dist_init; 
            std::uniform_int_distribution<int> dist_coord;   
            std::map<std::pair<int, int>, double> acceptance_probs; // Task 4

            // Functions.
            int SpinSum(int i, int j);
            void SetTemperature(double T_new, int index);
            std::pair<double, double> Observables();
            double AcceptanceProb(int i, int j);
            std::vector<std::pair<int, int>> GetSiteNeighbours(std::pair<int, int> coord);
            void AddBondsRecursively(const std::pair<int, int>& coord);
            void MetropolisFlip();
            void WolffFlip();
            void SwedsenWangFlip();
            int FindOriginalCluster(int k, int *M);
            void Sweep(int K);
            void MeasureForTemperature(int K, int I); 
            void RecordMRealtime(int i, int j);
            std::vector<double> CalculateDeltas(int max_l);
            void (IsingLattice::*FlipMethod)();                
    
        public:
            double kB = 1.0;  
            double J = 1.0;  
            // Variables for auto-correlation measurements.
            std::vector<double> deltas;
            double runtime;
            unsigned int steps;

            IsingLattice(int lattice_size, unsigned int rseed);
            
            void InitLatticeSpinUp(double p_spin_up);
            void PrintLattice();
            void PrintLattice(std::ofstream& file) const;
            void SetFlipMethod(int flip_method);
            void SimulateLattice(std::vector<double> T_list, int K);
            void SpinAutoCorrelation(int equilibration_sweeps, int max_l, double T);

            double GetTemperature() const;
            int GetLatticeSize() const;
       
            std::vector<double> E_measured; // Energy
            std::vector<double> Cv_measured; // Heat capacity
            std::vector<double> M_measured; // Magnetization
            std::vector<double> Chi_measured; // Susceptibility
            std::vector<double> U_L_measured; // Binder cumulant
    };

    void MC_speeds(int L, double T, int flip_method, int equilibration_sweeps, int l_max, unsigned int seed);
    void WriteObservablesToFile(const std::string& folder_name, const std::string& filename, const IsingLattice& Lattice, std::vector<double> T_list);
    void WriteAutocorrMeasToFile(const std::string& folder_name, const std::string& filename, const std::vector<IsingLattice>& Lattices);
} 

#endif