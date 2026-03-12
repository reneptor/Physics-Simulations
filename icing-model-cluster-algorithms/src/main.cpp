#include "ising_model.h"


void MC_speeds(int L, double T, int flip_method, int equilibration_sweeps, int l_max, unsigned int seed){
    IsingModel::IsingLattice Lattice = IsingModel::IsingLattice(L, seed);
    Lattice.SetFlipMethod(flip_method);
    Lattice.InitLatticeSpinUp(1.0);
    Lattice.SpinAutoCorrelation(equilibration_sweeps, l_max, T);
    std::cout << "L = " << L << std::endl;
    std::cout << "T = " << T << std::endl;
    std::cout << "method = " << flip_method << std::endl;
    std::cout << "l_max = " << l_max << std::endl;
    std::cout << "runtime = " << Lattice.runtime << std::endl;

    double tau = Lattice.deltas[l_max]/Lattice.deltas[0];
    int steps = std::pow(2, l_max+5);
    double mc_speed = steps / (Lattice.runtime * tau);

    std::cout << "mc_speed = " << mc_speed << std::endl;
    std::cout << std::endl;
}


int main(int argc, char *argv[]){
    // Usage: ./ISING (task_id) (folder_name)

    int task_id = std::stoi(argv[1]);
    const std::string folder_name = argv[2];
    
    int lattice_size = 32;
    double T_C = 2/std::log(1 + std::sqrt(2));
    unsigned int seed = 123456789;

    if (task_id == 1){ // Task 2
        std::vector<double> Temperatures;
        for (int i = 15; i <= 30; i++) {
            Temperatures.push_back(i * 0.1);
        }
        
        // Metropolis
        IsingModel::IsingLattice MetropolisLattice = IsingModel::IsingLattice(lattice_size, seed);
        MetropolisLattice.SetFlipMethod(1); 
        MetropolisLattice.InitLatticeSpinUp(1.0);
        MetropolisLattice.SimulateLattice(Temperatures, 500);
        IsingModel::WriteObservablesToFile(folder_name, "metropolis_sim_lattice", MetropolisLattice, Temperatures);

        // Wolff
        IsingModel::IsingLattice WolffLattice = IsingModel::IsingLattice(lattice_size, seed);
        WolffLattice.SetFlipMethod(2); 
        WolffLattice.InitLatticeSpinUp(1.0);
        WolffLattice.SimulateLattice(Temperatures, 1);
        IsingModel::WriteObservablesToFile(folder_name, "wolff_sim_lattice", WolffLattice, Temperatures);

        // Swedsen-Wang
        IsingModel::IsingLattice SwedsenWangLattice = IsingModel::IsingLattice(lattice_size, seed);
        SwedsenWangLattice.SetFlipMethod(3); 
        SwedsenWangLattice.InitLatticeSpinUp(1.0);
        SwedsenWangLattice.SimulateLattice(Temperatures, 1);
        IsingModel::WriteObservablesToFile(folder_name, "swedsen_wang_sim_lattice", SwedsenWangLattice, Temperatures);
    }
    else if (task_id == 2){ // Task 3.a - plots
        lattice_size = 16;

        std::vector<double> T_autocorr = {T_C-1, T_C-0.5, T_C, T_C+0.5, T_C+1};
        std::vector<IsingModel::IsingLattice> metropolis_lattices;
        std::vector<IsingModel::IsingLattice> wolf_lattices;
        std::vector<IsingModel::IsingLattice> swedsen_wang_lattices;

        // Metropolis
        for (auto T : T_autocorr){
            std::cout << "metropolis" << std::endl;
            IsingModel::IsingLattice MetropolisLattice = IsingModel::IsingLattice(lattice_size, seed);
            MetropolisLattice.SetFlipMethod(1); 
            MetropolisLattice.InitLatticeSpinUp(1.0);
            MetropolisLattice.SpinAutoCorrelation(100, 18, T); // 100 equilibration sweeps for Metropolis.
            metropolis_lattices.push_back(MetropolisLattice);
            std::cout << std::endl;
        }

        // Wolff
        for (auto T : T_autocorr){
            std::cout << "wolff" << std::endl;
            IsingModel::IsingLattice WolffLattice = IsingModel::IsingLattice(lattice_size, seed);
            WolffLattice.SetFlipMethod(2); 
            WolffLattice.InitLatticeSpinUp(1.0);
            WolffLattice.SpinAutoCorrelation(10, 13, T); // 10 equilibration sweeps for Wolff.
            wolf_lattices.push_back(WolffLattice);
            std::cout << std::endl;
        }

        // Swedsen Wang
        for (auto T : T_autocorr){
            std::cout << "swedsen wang" << std::endl;
            IsingModel::IsingLattice SwedsenWangLattice = IsingModel::IsingLattice(lattice_size, seed);
            SwedsenWangLattice.SetFlipMethod(3); 
            SwedsenWangLattice.InitLatticeSpinUp(1.0);
            SwedsenWangLattice.SpinAutoCorrelation(10, 9, T); // 10 equilibration sweeps for Swedsen Wang.
            swedsen_wang_lattices.push_back(SwedsenWangLattice);
            std::cout << std::endl;
        }

        WriteAutocorrMeasToFile(folder_name, "metropolis_autocorr", metropolis_lattices);
        WriteAutocorrMeasToFile(folder_name, "wolff_autocorr", wolf_lattices);
        WriteAutocorrMeasToFile(folder_name, "swedsen_wang_autocorr", swedsen_wang_lattices);
    }
    else if (task_id == 3){ // task 3.a - mc speeds

        // Metropolis
        MC_speeds(16, T_C-1, 1, 100, 13, seed);
        MC_speeds(16, T_C-0.5, 1, 100, 15, seed);
        MC_speeds(16, T_C, 1, 100, 18, seed);
        MC_speeds(16, T_C+0.5, 1, 100, 17, seed);
        MC_speeds(16, T_C+1, 1, 100, 17, seed);

        // Wolff
        MC_speeds(16, T_C-1, 2, 10, 10, seed);
        MC_speeds(16, T_C-0.5, 2, 10, 9, seed);
        MC_speeds(16, T_C, 2, 10, 11, seed);
        MC_speeds(16, T_C+0.5, 2, 10, 11, seed);
        MC_speeds(16, T_C+1, 2, 10, 12, seed);

        // Swedsen Wang
        MC_speeds(16, T_C-1, 3, 10, 9, seed);
        MC_speeds(16, T_C-0.5, 3, 10, 1, seed);
        MC_speeds(16, T_C, 3, 10, 7, seed);
        MC_speeds(16, T_C+0.5, 3, 10, 2, seed);
        MC_speeds(16, T_C+1, 3, 10, 4, seed);

    }
    else if (task_id == 4){ // Task 3.b - plots
        std::vector<IsingModel::IsingLattice> metropolis_lattices;
        std::vector<IsingModel::IsingLattice> wolf_lattices;
        std::vector<IsingModel::IsingLattice> swedsen_wang_lattices;
        std::vector<int> lattice_sizes = {10, 20, 30};

        // Metropolis
        for (auto L : lattice_sizes){
            std::cout << "metropolis, L = " << L << std::endl;
            IsingModel::IsingLattice MetropolisLattice = IsingModel::IsingLattice(L, seed);
            MetropolisLattice.SetFlipMethod(1); 
            MetropolisLattice.InitLatticeSpinUp(1.0);
            MetropolisLattice.SpinAutoCorrelation(100, 16, T_C); // 100 equilibration sweeps for Metropolis.
            metropolis_lattices.push_back(MetropolisLattice);
            std::cout << std::endl;
        }

        // Wolff
        for (auto L : lattice_sizes){
            std::cout << "wolff, L = " << L << std::endl;
            IsingModel::IsingLattice WolffLattice = IsingModel::IsingLattice(L, seed);
            WolffLattice.SetFlipMethod(2); 
            WolffLattice.InitLatticeSpinUp(1.0);
            WolffLattice.SpinAutoCorrelation(10, 8, T_C); // 10 equilibration sweeps for Wolff.
            wolf_lattices.push_back(WolffLattice);
            std::cout << std::endl;
        }

        // Wolff
        for (auto L : lattice_sizes){
            std::cout << "swedsen wang, L = " << L << std::endl;
            IsingModel::IsingLattice SwedsenWangLattice = IsingModel::IsingLattice(L, seed);
            SwedsenWangLattice.SetFlipMethod(3); 
            SwedsenWangLattice.InitLatticeSpinUp(1.0);
            SwedsenWangLattice.SpinAutoCorrelation(10, 6, T_C); // 10 equilibration sweeps for Wolff.
            swedsen_wang_lattices.push_back(SwedsenWangLattice);
            std::cout << std::endl;
        }

        WriteAutocorrMeasToFile(folder_name, "metropolis_autocorr", metropolis_lattices);
        WriteAutocorrMeasToFile(folder_name, "wolff_autocorr", wolf_lattices);
        WriteAutocorrMeasToFile(folder_name, "swedsen_wang_autocorr", swedsen_wang_lattices);
    }
    else if (task_id == 5){ // task 3.b - mc speeds

        // Metropolis
        MC_speeds(10, T_C, 1, 100, 16, seed);
        MC_speeds(20, T_C, 1, 100, 16, seed);
        MC_speeds(30, T_C, 1, 100, 16, seed);

        // Wolff
        MC_speeds(10, T_C, 2, 10, 8, seed);
        MC_speeds(20, T_C, 2, 10, 7, seed);
        MC_speeds(30, T_C, 2, 10, 7, seed);

        // Swedsen Wang
        MC_speeds(10, T_C, 3, 10, 3, seed);
        MC_speeds(20, T_C, 3, 10, 6, seed);
        MC_speeds(30, T_C, 3, 10, 1, seed);
    }
    else{
        return 0; // Unknown task id.
    }


    return 0;
}


