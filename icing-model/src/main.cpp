#include "icing_model.h"



int main(int argc, char *argv[]){
    // Usage: ./ICING (folder_name) (lattice size); 
    //        ./ICING (folder_name) (lattice size) (num_lattices)

    const std::string folder_name = argv[1];

    if(argc == 3){
        unsigned int seed = 123456789;
        int lattice_size = std::stoi(argv[2]);
        const std::string filename = std::string(argv[2]) + "x" + std::string(argv[2]) + "_lattice.csv";
    
        IcingModel::IcingLattice Lattice = IcingModel::IcingLattice(lattice_size, seed);
        
        std::vector<double> Temperatures;
        for (int i = 15; i <= 30; i++) {
            Temperatures.push_back(i * 0.1);
        }
    
        Lattice.Temperatures(Temperatures);
        Lattice.InitLatticeSpinUp();
        Lattice.Sweep(2000, 0); // Equilibrate to initial temperature.
     
        for (int T_index = 0; T_index < Temperatures.size(); T_index++){

            Lattice.Sweep(200, T_index);
            Lattice.MeasureForTemperature(100, T_index);
    
            // Print temperature and lattice to keep track of simulation progress.
            std::cout << "Temperature: " << Temperatures[T_index] << std::endl;
            Lattice.PrintLattice();
            std::cout << std::endl;
        }
    
        std::string folder_path = "/home/alex_unix/csp25ex/ex01/lattice_data/" + folder_name + "/";

        IcingModel::WriteDataToFile(folder_name, filename, Lattice, lattice_size);
    }
    else if(argc == 4){
        int lattice_size = std::stoi(argv[2]);
        int num_lattices = std::stoi(argv[3]);

        IcingModel::SimulateMultipleLattices(folder_name, lattice_size, num_lattices);
    }
    else {
        std::cout << "Too many arguments.";
    }

    return 0;
}



