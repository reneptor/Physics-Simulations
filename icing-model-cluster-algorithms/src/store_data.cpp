#include "ising_model.h"


namespace IsingModel{

    void WriteObservablesToFile(const std::string& folder_name, const std::string& filename, const IsingLattice& Lattice, std::vector<double> T_list){
        // Function to write lattice observables to csv file.

        std::cout << "writing to file: "<< filename << std::endl;
        std::filesystem::path base_path = std::filesystem::current_path();
        std::filesystem::path folder_path = base_path.parent_path() / "lattice_data/" / folder_name;
        std::filesystem::create_directories(folder_path);
        std::filesystem::path full_path = folder_path / (filename + ".csv");

        std::ofstream file(full_path);

        file << "# Stored observables for " << Lattice.GetLatticeSize() << "x" << Lattice.GetLatticeSize() << " lattice" << std::endl;
        file << "# CSV format with order:" << std::endl;
        file << "# Temperature" << std::endl;
        file << "# Energy" << std::endl;
        file << "# Heat capacity" << std::endl;
        file << "# Magnetization" << std::endl;
        file << "# Susceptibility" << std::endl;
        file << "# Binder cumulant" << std::endl;
        file << "# Lattice" << std::endl;

        for (auto& T : T_list){
            file << T << ',';
        }
        file << std::endl;
        for (auto& E : Lattice.E_measured){
            file << E << ',';
        }
        file << std::endl;
        for (auto& Cv : Lattice.Cv_measured){
            file << Cv << ',';
        }
        file << std::endl;
        for (auto& M : Lattice.M_measured){
            file << M << ',';
        }
        file << std::endl;
        for (auto& Chi : Lattice.Chi_measured){
            file << Chi << ',';
        }
        file << std::endl;
        for (auto& U_L : Lattice.U_L_measured){
            file << U_L << ',';
        }
        file << std::endl;

        Lattice.PrintLattice(file);

        file.close();
    }

    void WriteAutocorrMeasToFile(const std::string& folder_name, const std::string& filename, const std::vector<IsingLattice>& Lattices){
        // Function to write lattice autocorrelation measurements to csv file.

        std::cout << "writing to file: "<< filename << std::endl;
        std::filesystem::path base_path = std::filesystem::current_path();
        std::filesystem::path folder_path = base_path.parent_path() / "lattice_data/" / folder_name;
        std::filesystem::create_directories(folder_path);
        std::filesystem::path full_path = folder_path / (filename + ".csv");

        std::ofstream file(full_path);

        file << "# Stored autocorrelation measurements for " << Lattices.size() << " lattices." << std::endl;
        file << "# CSV format with order:" << std::endl;
        file << "# Number of lattices" << std::endl;
        file << "# Temperature" << std::endl;
        file << "# Deltas" << std::endl;
        file << "# Runtime" << std::endl;
        file << "# Steps" << std::endl;
        file << "# Lattice size" << std::endl;

        file << Lattices.size() << std::endl;
        for (auto& Lattice : Lattices){
            file << Lattice.GetTemperature() << std::endl;
            for (auto& delta : Lattice.deltas){
                file << delta << ',';
            }
            file << std::endl;
            file << Lattice.runtime << std::endl;
            file << Lattice.steps << std::endl;
            file << Lattice.GetLatticeSize() << std::endl;
        }
        

        file.close();
    }


}