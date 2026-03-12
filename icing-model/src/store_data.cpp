#include "icing_model.h"


namespace IcingModel{

    void WriteDataToFile(const std::string& folder_name, const std::string& filename, const IcingLattice& Lattice, int lattice_size){
        // Function to write Lattice instance to csv file.

        std::cout << "writing to file: "<< filename << std::endl;

        std::filesystem::path base_path = std::filesystem::current_path();

        std::filesystem::path folder_path = base_path.parent_path() / "lattice_data/" / folder_name;
        std::filesystem::create_directories(folder_path);
        std::filesystem::path full_path = folder_path / filename;

        std::ofstream file(full_path);

        file << "# Stored data for " << lattice_size << "x" << lattice_size << " lattice" << std::endl;
        file << "# CSV format with order:" << std::endl;
        file << "# Temperature" << std::endl;
        file << "# Energy" << std::endl;
        file << "# Heat capacity" << std::endl;
        file << "# Magnetization" << std::endl;
        file << "# Susceptibility" << std::endl;
        file << "# Lattice" << std::endl;

        for (auto& T : Lattice.Temperatures()){
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

        Lattice.PrintLattice(file);

        file.close();
    }


}