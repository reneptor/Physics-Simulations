#include "icing_model.h"



namespace IcingModel{


    IcingLattice::IcingLattice(int L, unsigned int seed) : 
    // Class for simulating single lattice of size LxL.

    rng(seed),
    dist_init(0, 1),
    dist_flip(0.0f, 1.0f),
    dist_coord(0, L-1)
    {
        lattice_size = L;
        N = lattice_size*lattice_size;
        lattice.assign(lattice_size, std::vector<int>(lattice_size, 1));
        InitLatticeSpinUp();
    }


    void IcingLattice::Temperatures(const std::vector<double> Temps){
        T_list = Temps;
    }
    

    std::vector<double> IcingLattice::Temperatures() const {
        return T_list;
    }


    void IcingLattice::InitLatticeRandomly(){
        lattice.assign(lattice_size, std::vector<int>(lattice_size, 0)); 

        for (auto& row : lattice) {
            for (auto& elem : row) {
                elem = (dist_init(rng) == 0) ? -1 : 1;
            }
        }
    }


    void IcingLattice::InitLatticeSpinUp(){
        lattice.assign(lattice_size, std::vector<int>(lattice_size, 1)); 
    }


    void IcingLattice::PrintLattice() {
        // Prints a 10x10 downsized version of the lattice.
        for (int i = 0; i < lattice_size; i += lattice_size / 10){  
            for (int j = 0; j < lattice_size; j += lattice_size / 10){  
                std::cout << lattice[i][j] << ' ';
            }
            std::cout << '\n';
        }
    }

    void IcingLattice::PrintLattice(std::ofstream& file) const {
        // Prints full lattice to given csv file.
        for (int i = 0; i < lattice_size; i += 1){  
            for (int j = 0; j < lattice_size; j += 1){  
                file << lattice[i][j] << ',';
            }
            file << '\n';
        }
    }
    

    int IcingLattice::SpinSum(int i, int j){
        int h = 0;
        int L = lattice_size;
        // Sum neighbours with periodic B.C.
        h += lattice[(i+1)%L][j];
        h += lattice[(L+i-1)%L][j];
        h += lattice[i][(j+1)%L];
        h += lattice[i][(L+j-1)%L];

        return lattice[i][j]*h;
    }


    double IcingLattice::AcceptanceProb(int i, int j, int T_index){
        int spins = SpinSum(i, j);
        double T = T_list[T_index];
        std::pair<int, int> key = std::make_pair(spins, T_index);

        if (acceptance_probs.find(key) == acceptance_probs.end()){
            double dE = 2*J*spins;
            double prob = std::min(1.0, exp(-dE/(kB*T)));
            acceptance_probs[key] = prob; // Store in map to avoid repeated computations.
        }
        return acceptance_probs[key];
    }


    void IcingLattice::Sweep(int K, int T_index){
 
        for (int k = 0; k < N*K; k++){ // Sweep N*K steps to decouple from lattice size.
            int i = dist_coord(rng);
            int j = dist_coord(rng);
            double acceptance_prob = AcceptanceProb(i, j, T_index);
            double sampled_prob = dist_flip(rng);
            if (sampled_prob < acceptance_prob){

                lattice[i][j] *= -1;
            }
        }
    }


    std::pair<double, double> IcingLattice::Observables(){
        // Measure energy and magnetization.
        double E = 0;
        double M = 0;
        for (int i = 0; i < lattice_size; i++){
            for (int j = 0; j < lattice_size; j++){
                M += lattice[i][j];
                E += -J * lattice[i][j] * 
                  (lattice[(i+1) % lattice_size][j] + lattice[i][(j+1) % lattice_size]); 
            }
        }
        return {E, M/N};
    }


    void IcingLattice::MeasureForTemperature(int K, int T_index){
        // Record observables based on measurements for E and M.
    
        double T = T_list[T_index];
        double E_mean = 0.0;
        double E2_mean = 0.0;
        double M_mean = 0.0;
        double M2_mean = 0.0;
        double M_abs_mean = 0.0;

        for (int k = 0; k < K; k++){
            Sweep(10, T_index);
            auto [E, M] = Observables();
            E_mean += E;
            E2_mean += E*E;
            M_mean += M;
            M2_mean += M*M;
            M_abs_mean += std::abs(M);
        }

        E_mean /= K;
        E2_mean /= K;
        M_mean /= K;
        M2_mean /= K;
        M_abs_mean /= K;

        double Cv = (E2_mean - E_mean*E_mean)/(kB*T*kB*T);
        double Chi = N*(M2_mean - M_mean*M_mean)/(kB*T);

        E_measured.push_back(E_mean);
        Cv_measured.push_back(Cv);
        M_measured.push_back(M_abs_mean);
        Chi_measured.push_back(Chi);
    }


}


