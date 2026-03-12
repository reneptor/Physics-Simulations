#include "ising_model.h"



namespace IsingModel{

    // Constructor
    IsingLattice::IsingLattice(int L, unsigned int seed) : 
    // Class for simulating single lattice of size LxL.

    rng(seed),
    dist_prob(0.0f, 1.0f),
    dist_init(-1, 1),
    dist_coord(0, L-1)
    {
        lattice_size = L;
        N = lattice_size*lattice_size;
    }


    void IsingLattice::RecordMRealtime(int i, int j){
        // Call after flip has been commited.
        M += 2*lattice[i][j];
    }


    void IsingLattice::InitLatticeSpinUp(double p_spin_up){
        // Init spins as up given probability p_spin_up.
        double sampled_prob = dist_prob(rng);
        int spin = (sampled_prob < p_spin_up) ? 1 : -1;
        lattice.assign(lattice_size, std::vector<int>(lattice_size, spin)); 

        auto [_, M_init] = Observables();
        M = M_init;
    }


    void IsingLattice::SetTemperature(double T_new, int index){
        T = T_new;
        T_index = index;
        bond_prob = 1-exp(-2*J/(kB*T_new));
    }


    double IsingLattice::GetTemperature() const{
        return T;
    }


    int IsingLattice::GetLatticeSize() const{
        return lattice_size;
    }


    void IsingLattice::PrintLattice() {
        // Prints a 10x10 downsized version of the lattice.
        for (int i = 0; i < lattice_size; i += lattice_size / 10){  
            for (int j = 0; j < lattice_size; j += lattice_size / 10){  
                std::cout << lattice[i][j] << ' ';
            }
            std::cout << '\n';
        }
    }


    void IsingLattice::PrintLattice(std::ofstream& file) const {
        // Prints full lattice to given csv file.
        for (int i = 0; i < lattice_size; i += 1){  
            for (int j = 0; j < lattice_size; j += 1){  
                file << lattice[i][j] << ',';
            }
            file << '\n';
        }
    }


    void IsingLattice::Sweep(int K){

        int steps = (flip_method == 1) ? N*K : lattice_size*K; // Reduce steps for cluster methods.
        for (int k = 0; k < steps; k++){ // Sweep N*K steps to decouple from lattice size.
            (this->*FlipMethod)();
        }
    }

    
    std::pair<double, double> IsingLattice::Observables(){
        // Measure energy and magnetization.
        double E = 0;
        int M = 0;
        int L = lattice_size;
        for (int i = 0; i < L; i++){
            for (int j = 0; j < L; j++){
                M += lattice[i][j];
                E += -J * lattice[i][j] * 
                  (lattice[(i+1) % L][j] + lattice[i][(j+1) % L]); 
            }
        }
        return {E, M};
    }


    std::vector<double> IsingLattice::CalculateDeltas(int max_l){
        // Assume lattice has already been equilibrated to temperature.

        std::vector<double> deltas;
        int M_0 = std::pow(2, max_l + 5); // + 5 such that minimum number of bins is 32.
        double *A_0 = new double[M_0];
        double A_0_mean = 0;
        double delta_0 = 0;

        for (int i = 0; i < M_0; i++){
            (this->*FlipMethod)();
            double m = (double)M/N;
            A_0[i] = m;
            A_0_mean += m;
            delta_0 += std::pow(m, 2);
        }
        A_0_mean /= M_0;
        delta_0 /= M_0;
        delta_0 -= std::pow(A_0_mean, 2);
        delta_0 /= M_0-1;
        delta_0 = std::sqrt(delta_0);
        deltas.push_back(delta_0);

        std::cout << "printing deltas" << std::endl;
        std::cout << "l="<< 0 << ", " << delta_0 << std::endl;

        double *A_l_prev = A_0;
        int M_l = M_0;

        for (int l = 1; l <= max_l; l++){

            double delta_l = 0; 
            M_l /= 2;
            double *A_l = new double[M_l];
            double A_l_mean = 0;
            for (int i = 0; i < M_l; i++){ // Calculate delta and new array in one go;
                A_l[i] = (A_l_prev[2*i] + A_l_prev[2*i + 1])/2;
                delta_l += std::pow(A_l[i], 2);
                A_l_mean += A_l[i];
            }
            A_l_mean /= M_l;
            delta_l /= M_l;
            delta_l -= std::pow(A_l_mean, 2);
            delta_l /= M_l-1;
            delta_l = std::sqrt(delta_l);

            std::cout << "l="<< l << ", " << delta_l << std::endl;
            deltas.push_back(delta_l);

            // Reassignment
            delete A_l_prev;
            A_l_prev = A_l;
        } 
        delete A_l_prev;
        return deltas;
    }   



    void IsingLattice::SpinAutoCorrelation(int equilibration_sweeps, int max_l, double T){
        // Measure runtime of deltas calculation and store results for this class instance.
        steps = std::pow(2, max_l + 5);
        SetTemperature(T, 0);
        std::cout << "equilibrating to T=" << T << std::endl;
        Sweep(equilibration_sweeps); 
        std::cout << "calculating deltas" << std::endl;
        std::cout << "M: " << M << std::endl;

        auto start = std::chrono::high_resolution_clock::now();
        deltas = CalculateDeltas(max_l);
        auto end = std::chrono::high_resolution_clock::now();

        std::chrono::duration<double> diff = end - start;
        runtime = diff.count();
        std::cout << "finnished in " << runtime << " seconds" << std::endl;
    } 


    void IsingLattice::MeasureForTemperature(int K, int I){
        // Record observables based on measurements for E and M.

        double E_mean = 0.0;
        double E2_mean = 0.0;
        double M_mean = 0.0;
        double M2_mean = 0.0;
        double M4_mean = 0.0;
        double M_abs_mean = 0.0;

        for (int i = 0; i < I; i++){
            Sweep(K);
            auto [E_obs, M_obs] = Observables();
            M_obs /= N;
            E_mean += E_obs;
            E2_mean += std::pow(E_obs, 2);
            M_mean += M_obs;
            M2_mean += std::pow(M_obs, 2);
            M4_mean += std::pow(M_obs, 4);
            M_abs_mean += std::abs(M_obs);
        }

        E_mean /= I;
        E2_mean /= I;
        M_mean /= I;
        M2_mean /= I;
        M4_mean /= I;
        M_abs_mean /= I;

        double Cv = (E2_mean - E_mean*E_mean)/(kB*T*kB*T);
        double Chi = N*(M2_mean - M_abs_mean*M_abs_mean)/(kB*T);
        double U_L = 1 - M4_mean/(3*std::pow(M2_mean, 2));

        E_measured.push_back(E_mean);
        Cv_measured.push_back(Cv);
        M_measured.push_back(M_abs_mean); // Normed to [0, 1]
        Chi_measured.push_back(Chi);
        U_L_measured.push_back(U_L);
    }


    void IsingLattice::SetFlipMethod(int flip_method){
        // Toggle between Metropolis and Wolff for spin updates.
        flip_method = flip_method;
        if (flip_method == 1){
            FlipMethod = &IsingLattice::MetropolisFlip;
        }
        else if (flip_method == 2){
            FlipMethod = &IsingLattice::WolffFlip;
        }
        else if (flip_method == 3){
            FlipMethod = &IsingLattice::SwedsenWangFlip;
        }
        else {
            return; // Unknown flip method index.
        }    
    }


    void IsingLattice::SimulateLattice(std::vector<double> T_list, int K){
        // Run full simulation given list of temperatures.
        SetTemperature(T_list[0], 0);

        Sweep(5*K); // Equilibrate to initial temperature.
     
        for (int T_index = 0; T_index < T_list.size(); T_index++){

            SetTemperature(T_list[T_index], T_index);
            Sweep(K); 
            MeasureForTemperature(K, 10);
    
            // Print temperature and lattice to keep track of simulation progress.
            std::cout << "Temperature: " << T_list[T_index] << std::endl;
            PrintLattice();
            std::cout << std::endl;
        }
    }


}


