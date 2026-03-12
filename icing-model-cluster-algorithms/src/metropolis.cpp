#include "ising_model.h"



namespace IsingModel{


    void IsingLattice::MetropolisFlip(){
        // Flip single spin based on M(RT)^2 method.
        int i = dist_coord(rng);
        int j = dist_coord(rng);
        double acceptance_prob = AcceptanceProb(i, j);
        double sampled_prob = dist_prob(rng);

        if (sampled_prob < acceptance_prob){
            lattice[i][j] *= -1;
            RecordMRealtime(i, j);
        }
    }


    double IsingLattice::AcceptanceProb(int i, int j){
        // Flip individual spin based on Metropolis method.
        int spins = SpinSum(i, j);
        std::pair<int, int> key = std::make_pair(spins, T_index);

        if (acceptance_probs.find(key) == acceptance_probs.end()){
            double dE = 2*J*spins;
            double prob = std::min(1.0, exp(-dE/(kB*T)));
            acceptance_probs[key] = prob; // Store in map to avoid repeated computations.
        }
        return acceptance_probs[key];
    }
    
    
    int IsingLattice::SpinSum(int i, int j){
        int h = 0;
        int L = lattice_size;
        // Sum neighbours with periodic B.C.
        h += lattice[(i+1)%L][j];
        h += lattice[(L+i-1)%L][j];
        h += lattice[i][(j+1)%L];
        h += lattice[i][(L+j-1)%L];

        return lattice[i][j]*h;
    }
}