#include "ising_model.h"


namespace IsingModel{


    void IsingLattice::WolffFlip(){
        // Flip spin cluster based on Wolff method.
        int i = dist_coord(rng);
        int j = dist_coord(rng);
        const std::pair<int, int> initial_coord(i, j);
        lattice[i][j] *= -1; // Always flip first cluster spin.
        RecordMRealtime(i, j);
    
        std::vector<std::pair<int, int>> initial_neighbours = GetSiteNeighbours(initial_coord);
        for (const auto& coord : initial_neighbours){ // FLip each spin in cluster.
            AddBondsRecursively(coord);
        }
    }


    void IsingLattice::AddBondsRecursively(const std::pair<int, int>& coord){
        // Form cluster by recursive dfs search of similar lattice sites,
        // based on bond probabilities.
        int i = coord.first; 
        int j = coord.second;
        double sampled_prob = dist_prob(rng);

        if(sampled_prob < bond_prob){
            lattice[i][j] *= -1;
            RecordMRealtime(i, j);
        }
        else{
            return;
        }

        std::vector<std::pair<int, int>> site_neighbors = GetSiteNeighbours(coord);
        
        for (const auto& neigh_coord : site_neighbors) {
            int ni = neigh_coord.first;
            int nj = neigh_coord.second;
            if (lattice[i][j] != lattice[ni][nj]) { // Expand towards spins that are anti-aligned. Boundary will take care of it self.
                AddBondsRecursively(neigh_coord);
            }
        }
    }


    std::vector<std::pair<int, int>> IsingLattice::GetSiteNeighbours(std::pair<int, int> coord) {
        // Yield neighbouring indices with periodic B.C.
        int L = lattice_size;
        int i = coord.first;
        int j = coord.second;
        std::vector<std::pair<int, int>> site_neighbors = {
            {(L+i-1)%L, j}, 
            {(i+1)%L, j}, 
            {i, (L+j-1)%L}, 
            {i, (j+1)%L} 
        };
        return site_neighbors;
    };
}












