#include "ising_model.h"
#include <cassert>


namespace IsingModel{


    void IsingLattice::SwedsenWangFlip(){
        // Did not figure out how to implement with periodic b.c.

        int L = lattice_size;
        // 4D mapping of bonds for efficient access.
        std::vector<std::vector<std::vector<std::vector<bool>>>> bond_mapping(
            L, std::vector<std::vector<std::vector<bool>>>(
                L, std::vector<std::vector<bool>>(
                    L, std::vector<bool>(L, false)
                )
            )
        );

        // Iterate over each bond.
        for (int i = 0; i < L; i++){
            for (int j = 0; j < L; j++){
                std::vector<std::pair<int, int>> neighbours = {
                    {i+1, j}, // right
                    {i, j+1} // bottom
                };
                for (auto &neighbour : neighbours){
                    if (neighbour.first == L || neighbour.second == L){ // Outside bounds.
                        continue;
                    } 
                    if (lattice[i][j] == lattice[neighbour.first][neighbour.second]){
                        double sampled_prob = dist_prob(rng);
                        if (sampled_prob < bond_prob){
                            bond_mapping[i][j][neighbour.first][neighbour.second] = true;
                            bond_mapping[neighbour.first][neighbour.second][i][j] = true;
                        }
                    }
                }

                
                
            }
        }

    
        // Run Hoshen-Kopelman.
        // std::vector<int> M;
        int *M = new int[N]();
        int k = 2;
        M[k] = 1;
        std::vector<std::vector<int>> labeled_lattice(L, std::vector<int>(L, 0));
        labeled_lattice[0][0] = k;
        for (int i = 0; i < L; i++){
            for (int j = 0; j < L; j++){

                if (i == 0 && j == 0){ // Already added first site.
                    continue;
                }

                bool top = (i > 0) ? bond_mapping[i][j][i-1][j] : false;
                bool left = (j > 0) ? bond_mapping[i][j][i][j-1] : false;

                if (left == false && top == false){
                    k += 1;
                    labeled_lattice[i][j] = k;
                    M[k] = 1;
                }
                else if (left == true && top == false){
                    int k_0 = FindOriginalCluster(labeled_lattice[i][j-1], M);
                    labeled_lattice[i][j] = k_0;
                    M[k_0] += 1;
                }
                else if (left == false && top == true){
                    int k_0 = FindOriginalCluster(labeled_lattice[i-1][j], M);
                    labeled_lattice[i][j] = k_0;
                    M[k_0] += 1;
                }
                else if (left == true && top == true){
                    int k_1 = FindOriginalCluster(labeled_lattice[i-1][j], M); // top
                    int k_2 = FindOriginalCluster(labeled_lattice[i][j-1], M); // left
                    if (k_1 != k_2){
                        labeled_lattice[i][j] = k_1;
                        M[k_1] += M[k_2] + 1; // Prioritize top.
                        M[k_2] = -k_1;
                    }
                    else {
                        labeled_lattice[i][j] = k_1;
                        M[k_1] += 1;
                    }
                }
            }
        }

        // Assign which clusters to flip.
        for (int l = 2; l < k; l++){ 
            if (M[l] > 0){
                double sampled_prob = dist_prob(rng);
                M[l] = (sampled_prob < 0.5) ? 0 : 1; // 0 means ignore, 1 means flip. 
            }
        }

        // Perform flips.
        k = 0;
        int flip = 0;
        for (int i = 0; i < L; i++){
            for (int j = 0; j < L; j++){
                if (labeled_lattice[i][j] != k){
                    k = labeled_lattice[i][j];
                    flip = M[FindOriginalCluster(labeled_lattice[i][j], M)];
                }
                if(flip == 1){ 
                    lattice[i][j] *= -1; 
                    RecordMRealtime(i, j);
                }
            }
        }
        delete M;   
    }


    int IsingLattice::FindOriginalCluster(int k, int *M){
        int iters = 0;
        while (M[k] < 0){
            k = -M[k];
            iters += 1;
            assert(iters < 100);
        }
        return k;
    }
}