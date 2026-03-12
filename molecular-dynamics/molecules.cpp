#include <array>
#include <tuple>
#include <cmath>
#include <random>
#include <iostream>
#include <unordered_set>
#include <fstream>
#include <numbers>


constexpr float time_span = 1.0f; // Time in seconds - dont change
constexpr float dt = 0.0001f; // Time-step size
constexpr int T = static_cast<int>(time_span/dt); // Number of time-steps
constexpr float m = 1.0f;
constexpr float sigma = 5.0f; // LJ potential strength
constexpr int N = 30; // Number of particles
constexpr float L = 10.0f; // Box size
constexpr float r_c2 = 2.5f*2.5f; // Cut-off radius squared
constexpr float d = 0.5f; // Intra-molecular distance between particles
constexpr float d2 = 0.25f; // Square of intra-molecular distance between particles
constexpr float pi = 3.14159265358f; 


using std::sqrt, std::pow, std::log, std::cos, std::sin, std::acos, std::round, std::fmod;


struct Particle{

    unsigned int particle_id; // For iterations 
    double x_new {}, y_new {}, z_new {}; // New time-step position
    double x_curr {}, y_curr {}, z_curr {}; // Current time-step position
    double x_prev {}, y_prev {}, z_prev {}; // Previous time-step position
    double vx {}, vy {}, vz {}; // Velocity
    double fx {}, fy {}, fz {}; // Inter molecular force
    double gx {}, gy {}, gz {}; // Intra molecular force
    float Ekin, Epot {}; // Individual particle energies
    Particle *twin; // Pointer to other particle in molecule
    
    Particle() {} // For particle array initialization
    Particle(unsigned int id) : particle_id(id) {} // Actual initialization


    void InitialPosition(double x_0, double y_0, double z_0){
        x_curr = x_0;
        y_curr = y_0;
        z_curr = z_0;
    }


    void InitialVelocity(double vx_0, double vy_0, double vz_0){
        x_prev = x_curr - vx_0*dt;
        y_prev = y_curr - vy_0*dt;
        z_prev = z_curr - vz_0*dt;
    }


    double EnergyKinetic(){
        Ekin = 0.5 * m * (vx*vx + vy*vy + vz*vz);
        return Ekin;
    }
};


struct RelativePosition{

    double x_rel_curr {};
    double y_rel_curr {};
    double z_rel_curr {};
    double x_rel_prev {};
    double y_rel_prev {};
    double z_rel_prev {};
    
    void CalculateRelativePosition(const Particle& p1, const Particle& p2){ // Calculations with periodic b.c.
        x_rel_prev = x_rel_curr;
        y_rel_prev = y_rel_curr;
        z_rel_prev = z_rel_curr;
        x_rel_curr = p2.x_curr - p1.x_curr - L*round((p2.x_curr - p1.x_curr)/L);
        y_rel_curr = p2.y_curr - p1.y_curr - L*round((p2.y_curr - p1.y_curr)/L);
        z_rel_curr = p2.z_curr - p1.z_curr - L*round((p2.z_curr - p1.z_curr)/L); 
    }

    double r2_curr() const {
        return x_rel_curr*x_rel_curr + y_rel_curr*y_rel_curr + z_rel_curr*z_rel_curr;
    }

    double r2_prev() const {
        return x_rel_prev*x_rel_prev + y_rel_prev*y_rel_prev + z_rel_prev*z_rel_prev;
    }

    double r_curr_prev () const {
        return x_rel_curr*x_rel_prev + y_rel_curr*y_rel_prev + z_rel_curr*z_rel_prev;
    }
};


struct Molecule{

    unsigned int molecule_id;
    Particle *p1 = nullptr; // Pointer to first particle
    Particle *p2 = nullptr; // Pointer to second particle
    double lambda {}; // Spring constant for bond
    float Ekin {}; // Kinetic energy
    float Epot {}; // Potential energy
    float Ebond {}; // Bond energy


    Molecule() {} // For molecule array initialization
    Molecule(unsigned int id, Particle *particle1, Particle *particle2) 
        : molecule_id(id), p1(particle1), p2(particle2) { // Actual initialization
        p1->twin = p2; // Set twin pointers
        p2->twin = p1; 
    }


    template<typename pos, typename vec> 
    void InitialPosition(pos initial_position, vec axis_orientation){
        auto [nx, ny, nz] = axis_orientation; 
        auto [x_0, y_0, z_0] = initial_position;
        p1->InitialPosition(x_0 + 0.5*nx*d, y_0 + 0.5*ny*d, z_0 + 0.5*nz*d);
        p2->InitialPosition(x_0 - 0.5*nx*d, y_0 - 0.5*ny*d, z_0 - 0.5*nz*d);
    }


    template<typename vel>
    void InitialVelocity(vel initial_velocity){
        auto [vx_0, vy_0, vz_0] = initial_velocity;
        p1->InitialVelocity(vx_0, vy_0, vz_0);
        p2->InitialVelocity(vx_0, vy_0, vz_0);
    }


    double MolecularBond(const RelativePosition& rp){
        double r2_c = rp.r2_curr();
        double r2_p = rp.r2_prev();
        double r_c_p = rp.r_curr_prev();    

        double lambda = (sqrt(r_c_p*r_c_p - r2_p*(r2_c - d2)) - r_c_p)/(2*dt*dt*r2_p/m); // Derived formula for diatomic mulecule
        Ebond = 0.5*lambda*pow(sqrt(r2_c) - d, 2); // Spring energy formula

        // Calculate intra molecular forces
        p1->gx = -lambda*rp.x_rel_curr;
        p1->gy = -lambda*rp.y_rel_curr;
        p1->gz = -lambda*rp.z_rel_curr;
        p2->gx = lambda*rp.x_rel_curr;
        p2->gy = lambda*rp.y_rel_curr;
        p2->gz = lambda*rp.z_rel_curr;

        return Ebond;
    }
        

    void ParticleEnergies(){
        Ekin = p1->Ekin + p2->Ekin;
        Epot = p1->Epot + p2->Epot;
    }
};


struct ParticleSimulation{

    float Ekin_total {};
    float Epot_total {}; 
    float Ebond_total {};
    float V_rc {};
    float V_grad_rc {};
    std::array<Particle, N> particles;
    std::array<Molecule, N/2> molecules;
    std::array<std::array<RelativePosition, N>, N> relative_positions;


    ParticleSimulation() {
        for (int i = 0; i < N; ++i){
            particles[i] = Particle(i);
        }
        for (int i = 0; i < N/2; ++i){
            molecules[i] = Molecule(i, &particles[2*i], &particles[2*i+1]);
        }
        for (int i = 0; i < N; ++i) {
            for (int j = i+1; j < N; ++j) {
                relative_positions[i][j] = RelativePosition(); // Only init upper triangle
            }
        }
        // Pre-calculate cut-off radius potential and gradient
        V_rc = 4*sigma*(pow(r_c2, -3) - 1) * pow(r_c2, -3); // Potential at cut-off radius
        V_grad_rc = 24*sigma*pow(r_c2, -4) * (1 - 2*pow(r_c2, -3)); // Gradient at cut-off radius
    }


    void AssignMoleculePositions(const std::array<std::tuple<int, int, int>, N/2>& positions, 
        const std::array<std::tuple<float, float, float>, N/2>& orientations){
        for (int i = 0; i < N/2; ++i) {
            molecules[i].InitialPosition(positions[i], orientations[i]);
            std::tuple<float, float, float> zero_velocity = std::make_tuple(0, 0, 0);
            molecules[i].InitialVelocity(zero_velocity); // In case we want to init particles in rest
        }
    }


    void AssignMoleculeVelocities(const std::array<std::tuple<float, float, float>, N/2>& velocities){
        for (int i = 0; i < N/2; ++i) {
            molecules[i].InitialVelocity(velocities[i]);
        }
    }


    void RunSimulation(std::array<std::tuple<float, float, float, float>, T>& logged_total_energies,
        std::array<std::array<std::tuple<float, float, float, float>, N/2>, T>& logged_molecule_energies,
        std::array<std::array<std::tuple<float, float, float>, N>, T>& logged_particle_positions){
        for (int t = 0; t < T; t++){
            IterTimestep(t, logged_total_energies, logged_molecule_energies, logged_particle_positions);
            // std::cout << "E_kin: " << Ekin_total << '\n';
            // std::cout << "E_pot: " << Epot_total << '\n';
            // std::cout << "E_tot: " << Ekin_total + Epot_total << '\n';
            if (t == 0){
                std::cout << "Initial energies\n";
                std::cout << "E_kin: " << Ekin_total << '\n';
                std::cout << "E_pot: " << Epot_total << '\n';
                std::cout << "E_bond: " << Ebond_total << '\n';
                std::cout << "E_tot: " << Ekin_total + Epot_total + Ebond_total << '\n';
            }
            if (t == T-1){
                std::cout << "Final energies\n";
                std::cout << "E_kin: " << Ekin_total << '\n';
                std::cout << "E_pot: " << Epot_total << '\n';
                std::cout << "E_bond: " << Ebond_total << '\n';
                std::cout << "E_tot: " << Ekin_total + Epot_total + Ebond_total << '\n';
            }
        }
    }


    void IterTimestep(int t, std::array<std::tuple<float, float, float, float>, T>& logged_total_energies,
        std::array<std::array<std::tuple<float, float, float, float>, N/2>, T>& logged_molecule_energies,
        std::array<std::array<std::tuple<float, float, float>, N>, T>& logged_particle_positions){
        // Reset energies
        Ekin_total = 0;
        Epot_total = 0;
        Ebond_total = 0;

        // Step 1: Calculate relative positions
        for (int i = 0; i < N; ++i) {
            for (int j = i+1; j < N; ++j) {
                relative_positions[i][j].CalculateRelativePosition(particles[i], particles[j]);
            }
        }

        // Step 2: Calculate intra-molecular forces
        for (int i = 0; i < N; ++i){
            for (int j = i+1; j < N; ++j) {
                double Eint = LJPotentialGrad(particles[i], particles[j], relative_positions[i][j]);
                particles[i].Epot += Eint/2;
                particles[j].Epot += Eint/2;
                Epot_total += Eint;
            }
        }

        // Step 3: Update particle positions based on intra-molecular forces
        for (auto &p : particles){

            // Calculate new positions using Verlet integration
            p.x_new = 2*p.x_curr - p.x_prev + dt*dt*p.fx/m;
            p.y_new = 2*p.y_curr - p.y_prev + dt*dt*p.fy/m;
            p.z_new = 2*p.z_curr - p.z_prev + dt*dt*p.fz/m;

            // Update positions with periodic boundary conditions
            p.x_new = fmod(L+p.x_new, L);
            p.y_new = fmod(L+p.y_new, L);
            p.z_new = fmod(L+p.z_new, L);

            // Reset forces
            p.fx = 0;
            p.fy = 0;
            p.fz = 0;
        }


        // Step 4: Recalculate relative positions
        for (int i = 0; i < N; ++i) {
            for (int j = i+1; j < N; ++j) {
                relative_positions[i][j].CalculateRelativePosition(particles[i], particles[j]);
            }
        }

        // Step 5: Calculate intra-molecular forces
        for (auto& m : molecules){
            double Ebond = m.MolecularBond(relative_positions[m.p1->particle_id][m.p2->particle_id]);
            Ebond_total += Ebond;
        }   

        // Step 6: Update particle positions based on inter-molecular forces
        for (auto &p : particles){

            // Calculate new positions using Verlet integration
            p.x_new += dt*dt*p.gx/m;
            p.y_new += dt*dt*p.gy/m;
            p.z_new += dt*dt*p.gz/m;

            // Update positions with periodic boundary conditions
            p.x_new = fmod(L+p.x_new, L);
            p.y_new = fmod(L+p.y_new, L);
            p.z_new = fmod(L+p.z_new, L);
        }

        // Step 7: Calculate velocities
        for (auto &p : particles){
            float delta_x = p.x_new - p.x_prev;
            float delta_y = p.y_new - p.y_prev;
            float delta_z = p.z_new - p.z_prev;

            delta_x = (delta_x > L / 2) ? delta_x - L : (delta_x < -L / 2) ? delta_x + L : delta_x;
            delta_y = (delta_y > L / 2) ? delta_y - L : (delta_y < -L / 2) ? delta_y + L : delta_y;
            delta_z = (delta_z > L / 2) ? delta_z - L : (delta_z < -L / 2) ? delta_z + L : delta_z;

            p.vx = 0.5*delta_x/dt;
            p.vy = 0.5*delta_y/dt;
            p.vz = 0.5*delta_z/dt;
            Ekin_total += p.EnergyKinetic();

            // Position reasignment
            p.x_prev = p.x_curr;
            p.y_prev = p.y_curr;
            p.z_prev = p.z_curr;
            p.x_curr = p.x_new;
            p.y_curr = p.y_new;
            p.z_curr = p.z_new; 
        }

        // Step 8: Log positions and energies
        logged_total_energies[t] = std::make_tuple(Ekin_total, Epot_total, Ebond_total, Ekin_total + Epot_total + Ebond_total);
        for (int i = 0; i < N/2; ++i){
            molecules[i].ParticleEnergies();
            logged_molecule_energies[t][i] = std::make_tuple(molecules[i].Ekin, molecules[i].Epot, molecules[i].Ebond, molecules[i].Ekin + molecules[i].Epot + molecules[i].Ebond);
            molecules[i].p1->Epot = 0; // Reset particle potential energies for next time-step
            molecules[i].p2->Epot = 0; 
        }
        for (int i = 0; i < N; ++i){
            logged_particle_positions[t][i] = std::make_tuple(particles[i].x_curr, particles[i].y_curr, particles[i].z_curr);
        }
    }


    double LJPotentialGrad(Particle& p1, Particle& p2, const RelativePosition& rp){

        if (p1.twin == &p2) { // Check if particles belong to the same molecule
            return 0; // No interaction
        }
    
        // Work directly with r2 for numerical stability, all exponents happen to be even anyway
        double r2 = rp.r2_curr();

        if (r2 > r_c2){ // Check cut-off radius criterion
            return 0; 
        }

        auto V_func = [](double r2) -> double { // Lambda function for potential
            return 4*sigma*(pow(r2, -3) - 1) * pow(r2, -3); 
        };
        auto V_grad_func = [](double r2) -> double { // Lambda function for gradient
            return 24*sigma*pow(r2, -4) * (1 - 2*pow(r2, -3)); 
        };

        // Potential and gradient shifted to be continious at r_c2
        double V = V_func(r2) - V_rc; 
        double V_grad = V_grad_func(r2) - V_grad_rc; 
    
        p1.fx += V_grad*rp.x_rel_curr;
        p1.fy += V_grad*rp.y_rel_curr;
        p1.fz += V_grad*rp.z_rel_curr;
        p2.fx -= V_grad*rp.x_rel_curr; // Newtons third law
        p2.fy -= V_grad*rp.y_rel_curr; // -||-
        p2.fz -= V_grad*rp.z_rel_curr; // -||-
        return V;
    }
};



template<typename pos> // Template instead of direct tuple for code readability
std::array<pos, N/2> InitMoleculePositions(unsigned int seed){
    // Init positions of N particles in a 3D box of size L with random, non-overlapping coordinates.

    int bounds = (int) L-1;
    std::mt19937 rng(seed); 
    std::uniform_int_distribution<int> dist_pos(1, bounds); 

    struct CoordHash {
        std::size_t operator()(const pos& coords) const {
            auto [x, y, z] = coords;
            return std::hash<int>{}(x) ^ (std::hash<int>{}(y) << 1) ^ (std::hash<int>{}(z) << 2);
        }
    };

    std::array<pos, N/2> initial_positions;
    std::unordered_set<pos, CoordHash> check_overlap;

    int count = 0;
    while (count < N/2){
        int x = dist_pos(rng);
        int y = dist_pos(rng);
        int z = dist_pos(rng);
        pos coord = {x, y, z};
        if (check_overlap.find(coord) == check_overlap.end()){
            initial_positions[count] = coord;
            check_overlap.insert(coord);
            count++;
        }
    }    
    return initial_positions;
} 


template<typename vec> // Template instead of direct tuple for code readability
std::array<vec, N/2> InitMoleculeOrientations(unsigned int seed){

    std::mt19937 rng(seed); 
    std::uniform_real_distribution<float> dist_uni_zero_one(0, 1); 
    std::array<vec, N/2> initial_orientations;

    for (int i = 0; i < N/2; ++i){
        float u = dist_uni_zero_one(rng);
        float v = dist_uni_zero_one(rng);

        // Sample spherical angles uniformly to ensure random directions
        float theta = 2*pi*u;
        float phi = acos(2*v - 1);

        float nx = sin(phi)*cos(theta);
        float ny = sin(phi)*sin(theta);
        float nz = cos(phi);

        vec orientation = {nx, ny, nz};
        initial_orientations[i] = orientation;
    }
    return initial_orientations;
}

template<typename vel>
std::array<vel, N/2> InitialVelocityBoltzmann(unsigned int seed, float particle_energy_mean){
    // Init velocities of N/2 molecules with kinetic energies distributed according to a Maxwell-Boltzmann distribution.

    std::mt19937 rng(seed); 
    std::uniform_real_distribution<float> dist_uni_zero_one(0, 1); 
    std::array<vel, N/2> initial_velocities;

    float vx_total = 0.0f, vy_total = 0.0f, vz_total = 0.0f;
    for (int i = 0; i < N/2; ++i){

        float u = dist_uni_zero_one(rng);
        float v = dist_uni_zero_one(rng);
        float w = dist_uni_zero_one(rng);

        // Sample spherical angles uniformly to ensure random directions
        float theta = 2*pi*u;
        float phi = acos(2*v - 1);

        // Sample Boltzmann distributed energy
        float particle_energy = -particle_energy_mean*log(1.0f - w);

        // Put results together and insert velocity to array
        float v_norm = sqrt(2*particle_energy/m);
        float vx = v_norm*sin(phi)*cos(theta);
        float vy = v_norm*sin(phi)*sin(theta);
        float vz = v_norm*cos(phi);

        vx_total += vx;
        vy_total += vy;
        vz_total += vz;

        vel velocity = {vx, vy, vz};
        initial_velocities[i] = velocity;
    }
    vx_total /= N/2;
    vy_total /= N/2;
    vz_total /= N/2;
    for (auto& v : initial_velocities) { // Normalize velocities to give molecule cloud zero net momentum
        std::get<0>(v) -= vx_total;
        std::get<1>(v) -= vy_total;
        std::get<2>(v) -= vz_total;
    }
    return initial_velocities;
}


void WritePositionsToCSV(const std::string& filename, const std::array<std::tuple<float, float, float, float>, T>& logged_total_energies,
    const std::array<std::array<std::tuple<float, float, float, float>, N/2>, T>& logged_molecule_energies,
    const std::array<std::array<std::tuple<float, float, float>, N>, T>& logged_particle_positions){
    // Function to write logged positions and energies to csv file.

    const std::string folder_name = "datafolder/";
    std::ofstream file(folder_name + filename);

    // Write header 
    file << "# Logged particle positions from simulation\n";
    file << "# N particles, T time-steps, dt single step time, L box side length 3 dimensions\n";
    file << "# Format:\n";
    file << "# T\n";
    file << "# dt\n";
    file << "# N\n";
    file << "# L\n";
    file << "# Ekin, Epot, Ebond, Etot (total energies)\n";
    file << "# Ekin1, Epot1, Ebond1, Etot1... EkinN, EpotN, EbondN, EtotN (molecule energies)\n";
    file << "# x1,y1,z1... xN, yN, zN (particle positions)\n";
    file << T << '\n';
    file << dt << '\n';
    file << N << '\n';
    file << L << '\n';

    // Write data
    for (int t = 0; t < T; ++t){ 
        auto total_energies = logged_total_energies[t];
        auto molecule_energies = logged_molecule_energies[t];
        auto positions = logged_particle_positions[t];
        file << std::get<0>(total_energies) << ',' << std::get<1>(total_energies) << ',' << 
        std::get<2>(total_energies) << ',' << std::get<3>(total_energies) << '\n'; // Total energies

        for (int i = 0; i < (N/2)-1; ++i) { // Particle energies
            file << std::get<0>(molecule_energies[i]);
            file << ',';
            file << std::get<1>(molecule_energies[i]);
            file << ',';
            file << std::get<2>(molecule_energies[i]);
            file << ',';
            file << std::get<3>(molecule_energies[i]);
            file << ',';
        }
        file << std::get<0>(molecule_energies[(N/2)-1]);
        file << ',';
        file << std::get<1>(molecule_energies[(N/2)-1]);
        file << ',';
        file << std::get<2>(molecule_energies[(N/2)-1]);
        file << ',';
        file << std::get<3>(molecule_energies[(N/2)-1]);
        file << '\n';

        for (int i = 0; i < N-1; ++i) { // Particle positions
            file << std::get<0>(positions[i]);
            file << ',';
            file << std::get<1>(positions[i]);
            file << ',';
            file << std::get<2>(positions[i]);
            file << ',';
        }
        file << std::get<0>(positions[N-1]);
        file << ',';
        file << std::get<1>(positions[N-1]);
        file << ',';
        file << std::get<2>(positions[N-1]);
        file << '\n';
    }
    file.close();
    std::cout << "Data written to " << filename << '\n';
}


int main(int argc, char *argv[]){

    std::string filename = "molecule_data.csv";
    if (argc >= 2) {
        int output_index = std::stoi(argv[1]);
        filename = "molecule_data_" + std::to_string(output_index) + ".csv";
    }
    
    unsigned int seed = 123456789;
    float particle_energy_mean = 10.0f; // Mean energy per particle

    std::array<std::tuple<float, float, float, float>, T> logged_total_energies;
    std::array<std::array<std::tuple<float, float, float, float>, N/2>, T> logged_molecule_energies;
    std::array<std::array<std::tuple<float, float, float>, N>, T> logged_particle_positions;
    
    auto initial_positions = InitMoleculePositions<std::tuple<int, int, int>>(seed);
    auto initial_orientations = InitMoleculeOrientations<std::tuple<float, float, float>>(seed);
    auto initial_velocities = InitialVelocityBoltzmann<std::tuple<float, float, float>>(seed, particle_energy_mean);

    ParticleSimulation sim = ParticleSimulation();
    sim.AssignMoleculePositions(initial_positions, initial_orientations);
    sim.AssignMoleculeVelocities(initial_velocities);
    sim.RunSimulation(logged_total_energies, logged_molecule_energies, logged_particle_positions);
    std::cout << "Simulation complete.\n";

    WritePositionsToCSV(filename, logged_total_energies, logged_molecule_energies, logged_particle_positions);

    return 0;
}




