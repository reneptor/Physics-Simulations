#include <array>
#include <tuple>
#include <cmath>
#include <random>
#include <iostream>
#include <unordered_set>
#include <fstream>
#include <numbers>


constexpr float time_span = 10; // Time in seconds - dont change
constexpr float dt = 0.001f; // Time-step size
constexpr int T = static_cast<int>(time_span/dt); // Number of time-steps
constexpr float m = 1.0f;
constexpr int N = 30; // Number of particles
constexpr float L = 10.0f; // Box size
constexpr float r_c2 = 2.5f*2.5f; // Cut-off radius squared
constexpr float pi = 3.14159265358f; 


using std::sqrt, std::pow, std::log, std::cos, std::sin, std::acos, std::round, std::fmod;


struct Particle{

    int particle_id; // For iterations 
    double x_curr {}, y_curr {}, z_curr {}; // Current time-step position
    double x_prev {}, y_prev {}, z_prev {}; // Previous time-step position
    double vx {}, vy {}, vz {}; // Velocity
    double fx {}, fy {}, fz {}; // Force
    float Ekin, Epot {}; // Individual particle energies

    
    Particle() {} // For particle array initialization
    Particle(int id) : particle_id(id) {} // Actual initialization


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

    double x_rel {};
    double y_rel {};
    double z_rel {};
    
    void CalculateRelativePosition(const Particle& p1, const Particle& p2){ // Calculations with periodic b.c.
        x_rel = p2.x_curr - p1.x_curr - L*round((p2.x_curr - p1.x_curr)/L);
        y_rel = p2.y_curr - p1.y_curr - L*round((p2.y_curr - p1.y_curr)/L);
        z_rel = p2.z_curr - p1.z_curr - L*round((p2.z_curr - p1.z_curr)/L); 
    }
};


struct ParticleSimulation{

    float Ekin_total {};
    float Epot_total {}; 
    float V_rc {};
    float V_grad_rc {};
    std::array<Particle, N> particles;
    std::array<std::array<RelativePosition, N>, N> relative_positions;


    ParticleSimulation() {
        for (int i = 0; i < N; ++i){
            particles[i] = Particle(i);
        }
        for (int i = 0; i < N; ++i) {
            for (int j = i+1; j < N; ++j) {
                relative_positions[i][j] = RelativePosition(); // Only init upper triangle
            }
        }
        // Pre-calculate cut-off radius potential and gradient
        V_rc = 4*(pow(r_c2, -3) - 1) * pow(r_c2, -3); // Potential at cut-off radius
        V_grad_rc = 24*pow(r_c2, -4) * (1 - 2*pow(r_c2, -3)); // Gradient at cut-off radius
    }


    void AssignParticlePositions(const std::array<std::tuple<int, int, int>, N>& positions){
        for (std::size_t i = 0; i < positions.size(); ++i) {
            const auto& [x, y, z] = positions[i];
            particles[i].InitialPosition(static_cast<double>(x), static_cast<double>(y), static_cast<double>(z));
            particles[i].InitialVelocity(0, 0, 0); // In case we want to init particles in rest
        }
    }


    void AssignParticleVelocities(const std::array<std::tuple<float, float, float>, N>& velocities){
        for (std::size_t i = 0; i < velocities.size(); ++i) {
            const auto& [vx, vy, vz] = velocities[i];
            particles[i].InitialVelocity(vx, vy, vz);
        }
    }


    void RunSimulation(std::array<std::tuple<float, float, float>, T>& logged_total_energies,
        std::array<std::array<std::tuple<float, float, float>, N>, T>& logged_particle_energies,
        std::array<std::array<std::tuple<float, float, float>, N>, T>& logged_positions){
        for (int t = 0; t < T; t++){
            IterTimestep(t, logged_total_energies, logged_particle_energies, logged_positions);
            // std::cout << "E_kin: " << Ekin_total << '\n';
            // std::cout << "E_pot: " << Epot_total << '\n';
            // std::cout << "E_tot: " << Ekin_total + Epot_total << '\n';
            if (t == 0){
                std::cout << "Initial energies\n";
                std::cout << "E_kin: " << Ekin_total << '\n';
                std::cout << "E_pot: " << Epot_total << '\n';
                std::cout << "E_tot: " << Ekin_total + Epot_total << '\n';
            }
            if (t == T-1){
                std::cout << "Final energies\n";
                std::cout << "E_kin: " << Ekin_total << '\n';
                std::cout << "E_pot: " << Epot_total << '\n';
                std::cout << "E_tot: " << Ekin_total + Epot_total << '\n';
            }
        }
    }


    void IterTimestep(int t, std::array<std::tuple<float, float, float>, T>& logged_total_energies,
        std::array<std::array<std::tuple<float, float, float>, N>, T>& logged_particle_energies,
        std::array<std::array<std::tuple<float, float, float>, N>, T>& logged_positions){
        // Reset energies
        Ekin_total = 0;
        Epot_total = 0;

        // Step 1: Calculate relative positions
        for (int i = 0; i < N; ++i) {
            for (int j = i+1; j < N; ++j) {
                relative_positions[i][j].CalculateRelativePosition(particles[i], particles[j]);
            }
        }

        // Step 2: Force calculations
        for (int i = 0; i < N; ++i){
            for (int j = i+1; j < N; ++j) {
                double Eint = LJPotentialGrad(particles[i], particles[j], relative_positions[i][j]);
                particles[i].Epot += Eint/2;
                particles[j].Epot += Eint/2;
                Epot_total += Eint;
            }
        }

        // Step 3: Update vectors
        for (auto &p : particles){

            // Calculate new positions using Verlet integration
            double x_new = 2*p.x_curr - p.x_prev + dt*dt*p.fx/m;
            double y_new = 2*p.y_curr - p.y_prev + dt*dt*p.fy/m;
            double z_new = 2*p.z_curr - p.z_prev + dt*dt*p.fz/m;

            // Update positions with periodic boundary conditions
            x_new = fmod(L+x_new, L);
            y_new = fmod(L+y_new, L);
            z_new = fmod(L+z_new, L);

            // Reset forces
            p.fx = 0;
            p.fy = 0;
            p.fz = 0;
            
            // Calculate velocities
            float delta_x = x_new - p.x_prev;
            float delta_y = y_new - p.y_prev;
            float delta_z = z_new - p.z_prev;

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
            p.x_curr = x_new;
            p.y_curr = y_new;
            p.z_curr = z_new; 
        }

        // Step 4: Log positions and energies
        logged_total_energies[t] = std::make_tuple(Ekin_total, Epot_total, Ekin_total + Epot_total);
        for (int i = 0; i < N; ++i){
            logged_particle_energies[t][i] = std::make_tuple(particles[i].Ekin, particles[i].Epot, particles[i].Ekin + particles[i].Epot);
            logged_positions[t][i] = std::make_tuple(particles[i].x_curr, particles[i].y_curr, particles[i].z_curr);
            particles[i].Epot = 0; // Reset particle energies for next time-step
            particles[i].Ekin = 0; // (Ekin redundant)
        }
    }


    double LJPotentialGrad(Particle& p1, Particle& p2, const RelativePosition& rp){
    
        // Work directly with r2 for numerical stability, all exponents happen to be even anyway
        double r2 = rp.x_rel*rp.x_rel + rp.y_rel*rp.y_rel + rp.z_rel*rp.z_rel; 

        auto V_func = [](double r2) -> double { // Lambda function for potential
            return 4*(pow(r2, -3) - 1) * pow(r2, -3); 
        };
        auto V_grad_func = [](double r2) -> double { // Lambda function for gradient
            return 24*pow(r2, -4) * (1 - 2*pow(r2, -3)); 
        };

        if (r2 > r_c2){ // Check cut-off radius criterion
            return 0; 
        }

        // Potential and gradient shifted to be continious at r_c2
        double V = V_func(r2) - V_rc; 
        double V_grad = V_grad_func(r2) - V_grad_rc; 
    
        p1.fx += V_grad*rp.x_rel;
        p1.fy += V_grad*rp.y_rel;
        p1.fz += V_grad*rp.z_rel;
        p2.fx -= V_grad*rp.x_rel; // Newtons third law
        p2.fy -= V_grad*rp.y_rel; // -||-
        p2.fz -= V_grad*rp.z_rel; // -||-
        return V;
    }
};



template<typename pos> // Template instead of direct tuple for code readability
std::array<pos, N> InitParticlePositions(unsigned int seed){
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

    std::array<pos, N> initial_positions;
    std::unordered_set<pos, CoordHash> check_overlap;

    int count = 0;
    while (count < N){
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


template<typename pos> 
std::array<pos, 2> InitTwoParticlePositions(){

    std::array<pos, 2> initial_positions;
    initial_positions[0] = {5, 5, 4};
    initial_positions[1] = {5, 5, 6};

    return initial_positions;
}


template<typename vel> 
std::array<vel, 2> InitTwoParticleVelocities(float eps_mean){

    std::array<vel, 2> initial_velocities;
    float v = std::sqrt(2*eps_mean/m);
    vel vel1 = {0, 0, v};
    vel vel2 = {0, 0, -v};
    initial_velocities[0] = vel1;
    initial_velocities[1] = vel2;

    return initial_velocities;
}


template<typename vel>
std::array<vel, N> InitialVelocityBoltzmann(unsigned int seed, float eps_mean){
    // Init velocities of N particles with kinetic energies distributed according to a Maxwell-Boltzmann distribution.

    std::mt19937 rng(seed); 
    std::uniform_real_distribution<float> dist_uni_zero_one(0, 1); 
    std::array<vel, N> initial_velocities;

    float vx_total = 0.0f, vy_total = 0.0f, vz_total = 0.0f;
    for (int i = 0; i < N; ++i){

        float u = dist_uni_zero_one(rng);
        float v = dist_uni_zero_one(rng);
        float w = dist_uni_zero_one(rng);

        // Sample spherical angles uniformly to ensure random directions
        float theta = 2*pi*u;
        float phi = acos(2*v - 1);

        // Sample Boltzmann distributed energy
        float eps = -eps_mean*log(1.0f - w);

        // Put results together and insert velocity to array
        float v_norm = sqrt(2*eps/m);
        float vx = v_norm*sin(phi)*cos(theta);
        float vy = v_norm*sin(phi)*sin(theta);
        float vz = v_norm*cos(phi);

        vx_total += vx;
        vy_total += vy;
        vz_total += vz;

        vel velocity = {vx, vy, vz};
        initial_velocities[i] = velocity;
    }
    vx_total /= N;
    vy_total /= N;
    vz_total /= N;
    for (auto& v : initial_velocities) { // Normalize velocities to give particle cloud zero net momentum
        std::get<0>(v) -= vx_total;
        std::get<1>(v) -= vy_total;
        std::get<2>(v) -= vz_total;
    }
    return initial_velocities;
}


void WritePositionsToCSV(const std::string& filename, const std::array<std::tuple<float, float, float>, T>& logged_total_energies,
    const std::array<std::array<std::tuple<float, float, float>, N>, T>& logged_particle_energies,
    const std::array<std::array<std::tuple<float, float, float>, N>, T>& logged_positions){
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
    file << "# Ekin, Epot, Etot (total energies)\n";
    file << "# Ekin1, Epot1, Etot1... EkinN, EpotN, EtotN (particle energies)\n";
    file << "# x1,y1,z1... xN, yN, zN (positions)\n";
    file << T << '\n';
    file << dt << '\n';
    file << N << '\n';
    file << L << '\n';

    // Write data
    for (int t = 0; t < T; ++t){ 
        auto total_energies = logged_total_energies[t];
        auto particle_energies = logged_particle_energies[t];
        auto positions = logged_positions[t];
        file << std::get<0>(total_energies) << ',' << std::get<1>(total_energies) << ',' << std::get<2>(total_energies) << '\n'; // Total energies

        for (size_t i = 0; i < N-1; ++i) { // Particle energies
            file << std::get<0>(particle_energies[i]);
            file << ',';
            file << std::get<1>(particle_energies[i]);
            file << ',';
            file << std::get<2>(particle_energies[i]);
            file << ',';
        }
        file << std::get<0>(particle_energies[N-1]);
        file << ',';
        file << std::get<1>(particle_energies[N-1]);
        file << ',';
        file << std::get<2>(particle_energies[N-1]);
        file << '\n';

        for (size_t i = 0; i < N-1; ++i) { // Particle positions
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

    std::string filename = "particle_data.csv";
    if (argc >= 2) {
        int output_index = std::stoi(argv[1]);
        filename = "particle_data_" + std::to_string(output_index) + ".csv";
    }
    
    unsigned int seed = 123456789;
    float eps_mean = 1.0f; // Mean energy per particle

    std::array<std::tuple<float, float, float>, T> logged_total_energies;
    std::array<std::array<std::tuple<float, float, float>, N>, T> logged_particle_energies;
    std::array<std::array<std::tuple<float, float, float>, N>, T> logged_positions;
    
    auto initial_positions = InitParticlePositions<std::tuple<int, int, int>>(seed);
    auto initial_velocities = InitialVelocityBoltzmann<std::tuple<float, float, float>>(seed, eps_mean);

    // auto initial_positions = InitTwoParticlePositions<std::tuple<int, int, int>>();
    // auto initial_velocities = InitTwoParticleVelocities<std::tuple<float, float, float>>(eps_mean);

    ParticleSimulation sim = ParticleSimulation();
    sim.AssignParticlePositions(initial_positions);
    sim.AssignParticleVelocities(initial_velocities);
    sim.RunSimulation(logged_total_energies, logged_particle_energies, logged_positions);
    std::cout << "Simulation complete.\n";

    WritePositionsToCSV(filename, logged_total_energies, logged_particle_energies, logged_positions);

    return 0;
}




