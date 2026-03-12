#include <array>
#include <tuple>
#include <cmath>
#include <random>
#include <iostream>
#include <unordered_set>
#include <fstream>
#include <numbers>


constexpr float time_span = 20; // Time in seconds 
constexpr float dt = 0.001f; // Time-step size
constexpr int T = static_cast<int>(time_span/dt); // Number of time-steps
constexpr int T_s = 4; // Selected time-steps
constexpr float m = 1.0f;
constexpr float sigma = 5.0f; // LJ potential strength
constexpr int N = 500; // Number of particles
constexpr float L = 20.0f; // Box size
constexpr float r_c2 = 2.5f*2.5f; // Cut-off radius squared
constexpr float pi = 3.14159265358f; 
constexpr float k_B = 1.0f; // Boltzmann constant


using std::sqrt, std::pow, std::log, std::cos, std::sin, std::acos, std::round, std::fmod;


struct Particle{

    int particle_id; // For iterations 
    double x_curr {}, y_curr {}, z_curr {}; // Current time-step position
    // double x_new {}, y_new {}, z_new {}; // Next time-step position
    double vx_curr {}, vy_curr {}, vz_curr {}; // Current time-step velocity
    double vx_half {}, vy_half {}, vz_half {}; // Half time-step velocity
    // double vx_new {}, vy_new {}, vz_new {}; // Next time-step velocity
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
        vx_curr = vx_0;
        vy_curr = vy_0;
        vz_curr = vz_0;
    }


    double EnergyKinetic(){ // Kinetic energy as used in thermostat
        Ekin = 0.5 * m * ((vx_curr*vx_curr + vy_curr*vy_curr + vz_curr*vz_curr) 
        + (vx_half*vx_half + vy_half*vy_half + vz_half*vz_half))/2;
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

    RelativePosition() {};
    RelativePosition(const Particle& p1, const Particle& p2){
        CalculateRelativePosition(p1, p2);
    }
};


struct NoseHooverThermostat{

    float eps {}; // Friction coeff
    float Q {}; // Coupling strength
    float Temp {}; // Temperature

    NoseHooverThermostat() {} // For reference initialization
    NoseHooverThermostat(float coupling, float temp) : Q(coupling), Temp(temp) {}

    void IterThermostat(float Ekin_total) {
        float d_eps_dt = (2*Ekin_total - (3*N + 1)*k_B*Temp)/Q;
        eps += d_eps_dt*dt; 
    }
};


struct ParticleSimulation{

    float Ekin_total {};
    float Epot_total {}; 
    float V_rc {};
    float V_grad_rc {};
    unsigned int log_count {};
    std::array<Particle, N> particles;
    std::array<std::array<RelativePosition, N>, N> relative_positions;
    NoseHooverThermostat thermostat;


    ParticleSimulation(float Q, float Temp, const std::array<std::tuple<int, int, int>, N>& positions,
    const std::array<std::tuple<float, float, float>, N>& velocities) {
        for (int i = 0; i < N; ++i){
            particles[i] = Particle(i);
        }

        AssignParticlePositions(positions);
        AssignParticleVelocities(velocities);

        for (int i = 0; i < N; ++i) { // Calculate initial forces
            for (int j = i+1; j < N; ++j) {
                relative_positions[i][j] = RelativePosition(particles[i], particles[j]); // Only init upper triangle
                double Eint = LJPotentialGrad(particles[i], particles[j], relative_positions[i][j]);
                particles[i].Epot += Eint/2;
                particles[j].Epot += Eint/2;
                Epot_total += Eint;
            }
        }

        // Pre-calculate cut-off radius potential and gradient
        V_rc = 4*sigma*(pow(r_c2, -3) - 1) * pow(r_c2, -3); // Potential at cut-off radius
        V_grad_rc = 24*sigma*pow(r_c2, -4) * (1 - 2*pow(r_c2, -3)); // Gradient at cut-off radius

        thermostat = NoseHooverThermostat(Q, Temp);
    }


    void AssignParticlePositions(const std::array<std::tuple<int, int, int>, N>& positions){
        for (int i = 0; i < N; ++i) {
            auto [x0, y0, z0] = positions[i]; 
            particles[i].InitialPosition(static_cast<double>(x0), static_cast<double>(y0), static_cast<double>(z0));
            particles[i].InitialVelocity(0, 0, 0); // In case we want to init particles in rest
        }
    }


    void AssignParticleVelocities(const std::array<std::tuple<float, float, float>, N>& velocities){
        for (int i = 0; i < N; ++i) {
            auto [vx0, vy0, vz0] = velocities[i];
            particles[i].InitialVelocity(vx0, vy0, vz0);
        }
    }


    float TempInstantaneous(){
        return 2.0f * Ekin_total / (3*(N-1) * k_B);
    }


    void RunSimulation(
        std::array<int, T_s>& selected_timesteps,
        std::array<float, T_s>& logged_temperatures,
        std::array<std::tuple<float, float, float>, T_s>& logged_total_energies,
        std::array<std::array<std::tuple<float, float, float>, N>, T_s>& logged_particle_energies,
        std::array<std::array<std::tuple<float, float, float>, N>, T_s>& logged_positions){

        IterTimestep(0, selected_timesteps, logged_temperatures, logged_total_energies, logged_particle_energies, logged_positions);
        std::cout << "Initial energies and temperatures\n";
        std::cout << "Temperature: " << TempInstantaneous() << '\n';
        std::cout << "E_kin: " << Ekin_total << '\n';
        std::cout << "E_pot: " << Epot_total << '\n';
        std::cout << "E_tot: " << Ekin_total + Epot_total << '\n';
        for (int t = 1; t < T; t++){
            IterTimestep(t, selected_timesteps, logged_temperatures, logged_total_energies, logged_particle_energies, logged_positions);
        }
        std::cout << "Final energies and temperatures\n";
        std::cout << "Temperature: " << TempInstantaneous() << '\n';
        std::cout << "E_kin: " << Ekin_total << '\n';
        std::cout << "E_pot: " << Epot_total << '\n';
        std::cout << "E_tot: " << Ekin_total + Epot_total << '\n';
    }


    void IterTimestep(int t, 
        std::array<int, T_s>& selected_timesteps,
        std::array<float, T_s>& logged_temperatures,
        std::array<std::tuple<float, float, float>, T_s>& logged_total_energies,
        std::array<std::array<std::tuple<float, float, float>, N>, T_s>& logged_particle_energies,
        std::array<std::array<std::tuple<float, float, float>, N>, T_s>& logged_positions){
        // Reset energies
        Ekin_total = 0;
        Epot_total = 0;
        float eps = thermostat.eps;

        // Step 1: Update positions
        for (auto &p : particles){

            // Calculate new positions and half time-step velocities using Verlet integration
            p.x_curr += dt*p.vx_curr + 0.5*dt*dt*(p.fx/m - eps*p.vx_curr);
            p.y_curr += dt*p.vy_curr + 0.5*dt*dt*(p.fy/m - eps*p.vy_curr);
            p.z_curr += dt*p.vz_curr + 0.5*dt*dt*(p.fz/m - eps*p.vz_curr);

            p.vx_half = p.vx_curr + 0.5*dt*(p.fx/m - eps*p.vx_curr);
            p.vy_half = p.vy_curr + 0.5*dt*(p.fy/m - eps*p.vy_curr);
            p.vz_half = p.vz_curr + 0.5*dt*(p.fz/m - eps*p.vz_curr);

            // Record total kinetic energy to be used in thermostat calculation
            Ekin_total += p.EnergyKinetic();

            // Update positions with periodic boundary conditions
            p.x_curr = fmod(L+p.x_curr, L);
            p.y_curr = fmod(L+p.y_curr, L);
            p.z_curr = fmod(L+p.z_curr, L);

            // Reset forces
            p.fx = 0;
            p.fy = 0;
            p.fz = 0;
        }

        // Step 2: Force calculations
        for (int i = 0; i < N; ++i){
            for (int j = i+1; j < N; ++j) {
                relative_positions[i][j].CalculateRelativePosition(particles[i], particles[j]);
                double Eint = LJPotentialGrad(particles[i], particles[j], relative_positions[i][j]);
                particles[i].Epot += Eint/2;
                particles[j].Epot += Eint/2;
                Epot_total += Eint;
            }
        }

        // Step 3: Update Nose-Hoover thermostat
        thermostat.IterThermostat(Ekin_total); 

        // Step 4: Update velocities
        eps = thermostat.eps;
        for (auto &p : particles){

            p.vx_curr = (p.vx_half + 0.5*dt*p.fx/m)/(1 + 0.5*dt*eps);
            p.vy_curr = (p.vy_half + 0.5*dt*p.fy/m)/(1 + 0.5*dt*eps);
            p.vz_curr = (p.vz_half + 0.5*dt*p.fz/m)/(1 + 0.5*dt*eps);
        }

        // Step 5: Log positions and energies
        if (t == selected_timesteps[log_count]){
            // std::cout << t << '\n';
            logged_temperatures[log_count] = TempInstantaneous();
            logged_total_energies[log_count] = std::make_tuple(Ekin_total, Epot_total, Ekin_total + Epot_total);
            for (int i = 0; i < N; ++i){
                logged_particle_energies[log_count][i] = std::make_tuple(particles[i].Ekin, particles[i].Epot, particles[i].Ekin + particles[i].Epot);
                logged_positions[log_count][i] = std::make_tuple(particles[i].x_curr, particles[i].y_curr, particles[i].z_curr);
            }
            log_count += 1;
        }

        for (int i = 0; i < N; ++i){
            particles[i].Epot = 0; // Reset particle potential energies for next time-step
        }
    }


    double LJPotentialGrad(Particle& p1, Particle& p2, const RelativePosition& rp){
    
        // Work directly with r2 for numerical stability, all exponents happen to be even anyway
        double r2 = rp.x_rel*rp.x_rel + rp.y_rel*rp.y_rel + rp.z_rel*rp.z_rel; 

        auto V_func = [](double r2) -> double { // Lambda function for potential
            return 4*sigma*(pow(r2, -3) - 1) * pow(r2, -3); 
        };
        auto V_grad_func = [](double r2) -> double { // Lambda function for gradient
            return 24*sigma*pow(r2, -4) * (1 - 2*pow(r2, -3)); 
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


template<typename vel>
std::array<vel, N> InitialVelocityBoltzmann(unsigned int seed, float particle_energy_mean){
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


template<typename vel>
std::array<vel, N> InitialVelocityUniform(unsigned int seed, float particle_energy_mean){
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

        // Put results together and insert velocity to array
        float v_norm = sqrt(2*particle_energy_mean/m);
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


void WritePositionsToCSV(const std::string& filename, 
    const std::array<int, T_s> selected_timesteps, 
    const std::array<float, T_s> logged_temperatures, 
    const std::array<std::tuple<float, float, float>, T_s>& logged_total_energies,
    const std::array<std::array<std::tuple<float, float, float>, N>, T_s>& logged_particle_energies,
    const std::array<std::array<std::tuple<float, float, float>, N>, T_s>& logged_positions){
    // Function to write logged positions and energies to csv file.

    const std::string folder_name = "datafolder/";
    std::ofstream file(folder_name + filename);

    // Write header 
    file << "# Logged particle data from simulation at selected time-steps\n";
    file << "# N particles, T_s number of selected time-steps, dt single step time, L box side length 3 dimensions\n";
    file << "# Format:\n";
    file << "# --------------------- General parameters\n";
    file << "# T\n";
    file << "# T_s\n";
    file << "# dt\n";
    file << "# N\n";
    file << "# L\n";
    file << "# --------------------- Time-step specific parameters\n";
    file << "# t (time-step)\n";
    file << "# T (temperature)\n";
    file << "# Ekin, Epot, Etot (total energies)\n";
    file << "# Ekin1, Epot1, Etot1... EkinN, EpotN, EtotN (particle energies)\n";
    file << "# x1,y1,z1... xN, yN, zN (positions)\n";
    file << T << '\n';
    file << T_s << '\n';
    file << dt << '\n';
    file << N << '\n';
    file << L << '\n';

    // Write data
    for (int t = 0; t < T_s; ++t){ 
        auto time_step = selected_timesteps[t];
        auto temperature = logged_temperatures[t];
        auto total_energies = logged_total_energies[t];
        auto particle_energies = logged_particle_energies[t];
        auto positions = logged_positions[t];

        file << time_step << '\n'; // Time-step
        file << temperature << '\n'; // Temperature
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

    float Q = 1.0f;

    std::string filename = "particle_data.csv";
    if (argc >= 2) {
        int output_index = std::stoi(argv[1]);
        Q = pow(10, std::stoi(argv[2]));
        filename = "particle_data_" + std::to_string(output_index) + ".csv";
    }
    
    unsigned int seed = 123456789;
    float Temp = 10.0f;
    float particle_energy_mean = 1.5*k_B*Temp; // Mean energy per particle (3 d.o.f.)

    std::array<int, T_s> selected_timesteps = {10, 500, 2500, 19999};
    std::array<float, T_s> logged_temperatures;
    std::array<std::tuple<float, float, float>, T_s> logged_total_energies;
    std::array<std::array<std::tuple<float, float, float>, N>, T_s> logged_particle_energies;
    std::array<std::array<std::tuple<float, float, float>, N>, T_s> logged_positions;
    
    auto initial_positions = InitParticlePositions<std::tuple<int, int, int>>(seed);
    auto initial_velocities = InitialVelocityUniform<std::tuple<float, float, float>>(seed, particle_energy_mean);
    // auto initial_velocities = InitialVelocityBoltzmann<std::tuple<float, float, float>>(seed, particle_energy_mean);

    ParticleSimulation sim = ParticleSimulation(Q, Temp, initial_positions, initial_velocities);
    sim.RunSimulation(selected_timesteps, logged_temperatures, logged_total_energies, logged_particle_energies, logged_positions);
    std::cout << "Simulation complete.\n";

    WritePositionsToCSV(filename, selected_timesteps, logged_temperatures, logged_total_energies, logged_particle_energies, logged_positions);

    return 0;
}




