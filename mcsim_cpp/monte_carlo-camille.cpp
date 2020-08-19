#include <iostream>
#include <fstream>
#include <vector>
#include <array>
#include <string>
#include <cmath>
#include <algorithm>

#include <random> // for random numbers
#include <chrono> // for generating the seed



typedef std::array<double, 3> AtomCoord;
typedef std::vector<AtomCoord> Coordinates;

// A Global! Probably shouldn't be used in real code
std::default_random_engine re;

/*! Generate a random double within a given range */
double random_double(double lower_bound, double upper_bound)
{
   std::uniform_real_distribution<double> dist(lower_bound, upper_bound);
   return dist(re);
}

/*! Generate a random integer within a given range
    The generated integer will be on the range [a,b)
*/
double random_integer(int lower_bound, int upper_bound)
{           
   //dist will return [a,b] but we want [a,b)
   std::uniform_int_distribution<int> dist(lower_bound, upper_bound-1);
   return dist(re);
} 

double calculate_LJ(double r_ij) {

    double inv_r_ij = 1.0 / r_ij;

    double r6_term = pow(inv_r_ij, 6.0);

    double r12_term = r6_term * r6_term;

    double pairwise_energy = 4.0 * (r12_term - r6_term);

    return pairwise_energy;
}

double calculate_distance(AtomCoord coord1, AtomCoord coord2, double box_length=-1.0) {

    double distance = 0;

    for (int i = 0; i < 3; i++) {
        double dim_dist = (coord1[i] - coord2[i]);

        if (box_length != -1.0) {
            dim_dist = dim_dist - box_length * round(dim_dist / box_length);
        }
        dim_dist = dim_dist * dim_dist;
        distance += dim_dist;
    }

    distance = sqrt(distance);

    return distance;
}

double calculate_total_energy(Coordinates coordinates, double box_length, double cutoff) {

    double total_energy = 0.0;

    int num_atoms = coordinates.size();

    for (int i = 0; i < num_atoms; i++) {
        for (int j = i + 1; j < num_atoms; j++) {
            double dist_ij = calculate_distance(coordinates[i], coordinates[j], box_length);

            if (dist_ij < cutoff) {
                double interaction_energy = calculate_LJ(dist_ij);
                total_energy += interaction_energy;
            }
        }
    }

    return total_energy;
}

double calculate_tail_correction(int num_particles, double box_length, double cutoff) {
    
    double const1 = (8.0 * M_PI * num_particles * num_particles) / (3 * pow(box_length, 3));
    double const2 = (1.0/3.0) * pow((1.0 / cutoff), 9) - pow((1.0 / cutoff) , 3);
    return const1 * const2;
}

double calculate_pair_energy (Coordinates coordinates, int i_particle, double box_length, double cutoff) {
    
    double e_total = 0;

    int num_atoms = coordinates.size();

    for (int j = 0; j < num_atoms; j++) {
        if (j != i_particle) {
            double dist_ij = calculate_distance(coordinates[i_particle], coordinates[j], box_length);

            if (dist_ij < cutoff) {
                double interaction_energy = calculate_LJ(dist_ij);
                e_total += interaction_energy;
            }
        }
    }
    
    return e_total;
}

bool accept_or_reject(double delta_e, double beta) {

    bool accept = false;
    if (delta_e <= 0) {
        accept = true;
    }
    else {
        double random_number = random_double(0, 1);
        double p_acc = exp(-beta * delta_e);
        if (random_number < p_acc) {
            accept = true;
        }
        else {
            accept = false;
        }
    }

    return accept; 
}

std::pair<Coordinates, double> generate_config(int num_particles, double box_length) {
    Coordinates coords;

    for (int i = 0; i < num_particles; i++) {
        double x_rand = random_double(0, box_length);
        double y_rand = random_double(0, box_length);
        double z_rand = random_double(0, box_length);

        AtomCoord coordinate = {x_rand, y_rand, z_rand};

        while (std::count(coords.begin(), coords.end(), coordinate) != 0) {
            x_rand = random_double(0, box_length);
            y_rand = random_double(0, box_length);
            z_rand = random_double(0, box_length);

            coordinate = {x_rand, y_rand, z_rand};
        }

        coords.push_back(coordinate);

    }

    return std::make_pair(coords, box_length);
}

std::pair<Coordinates, double> read_xyz(std::string file_path) {
    // Opens up a file stream for input
    std::ifstream infile(file_path);

    // Check that it was successfully opened
    if(!infile.is_open()) {   
        throw std::runtime_error("File path in read_xyz does not exist!");
    }
    
    double dummy; // Data that is thrown away (box length, atom indices)
    double box_length;
    int num_atoms;
    
    // Grab box_length from first number, throw the rest away
    infile >> box_length >> dummy >> dummy;
    
    // now the number of atoms
    infile >> num_atoms;
    
    // Uncomment to help troubleshoot
    //std::cout << "Box length: " << box_length << " natoms: " << num_atoms << std::endl;
    
    // Stores the atomic coordinates
    // Remember, this is a vector of arrays
    Coordinates coords;
    
    for(int i = 0; i < num_atoms; i++) {   
        AtomCoord coord;
        
        // Throws away the atom index
        infile >> dummy >> coord[0] >> coord[1] >> coord[2];
        
        // Add to the vector
        coords.push_back(coord);
    }

    // Makes an appropriate pair object
    return std::make_pair(coords, box_length);
}

Coordinates run_simulation(Coordinates coordinates, const double box_length, const double cutoff, 
    const double reduced_temperature, const int num_steps, const double max_displacement=0.1, const int freq=1000) {
    
    std::vector<int> steps;
    std::vector<double> energies;
    std::vector<Coordinates> all_coordinates;

    double beta = 1.0 / reduced_temperature;

    int num_particles = coordinates.size();

    double total_energy = calculate_total_energy(coordinates, box_length, cutoff);
    total_energy += calculate_tail_correction(num_particles, box_length, cutoff);

    for (int i = 0; i < num_steps; i++) {

        int random_particle = random_integer(0, num_particles);

        double current_energy = calculate_pair_energy(coordinates, random_particle, box_length, cutoff);

        double x_rand = random_double(-max_displacement, max_displacement);
        double y_rand = random_double(-max_displacement, max_displacement);
        double z_rand = random_double(-max_displacement, max_displacement);

        coordinates[random_particle][0] += x_rand;
        coordinates[random_particle][1] += y_rand;
        coordinates[random_particle][2] += z_rand;

        double proposed_energy = calculate_pair_energy(coordinates, random_particle, box_length, cutoff);

        double delta_energy = proposed_energy - current_energy;

        bool accept = accept_or_reject(delta_energy, beta);

        if (accept) {
            total_energy += delta_energy;
        }
        else {
            coordinates[random_particle][0] -= x_rand;
            coordinates[random_particle][1] -= y_rand;
            coordinates[random_particle][2] -= z_rand;
        }

        if (i % freq == 0) {
            std::cout << i << ' ' << total_energy / num_particles << std::endl;
            steps.push_back(i);
            energies.push_back(total_energy / num_particles);
            all_coordinates.push_back(coordinates);
        }

    }

    return coordinates;
}

int main(void) {

    re.seed(std::chrono::system_clock::now().time_since_epoch().count());

    std::cout << "calculate_LJ tests: " << std::endl;
    std::cout << calculate_LJ(1) << std::endl;

    std::cout << "accept_or_reject tests: " << std::endl;
    std::cout << accept_or_reject(1, 1) << std::endl;
    std::cout << accept_or_reject(-100, 1) << std::endl;

    std::cout << "calculate_distance tests: " << std::endl;

    AtomCoord coord1 = {0.0, 0.0, 0.0};
    AtomCoord coord2 = {0.0, 0.0, 1.0};
    std::cout << calculate_distance(coord1, coord2, 10.0) << std::endl;
    std::cout << calculate_distance({0.0, 0.0, 0.0}, {0.0, 0.0, 8.0}, 10.0) << std::endl;

    std::pair<Coordinates, double> xyz_info = read_xyz("../lj_sample_configurations/lj_sample_config_periodic1.txt");
    Coordinates coords = xyz_info.first;
    double box_length = xyz_info.second;

    std::cout << "Test calculate_total_energy: " << std::endl;
    std::cout << calculate_total_energy(coords, box_length, 3) << std::endl;

    std::cout << "Test calculate_tail_correction: " << std::endl;
    std::cout << calculate_tail_correction(coords.size(), box_length, 3) << std::endl;

    std::cout << "Test calculate_pair_energy: " << std::endl;
    Coordinates coords2;
    coords2.push_back({0.0, 0.0, 0.0});
    coords2.push_back({0.0, 0.0, pow(2.0, 1.0 / 6.0)});
    coords2.push_back({0.0, 0.0, 2.0 * pow(2.0, 1.0 /6.0)});
    std::cout << calculate_pair_energy(coords2, 1, 10.0, 3) << std::endl;

    std::cout << "Test generate_config: " << std::endl;
    std::pair<Coordinates, double> randsys = generate_config(50, 10);
    Coordinates randcoords = randsys.first;
    double rand_box_length = randsys.second;
    std::cout << randcoords.size() << std::endl;
    for (size_t i = 0; i < randcoords.size(); i++) {
        for (int j = 0; j < 3; j++) {
            std::cout << randcoords[i][j] << ' ';
        }
        std::cout << std::endl;
    }
    std::cout << rand_box_length << std::endl;

    std::cout << "Running simulation: " << std::endl;
    Coordinates newcoords = run_simulation(coords, box_length, 3, 1.4, 10000);

    std::cout << newcoords.size() << std::endl;

    return 0;
}