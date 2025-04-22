#include <random>
#include <vector>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <filesystem>

// Define constants
constexpr int DEFAULT_STEPS_PER_LATTICEPOINT = 1024; // Height of the 1D lattice (number of puyos)
constexpr int DEFAULT_N_SPECIES = 6; // Number of species

std::random_device rd;
std::mt19937 gen(rd());

void placePuyo(std::vector<int> &lattice, std::uniform_int_distribution<> &dis_species, int N_PUYOS)
{
    // Select a random species
    int species = dis_species(gen);

    // Place the puyo at the top of the lattice
    lattice.push_back(species);

    // Check the cell below for elimination
    if (lattice.size() > 1 && lattice[lattice.size() - 2] == species)
    {
        // Eliminate both the current and the below cell
        lattice.pop_back(); // Remove the current puyo
        lattice.pop_back(); // Remove the below puyo
    }
}

void run(std::ofstream &file, int N_PUYOS, int N_SPECIES)
{
    // Define distributions
    std::uniform_int_distribution<> dis_species(1, N_SPECIES);

    // Initialize the 1D lattice
    std::vector<int> lattice;

    for (int step = 0; step < N_PUYOS; ++step)
    {
        // Drop a puyo into the lattice
        placePuyo(lattice, dis_species, N_PUYOS);

        // Record the number of filled cells
        int filled_cells = lattice.size();

        // Write the results to the file
        file << step << "\t" << filled_cells << "\n";

        // Print progress
        std::cout << "Progress: " << std::fixed << std::setprecision(2)
                  << static_cast<double>(step) / N_PUYOS * 100 << "%\r" << std::flush;
    }
}

int main(int argc, char *argv[])
{
    int N_SPECIES = DEFAULT_N_SPECIES;
    int STEPS_PER_LATTICEPOINT = DEFAULT_STEPS_PER_LATTICEPOINT;
    if (argc > 1)
        N_SPECIES = std::stoi(argv[1]);
    if (argc > 2)
        STEPS_PER_LATTICEPOINT = std::stoi(argv[2]);

    std::string exePath = argv[0];
    std::string exeDir = std::filesystem::path(exePath).parent_path().string();
    std::ostringstream filePathStream;
    filePathStream << exeDir << "\\outputs\\gravity1D\\N_" << N_SPECIES << "_steps_" << STEPS_PER_LATTICEPOINT << ".tsv";
    std::string filePath = filePathStream.str();

    std::ofstream file;
    file.open(filePath);
    file << "step\tmass\n";

    run(file, STEPS_PER_LATTICEPOINT, N_SPECIES);

    file.close();

    return 0;
}