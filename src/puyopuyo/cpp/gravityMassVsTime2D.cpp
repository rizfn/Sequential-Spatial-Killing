#include <random>
#include <vector>
#include <unordered_map>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <filesystem>
#include <queue>
#include <set>

// #pragma GCC optimize("Ofast","inline","fast-math","unroll-loops","no-stack-protector")
#pragma GCC optimize("inline", "unroll-loops", "no-stack-protector")
#pragma GCC target("sse,sse2,sse3,ssse3,sse4,popcnt,abm,mmx,avx,avx2,tune=native", "f16c")

static auto _ = []()
{std::ios_base::sync_with_stdio(false);std::cin.tie(nullptr);std::cout.tie(nullptr);return 0; }();

std::random_device rd;
std::mt19937 gen(rd());

// Define constants
constexpr int DEFAULT_L = 128;       // side length of the square lattice
constexpr int DEFAULT_N_SPECIES = 6; // number of species
constexpr int DEFAULT_STEPS_PER_LATTICEPOINT = 128;

void placePuyo(std::vector<std::vector<int>> &lattice, std::vector<std::vector<bool>> &movedSites,
               std::uniform_int_distribution<> &dis_l, std::uniform_int_distribution<> &dis_species, int L, int H)
{
    // Select a column and species
    int col = dis_l(gen);
    int species = dis_species(gen);

    for (int row = H - 1; row >= 0; --row)
    {
        if (lattice[row][col] == 0)
        {
            lattice[row][col] = species;
            movedSites[row][col] = true;
            return;
        }
    }
}

void annihilatePuyo(std::vector<std::vector<int>> &lattice, std::vector<std::vector<bool>> &movedSites, int L, int H)
{
    std::vector<std::vector<bool>> visited(H, std::vector<bool>(L, false));

    // Directions for neighbors (up, down, left, right)
    std::vector<std::pair<int, int>> directions = {{-1, 0}, {1, 0}, {0, -1}, {0, 1}};

    for (int row = 0; row < H; ++row)
    {
        for (int col = 0; col < L; ++col)
        {
            if (movedSites[row][col] && !visited[row][col] && lattice[row][col] != 0)
            {
                // Perform flood-fill to find the cluster
                std::vector<std::pair<int, int>> cluster;
                std::queue<std::pair<int, int>> q;
                q.push({row, col});
                visited[row][col] = true;

                while (!q.empty())
                {
                    auto [r, c] = q.front();
                    q.pop();
                    cluster.push_back({r, c});

                    // Check neighbors
                    for (auto [dr, dc] : directions)
                    {
                        int nr = r + dr, nc = c + dc;
                        if (nr >= 0 && nr < H && nc >= 0 && nc < L && !visited[nr][nc] && lattice[nr][nc] == lattice[row][col])
                        {
                            q.push({nr, nc});
                            visited[nr][nc] = true;
                        }
                    }
                }

                // If the cluster size is greater than 1, remove it
                if (cluster.size() > 1)
                {
                    for (auto [r, c] : cluster)
                    {
                        lattice[r][c] = 0;
                    }
                }
            }
        }
    }
}

void fall(std::vector<std::vector<int>> &lattice, std::vector<std::vector<bool>> &movedSites, int L, int H)
{
    for (int col = 0; col < L; ++col)
    {
        int writeRow = H - 1; // Start from the bottom of the column
        for (int row = H - 1; row >= 0; --row)
        {
            if (lattice[row][col] != 0)
            {
                if (writeRow != row)
                {
                    lattice[writeRow][col] = lattice[row][col];
                    lattice[row][col] = 0;
                    movedSites[writeRow][col] = true;
                }
                writeRow--;
            }
        }
        // Clear any remaining cells above the last written row
        for (int row = writeRow; row >= 0; --row)
        {
            movedSites[row][col] = false;
        }
    }
}

void run(std::ofstream &file, int L, int N_SPECIES, int STEPS_PER_LATTICEPOINT)
{
    // Define distributions
    std::uniform_int_distribution<> dis_species(1, N_SPECIES);
    std::uniform_int_distribution<> dis_l(0, L - 1);

    int H = STEPS_PER_LATTICEPOINT; // Height of the lattice
    std::vector<std::vector<int>> lattice(H, std::vector<int>(L, 0));
    std::vector<std::vector<bool>> movedSites(H, std::vector<bool>(L, false));

    for (int step = 0; step <= STEPS_PER_LATTICEPOINT; ++step)
    {
        // Add L random puyos to random columns
        for (int i = 0; i < L; ++i)
        {
            placePuyo(lattice, movedSites, dis_l, dis_species, L, H);

            // Annihilation-fall cycle
            while (true)
            {
                annihilatePuyo(lattice, movedSites, L, H);
                std::vector<std::vector<bool>> newMovedSites(H, std::vector<bool>(L, false));
                fall(lattice, newMovedSites, L, H);

                if (newMovedSites == movedSites)
                    break;

                movedSites = newMovedSites;
            }
        }

        // Record both the number of filled cells and the max height of any column
        int filled_cells = 0;
        int max_height = 0;
        for (int col = 0; col < L; ++col)
        {
            for (int row = 0; row < H; ++row)
            {
                if (lattice[row][col] != 0)
                {
                    filled_cells++;
                    max_height = std::max(max_height, H - row);
                }
            }
        }

        // Record the number of filled cells and the max height
        file << step << "\t" << filled_cells << "\t" << max_height << "\n";

        // Print progress
        std::cout << "Progress: " << std::fixed << std::setprecision(2)
                  << static_cast<double>(step) / STEPS_PER_LATTICEPOINT * 100 << "%\r" << std::flush;
    }
}

int main(int argc, char *argv[])
{
    int L = DEFAULT_L;
    int N_SPECIES = DEFAULT_N_SPECIES;
    int STEPS_PER_LATTICEPOINT = DEFAULT_STEPS_PER_LATTICEPOINT;
    if (argc > 1)
        L = std::stoi(argv[1]);
    if (argc > 2)
        N_SPECIES = std::stoi(argv[2]);
    if (argc > 3)
        STEPS_PER_LATTICEPOINT = std::stoi(argv[3]);

    std::string exePath = argv[0];
    std::string exeDir = std::filesystem::path(exePath).parent_path().string();
    std::ostringstream filePathStream;
    filePathStream << exeDir << "\\outputs\\gravity2D\\L_" << L << "_N_" << N_SPECIES << "_steps_" << STEPS_PER_LATTICEPOINT << ".tsv";
    std::string filePath = filePathStream.str();

    std::ofstream file;
    file.open(filePath);
    file << "step\tmass\theight\n";

    run(file, L, N_SPECIES, STEPS_PER_LATTICEPOINT);

    file.close();

    return 0;
}