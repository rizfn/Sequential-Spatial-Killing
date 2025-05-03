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
constexpr int DEFAULT_L = 24;       // side length of the square lattice
constexpr int DEFAULT_N_SPECIES = 6; // number of species
constexpr int DEFAULT_STEPS_PER_LATTICEPOINT = 128;

void placePuyo(std::vector<std::vector<std::vector<int>>> &lattice, std::vector<std::vector<std::vector<bool>>> &movedSites,
               std::uniform_int_distribution<> &dis_x, std::uniform_int_distribution<> &dis_y,
               std::uniform_int_distribution<> &dis_species, int L, int H)
{
    // Select a column (x, y) and species
    int x = dis_x(gen);
    int y = dis_y(gen);
    int species = dis_species(gen);

    for (int z = H - 1; z >= 0; --z)
    {
        if (lattice[z][x][y] == 0)
        {
            lattice[z][x][y] = species;
            movedSites[z][x][y] = true;
            return;
        }
    }
}

void annihilatePuyo(std::vector<std::vector<std::vector<int>>> &lattice, std::vector<std::vector<std::vector<bool>>> &movedSites, int L, int H)
{
    std::vector<std::vector<std::vector<bool>>> visited(H, std::vector<std::vector<bool>>(L, std::vector<bool>(L, false)));

    // Directions for neighbors (up, down, left, right, above, below)
    std::vector<std::tuple<int, int, int>> directions = {
        {-1, 0, 0}, {1, 0, 0}, // vertical (z)
        {0, -1, 0}, {0, 1, 0}, // x
        {0, 0, -1}, {0, 0, 1}  // y
    };

    for (int z = 0; z < H; ++z)
    {
        for (int x = 0; x < L; ++x)
        {
            for (int y = 0; y < L; ++y)
            {
                if (movedSites[z][x][y] && !visited[z][x][y] && lattice[z][x][y] != 0)
                {
                    // Perform flood-fill to find the cluster
                    std::vector<std::tuple<int, int, int>> cluster;
                    std::queue<std::tuple<int, int, int>> q;
                    q.push({z, x, y});
                    visited[z][x][y] = true;

                    while (!q.empty())
                    {
                        auto [cz, cx, cy] = q.front();
                        q.pop();
                        cluster.push_back({cz, cx, cy});

                        // Check neighbors with periodic boundary in x and y (but not z)
                        for (auto [dz, dx, dy] : directions)
                        {
                            int nz = cz + dz;
                            int nx = (cx + dx + L) % L; // periodic in x
                            int ny = (cy + dy + L) % L; // periodic in y
                            // Only wrap x and y, not z
                            if (nz >= 0 && nz < H && !visited[nz][nx][ny] && lattice[nz][nx][ny] == lattice[z][x][y])
                            {
                                q.push({nz, nx, ny});
                                visited[nz][nx][ny] = true;
                            }
                        }
                    }

                    // If the cluster size is greater than 1, remove it
                    if (cluster.size() > 1)
                    {
                        for (auto [cz, cx, cy] : cluster)
                        {
                            lattice[cz][cx][cy] = 0;
                        }
                    }
                }
            }
        }
    }
}

void fall(std::vector<std::vector<std::vector<int>>> &lattice, std::vector<std::vector<std::vector<bool>>> &movedSites, int L, int H)
{
    for (int x = 0; x < L; ++x)
    {
        for (int y = 0; y < L; ++y)
        {
            int writeZ = H - 1; // Start from the bottom of the column
            for (int z = H - 1; z >= 0; --z)
            {
                if (lattice[z][x][y] != 0)
                {
                    if (writeZ != z)
                    {
                        lattice[writeZ][x][y] = lattice[z][x][y];
                        lattice[z][x][y] = 0;
                        movedSites[writeZ][x][y] = true;
                    }
                    writeZ--;
                }
            }
            // Clear any remaining cells above the last written row
            for (int z = writeZ; z >= 0; --z)
            {
                movedSites[z][x][y] = false;
            }
        }
    }
}

void run(std::ofstream &file, int L, int N_SPECIES, int STEPS_PER_LATTICEPOINT)
{
    // Define distributions
    std::uniform_int_distribution<> dis_species(1, N_SPECIES);
    std::uniform_int_distribution<> dis_x(0, L - 1);
    std::uniform_int_distribution<> dis_y(0, L - 1);

    int H = STEPS_PER_LATTICEPOINT; // Height of the lattice
    std::vector<std::vector<std::vector<int>>> lattice(H, std::vector<std::vector<int>>(L, std::vector<int>(L, 0)));
    std::vector<std::vector<std::vector<bool>>> movedSites(H, std::vector<std::vector<bool>>(L, std::vector<bool>(L, false)));

    for (int step = 0; step <= STEPS_PER_LATTICEPOINT; ++step)
    {
        // Add L * L random puyos to random columns
        for (int i = 0; i < L * L; ++i)
        {
            placePuyo(lattice, movedSites, dis_x, dis_y, dis_species, L, H);

            // Annihilation-fall cycle
            while (true)
            {
                annihilatePuyo(lattice, movedSites, L, H);
                std::vector<std::vector<std::vector<bool>>> newMovedSites(H, std::vector<std::vector<bool>>(L, std::vector<bool>(L, false)));
                fall(lattice, newMovedSites, L, H);

                if (newMovedSites == movedSites)
                    break;

                movedSites = newMovedSites;
            }
        }

        // Record both the number of filled cells and the max height of any column
        int filled_cells = 0;
        int max_height = 0;
        for (int x = 0; x < L; ++x)
        {
            for (int y = 0; y < L; ++y)
            {
                for (int z = 0; z < H; ++z)
                {
                    if (lattice[z][x][y] != 0)
                    {
                        filled_cells++;
                        max_height = std::max(max_height, H - z);
                    }
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
    filePathStream << exeDir << "\\outputs\\gravity3D\\L_" << L << "_N_" << N_SPECIES << "_steps_" << STEPS_PER_LATTICEPOINT << ".tsv";
    std::string filePath = filePathStream.str();

    std::ofstream file;
    file.open(filePath);
    file << "step\tmass\theight\n";

    run(file, L, N_SPECIES, STEPS_PER_LATTICEPOINT);

    file.close();

    return 0;
}