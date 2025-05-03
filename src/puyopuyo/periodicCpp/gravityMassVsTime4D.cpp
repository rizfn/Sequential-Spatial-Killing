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

#pragma GCC optimize("inline", "unroll-loops", "no-stack-protector")
#pragma GCC target("sse,sse2,sse3,ssse3,sse4,popcnt,abm,mmx,avx,avx2,tune=native", "f16c")

static auto _ = []()
{std::ios_base::sync_with_stdio(false);std::cin.tie(nullptr);std::cout.tie(nullptr);return 0; }();

std::random_device rd;
std::mt19937 gen(rd());

// Define constants
constexpr int DEFAULT_L = 8;         // side length of the hypercube lattice
constexpr int DEFAULT_N_SPECIES = 6; // number of species
constexpr int DEFAULT_STEPS_PER_LATTICEPOINT = 128;

void placePuyo(std::vector<std::vector<std::vector<std::vector<int>>>> &lattice,
               std::vector<std::vector<std::vector<std::vector<bool>>>> &movedSites,
               std::uniform_int_distribution<> &dis_x, std::uniform_int_distribution<> &dis_y,
               std::uniform_int_distribution<> &dis_w, std::uniform_int_distribution<> &dis_species, int L, int H)
{
    // Select a column (x, y, w) and species
    int x = dis_x(gen);
    int y = dis_y(gen);
    int w = dis_w(gen);
    int species = dis_species(gen);

    for (int z = H - 1; z >= 0; --z)
    {
        if (lattice[z][x][y][w] == 0)
        {
            lattice[z][x][y][w] = species;
            movedSites[z][x][y][w] = true;
            return;
        }
    }
}

void annihilatePuyo(std::vector<std::vector<std::vector<std::vector<int>>>> &lattice,
                    std::vector<std::vector<std::vector<std::vector<bool>>>> &movedSites, int L, int H)
{
    std::vector<std::vector<std::vector<std::vector<bool>>>> visited(
        H, std::vector<std::vector<std::vector<bool>>>(L, std::vector<std::vector<bool>>(L, std::vector<bool>(L, false))));

    // Directions for neighbors (up, down, left, right, above, below, w+1, w-1)
    std::vector<std::tuple<int, int, int, int>> directions = {
        {-1, 0, 0, 0}, {1, 0, 0, 0}, // z (vertical, no wrap)
        {0, -1, 0, 0}, {0, 1, 0, 0}, // x
        {0, 0, -1, 0}, {0, 0, 1, 0}, // y
        {0, 0, 0, -1}, {0, 0, 0, 1} // w
    };

    for (int z = 0; z < H; ++z)
    {
        for (int x = 0; x < L; ++x)
        {
            for (int y = 0; y < L; ++y)
            {
                for (int w = 0; w < L; ++w)
                {
                    if (movedSites[z][x][y][w] && !visited[z][x][y][w] && lattice[z][x][y][w] != 0)
                    {
                        // Perform flood-fill to find the cluster
                        std::vector<std::tuple<int, int, int, int>> cluster;
                        std::queue<std::tuple<int, int, int, int>> q;
                        q.push({z, x, y, w});
                        visited[z][x][y][w] = true;

                        while (!q.empty())
                        {
                            auto [cz, cx, cy, cw] = q.front();
                            q.pop();
                            cluster.push_back({cz, cx, cy, cw});

                            // Check neighbors with periodic boundary in x, y, w (but not z)
                            for (auto [dz, dx, dy, dw] : directions)
                            {
                                int nz = cz + dz;
                                int nx = (cx + dx + L) % L; // periodic in x
                                int ny = (cy + dy + L) % L; // periodic in y
                                int nw = (cw + dw + L) % L; // periodic in w
                                // Only wrap x, y, w; not z
                                if (nz >= 0 && nz < H && !visited[nz][nx][ny][nw] && lattice[nz][nx][ny][nw] == lattice[z][x][y][w])
                                {
                                    q.push({nz, nx, ny, nw});
                                    visited[nz][nx][ny][nw] = true;
                                }
                            }
                        }

                        // If the cluster size is greater than 1, remove it
                        if (cluster.size() > 1)
                        {
                            for (auto [cz, cx, cy, cw] : cluster)
                            {
                                lattice[cz][cx][cy][cw] = 0;
                            }
                        }
                    }
                }
            }
        }
    }
}

void fall(std::vector<std::vector<std::vector<std::vector<int>>>> &lattice,
          std::vector<std::vector<std::vector<std::vector<bool>>>> &movedSites, int L, int H)
{
    for (int x = 0; x < L; ++x)
    {
        for (int y = 0; y < L; ++y)
        {
            for (int w = 0; w < L; ++w)
            {
                int writeZ = H - 1; // Start from the bottom of the column
                for (int z = H - 1; z >= 0; --z)
                {
                    if (lattice[z][x][y][w] != 0)
                    {
                        if (writeZ != z)
                        {
                            lattice[writeZ][x][y][w] = lattice[z][x][y][w];
                            lattice[z][x][y][w] = 0;
                            movedSites[writeZ][x][y][w] = true;
                        }
                        writeZ--;
                    }
                }
                // Clear any remaining cells above the last written row
                for (int z = writeZ; z >= 0; --z)
                {
                    movedSites[z][x][y][w] = false;
                }
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
    std::uniform_int_distribution<> dis_w(0, L - 1);

    int H = STEPS_PER_LATTICEPOINT; // Height of the lattice
    std::vector<std::vector<std::vector<std::vector<int>>>> lattice(
        H, std::vector<std::vector<std::vector<int>>>(L, std::vector<std::vector<int>>(L, std::vector<int>(L, 0))));
    std::vector<std::vector<std::vector<std::vector<bool>>>> movedSites(
        H, std::vector<std::vector<std::vector<bool>>>(L, std::vector<std::vector<bool>>(L, std::vector<bool>(L, false))));

    for (int step = 0; step <= STEPS_PER_LATTICEPOINT; ++step)
    {
        // Add L * L * L random puyos to random columns
        for (int i = 0; i < L * L * L; ++i)
        {
            placePuyo(lattice, movedSites, dis_x, dis_y, dis_w, dis_species, L, H);

            // Annihilation-fall cycle
            while (true)
            {
                annihilatePuyo(lattice, movedSites, L, H);
                std::vector<std::vector<std::vector<std::vector<bool>>>> newMovedSites(
                    H, std::vector<std::vector<std::vector<bool>>>(L, std::vector<std::vector<bool>>(L, std::vector<bool>(L, false))));
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
                for (int w = 0; w < L; ++w)
                {
                    for (int z = 0; z < H; ++z)
                    {
                        if (lattice[z][x][y][w] != 0)
                        {
                            filled_cells++;
                            max_height = std::max(max_height, H - z);
                        }
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
    filePathStream << exeDir << "\\outputs\\gravity4D\\L_" << L << "_N_" << N_SPECIES << "_steps_" << STEPS_PER_LATTICEPOINT << ".tsv";
    std::string filePath = filePathStream.str();

    std::ofstream file;
    file.open(filePath);
    file << "step\tmass\theight\n";

    run(file, L, N_SPECIES, STEPS_PER_LATTICEPOINT);

    file.close();

    return 0;
}