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
#include <cmath> // For std::floor
#include <memory>

// #pragma GCC optimize("Ofast","inline","fast-math","unroll-loops","no-stack-protector")
#pragma GCC optimize("inline", "unroll-loops", "no-stack-protector")
#pragma GCC target("sse,sse2,sse3,ssse3,sse4,popcnt,abm,mmx,avx,avx2,tune=native", "f16c")

static auto _ = []()
{std::ios_base::sync_with_stdio(false);std::cin.tie(nullptr);std::cout.tie(nullptr);return 0; }();

std::random_device rd;
std::mt19937 gen(rd());

// Define constants
constexpr int DEFAULT_L = 128;       // side length of the square lattice
constexpr double DEFAULT_N_SPECIES = 6.0; // number of species (now double)
constexpr int DEFAULT_STEPS_PER_LATTICEPOINT = 128;

// Helper to create a discrete distribution for floating-point N_SPECIES
std::discrete_distribution<> createSpeciesDistribution(double N_SPECIES) {
    int n_int = static_cast<int>(std::floor(N_SPECIES));
    double frac = N_SPECIES - n_int;
    std::vector<double> weights;
    for (int i = 0; i < n_int; ++i) {
        weights.push_back(1.0);
    }
    if (frac > 0) {
        weights.push_back(frac);
    }
    return std::discrete_distribution<>(weights.begin(), weights.end());
}

inline int index2D(int row, int col, int L)
{
    return row * L + col;
}

bool placePuyo(std::vector<int> &lattice, bool *movedSites,
               std::uniform_int_distribution<> &dis_l, std::discrete_distribution<> &species_dist, int L, int H)
{
    int col = dis_l(gen);
    if (col < 0 || col >= L)
        return false;

    int species = species_dist(gen) + 1;

    for (int row = H - 1; row >= 0; --row)
    {
        int pos = index2D(row, col, L);
        if (lattice[pos] == 0)
        {
            lattice[pos] = species;
            movedSites[pos] = true;
            return true;
        }
    }

    return false;
}

std::vector<std::pair<int, int>> getNeighbors(int row, int col, int L, int H)
{
    std::vector<std::pair<int, int>> neighbors;

    for (int dr : {-1, 1})
    {
        int nr = row + dr;
        if (nr >= 0 && nr < H)
        {
            neighbors.push_back({nr, col});
        }
    }

    std::vector<int> sideRows;
    if (col % 2 == 0)
        sideRows = {row, row + 1};
    else
        sideRows = {row - 1, row};

    for (int dc : {-1, 1})
    {
        int nc = (col + dc + L) % L;
        for (int nr : sideRows)
        {
            if (nr >= 0 && nr < H)
            {
                neighbors.push_back({nr, nc});
            }
        }
    }

    return neighbors;
}

void annihilatePuyo(std::vector<int> &lattice, bool *movedSites, int L, int H)
{
    std::vector<bool> visited(H * L, false);

    for (int row = 0; row < H; ++row)
    {
        for (int col = 0; col < L; ++col)
        {
            int start = index2D(row, col, L);
            if (movedSites[start] && !visited[start] && lattice[start] != 0)
            {
                std::vector<std::pair<int, int>> cluster;
                std::queue<std::pair<int, int>> q;
                q.push({row, col});
                visited[start] = true;
                int species = lattice[start];

                while (!q.empty())
                {
                    auto [r, c] = q.front();
                    q.pop();
                    cluster.push_back({r, c});

                    for (auto [nr, nc] : getNeighbors(r, c, L, H))
                    {
                        int next = index2D(nr, nc, L);
                        if (!visited[next] && lattice[next] == species)
                        {
                            q.push({nr, nc});
                            visited[next] = true;
                        }
                    }
                }

                if (cluster.size() > 1)
                {
                    for (auto [r, c] : cluster)
                    {
                        int pos = index2D(r, c, L);
                        lattice[pos] = 0;
                        movedSites[pos] = false;
                    }
                }
            }
        }
    }
}

void fall(std::vector<int> &lattice, bool *movedSites, int L, int H)
{
    std::unique_ptr<bool[]> newMovedSites(new bool[H * L]());

    for (int col = 0; col < L; ++col)
    {
        int writeRow = H - 1;
        for (int row = H - 1; row >= 0; --row)
        {
            int src = index2D(row, col, L);
            if (lattice[src] != 0)
            {
                if (writeRow != row)
                {
                    int dst = index2D(writeRow, col, L);
                    lattice[dst] = lattice[src];
                    lattice[src] = 0;
                    newMovedSites[dst] = true;
                }
                writeRow--;
            }
        }
    }

    for (int i = 0; i < H * L; ++i)
        movedSites[i] = newMovedSites[i];
}

void run(std::ofstream &file, int L, double N_SPECIES, int STEPS_PER_LATTICEPOINT)
{
    std::discrete_distribution<> species_dist = createSpeciesDistribution(N_SPECIES);
    std::uniform_int_distribution<> dis_l(0, L - 1);

    int H = STEPS_PER_LATTICEPOINT;
    int total_drops = (STEPS_PER_LATTICEPOINT + 1) * L;
    std::unique_ptr<int[]> random_columns(new int[total_drops]);
    std::unique_ptr<int[]> random_species(new int[total_drops]);
    for (int drop = 0; drop < total_drops; ++drop)
    {
        random_columns[drop] = dis_l(gen);
        random_species[drop] = species_dist(gen) + 1;
    }

    std::vector<int> lattice(H * L, 0);
    std::unique_ptr<bool[]> movedSites(new bool[H * L]());

    int drop = 0;
    for (int step = 0; step <= STEPS_PER_LATTICEPOINT; ++step)
    {
        for (int i = 0; i < L; ++i)
        {
            int col = random_columns[drop];
            int species = random_species[drop];
            ++drop;

            bool placed = false;
            for (int row = H - 1; row >= 0; --row)
            {
                int pos = index2D(row, col, L);
                if (lattice[pos] == 0)
                {
                    lattice[pos] = species;
                    movedSites[pos] = true;
                    placed = true;
                    break;
                }
            }
            if (!placed)
                continue;

            while (true)
            {
                annihilatePuyo(lattice, movedSites.get(), L, H);
                std::unique_ptr<bool[]> newMovedSites(new bool[H * L]());
                fall(lattice, newMovedSites.get(), L, H);

                bool same = true;
                for (int idx = 0; idx < H * L; ++idx)
                {
                    if (newMovedSites[idx] != movedSites[idx])
                    {
                        same = false;
                        break;
                    }
                }

                if (same)
                    break;

                for (int idx = 0; idx < H * L; ++idx)
                    movedSites[idx] = newMovedSites[idx];
            }
        }

        int filled_cells = 0;
        int max_height = 0;
        for (int col = 0; col < L; ++col)
        {
            for (int row = 0; row < H; ++row)
            {
                if (lattice[index2D(row, col, L)] != 0)
                {
                    filled_cells++;
                    max_height = std::max(max_height, H - row);
                }
            }
        }

        file << step << "\t" << filled_cells << "\t" << max_height << "\n";

        std::cout << "Progress: " << std::fixed << std::setprecision(2)
                  << static_cast<double>(step) / STEPS_PER_LATTICEPOINT * 100 << "%\r" << std::flush;
    }
}

int main(int argc, char *argv[])
{
    int L = DEFAULT_L;
    double N_SPECIES = DEFAULT_N_SPECIES;
    int STEPS_PER_LATTICEPOINT = DEFAULT_STEPS_PER_LATTICEPOINT;
    int SIM_NO = 0;
    if (argc > 1)
        L = std::stoi(argv[1]);
    if (argc > 2)
        N_SPECIES = std::stod(argv[2]);
    if (argc > 3)
        STEPS_PER_LATTICEPOINT = std::stoi(argv[3]);
    if (argc > 4)
        SIM_NO = std::stoi(argv[4]);

    std::filesystem::path exeDir = std::filesystem::path(argv[0]).parent_path();
    std::filesystem::path filePath = exeDir / "outputs" / "massVsTime2D" / ("L_" + std::to_string(L) + "_N_" + std::to_string(N_SPECIES) + "_steps_" + std::to_string(STEPS_PER_LATTICEPOINT) + "_sim_" + std::to_string(SIM_NO) + ".tsv");

    std::filesystem::create_directories(filePath.parent_path());

    std::ofstream file;
    file.open(filePath);
    if (!file.is_open())
    {
        std::cerr << "Failed to open output file: " << filePath.string() << "\n";
        return 1;
    }
    file << "step\tmass\theight\n";

    run(file, L, N_SPECIES, STEPS_PER_LATTICEPOINT);

    file.close();

    return 0;
}
