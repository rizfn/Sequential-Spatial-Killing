#include <algorithm>
#include <cmath>
#include <fstream>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <memory>
#include <queue>
#include <random>
#include <sstream>
#include <utility>
#include <vector>

// #pragma GCC optimize("Ofast","inline","fast-math","unroll-loops","no-stack-protector")
#pragma GCC optimize("inline", "unroll-loops", "no-stack-protector")
#pragma GCC target("sse,sse2,sse3,ssse3,sse4,popcnt,abm,mmx,avx,avx2,tune=native", "f16c")

static auto _ = []()
{std::ios_base::sync_with_stdio(false);std::cin.tie(nullptr);std::cout.tie(nullptr);return 0; }();

std::random_device rd;
std::mt19937 gen(rd());

constexpr int DEFAULT_L = 128;
constexpr int DEFAULT_N_COLORS = 16;
constexpr int DEFAULT_N_STEPS = 4096;
constexpr double DEFAULT_DENSITY = 0.25;
constexpr int DEFAULT_SIM_NO = 0;

inline int index2D(int row, int col, int L)
{
    return row * L + col;
}

std::vector<std::vector<bool>> createInteractionMatrix(double density, int nColors)
{
    density = std::clamp(density, 0.0, 1.0);

    std::vector<std::pair<int, int>> upperEntries;
    upperEntries.reserve(nColors * (nColors + 1) / 2);
    for (int i = 1; i <= nColors; ++i)
    {
        for (int j = i; j <= nColors; ++j)
            upperEntries.emplace_back(i, j);
    }

    std::shuffle(upperEntries.begin(), upperEntries.end(), gen);

    int onesTarget = static_cast<int>(std::llround(density * static_cast<double>(upperEntries.size())));
    onesTarget = std::clamp(onesTarget, 0, static_cast<int>(upperEntries.size()));

    std::vector<std::vector<bool>> J(nColors + 1, std::vector<bool>(nColors + 1, false));
    for (int idx = 0; idx < onesTarget; ++idx)
    {
        auto [i, j] = upperEntries[idx];
        J[i][j] = true;
        J[j][i] = true;
    }

    // for (int i = 1; i <= nColors; ++i)
    // {
    //     bool rowHasOne = false;
    //     for (int j = 1; j <= nColors; ++j)
    //     {
    //         if (J[i][j])
    //         {
    //             rowHasOne = true;
    //             break;
    //         }
    //     }
    //     if (!rowHasOne)
    //         J[i][i] = true;
    // }

    return J;
}

std::vector<std::pair<int, int>> getNeighbors(int row, int col, int L, int H)
{
    std::vector<std::pair<int, int>> neighbors;

    if (row > 0)
        neighbors.push_back({row - 1, col});
    if (row + 1 < H)
        neighbors.push_back({row + 1, col});

    neighbors.push_back({row, (col - 1 + L) % L});
    neighbors.push_back({row, (col + 1) % L});

    return neighbors;
}

int annihilatePuyo(std::vector<int> &lattice, bool *movedSites,
                   const std::vector<std::vector<bool>> &J, int L, int H)
{
    std::vector<bool> seedProcessed(H * L, false);
    std::vector<bool> toEliminate(H * L, false);
    int totalEliminated = 0;

    for (int row = 0; row < H; ++row)
    {
        for (int col = 0; col < L; ++col)
        {
            int start = index2D(row, col, L);
            if (!movedSites[start] || seedProcessed[start] || lattice[start] == 0)
                continue;

            std::vector<bool> visited(H * L, false);
            std::vector<int> component;
            std::queue<int> q;
            q.push(start);
            visited[start] = true;
            seedProcessed[start] = true;

            bool hasReaction = false;

            while (!q.empty())
            {
                int current = q.front();
                q.pop();
                component.push_back(current);

                int r = current / L;
                int c = current % L;
                int speciesCurrent = lattice[current];

                for (auto [nr, nc] : getNeighbors(r, c, L, H))
                {
                    int next = index2D(nr, nc, L);
                    int speciesNext = lattice[next];
                    if (speciesNext == 0)
                        continue;

                    if (!J[speciesCurrent][speciesNext])
                        continue;

                    hasReaction = true;
                    if (!visited[next])
                    {
                        visited[next] = true;
                        q.push(next);
                        if (movedSites[next])
                            seedProcessed[next] = true;
                    }
                }
            }

            if (hasReaction)
            {
                for (int idx : component)
                {
                    if (!toEliminate[idx])
                    {
                        toEliminate[idx] = true;
                        ++totalEliminated;
                    }
                }
            }
        }
    }

    for (int idx = 0; idx < H * L; ++idx)
    {
        if (toEliminate[idx])
        {
            lattice[idx] = 0;
            movedSites[idx] = false;
        }
    }

    return totalEliminated;
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

int runSimulation(int L, int nColors, int nSteps, double density)
{
    std::uniform_int_distribution<> dis_l(0, L - 1);
    std::uniform_int_distribution<> dis_species(1, nColors);
    std::vector<std::vector<bool>> J = createInteractionMatrix(density, nColors);

    int H = nSteps;
    int total_drops = (nSteps + 1) * L;
    std::unique_ptr<int[]> random_columns(new int[total_drops]);
    std::unique_ptr<int[]> random_species(new int[total_drops]);
    for (int drop = 0; drop < total_drops; ++drop)
    {
        random_columns[drop] = dis_l(gen);
        random_species[drop] = dis_species(gen);
    }

    std::vector<int> lattice(H * L, 0);
    std::unique_ptr<bool[]> movedSites(new bool[H * L]());

    int drop = 0;
    for (int step = 0; step <= nSteps; ++step)
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
                annihilatePuyo(lattice, movedSites.get(), J, L, H);
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
    }

    int final_mass = 0;
    for (int col = 0; col < L; ++col)
    {
        for (int row = 0; row < H; ++row)
        {
            if (lattice[index2D(row, col, L)] != 0)
                final_mass++;
        }
    }

    return final_mass;
}

int main(int argc, char *argv[])
{
    int L = DEFAULT_L;
    int N_COLORS = DEFAULT_N_COLORS;
    int N_STEPS = DEFAULT_N_STEPS;
    double density = DEFAULT_DENSITY;
    int SIM_NO = DEFAULT_SIM_NO;
    
    if (argc > 1)
        L = std::stoi(argv[1]);
    if (argc > 2)
        N_COLORS = std::stoi(argv[2]);
    if (argc > 3)
        N_STEPS = std::stoi(argv[3]);
    if (argc > 4)
        density = std::stod(argv[4]);
    if (argc > 5)
        SIM_NO = std::stoi(argv[5]);

    std::filesystem::path exeDir = std::filesystem::path(argv[0]).parent_path();
    std::filesystem::path filePath = exeDir / "outputs" / "finalMass" /
                                      ("L_" + std::to_string(L) +
                                       "_N_" + std::to_string(N_COLORS) + "/" +
                                       "rho_" + std::to_string(density) +
                                       "_steps_" + std::to_string(N_STEPS) +
                                       "_sim_" + std::to_string(SIM_NO) + ".txt");

    std::filesystem::create_directories(filePath.parent_path());

    std::ofstream file;
    file.open(filePath);
    if (!file.is_open())
    {
        std::cerr << "Failed to open output file: " << filePath.string() << "\n";
        return 1;
    }

    int finalMass = runSimulation(L, N_COLORS, N_STEPS, density);
    file << finalMass << "\n";
    file.close();

    return 0;
}