#include <random>
#include <vector>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <filesystem>
#include <queue>
#include <memory>

// #pragma GCC optimize("Ofast","inline","fast-math","unroll-loops","no-stack-protector")
#pragma GCC optimize("inline", "unroll-loops", "no-stack-protector")
#pragma GCC target("sse,sse2,sse3,ssse3,sse4,popcnt,abm,mmx,avx,avx2,tune=native", "f16c")

static auto _ = []()
{std::ios_base::sync_with_stdio(false);std::cin.tie(nullptr);std::cout.tie(nullptr);return 0; }();

std::mt19937 gen;  // seeded deterministically from SIM_NO in main (reproducible)

// Define constants
constexpr int DEFAULT_L = 128;       // side length of the square lattice
constexpr int DEFAULT_N_SPECIES = 6;
constexpr int DEFAULT_STEPS_PER_LATTICEPOINT = 128;

inline int index2D(int row, int col, int L)
{
    return row * L + col;
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

int annihilatePuyo(std::vector<int> &lattice, bool *movedSites, int L, int H)
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

                    if (speciesNext != speciesCurrent)
                        continue;

                    if (!visited[next])
                    {
                        visited[next] = true;
                        q.push(next);
                        if (movedSites[next])
                            seedProcessed[next] = true;
                    }
                }
            }

            if (component.size() > 1)
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

void run(std::ofstream &file, int L, int N_SPECIES, int STEPS_PER_LATTICEPOINT, int H_BOX)
{
    std::uniform_int_distribution<> dis_l(0, L - 1);
    std::uniform_int_distribution<> dis_species(1, N_SPECIES);

    file << "step\tfirst_col_height\tslope_distribution\n";

    int H = H_BOX;
    long long ceiling_hits = 0;
    int total_drops = (STEPS_PER_LATTICEPOINT + 1) * L;
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
            {
                ++ceiling_hits;
                continue;
            }

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

        std::vector<int> heights(L, 0);
        for (int c = 0; c < L; ++c)
        {
            for (int r = 0; r < H; ++r)
            {
                if (lattice[index2D(r, c, L)] != 0)
                {
                    heights[c] = H - r;
                    break;
                }
            }
        }

        file << std::fixed << std::setprecision(6) << static_cast<double>(step) << "\t" << heights[0] << "\t";
        for (int c = 0; c < L; ++c)
        {
            int next = (c + 1) % L;
            int slope = heights[next] - heights[c];
            file << slope;
            if (c < L - 1)
                file << ",";
        }
        file << "\n";

        std::cout << "Progress: " << std::fixed << std::setprecision(2)
                  << static_cast<double>(step) / STEPS_PER_LATTICEPOINT * 100 << "%\r" << std::flush;
    }
    if (ceiling_hits > 0)
        std::cerr << "WARNING: " << ceiling_hits << " deposition(s) hit the box ceiling H=" << H
                  << " (L=" << L << ", N=" << N_SPECIES << "). Increase H.\n";
}

int main(int argc, char *argv[])
{
    int L = DEFAULT_L;
    int N_SPECIES = DEFAULT_N_SPECIES;
    int STEPS_PER_LATTICEPOINT = DEFAULT_STEPS_PER_LATTICEPOINT;
    int SIM_NO = 0;
    int H_BOX = -1; // -1 => default to STEPS_PER_LATTICEPOINT (original behaviour)
    if (argc > 1)
        L = std::stoi(argv[1]);
    if (argc > 2)
        N_SPECIES = std::stoi(argv[2]);
    if (argc > 3)
        STEPS_PER_LATTICEPOINT = std::stoi(argv[3]);
    if (argc > 4)
        SIM_NO = std::stoi(argv[4]);
    if (argc > 5)
        H_BOX = std::stoi(argv[5]);
    if (H_BOX <= 0)
        H_BOX = STEPS_PER_LATTICEPOINT;

    gen.seed(2654435761u * static_cast<unsigned>(SIM_NO + 1));  // reproducible per sim

    std::filesystem::path exeDir = std::filesystem::path(argv[0]).parent_path();
    std::filesystem::path filePath = exeDir / "outputs" / "slopeDist" / ("L_" + std::to_string(L) + "_N_" + std::to_string(N_SPECIES) + "_steps_" + std::to_string(STEPS_PER_LATTICEPOINT) + "_sim_" + std::to_string(SIM_NO) + ".tsv");

    std::filesystem::create_directories(filePath.parent_path());

    std::ofstream file;
    file.open(filePath);
    if (!file.is_open())
    {
        std::cerr << "Failed to open output file: " << filePath.string() << "\n";
        return 1;
    }
    run(file, L, N_SPECIES, STEPS_PER_LATTICEPOINT, H_BOX);

    file.close();

    return 0;
}
