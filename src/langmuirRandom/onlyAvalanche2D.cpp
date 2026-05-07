#include <random>
#include <vector>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <filesystem>
#include <queue>
#include <algorithm>
#include <memory>

#pragma GCC optimize("inline", "unroll-loops", "no-stack-protector")
#pragma GCC target("sse,sse2,sse3,ssse3,sse4,popcnt,abm,mmx,avx,avx2,tune=native", "f16c")

static auto _ = []()
{std::ios_base::sync_with_stdio(false);std::cin.tie(nullptr);std::cout.tie(nullptr);return 0; }();

std::random_device rd;
std::mt19937 gen(rd());

constexpr int DEFAULT_L = 128;
constexpr int DEFAULT_N_SPECIES = 6;
constexpr int DEFAULT_STEPS_PER_LATTICEPOINT = 128;

inline int index2D(int row, int col, int L)
{
    return row * L + col;
}

std::vector<std::vector<bool>> createInteractionMatrix(int nSpecies)
{
    std::vector<std::vector<bool>> J(nSpecies + 1, std::vector<bool>(nSpecies + 1, false));
    std::bernoulli_distribution dis_interaction(0.5);

    for (int i = 1; i <= nSpecies; ++i)
    {
        for (int j = i; j <= nSpecies; ++j)
        {
            bool canAnnihilate = dis_interaction(gen);
            J[i][j] = canAnnihilate;
            J[j][i] = canAnnihilate;
        }
    }

    for (int i = 1; i <= nSpecies; ++i)
    {
        bool rowHasOne = false;
        for (int j = 1; j <= nSpecies; ++j)
        {
            if (J[i][j])
            {
                rowHasOne = true;
                break;
            }
        }
        if (!rowHasOne)
            J[i][i] = true;
    }

    return J;
}

void writeInteractionMatrix(std::ofstream &file, const std::vector<std::vector<bool>> &J)
{
    int nSpecies = static_cast<int>(J.size()) - 1;
    file << "# J matrix (symmetric, 1-indexed)\n";
    for (int i = 1; i <= nSpecies; ++i)
    {
        file << "#";
        for (int j = 1; j <= nSpecies; ++j)
            file << (j == 1 ? " " : "\t") << (J[i][j] ? 1 : 0);
        file << "\n";
    }
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
                   const std::vector<std::vector<bool>> &J, int L, int H, int &avalancheCount)
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
                ++avalancheCount;
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

void run(std::ofstream &file, int L, int N_SPECIES, int STEPS_PER_LATTICEPOINT)
{
    std::uniform_int_distribution<> dis_l(0, L - 1);
    std::uniform_int_distribution<> dis_species(1, N_SPECIES);
    std::vector<std::vector<bool>> J = createInteractionMatrix(N_SPECIES);

    writeInteractionMatrix(file, J);
    file << "# step\tavalanches\ttotal_eliminated\n";

    int H = STEPS_PER_LATTICEPOINT;
    int total_drops = L * STEPS_PER_LATTICEPOINT;
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
    for (int i = 0; i < L * STEPS_PER_LATTICEPOINT; ++i)
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

        int avalancheCount = 0;
        int totalEliminated = 0;

        while (true)
        {
            int eliminated = annihilatePuyo(lattice, movedSites.get(), J, L, H, avalancheCount);
            totalEliminated += eliminated;

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

        // Only write if there was an avalanche
        if (avalancheCount > 0) {
            file << std::fixed << std::setprecision(6)
                 << static_cast<double>(i) / L << "\t"
                 << avalancheCount << "\t"
                 << totalEliminated << "\n";
        }
    }
}

int main(int argc, char *argv[])
{
    int L = DEFAULT_L;
    int N_SPECIES = DEFAULT_N_SPECIES;
    int STEPS_PER_LATTICEPOINT = DEFAULT_STEPS_PER_LATTICEPOINT;
    int SIM_NO = 0;
    if (argc > 1)
        L = std::stoi(argv[1]);
    if (argc > 2)
        N_SPECIES = std::stoi(argv[2]);
    if (argc > 3)
        STEPS_PER_LATTICEPOINT = std::stoi(argv[3]);
    if (argc > 4)
        SIM_NO = std::stoi(argv[4]);

    std::filesystem::path exeDir = std::filesystem::path(argv[0]).parent_path();
    std::ostringstream fileNameStream;
    fileNameStream << "L_" << L << "_N_" << N_SPECIES << "_" << SIM_NO << ".tsv";
    std::filesystem::path filePath = exeDir / "outputs" / "avalanche2D" / fileNameStream.str();

    std::filesystem::create_directories(filePath.parent_path());

    std::ofstream file;
    file.open(filePath);

    run(file, L, N_SPECIES, STEPS_PER_LATTICEPOINT);

    file.close();

    return 0;
}