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
#include <cmath>
#include <algorithm>

#pragma GCC optimize("inline", "unroll-loops", "no-stack-protector")
#pragma GCC target("sse,sse2,sse3,ssse3,sse4,popcnt,abm,mmx,avx,avx2,tune=native", "f16c")

static auto _ = []()
{std::ios_base::sync_with_stdio(false);std::cin.tie(nullptr);std::cout.tie(nullptr);return 0; }();

std::random_device rd;
std::mt19937 gen(rd());

constexpr int DEFAULT_L = 128;
constexpr int DEFAULT_N_SPECIES = 6;
constexpr int DEFAULT_STEPS_PER_LATTICEPOINT = 128;

// Generate a sorted random probability vector and return a discrete_distribution
std::discrete_distribution<> createRandomSpeciesDistribution(std::vector<double> &sorted_probs, int N_SPECIES)
{
    std::vector<double> weights(N_SPECIES);
    std::uniform_real_distribution<> dis(0.0, 1.0);
    double sum = 0.0;
    for (int i = 0; i < N_SPECIES; ++i)
    {
        weights[i] = dis(gen);
        sum += weights[i];
    }
    for (int i = 0; i < N_SPECIES; ++i)
        weights[i] /= sum;
    sorted_probs = weights;
    std::sort(sorted_probs.begin(), sorted_probs.end(), std::greater<>());
    return std::discrete_distribution<>(sorted_probs.begin(), sorted_probs.end());
}

void placePuyo(std::vector<std::vector<int>> &lattice, std::vector<std::vector<bool>> &movedSites,
               std::uniform_int_distribution<> &dis_l, std::discrete_distribution<> &species_dist, int L, int H)
{
    int col = dis_l(gen);
    int species = species_dist(gen) + 1; // 1-indexed

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

int annihilatePuyo(std::vector<std::vector<int>> &lattice,
                   std::vector<std::vector<bool>> &movedSites,
                   int L, int H, int &avalancheCount)
{
    std::vector<std::vector<bool>> visited(H, std::vector<bool>(L, false));
    int totalEliminated = 0;
    std::vector<std::pair<int, int>> directions = {{-1, 0}, {1, 0}, {0, -1}, {0, 1}};

    for (int row = 0; row < H; ++row)
    {
        for (int col = 0; col < L; ++col)
        {
            if (!movedSites[row][col] || visited[row][col] || lattice[row][col] == 0)
                continue;
            std::queue<std::pair<int, int>> q;
            std::vector<std::pair<int, int>> cluster;
            q.push({row, col});
            visited[row][col] = true;
            int species = lattice[row][col];
            while (!q.empty())
            {
                auto [r, c] = q.front();
                q.pop();
                cluster.emplace_back(r, c);
                for (auto [dr, dc] : directions)
                {
                    int nr = r + dr;
                    int nc = (c + dc + L) % L;
                    if (nr >= 0 && nr < H && !visited[nr][nc] && lattice[nr][nc] == species)
                    {
                        visited[nr][nc] = true;
                        q.push({nr, nc});
                    }
                }
            }
            if (cluster.size() > 1)
            {
                for (auto [r, c] : cluster)
                {
                    lattice[r][c] = 0;
                }
                totalEliminated += cluster.size();
                ++avalancheCount;
            }
        }
    }
    return totalEliminated;
}

void fall(std::vector<std::vector<int>> &lattice, std::vector<std::vector<bool>> &movedSites, int L, int H)
{
    for (int col = 0; col < L; ++col)
    {
        int writeRow = H - 1;
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
        for (int row = writeRow; row >= 0; --row)
        {
            movedSites[row][col] = false;
        }
    }
}

void run(std::ofstream &file, int L, std::discrete_distribution<> &species_dist, int STEPS_PER_LATTICEPOINT)
{
    std::uniform_int_distribution<> dis_l(0, L - 1);

    int H = STEPS_PER_LATTICEPOINT;
    std::vector<std::vector<int>> lattice(H, std::vector<int>(L, 0));
    std::vector<std::vector<bool>> movedSites(H, std::vector<bool>(L, false));

    double step = 0.0;

    for (int i = 0; i < L * STEPS_PER_LATTICEPOINT; ++i)
    {
        placePuyo(lattice, movedSites, dis_l, species_dist, L, H);

        int avalancheCount = 0;
        int totalEliminated = 0;

        while (true)
        {
            int eliminated = annihilatePuyo(lattice, movedSites, L, H, avalancheCount);
            totalEliminated += eliminated;

            std::vector<std::vector<bool>> newMovedSites(H, std::vector<bool>(L, false));
            fall(lattice, newMovedSites, L, H);

            if (newMovedSites == movedSites)
                break;

            movedSites = newMovedSites;
        }

        // Only write if there was an avalanche
        if (avalancheCount > 0) {
            file << std::fixed << std::setprecision(6)
                 << step << "\t"
                 << avalancheCount << "\t"
                 << totalEliminated << "\n";
        }

        step += 1.0 / L;
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

    // Generate sorted random probabilities
    std::vector<double> sorted_probs;
    std::discrete_distribution<> species_dist = createRandomSpeciesDistribution(sorted_probs, N_SPECIES);

    // Build probability string for filename
    std::ostringstream probStream;
    probStream << std::fixed << std::setprecision(4);
    for (size_t i = 0; i < sorted_probs.size(); ++i) {
        if (i > 0) probStream << "-";
        probStream << sorted_probs[i];
    }
    std::string probStr = probStream.str();

    std::string exePath = argv[0];
    std::string exeDir = std::filesystem::path(exePath).parent_path().string();
    std::ostringstream filePathStream;
    filePathStream << exeDir << "\\outputs\\avalanche2D\\onlyAvalancheRandomProbs\\L_" << L << "_P_" << probStr << ".tsv";
    std::string filePath = filePathStream.str();

    std::ofstream file;
    file.open(filePath);
    file << "step\tavalanches\ttotal_eliminated\n";

    run(file, L, species_dist, STEPS_PER_LATTICEPOINT);

    file.close();

    return 0;
}