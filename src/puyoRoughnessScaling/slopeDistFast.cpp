// Optimized Puyo deposition simulation.
//
// Produces BIT-IDENTICAL output to slopeDist.cpp (same deterministic seeding,
// same dynamics) but with per-drop cost O(active region) instead of O(H*L), so
// a tall box costs only memory, never compute.  This is what makes large N and
// long times affordable.  Correctness is verified by diffing against
// slopeDist.cpp over many (L, N, steps, seed) combinations.
//
// Key data structures (vs the reference's full-lattice scans):
//   colH[c]          : number of blocks in column c   -> O(1) placement, O(L) heights
//   movedList        : sites that just moved          -> annihilation seeds (no H*L scan)
//   visitedGen[]     : generation-stamped BFS visited  -> no per-seed alloc/clear
//   dirty columns    : only these are re-compacted by gravity
//
// CLI:  L  N_SPECIES  STEPS  SIM_NO  H_BOX  [REC_INTERVAL]
//   REC_INTERVAL (optional, default 1): record every step for step<256, then
//   only every REC_INTERVAL-th step (plus the final step).
#include <random>
#include <vector>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <filesystem>
#include <memory>

#pragma GCC optimize("inline", "unroll-loops", "no-stack-protector")
#pragma GCC target("sse,sse2,sse3,ssse3,sse4,popcnt,abm,mmx,avx,avx2,tune=native", "f16c")

static auto _ = []()
{ std::ios_base::sync_with_stdio(false); std::cin.tie(nullptr); std::cout.tie(nullptr); return 0; }();

std::mt19937 gen;  // seeded deterministically from SIM_NO in main

constexpr int DEFAULT_L = 128;
constexpr int DEFAULT_N_SPECIES = 6;
constexpr int DEFAULT_STEPS_PER_LATTICEPOINT = 128;

void run(std::ofstream &file, int L, int N_SPECIES, int STEPS_PER_LATTICEPOINT, int H, int REC_INTERVAL)
{
    std::uniform_int_distribution<> dis_l(0, L - 1);
    std::uniform_int_distribution<> dis_species(1, N_SPECIES);

    file << "step\tfirst_col_height\tslope_distribution\n";

    long long ceiling_hits = 0;

    // column c occupies rows [H - colH[c], H - 1]; row H-1 is the floor
    std::vector<int> lat(static_cast<size_t>(H) * L, 0);
    std::vector<int> colH(L, 0);

    std::vector<int> visitedGen(static_cast<size_t>(H) * L, 0);  // BFS stamp
    std::vector<int> colDirtyGen(L, 0);                          // dirty-column stamp
    std::vector<int> lowestElim(L, 0);                           // per dirty col
    int stamp = 0;

    std::vector<int> movedList, newMovedList, dirtyCols, component, bfs;
    movedList.reserve(1024); newMovedList.reserve(1024);
    dirtyCols.reserve(256); component.reserve(1024); bfs.reserve(1024);

    auto annihilate = [&]()
    {
        // eliminates every maximal same-species connected component that contains
        // >=1 moved seed and has size > 1; records dirty columns + lowest elim row
        int g = ++stamp;
        dirtyCols.clear();
        for (int s : movedList)
        {
            if (lat[s] == 0 || visitedGen[s] == g)
                continue;
            int sp = lat[s];
            component.clear();
            bfs.clear();
            bfs.push_back(s);
            visitedGen[s] = g;
            size_t head = 0;
            while (head < bfs.size())
            {
                int cur = bfs[head++];
                component.push_back(cur);
                int r = cur / L, c = cur % L;
                if (r > 0) { int nb = cur - L; if (lat[nb] == sp && visitedGen[nb] != g) { visitedGen[nb] = g; bfs.push_back(nb); } }
                if (r + 1 < H) { int nb = cur + L; if (lat[nb] == sp && visitedGen[nb] != g) { visitedGen[nb] = g; bfs.push_back(nb); } }
                { int nc = (c - 1 + L) % L; int nb = r * L + nc; if (lat[nb] == sp && visitedGen[nb] != g) { visitedGen[nb] = g; bfs.push_back(nb); } }
                { int nc = (c + 1) % L; int nb = r * L + nc; if (lat[nb] == sp && visitedGen[nb] != g) { visitedGen[nb] = g; bfs.push_back(nb); } }
            }
            if (component.size() > 1)
            {
                for (int id : component)
                {
                    lat[id] = 0;
                    int c = id % L, r = id / L;
                    if (colDirtyGen[c] != g) { colDirtyGen[c] = g; lowestElim[c] = r; dirtyCols.push_back(c); }
                    else if (r > lowestElim[c]) lowestElim[c] = r;
                }
            }
        }
    };

    auto fallDirty = [&]()
    {
        // compact each dirty column above its lowest gap; blocks below are inert.
        newMovedList.clear();
        for (int c : dirtyCols)
        {
            int Le = lowestElim[c];        // lowest (largest-index) eliminated row -> now empty
            int topRow = H - colH[c];      // pre-elimination top (upper bound on filled rows)
            int write = Le;
            for (int r = Le - 1; r >= topRow; --r)
            {
                int id = r * L + c;
                if (lat[id] != 0)
                {
                    int dst = write * L + c;   // write < r always here, so it moved
                    lat[dst] = lat[id];
                    lat[id] = 0;
                    newMovedList.push_back(dst);
                    --write;
                }
            }
            colH[c] = H - 1 - write;
        }
    };

    for (int step = 0; step <= STEPS_PER_LATTICEPOINT; ++step)
    {
        for (int i = 0; i < L; ++i)
        {
            int col = dis_l(gen);        // same draw order as before -> identical output
            int species = dis_species(gen);

            if (colH[col] >= H)          // column full -> ceiling hit
            {
                ++ceiling_hits;
                continue;
            }
            int pos = (H - 1 - colH[col]) * L + col;
            lat[pos] = species;
            ++colH[col];

            movedList.clear();
            movedList.push_back(pos);
            while (true)
            {
                annihilate();
                fallDirty();
                if (newMovedList.empty())
                    break;
                movedList.swap(newMovedList);
            }
        }

        if (step < 256 || step % REC_INTERVAL == 0 || step == STEPS_PER_LATTICEPOINT)
        {
            file << std::fixed << std::setprecision(6) << static_cast<double>(step) << "\t" << colH[0] << "\t";
            for (int c = 0; c < L; ++c)
            {
                int next = (c + 1) % L;
                int slope = colH[next] - colH[c];
                file << slope;
                if (c < L - 1)
                    file << ",";
            }
            file << "\n";
        }

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
    int H_BOX = -1;
    int REC_INTERVAL = 1;
    if (argc > 1) L = std::stoi(argv[1]);
    if (argc > 2) N_SPECIES = std::stoi(argv[2]);
    if (argc > 3) STEPS_PER_LATTICEPOINT = std::stoi(argv[3]);
    if (argc > 4) SIM_NO = std::stoi(argv[4]);
    if (argc > 5) H_BOX = std::stoi(argv[5]);
    if (argc > 6) REC_INTERVAL = std::stoi(argv[6]);
    if (H_BOX <= 0) H_BOX = STEPS_PER_LATTICEPOINT;
    if (REC_INTERVAL <= 0) REC_INTERVAL = 1;

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
    run(file, L, N_SPECIES, STEPS_PER_LATTICEPOINT, H_BOX, REC_INTERVAL);
    file.close();
    return 0;
}
