// Finite-size scaling of the growth-arrest transition at N_c.
//
// The question this exists to answer: is the sublinear growth <h> ~ t^0.32 at
// N_c genuine critical scaling, or a slow crossover to linear growth at a small
// positive velocity?  Two decades at a single L cannot tell the difference.
//
// The discriminating observable is the velocity's finite-size scaling.  If the
// transition is a real continuous one with correlation length xi ~ |N-N_c|^-nu,
// then standard FSS says
//
//     v(N, L) = L^{-theta/nu} G( (N-N_c) L^{1/nu} )
//
// so AT N_c the velocity must vanish as a power of L, v(N_c, L) ~ L^{-theta/nu},
// and h(t) ~ t^{1-theta/(nu z)} only for t << L^z before crossing over to linear
// growth at that L-dependent velocity.  If instead the measured N_c sits slightly
// above the true one, v(N_c, L) converges to a nonzero constant as L grows and
// the whole t^0.32 regime is a transient.  These two predictions differ by a
// power of L, which is not a subtle distinction.
//
// So the primary output here is h(t) sampled at log-spaced times, for many L at
// several N bracketing N_c.  v is extracted from the late-time slope, and W(t)
// (interface width) comes along for free since it costs one extra O(L) pass.
//
// DYNAMICS ARE IDENTICAL to ../avalancheScaling/avalancheDist.cpp (which in turn
// carries ../slopeDistFast.cpp's optimizations): continuous N via a fractional
// last species, BFS annihilation of same-species components >= 2 seeded from
// moved sites, gravity on dirty columns only.  Nothing about the model is
// changed here -- only what is recorded.  Verified bit-identical against
// avalancheDist.cpp's moments output (see readme).
//
// CLI:  L  N_SPECIES  STEPS  SIM_NO  H_BOX
#include <random>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <filesystem>
#include <cstdint>
#include <algorithm>
#include <cmath>

#pragma GCC optimize("inline", "unroll-loops", "no-stack-protector")
#pragma GCC target("sse,sse2,sse3,ssse3,sse4,popcnt,abm,mmx,avx,avx2,tune=native", "f16c")

static auto _ = []()
{ std::ios_base::sync_with_stdio(false); std::cin.tie(nullptr); std::cout.tie(nullptr); return 0; }();

std::mt19937 gen;

constexpr int NSAMPLE = 400;   // log-spaced time samples of h(t), W(t)

std::discrete_distribution<> createSpeciesDistribution(double N_SPECIES)
{
    int n_int = static_cast<int>(std::floor(N_SPECIES));
    double frac = N_SPECIES - n_int;
    std::vector<double> w(n_int, 1.0);
    if (frac > 0)
        w.push_back(frac);
    return std::discrete_distribution<>(w.begin(), w.end());
}

std::string fmtN(double N)
{
    std::ostringstream s;
    s << std::fixed << std::setprecision(4) << N;
    return s.str();
}

void run(std::ofstream &file, int L, double N_SPECIES, int STEPS, int H)
{
    std::uniform_int_distribution<> dis_l(0, L - 1);
    std::discrete_distribution<> dis_species = createSpeciesDistribution(N_SPECIES);

    // log-spaced sample steps, deduplicated (small t has fewer distinct integers
    // than samples, so the early part of the schedule collapses onto every step)
    std::vector<int> sampleAt;
    for (int k = 0; k < NSAMPLE; ++k)
    {
        int s = static_cast<int>(std::llround(std::pow(static_cast<double>(STEPS),
                                                       static_cast<double>(k) / (NSAMPLE - 1))));
        if (sampleAt.empty() || s > sampleAt.back())
            sampleAt.push_back(s);
    }
    size_t nextSample = 0;
    std::vector<int> outStep;
    std::vector<double> outH, outW, outActive, outMeanS;

    long long ceiling_hits = 0;

    std::vector<uint8_t> lat(static_cast<size_t>(H) * L, 0);
    std::vector<int> colH(L, 0);

    std::vector<int> visitedGen(static_cast<size_t>(H) * L, 0);
    std::vector<int> colDirtyGen(L, 0);
    std::vector<int> lowestElim(L, 0);
    int stamp = 0;

    std::vector<int> movedList, newMovedList, dirtyCols, component, bfs;
    movedList.reserve(1024); newMovedList.reserve(1024);
    dirtyCols.reserve(256); component.reserve(1024); bfs.reserve(1024);

    long long casMass = 0;
    // activity accumulators, reset after each sample so each reported value is
    // local in time rather than a running average over the whole transient
    long long segDrops = 0, segActive = 0, segMass = 0;

    auto annihilate = [&]()
    {
        int g = ++stamp;
        dirtyCols.clear();
        long long genMass = 0;
        for (int s : movedList)
        {
            if (lat[s] == 0 || visitedGen[s] == g)
                continue;
            uint8_t sp = lat[s];
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
                genMass += static_cast<long long>(component.size());
            }
        }
        if (genMass > 0)
            casMass += genMass;
    };

    auto fallDirty = [&]()
    {
        newMovedList.clear();
        for (int c : dirtyCols)
        {
            int Le = lowestElim[c];
            int topRow = H - colH[c];
            int write = Le;
            for (int r = Le - 1; r >= topRow; --r)
            {
                int id = r * L + c;
                if (lat[id] != 0)
                {
                    int dst = write * L + c;
                    lat[dst] = lat[id];
                    lat[id] = 0;
                    newMovedList.push_back(dst);
                    --write;
                }
            }
            colH[c] = H - 1 - write;
        }
    };

    for (int step = 0; step <= STEPS; ++step)
    {
        for (int i = 0; i < L; ++i)
        {
            int col = dis_l(gen);
            int species = dis_species(gen) + 1;

            if (colH[col] >= H)
            {
                ++ceiling_hits;
                continue;
            }

            int pos = (H - 1 - colH[col]) * L + col;
            lat[pos] = static_cast<uint8_t>(species);
            ++colH[col];

            casMass = 0;
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

            ++segDrops;
            if (casMass > 0) { ++segActive; segMass += casMass; }
        }

        if (nextSample < sampleAt.size() && step == sampleAt[nextSample])
        {
            double sum = 0;
            for (int c = 0; c < L; ++c) sum += colH[c];
            double mean = sum / L;
            double var = 0;
            for (int c = 0; c < L; ++c) { double d = colH[c] - mean; var += d * d; }
            outStep.push_back(step);
            outH.push_back(mean);
            outW.push_back(std::sqrt(var / L));
            outActive.push_back(segDrops ? static_cast<double>(segActive) / segDrops : 0.0);
            outMeanS.push_back(segActive ? static_cast<double>(segMass) / segActive : 0.0);
            segDrops = segActive = segMass = 0;
            ++nextSample;
        }
    }

    file << "# L=" << L << " N=" << std::setprecision(6) << N_SPECIES
         << " steps=" << STEPS << " H=" << H << " ceiling_hits=" << ceiling_hits << "\n";
    file << "step\tmean_h\twidth\tactive_frac\tmean_s\n";
    for (size_t i = 0; i < outStep.size(); ++i)
        file << outStep[i] << "\t" << std::fixed << std::setprecision(6)
             << outH[i] << "\t" << outW[i] << "\t" << outActive[i] << "\t" << outMeanS[i] << "\n";

    if (ceiling_hits > 0)
        std::cerr << "WARNING: " << ceiling_hits << " deposition(s) hit ceiling H=" << H
                  << " (L=" << L << ", N=" << N_SPECIES << "). Increase H.\n";
}

int main(int argc, char *argv[])
{
    int L = 128;
    double N_SPECIES = 6.0;
    int STEPS = 128;
    int SIM_NO = 0;
    int H_BOX = -1;
    if (argc > 1) L = std::stoi(argv[1]);
    if (argc > 2) N_SPECIES = std::stod(argv[2]);
    if (argc > 3) STEPS = std::stoi(argv[3]);
    if (argc > 4) SIM_NO = std::stoi(argv[4]);
    if (argc > 5) H_BOX = std::stoi(argv[5]);
    if (H_BOX <= 0) H_BOX = STEPS;

    gen.seed(2654435761u * static_cast<unsigned>(SIM_NO + 1));

    std::string tag = "L_" + std::to_string(L) + "_N_" + fmtN(N_SPECIES) +
                      "_steps_" + std::to_string(STEPS) +
                      "_sim_" + std::to_string(SIM_NO) + ".tsv";
    std::filesystem::path exeDir = std::filesystem::path(argv[0]).parent_path();
    std::filesystem::path filePath = exeDir / "outputs" / "growth" / tag;
    std::filesystem::create_directories(filePath.parent_path());

    std::ofstream file(filePath);
    if (!file.is_open())
    {
        std::cerr << "Failed to open " << filePath.string() << "\n";
        return 1;
    }
    run(file, L, N_SPECIES, STEPS, H_BOX);
    file.close();
    return 0;
}
