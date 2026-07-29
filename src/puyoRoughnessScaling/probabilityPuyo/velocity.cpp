// Growth velocity v(N) across the fractional-species interpolation.
//
// WHY: for INTEGER N the species are uniform, so the probability that two
// independently drawn blocks carry the same species is
//
//     sum_i p_i^2 = N * (1/N)^2 = 1/N
//
// and empirically v = 1 - N_c/N is a straight line in 1/N.  The claim to test is
// that 1/N was never the natural variable -- the COLLISION PROBABILITY was, and
// the two coincide only for uniform weights.  For N = n + f the weights are
// (1,...,1,f)/N, so
//
//     sum_i p_i^2 = (n + f^2) / N^2   !=   1/N
//
// which is why v vs non-integer N is not straight.  sum_i p_i^2 is the Simpson
// index = exp(-H_2), H_2 the Renyi entropy of order 2 (the collision entropy),
// and 1/sum p_i^2 is the Hill number of order 2 (the "effective number of
// species").  So the prediction is v = 1 - c * exp(-H_2), linear in exp(-H_2)
// for ALL N, integer or not.
//
// The confound this must also measure: the pile's composition is NOT the
// deposition composition.  A rare fractional species is deposited seldom but
// removable only in pairs, so it accumulates (~5.6x enrichment at N=6.05, see
// avalancheScaling/readme.md).  The collision probability that the DYNAMICS
// actually sees is therefore sum_i rho_i^2 over the PILE densities rho_i, not
// sum_i p_i^2 over the deposition weights.  Both are written out so the analysis
// can decide which one linearises v.
//
// DYNAMICS IDENTICAL to ../criticalScaling/criticalScaling.cpp and
// ../avalancheScaling/avalancheDist.cpp; only the recording differs.
//
// CLI:  L  N_SPECIES  STEPS  SIM_NO  H_BOX  [WARMUP]
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

constexpr int NSAMPLE = 200;   // LINEARLY spaced samples: v is a straight-line fit of h vs t

std::vector<double> speciesWeights(double N_SPECIES)
{
    int n_int = static_cast<int>(std::floor(N_SPECIES));
    double frac = N_SPECIES - n_int;
    std::vector<double> w(n_int, 1.0);
    if (frac > 0)
        w.push_back(frac);
    return w;
}

std::string fmtN(double N)
{
    std::ostringstream s;
    s << std::fixed << std::setprecision(4) << N;
    return s.str();
}

void run(std::ofstream &file, int L, double N_SPECIES, int STEPS, int H, int WARMUP)
{
    std::uniform_int_distribution<> dis_l(0, L - 1);
    std::vector<double> w = speciesWeights(N_SPECIES);
    std::discrete_distribution<> dis_species(w.begin(), w.end());

    long long ceiling_hits = 0;
    // r_i: blocks of species i eliminated, counted only after warmup.  This is
    // the quantity the steady-state balance p_i - r_i = v*rho_i needs, and
    // measuring it directly beats guessing a closure for it.
    std::vector<long long> elimSp(w.size() + 1, 0);
    long long elimDrops = 0;
    bool recording = false;
    std::vector<uint8_t> lat(static_cast<size_t>(H) * L, 0);
    std::vector<int> colH(L, 0);
    std::vector<int> visitedGen(static_cast<size_t>(H) * L, 0);
    std::vector<int> colDirtyGen(L, 0);
    std::vector<int> lowestElim(L, 0);
    int stamp = 0;

    std::vector<int> movedList, newMovedList, dirtyCols, component, bfs;
    movedList.reserve(1024); newMovedList.reserve(1024);
    dirtyCols.reserve(256); component.reserve(1024); bfs.reserve(1024);

    std::vector<int> outStep;
    std::vector<double> outH;

    auto annihilate = [&]()
    {
        int g = ++stamp;
        dirtyCols.clear();
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
                if (recording) elimSp[sp] += static_cast<long long>(component.size());
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

    int sampleEvery = std::max(1, (STEPS - WARMUP) / NSAMPLE);

    for (int step = 0; step <= STEPS; ++step)
    {
        recording = (step >= WARMUP);
        for (int i = 0; i < L; ++i)
        {
            if (recording) ++elimDrops;
            int col = dis_l(gen);
            int species = dis_species(gen) + 1;

            if (colH[col] >= H) { ++ceiling_hits; continue; }

            int pos = (H - 1 - colH[col]) * L + col;
            lat[pos] = static_cast<uint8_t>(species);
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

        if (step >= WARMUP && (step - WARMUP) % sampleEvery == 0)
        {
            double sum = 0;
            for (int c = 0; c < L; ++c) sum += colH[c];
            outStep.push_back(step);
            outH.push_back(sum / L);   // = mass/L exactly: gravity leaves no holes
        }
    }

    // Pile composition at the end: the dynamics collides PILE blocks, not
    // freshly deposited ones, so sum_i rho_i^2 may be the operative collision
    // probability rather than sum_i p_i^2.
    int nsp = static_cast<int>(w.size());
    std::vector<long long> pileSp(nsp + 1, 0);
    for (size_t i = 0; i < lat.size(); ++i)
        if (lat[i]) ++pileSp[lat[i]];
    long long pileTot = 0;
    for (int i = 1; i <= nsp; ++i) pileTot += pileSp[i];

    file << "# L=" << L << " N=" << std::setprecision(6) << N_SPECIES
         << " steps=" << STEPS << " H=" << H << " warmup=" << WARMUP
         << " ceiling_hits=" << ceiling_hits << "\n";
    file << "# pile_total=" << pileTot << " pile_by_species=";
    for (int i = 1; i <= nsp; ++i)
        file << pileSp[i] << (i < nsp ? "," : "");
    file << "\n";
    file << "# elim_drops=" << elimDrops << " elim_by_species=";
    for (size_t i = 1; i <= w.size(); ++i)
        file << elimSp[i] << (i < w.size() ? "," : "");
    file << "\n";
    file << "# deposition_weights=";
    for (int i = 0; i < nsp; ++i)
        file << w[i] << (i + 1 < nsp ? "," : "");
    file << "\n";
    file << "step\tmean_h\n";
    for (size_t i = 0; i < outStep.size(); ++i)
        file << outStep[i] << "\t" << std::fixed << std::setprecision(6) << outH[i] << "\n";

    if (ceiling_hits > 0)
        std::cerr << "WARNING: " << ceiling_hits << " deposition(s) hit ceiling H=" << H
                  << " (L=" << L << ", N=" << N_SPECIES << "). Increase H.\n";
}

int main(int argc, char *argv[])
{
    int L = 512;
    double N_SPECIES = 6.0;
    int STEPS = 30000;
    int SIM_NO = 0;
    int H_BOX = -1, WARMUP = -1;
    if (argc > 1) L = std::stoi(argv[1]);
    if (argc > 2) N_SPECIES = std::stod(argv[2]);
    if (argc > 3) STEPS = std::stoi(argv[3]);
    if (argc > 4) SIM_NO = std::stoi(argv[4]);
    if (argc > 5) H_BOX = std::stoi(argv[5]);
    if (argc > 6) WARMUP = std::stoi(argv[6]);
    if (H_BOX <= 0) H_BOX = STEPS;
    if (WARMUP < 0) WARMUP = STEPS / 4;

    gen.seed(2654435761u * static_cast<unsigned>(SIM_NO + 1));

    std::string tag = "L_" + std::to_string(L) + "_N_" + fmtN(N_SPECIES) +
                      "_steps_" + std::to_string(STEPS) + "_sim_" + std::to_string(SIM_NO) + ".tsv";
    std::filesystem::path exeDir = std::filesystem::path(argv[0]).parent_path();
    std::filesystem::path filePath = exeDir / "outputs" / "velocity" / tag;
    std::filesystem::create_directories(filePath.parent_path());

    std::ofstream file(filePath);
    if (!file.is_open()) { std::cerr << "Failed to open " << filePath.string() << "\n"; return 1; }
    run(file, L, N_SPECIES, STEPS, H_BOX, WARMUP);
    file.close();
    return 0;
}
