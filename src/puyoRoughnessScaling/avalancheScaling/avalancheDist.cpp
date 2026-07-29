// Avalanche-size distribution for Puyo deposition, with the O(active region)
// optimizations of ../slopeDistFast.cpp so that large L is affordable.
//
// One "avalanche" = the entire cascade triggered by a single deposition:
//   mass     s : total number of blocks eliminated before the pile is quiet again
//   clusters n : number of maximal same-species components eliminated
//   duration d : number of chain generations that eliminated at least one block
// Depositions that eliminate nothing (s = 0) are counted but not histogrammed.
//
// Histograms are accumulated in RAM and dumped once at the end, so the per-drop
// cost is pure dynamics and the output is O(s_max) rather than O(L * steps).
//
// Two files are written:
//   outputs/avalancheDist/... : the s/n/d histograms over t >= WARMUP
//   outputs/moments/...       : per-time-window moments over ALL t, including
//     the transient.  These exist to *verify* steady state rather than assume
//     it: <s^2>/<s> vs t must plateau well before WARMUP.  Windows are log
//     spaced, since the approach to steady state is a power law in t.
//
// Data structures inherited from slopeDistFast.cpp:
//   colH[c]     : blocks in column c        -> O(1) placement, no lattice scan
//   movedList   : sites that just moved     -> annihilation seeds (no H*L scan)
//   visitedGen  : generation-stamped BFS    -> no per-seed alloc/clear
//   dirtyCols   : only these are compacted by gravity
// lat[] is uint8_t (species <= 255) to keep tall boxes cheap in memory.
//
// CLI:  L  N_SPECIES  STEPS  SIM_NO  H_BOX  [WARMUP]
//   WARMUP (default STEPS/4): steps discarded before histogramming, so that
//   only the statistically steady surface contributes.
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

std::mt19937 gen;  // seeded deterministically from SIM_NO in main

constexpr int DEFAULT_L = 128;
constexpr double DEFAULT_N_SPECIES = 6.0;
constexpr int DEFAULT_STEPS_PER_LATTICEPOINT = 128;

// N_SPECIES is CONTINUOUS (as in probabilityPuyoPuyo/onlyAvalanche2D.cpp): the
// first floor(N) species get weight 1 and a final partial species gets the
// fractional part, so total weight = N and each full species has probability
// 1/N.  This is what lets us tune through the growth-arrest transition at
// N_c in (5,6), where the interesting physics is.
std::discrete_distribution<> createSpeciesDistribution(double N_SPECIES)
{
    int n_int = static_cast<int>(std::floor(N_SPECIES));
    double frac = N_SPECIES - n_int;
    std::vector<double> w(n_int, 1.0);
    if (frac > 0)
        w.push_back(frac);
    return std::discrete_distribution<>(w.begin(), w.end());
}

// N appears in filenames; fix the format so globs are predictable
std::string fmtN(double N)
{
    std::ostringstream s;
    s << std::fixed << std::setprecision(3) << N;
    return s.str();
}

// grow-on-demand histogram: bump(h, v) increments the bin for value v
static inline void bump(std::vector<long long> &h, size_t v)
{
    if (v >= h.size())
        h.resize(v + 1, 0);
    ++h[v];
}

constexpr int NWIN = 64;   // log-spaced time windows for the steady-state check
constexpr int MRANGE = 48; // local slope m is binned over [-MRANGE, MRANGE]

void run(std::ofstream &file, std::ofstream &mfile, std::ofstream &sfile,
         int L, double N_SPECIES, int STEPS_PER_LATTICEPOINT, int H, int WARMUP)
{
    // Slope-resolved cascade mass: <s|m> vs m, a DIRECT test of the mechanism
    // P(m) ~ exp(-lambda m) & s ~ m^d  =>  P(s) is Weibull with shape 1/d.
    // m is the local slope h[c+1]-h[c] at the deposition column, read BEFORE
    // the block lands.
    std::vector<long long> mCount(2 * MRANGE + 1, 0), mSum(2 * MRANGE + 1, 0),
        mSum2(2 * MRANGE + 1, 0), mAct(2 * MRANGE + 1, 0);

    // Spatial extent w = number of distinct columns the cascade eliminates in.
    // If P(w) ~ exp(-w/w0) and s ~ w^d, then s is Weibull with shape 1/d -- the
    // same argument as for the slope, but with a variable that actually
    // correlates with mass.  wSum[w] gives <s|w> directly.
    std::vector<int> colTouch(L, 0);
    int casStamp = 0, casWidth = 0;
    std::vector<long long> wCount, wSum, wSum2;
    std::uniform_int_distribution<> dis_l(0, L - 1);
    std::discrete_distribution<> dis_species = createSpeciesDistribution(N_SPECIES);

    // step -> window index, precomputed so the hot loop does one array lookup
    std::vector<int> winOf(STEPS_PER_LATTICEPOINT + 1, 0);
    std::vector<double> winEdge(NWIN + 1);
    for (int w = 0; w <= NWIN; ++w)
        winEdge[w] = std::pow(static_cast<double>(STEPS_PER_LATTICEPOINT + 1),
                              static_cast<double>(w) / NWIN);
    for (int s = 0; s <= STEPS_PER_LATTICEPOINT; ++s)
    {
        int w = 0;
        while (w < NWIN - 1 && s + 1 > winEdge[w + 1]) ++w;
        winOf[s] = w;
    }
    std::vector<long long> wDrops(NWIN, 0), wActive(NWIN, 0), wS(NWIN, 0), wS2(NWIN, 0), wMax(NWIN, 0);
    // mean pile height per window -> growth velocity v = d<h>/dt, the order
    // parameter of the growth-arrest transition (v -> 0 at N_c)
    std::vector<double> wH(NWIN, 0.0);
    std::vector<long long> wHn(NWIN, 0);

    long long ceiling_hits = 0;
    long long drops_counted = 0;   // depositions after warmup (normalization)
    long long drops_active = 0;    // of those, ones that triggered an elimination

    // column c occupies rows [H - colH[c], H - 1]; row H-1 is the floor
    std::vector<uint8_t> lat(static_cast<size_t>(H) * L, 0);
    std::vector<int> colH(L, 0);

    std::vector<int> visitedGen(static_cast<size_t>(H) * L, 0);  // BFS stamp
    std::vector<int> colDirtyGen(L, 0);                          // dirty-column stamp
    std::vector<int> lowestElim(L, 0);                           // per dirty col
    int stamp = 0;

    std::vector<int> movedList, newMovedList, dirtyCols, component, bfs;
    movedList.reserve(1024); newMovedList.reserve(1024);
    dirtyCols.reserve(256); component.reserve(1024); bfs.reserve(1024);

    std::vector<long long> histMass, histClusters, histDuration;

    // per-cascade accumulators, reset at each deposition
    long long casMass = 0;
    int casClusters = 0, casDuration = 0;

    auto annihilate = [&]()
    {
        // eliminates every maximal same-species connected component that contains
        // >=1 moved seed and has size > 1; records dirty columns + lowest elim row
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
                    // distinct columns touched across the WHOLE cascade
                    if (colTouch[c] != casStamp) { colTouch[c] = casStamp; ++casWidth; }
                }
                genMass += static_cast<long long>(component.size());
                ++casClusters;
            }
        }
        if (genMass > 0)
        {
            casMass += genMass;
            ++casDuration;
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
        bool record = (step >= WARMUP);
        int w = winOf[step];
        for (int i = 0; i < L; ++i)
        {
            int col = dis_l(gen);
            int species = dis_species(gen) + 1;   // 1-indexed; 0 means empty

            if (colH[col] >= H)          // column full -> ceiling hit
            {
                ++ceiling_hits;
                continue;
            }
            int mloc = colH[(col + 1) % L] - colH[col];   // local slope, pre-landing

            int pos = (H - 1 - colH[col]) * L + col;
            lat[pos] = static_cast<uint8_t>(species);
            ++colH[col];

            casMass = 0; casClusters = 0; casDuration = 0;
            casWidth = 0; ++casStamp;   // stamp is per-cascade, distinct from `stamp`

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

            // steady-state check: accumulated over ALL t, transient included
            ++wDrops[w];
            if (casMass > 0)
            {
                ++wActive[w];
                wS[w] += casMass;
                wS2[w] += casMass * casMass;
                if (casMass > wMax[w]) wMax[w] = casMass;
            }

            if (record)
            {
                ++drops_counted;
                int mi = std::clamp(mloc, -MRANGE, MRANGE) + MRANGE;
                ++mCount[mi];
                mSum[mi] += casMass;
                mSum2[mi] += casMass * casMass;
                if (casMass > 0)
                {
                    ++mAct[mi];
                    ++drops_active;
                    bump(histMass, static_cast<size_t>(casMass));
                    bump(histClusters, static_cast<size_t>(casClusters));
                    bump(histDuration, static_cast<size_t>(casDuration));
                    size_t wi = static_cast<size_t>(casWidth);
                    if (wi >= wCount.size())
                    {
                        wCount.resize(wi + 1, 0); wSum.resize(wi + 1, 0); wSum2.resize(wi + 1, 0);
                    }
                    ++wCount[wi]; wSum[wi] += casMass; wSum2[wi] += casMass * casMass;
                }
            }
        }

        {   // O(L) per step, negligible against the L depositions above
            double sum = 0;
            for (int c = 0; c < L; ++c) sum += colH[c];
            wH[w] += sum / L;
            ++wHn[w];
        }

        if (step % std::max(1, STEPS_PER_LATTICEPOINT / 100) == 0)
            std::cout << "Progress: " << static_cast<double>(step) / STEPS_PER_LATTICEPOINT * 100 << "%\r" << std::flush;
    }

    // metadata as comments, then one row per value; a value's bins are 0 where
    // that observable never took it (the three histograms share a value axis).
    file << "# L=" << L << " N=" << N_SPECIES << " steps=" << STEPS_PER_LATTICEPOINT
         << " H=" << H << " warmup=" << WARMUP << "\n";
    // steps_recorded is the meaningful measure of run length; drops_counted is
    // only here because P(s) must be normalized per deposition.
    file << "# steps_recorded=" << (STEPS_PER_LATTICEPOINT - WARMUP + 1)
         << " drops_counted=" << drops_counted << " drops_active=" << drops_active
         << " ceiling_hits=" << ceiling_hits << "\n";
    file << "value\tmass\tclusters\tduration\n";
    size_t vmax = std::max({histMass.size(), histClusters.size(), histDuration.size()});
    for (size_t v = 1; v < vmax; ++v)
    {
        long long m = v < histMass.size() ? histMass[v] : 0;
        long long c = v < histClusters.size() ? histClusters[v] : 0;
        long long d = v < histDuration.size() ? histDuration[v] : 0;
        if (m == 0 && c == 0 && d == 0)
            continue;
        file << v << "\t" << m << "\t" << c << "\t" << d << "\n";
    }

    // moments per time window: lets the analysis SHOW that <s^2>/<s> plateaus
    // long before WARMUP, instead of assuming it.
    mfile << "# L=" << L << " N=" << N_SPECIES << " steps=" << STEPS_PER_LATTICEPOINT
          << " H=" << H << " warmup=" << WARMUP << "\n";
    mfile << "step_lo\tstep_hi\tdrops\tactive\tsum_s\tsum_s2\ts_max\tmean_h\n";
    for (int w = 0; w < NWIN; ++w)
    {
        if (wDrops[w] == 0)
            continue;
        long long lo = static_cast<long long>(std::floor(winEdge[w])) - 1;
        long long hi = static_cast<long long>(std::floor(winEdge[w + 1])) - 1;
        if (w == NWIN - 1) hi = STEPS_PER_LATTICEPOINT;
        mfile << lo << "\t" << hi << "\t" << wDrops[w] << "\t" << wActive[w] << "\t"
              << wS[w] << "\t" << wS2[w] << "\t" << wMax[w] << "\t"
              << std::fixed << std::setprecision(4)
              << (wHn[w] ? wH[w] / wHn[w] : 0.0) << "\n";
    }

    // Slope-resolved cascade mass, and the final species composition of the
    // pile.  The composition tests the sawtooth mechanism: if a fractional
    // ("impurity") species is enriched in the pile relative to how often it is
    // deposited, it is a frozen defect that fragments clusters.
    int nsp = static_cast<int>(std::floor(N_SPECIES)) + (N_SPECIES > std::floor(N_SPECIES) ? 1 : 0);
    std::vector<long long> pileSp(nsp + 1, 0);
    for (size_t i = 0; i < lat.size(); ++i)
        if (lat[i]) ++pileSp[lat[i]];
    long long pileTot = 0;
    for (int i = 1; i <= nsp; ++i) pileTot += pileSp[i];

    sfile << "# L=" << L << " N=" << N_SPECIES << " steps=" << STEPS_PER_LATTICEPOINT
          << " H=" << H << " warmup=" << WARMUP << "\n";
    sfile << "# pile_total=" << pileTot << " pile_by_species=";
    for (int i = 1; i <= nsp; ++i)
        sfile << pileSp[i] << (i < nsp ? "," : "");
    sfile << "\n";
    sfile << "m\tdrops\tactive\tsum_s\tsum_s2\n";
    for (int i = 0; i <= 2 * MRANGE; ++i)
    {
        if (mCount[i] == 0)
            continue;
        sfile << (i - MRANGE) << "\t" << mCount[i] << "\t" << mAct[i] << "\t"
              << mSum[i] << "\t" << mSum2[i] << "\n";
    }
    sfile << "# extent\n";
    sfile << "w\tcount\tsum_s\tsum_s2\n";
    for (size_t i = 1; i < wCount.size(); ++i)
    {
        if (wCount[i] == 0)
            continue;
        sfile << i << "\t" << wCount[i] << "\t" << wSum[i] << "\t" << wSum2[i] << "\n";
    }

    if (ceiling_hits > 0)
        std::cerr << "WARNING: " << ceiling_hits << " deposition(s) hit the box ceiling H=" << H
                  << " (L=" << L << ", N=" << N_SPECIES << "). Increase H.\n";
}

int main(int argc, char *argv[])
{
    int L = DEFAULT_L;
    double N_SPECIES = DEFAULT_N_SPECIES;
    int STEPS_PER_LATTICEPOINT = DEFAULT_STEPS_PER_LATTICEPOINT;
    int SIM_NO = 0;
    int H_BOX = -1;
    int WARMUP = -1;
    if (argc > 1) L = std::stoi(argv[1]);
    if (argc > 2) N_SPECIES = std::stod(argv[2]);
    if (argc > 3) STEPS_PER_LATTICEPOINT = std::stoi(argv[3]);
    if (argc > 4) SIM_NO = std::stoi(argv[4]);
    if (argc > 5) H_BOX = std::stoi(argv[5]);
    if (argc > 6) WARMUP = std::stoi(argv[6]);
    if (H_BOX <= 0) H_BOX = STEPS_PER_LATTICEPOINT;
    if (WARMUP < 0) WARMUP = STEPS_PER_LATTICEPOINT / 4;

    gen.seed(2654435761u * static_cast<unsigned>(SIM_NO + 1));  // reproducible per sim

    std::string tag = "L_" + std::to_string(L) + "_N_" + fmtN(N_SPECIES) +
                      "_steps_" + std::to_string(STEPS_PER_LATTICEPOINT) +
                      "_sim_" + std::to_string(SIM_NO) + ".tsv";
    std::filesystem::path exeDir = std::filesystem::path(argv[0]).parent_path();
    std::filesystem::path filePath = exeDir / "outputs" / "avalancheDist" / tag;
    std::filesystem::path momPath = exeDir / "outputs" / "moments" / tag;
    std::filesystem::path slopePath = exeDir / "outputs" / "slopeResolved" / tag;

    std::filesystem::create_directories(filePath.parent_path());
    std::filesystem::create_directories(momPath.parent_path());
    std::filesystem::create_directories(slopePath.parent_path());

    std::ofstream file(filePath), mfile(momPath), sfile(slopePath);
    if (!file.is_open() || !mfile.is_open() || !sfile.is_open())
    {
        std::cerr << "Failed to open output files under " << exeDir.string() << "\n";
        return 1;
    }
    run(file, mfile, sfile, L, N_SPECIES, STEPS_PER_LATTICEPOINT, H_BOX, WARMUP);
    file.close();
    mfile.close();
    sfile.close();
    return 0;
}
