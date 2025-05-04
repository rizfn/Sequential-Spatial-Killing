#include <random>
#include <vector>
#include <queue>
#include <iostream>
#include <fstream>
#include <sstream>
#include <filesystem>
#include <iomanip>
#include <algorithm>

#pragma GCC optimize("inline", "unroll-loops", "no-stack-protector")
#pragma GCC target("sse,sse2,sse3,ssse3,sse4,popcnt,abm,mmx,avx,avx2,tune=native", "f16c")

static auto _ = []()
{std::ios_base::sync_with_stdio(false);std::cin.tie(nullptr);std::cout.tie(nullptr);return 0; }();

void floodAnnihilate(std::vector<std::vector<int>> &grid,
                     int r0, int c0,
                     std::vector<std::pair<int, int>> &cluster)
{
    int L = grid.size();
    int species = grid[r0][c0];
    cluster.clear();
    std::vector<std::vector<bool>> vis(L, std::vector<bool>(L, false));
    std::queue<std::pair<int, int>> q;
    q.push({r0, c0});
    vis[r0][c0] = true;
    while (!q.empty())
    {
        auto [r, c] = q.front();
        q.pop();
        cluster.emplace_back(r, c);
        for (auto d : std::vector<std::pair<int, int>>{{1, 0}, {-1, 0}, {0, 1}, {0, -1}})
        {
            int nr = r + d.first, nc = c + d.second;
            if (nr >= 0 && nr < L && nc >= 0 && nc < L && !vis[nr][nc] && grid[nr][nc] == species)
            {
                vis[nr][nc] = true;
                q.push({nr, nc});
            }
        }
    }
    if (cluster.size() > 1)
    {
        for (auto &p : cluster)
            grid[p.first][p.second] = 0;
    }
    else
    {
        cluster.clear();
    }
}

inline bool hasFilledNeighbor(const std::vector<std::vector<int>> &grid, int r, int c)
{
    int L = grid.size();
    static const int D[4][2] = {{1, 0}, {-1, 0}, {0, 1}, {0, -1}};
    for (auto &d : D)
    {
        int nr = r + d[0], nc = c + d[1];
        if (nr >= 0 && nr < L && nc >= 0 && nc < L && grid[nr][nc] > 0)
            return true;
    }
    return false;
}

void run(std::ofstream &out, int L, int N_SPECIES, int STEPS)
{
    std::mt19937 gen(std::random_device{}());
    std::uniform_int_distribution<> dis_sp(1, N_SPECIES);

    // grid and boundary structures
    std::vector<std::vector<int>> grid(L, std::vector<int>(L, 0));
    std::vector<int> boundary;              // stores linear indices r*L + c
    std::vector<int> pos(L * L, -1);        // position in boundary vector or -1
    std::vector<char> isBoundary(L * L, 0); // mask

    // initialize center
    int rc = L / 2, cc = L / 2;
    grid[rc][cc] = dis_sp(gen);
    // add its empty neighbors to boundary
    for (auto d : std::vector<std::pair<int, int>>{{1, 0}, {-1, 0}, {0, 1}, {0, -1}})
    {
        int nr = rc + d.first, nc = cc + d.second;
        if (nr >= 0 && nr < L && nc >= 0 && nc < L)
        {
            int idx = nr * L + nc;
            isBoundary[idx] = 1;
            pos[idx] = boundary.size();
            boundary.push_back(idx);
        }
    }

    std::vector<std::pair<int, int>> cluster;
    for (int t = 0; t < STEPS; ++t)
    {
        if (boundary.empty())
            break;
        // pick random boundary site
        std::uniform_int_distribution<> dis_b(0, (int)boundary.size() - 1);
        int bi = dis_b(gen);
        int idx = boundary[bi];
        int r = idx / L, c = idx % L;
        // remove from boundary
        int last = boundary.back();
        std::swap(boundary[bi], boundary.back());
        pos[last] = bi;
        boundary.pop_back();
        pos[idx] = -1;
        isBoundary[idx] = 0;

        // occupy
        grid[r][c] = dis_sp(gen);

        // add new empty neighbors
        for (auto &d : std::vector<std::pair<int, int>>{{1, 0}, {-1, 0}, {0, 1}, {0, -1}})
        {
            int nr = r + d.first, nc = c + d.second;
            if (nr >= 0 && nr < L && nc >= 0 && nc < L && grid[nr][nc] == 0)
            {
                int nidx = nr * L + nc;
                if (!isBoundary[nidx])
                {
                    isBoundary[nidx] = 1;
                    pos[nidx] = boundary.size();
                    boundary.push_back(nidx);
                }
            }
        }

        // annihilate clusters
        floodAnnihilate(grid, r, c, cluster);
        // update boundary around annihilated
        for (auto &p : cluster)
        {
            int pr = p.first, pc = p.second, pidx = pr * L + pc;
            // this site is empty and may become boundary
            if (hasFilledNeighbor(grid, pr, pc) && !isBoundary[pidx])
            {
                isBoundary[pidx] = 1;
                pos[pidx] = boundary.size();
                boundary.push_back(pidx);
            }
            // its neighbors may lose their boundary status
            for (auto &d : std::vector<std::pair<int, int>>{{1, 0}, {-1, 0}, {0, 1}, {0, -1}})
            {
                int nr = pr + d.first, nc = pc + d.second;
                if (nr >= 0 && nr < L && nc >= 0 && nc < L && grid[nr][nc] == 0)
                {
                    int nidx = nr * L + nc;
                    if (isBoundary[nidx] && !hasFilledNeighbor(grid, nr, nc))
                    {
                        int ppos = pos[nidx], lidx = boundary.back();
                        std::swap(boundary[ppos], boundary.back());
                        pos[lidx] = ppos;
                        boundary.pop_back();
                        pos[nidx] = -1;
                        isBoundary[nidx] = 0;
                    }
                }
            }
        }

        // record stats: output step and the entire lattice as a flattened, comma-separated string
        out << t << "\t";
        for (int r = 0; r < L; ++r)
        {
            for (int c = 0; c < L; ++c)
            {
                out << grid[r][c];
                if (r != L - 1 || c != L - 1)
                    out << ",";
            }
        }
        out << "\n";

        std::cout << "\rProgress: " << std::fixed
                  << std::setprecision(2)
                  << (100.0 * (t + 1) / STEPS) << "% " << std::flush;
    }
    std::cout << "\nDone.\n";
}

int main(int argc, char *argv[])
{
    int L = 32, N = 3, N_STEPS = 1024;
    if (argc > 1)
        L = std::stoi(argv[1]);
    if (argc > 2)
        N = std::stoi(argv[2]);
    if (argc > 3)
        N_STEPS = std::stoi(argv[3]);

    auto exe = std::filesystem::path(argv[0]).parent_path();
    std::ostringstream fp;
    fp << exe.string() << "\\outputs\\lattice2D\\L_" << L
       << "_N_" << N << "_steps_" << N_STEPS << ".tsv";
    std::ofstream out(fp.str());
    out << "step\tlattice\n";

    run(out, L, N, N_STEPS);
    return 0;
}