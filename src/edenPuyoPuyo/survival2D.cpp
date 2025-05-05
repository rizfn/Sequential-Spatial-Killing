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

int run_survival(int L, int N_SPECIES, int STEPS) {
    std::mt19937 gen(std::random_device{}());
    std::uniform_int_distribution<> dis_sp(1, N_SPECIES);

    std::vector<std::vector<int>> grid(L, std::vector<int>(L, 0));
    std::vector<int> boundary;
    std::vector<int> pos(L * L, -1);
    std::vector<char> isBoundary(L * L, 0);

    int rc = L / 2, cc = L / 2;
    grid[rc][cc] = dis_sp(gen);
    for (auto d : std::vector<std::pair<int, int>>{{1, 0}, {-1, 0}, {0, 1}, {0, -1}}) {
        int nr = rc + d.first, nc = cc + d.second;
        if (nr >= 0 && nr < L && nc >= 0 && nc < L) {
            int idx = nr * L + nc;
            isBoundary[idx] = 1;
            pos[idx] = boundary.size();
            boundary.push_back(idx);
        }
    }

    std::vector<std::pair<int, int>> cluster;
    int t_dead = -1;
    for (int t = 0; t < STEPS; ++t) {
        if (boundary.empty())
            break;
        std::uniform_int_distribution<> dis_b(0, (int)boundary.size() - 1);
        int bi = dis_b(gen);
        int idx = boundary[bi];
        int r = idx / L, c = idx % L;
        int last = boundary.back();
        std::swap(boundary[bi], boundary.back());
        pos[last] = bi;
        boundary.pop_back();
        pos[idx] = -1;
        isBoundary[idx] = 0;

        grid[r][c] = dis_sp(gen);

        for (auto &d : std::vector<std::pair<int, int>>{{1, 0}, {-1, 0}, {0, 1}, {0, -1}}) {
            int nr = r + d.first, nc = c + d.second;
            if (nr >= 0 && nr < L && nc >= 0 && nc < L && grid[nr][nc] == 0) {
                int nidx = nr * L + nc;
                if (!isBoundary[nidx]) {
                    isBoundary[nidx] = 1;
                    pos[nidx] = boundary.size();
                    boundary.push_back(nidx);
                }
            }
        }

        floodAnnihilate(grid, r, c, cluster);
        for (auto &p : cluster) {
            int pr = p.first, pc = p.second, pidx = pr * L + pc;
            if (hasFilledNeighbor(grid, pr, pc) && !isBoundary[pidx]) {
                isBoundary[pidx] = 1;
                pos[pidx] = boundary.size();
                boundary.push_back(pidx);
            }
            for (auto &d : std::vector<std::pair<int, int>>{{1, 0}, {-1, 0}, {0, 1}, {0, -1}}) {
                int nr = pr + d.first, nc = pc + d.second;
                if (nr >= 0 && nr < L && nc >= 0 && nc < L && grid[nr][nc] == 0) {
                    int nidx = nr * L + nc;
                    if (isBoundary[nidx] && !hasFilledNeighbor(grid, nr, nc)) {
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
        t_dead = t;
    }

    for (int i = 0; i < L; ++i) {
        for (int j = 0; j < L; ++j) {
            if (grid[i][j] > 0) {
                t_dead = -1;
                break;
            }
        }
    }
    
    return t_dead;
}

int main(int argc, char *argv[])
{
    int L = 128, N = 3, N_STEPS = 1024*4, N_sims = 1000;
    if (argc > 1)
        L = std::stoi(argv[1]);
    if (argc > 2)
        N = std::stoi(argv[2]);
    if (argc > 3)
        N_STEPS = std::stoi(argv[3]);
    if (argc > 4)
        N_sims = std::stoi(argv[4]);

    auto exe = std::filesystem::path(argv[0]).parent_path();
    std::ostringstream fp;
    fp << exe.string() << "\\outputs\\survival\\2D\\L_" << L
       << "_N_" << N << "_steps_" << N_STEPS << ".tsv";
    std::ofstream out(fp.str());
    out << "sim\tt_dead\n";

    for (int sim = 0; sim < N_sims; ++sim) {
        int t_dead = run_survival(L, N, N_STEPS);
        out << sim << "\t" << t_dead << "\n";
        std::cout << "\rSim " << (sim + 1) << "/" << N_sims << std::flush;
    }
    std::cout << "\nDone.\n";
    return 0;
}
