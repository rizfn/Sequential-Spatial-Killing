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

void floodAnnihilate3D(std::vector<std::vector<std::vector<int>>> &grid,
                       int x0, int y0, int z0,
                       std::vector<std::tuple<int, int, int>> &cluster)
{
    int L = grid.size();
    int species = grid[x0][y0][z0];
    cluster.clear();
    std::vector<std::vector<std::vector<bool>>> vis(L, std::vector<std::vector<bool>>(L, std::vector<bool>(L, false)));
    std::queue<std::tuple<int, int, int>> q;
    q.push({x0, y0, z0});
    vis[x0][y0][z0] = true;
    while (!q.empty())
    {
        auto [x, y, z] = q.front();
        q.pop();
        cluster.emplace_back(x, y, z);
        for (auto d : std::vector<std::tuple<int, int, int>>{
                 {1, 0, 0}, {-1, 0, 0},
                 {0, 1, 0}, {0, -1, 0},
                 {0, 0, 1}, {0, 0, -1}})
        {
            int nx = x + std::get<0>(d), ny = y + std::get<1>(d), nz = z + std::get<2>(d);
            if (nx >= 0 && nx < L && ny >= 0 && ny < L && nz >= 0 && nz < L &&
                !vis[nx][ny][nz] && grid[nx][ny][nz] == species)
            {
                vis[nx][ny][nz] = true;
                q.push({nx, ny, nz});
            }
        }
    }
    if (cluster.size() > 1)
    {
        for (auto &p : cluster)
        {
            auto [x, y, z] = p;
            grid[x][y][z] = 0;
        }
    }
    else
    {
        cluster.clear();
    }
}

inline bool hasFilledNeighbor3D(const std::vector<std::vector<std::vector<int>>> &grid, int x, int y, int z)
{
    int L = grid.size();
    static const int D[6][3] = {
        {1, 0, 0}, {-1, 0, 0},
        {0, 1, 0}, {0, -1, 0},
        {0, 0, 1}, {0, 0, -1}};
    for (auto &d : D)
    {
        int nx = x + d[0], ny = y + d[1], nz = z + d[2];
        if (nx >= 0 && nx < L && ny >= 0 && ny < L && nz >= 0 && nz < L && grid[nx][ny][nz] > 0)
            return true;
    }
    return false;
}

int run_survival3D(int L, int N_SPECIES, int STEPS)
{
    std::mt19937 gen(std::random_device{}());
    std::uniform_int_distribution<> dis_sp(1, N_SPECIES);

    std::vector<std::vector<std::vector<int>>> grid(L, std::vector<std::vector<int>>(L, std::vector<int>(L, 0)));
    std::vector<int> boundary;
    std::vector<int> pos(L * L * L, -1);
    std::vector<char> isBoundary(L * L * L, 0);

    int cx = L / 2, cy = L / 2, cz = L / 2;
    grid[cx][cy][cz] = dis_sp(gen);
    for (auto d : std::vector<std::tuple<int, int, int>>{
             {1, 0, 0}, {-1, 0, 0},
             {0, 1, 0}, {0, -1, 0},
             {0, 0, 1}, {0, 0, -1}})
    {
        int nx = cx + std::get<0>(d), ny = cy + std::get<1>(d), nz = cz + std::get<2>(d);
        if (nx >= 0 && nx < L && ny >= 0 && ny < L && nz >= 0 && nz < L)
        {
            int idx = nx * L * L + ny * L + nz;
            isBoundary[idx] = 1;
            pos[idx] = boundary.size();
            boundary.push_back(idx);
        }
    }

    std::vector<std::tuple<int, int, int>> cluster;
    int t_dead = -1;
    for (int t = 0; t < STEPS; ++t)
    {
        if (boundary.empty())
            break;
        std::uniform_int_distribution<> dis_b(0, (int)boundary.size() - 1);
        int bi = dis_b(gen);
        int idx = boundary[bi];
        int x = idx / (L * L), y = (idx / L) % L, z = idx % L;
        int last = boundary.back();
        std::swap(boundary[bi], boundary.back());
        pos[last] = bi;
        boundary.pop_back();
        pos[idx] = -1;
        isBoundary[idx] = 0;

        grid[x][y][z] = dis_sp(gen);

        for (auto &d : std::vector<std::tuple<int, int, int>>{
                 {1, 0, 0}, {-1, 0, 0},
                 {0, 1, 0}, {0, -1, 0},
                 {0, 0, 1}, {0, 0, -1}})
        {
            int nx = x + std::get<0>(d), ny = y + std::get<1>(d), nz = z + std::get<2>(d);
            if (nx >= 0 && nx < L && ny >= 0 && ny < L && nz >= 0 && nz < L && grid[nx][ny][nz] == 0)
            {
                int nidx = nx * L * L + ny * L + nz;
                if (!isBoundary[nidx])
                {
                    isBoundary[nidx] = 1;
                    pos[nidx] = boundary.size();
                    boundary.push_back(nidx);
                }
            }
        }

        floodAnnihilate3D(grid, x, y, z, cluster);
        for (auto &p : cluster)
        {
            int px, py, pz;
            std::tie(px, py, pz) = p;
            int pidx = px * L * L + py * L + pz;
            if (hasFilledNeighbor3D(grid, px, py, pz) && !isBoundary[pidx])
            {
                isBoundary[pidx] = 1;
                pos[pidx] = boundary.size();
                boundary.push_back(pidx);
            }
            for (auto &d : std::vector<std::tuple<int, int, int>>{
                     {1, 0, 0}, {-1, 0, 0},
                     {0, 1, 0}, {0, -1, 0},
                     {0, 0, 1}, {0, 0, -1}})
            {
                int nx = px + std::get<0>(d), ny = py + std::get<1>(d), nz = pz + std::get<2>(d);
                if (nx >= 0 && nx < L && ny >= 0 && ny < L && nz >= 0 && nz < L && grid[nx][ny][nz] == 0)
                {
                    int nidx = nx * L * L + ny * L + nz;
                    if (isBoundary[nidx] && !hasFilledNeighbor3D(grid, nx, ny, nz))
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
        t_dead = t;
    }

    for (int i = 0; i < L; ++i)
        for (int j = 0; j < L; ++j)
            for (int k = 0; k < L; ++k)
                if (grid[i][j][k] > 0)
                {
                    t_dead = -1;
                    break;
                }

    return t_dead;
}

int main(int argc, char *argv[])
{
    int L = 16, N = 3, N_STEPS = 1024, N_sims = 1000;
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
    fp << exe.string() << "\\outputs\\survival\\3D\\L_" << L
       << "_N_" << N << "_steps_" << N_STEPS << ".tsv";
    std::ofstream out(fp.str());
    out << "sim\tt_dead\n";

    for (int sim = 0; sim < N_sims; ++sim)
    {
        int t_dead = run_survival3D(L, N, N_STEPS);
        out << sim << "\t" << t_dead << "\n";
        std::cout << "\rSim " << (sim + 1) << "/" << N_sims << std::flush;
    }
    std::cout << "\nDone.\n";
    return 0;
}