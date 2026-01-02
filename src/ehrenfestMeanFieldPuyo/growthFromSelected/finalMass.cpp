#include <random>
#include <vector>
#include <unordered_map>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <filesystem>
#include <numeric>

// #pragma GCC optimize("Ofast","inline","fast-math","unroll-loops","no-stack-protector")
#pragma GCC optimize("inline", "unroll-loops", "no-stack-protector")
#pragma GCC target("sse,sse2,sse3,ssse3,sse4,popcnt,abm,mmx,avx,avx2,tune=native", "f16c")

static auto _ = []()
{std::ios_base::sync_with_stdio(false);std::cin.tie(nullptr);std::cout.tie(nullptr);return 0; }();

std::random_device rd;
std::mt19937 gen(rd());

// Define constants
constexpr int DEFAULT_N_COLORS = 500;
constexpr int DEFAULT_K_SELECTIONS = 15;
constexpr int DEFAULT_N_STEPS = 10000;
constexpr int DEFAULT_INITIAL_TOTAL_BALLS = 1000;

// Efficient sampling without replacement using Fisher-Yates shuffle
void sampleWithoutReplacement(std::vector<int>& all_balls, int K, std::vector<int>& drawn_balls)
{
    int total = all_balls.size();
    drawn_balls.clear();
    
    for (int i = 0; i < K; ++i)
    {
        std::uniform_int_distribution<> dis(i, total - 1);
        int j = dis(gen);
        std::swap(all_balls[i], all_balls[j]);
        drawn_balls.push_back(all_balls[i]);
    }
}

int simulate_urn(int N_colors, int K_selections, int N_steps, int initial_total_balls)
{
    // Initialize urn with equal number of balls of each color
    int initial_balls_per_color = initial_total_balls / N_colors;
    std::vector<int> urn(N_colors, initial_balls_per_color);
    
    // Pre-allocate vectors to avoid repeated allocations
    std::vector<int> all_balls;
    all_balls.reserve(initial_total_balls * 2); // Reserve extra space for growth
    std::vector<int> drawn_balls;
    drawn_balls.reserve(K_selections);
    
    for (int step = 0; step < N_steps; ++step)
    {
        int total_balls = std::accumulate(urn.begin(), urn.end(), 0);
        
        // Can't draw K balls if we don't have enough
        if (total_balls < K_selections)
            break;
        
        // Create list of all balls by expanding urn counts
        all_balls.clear();
        for (int color = 0; color < N_colors; ++color)
        {
            for (int i = 0; i < urn[color]; ++i)
            {
                all_balls.push_back(color);
            }
        }
        
        // Randomly draw K balls without replacement
        sampleWithoutReplacement(all_balls, K_selections, drawn_balls);
        
        // Count occurrences of each color in drawn balls
        std::unordered_map<int, int> color_counts;
        for (int color : drawn_balls)
        {
            color_counts[color]++;
        }
        
        // Check if all drawn balls are unique (no duplicates)
        if (color_counts.size() == static_cast<size_t>(K_selections))
        {
            // All K balls are different colors - add one to a random color from drawn
            std::uniform_int_distribution<> dis(0, drawn_balls.size() - 1);
            int random_idx = dis(gen);
            int random_color = drawn_balls[random_idx];
            urn[random_color]++;
        }
        else
        {
            // There are duplicates - remove all balls of colors that appeared more than once
            for (const auto& [color, count] : color_counts)
            {
                if (count > 1)
                {
                    urn[color] -= count;
                }
            }
        }
        
        // Print progress every 1000 steps
        if (step % 1000 == 0)
        {
            // std::cout << "Progress: " << std::fixed << std::setprecision(2)
            //           << static_cast<double>(step) / N_steps * 100 << "%\r" << std::flush;
        }
    }
    
    // std::cout << "Progress: 100.00%  \n" << std::flush;
    
    // Return final total
    return std::accumulate(urn.begin(), urn.end(), 0);
}

int main(int argc, char *argv[])
{
    int N_colors = DEFAULT_N_COLORS;
    int K_selections = DEFAULT_K_SELECTIONS;
    int N_steps = DEFAULT_N_STEPS;
    int initial_total_balls = DEFAULT_INITIAL_TOTAL_BALLS;
    
    if (argc > 1)
        N_colors = std::stoi(argv[1]);
    if (argc > 2)
        K_selections = std::stoi(argv[2]);
    if (argc > 3)
        N_steps = std::stoi(argv[3]);
    if (argc > 4)
        initial_total_balls = std::stoi(argv[4]);
    
    int final_size = simulate_urn(N_colors, K_selections, N_steps, initial_total_balls);
    
    std::string exePath = argv[0];
    std::string exeDir = std::filesystem::path(exePath).parent_path().string();
    std::ostringstream filePathStream;
    filePathStream << exeDir << "/outputs/finalMass/N" << N_colors 
                   << "_K" << K_selections 
                   << "_steps" << N_steps 
                   << "_init" << initial_total_balls << ".txt";
    std::string filePath = filePathStream.str();
    
    // Create output directory if it doesn't exist
    std::filesystem::create_directories(exeDir + "/outputs/finalMass");
    
    std::ofstream file(filePath);
    file << final_size << "\n";
    file.close();
        
    return 0;
}
