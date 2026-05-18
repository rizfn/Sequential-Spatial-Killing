import numpy as np
import matplotlib.pyplot as plt
import os
import multiprocessing as mp

def run_single_sim(seed, L, T, dt, N, a, c, dx, noise_strength, tail_fraction=0.1):
    """
    Integrates the PDE for a single realization:
    dh/dt = (a^2 p^2 / 2) \\nabla^2 h - (c/N) (\\nabla h)^2 - a p^2 |\\nabla h| + \\eta(x, t)
    where p = 1/N
    """
    np.random.seed(seed)
    p = 1.0 / N
    
    # Coefficients
    c_diff = (a**2 * p**2) / 2.0
    c_nonlin2 = c / N
    c_nonlin1 = a * p**2
    
    # Initialize h
    h = np.zeros(L)
    
    # Track roughness over time
    steps = int(T / dt)
    roughness = np.zeros(steps)
    
    tail_steps = max(1, int(steps * tail_fraction))
    slopes_tail = []
    
    for i in range(steps):
        # Calculate spatial derivatives with periodic boundary conditions
        h_left = np.roll(h, 1)
        h_right = np.roll(h, -1)
        
        # Central difference for first derivative
        grad_h = (h_right - h_left) / (2.0 * dx)
        
        # Central difference for second derivative
        laplace_h = (h_right - 2.0 * h + h_left) / (dx**2)
        
        # Deterministic part of dh/dt
        dhdt = (c_diff * laplace_h - 
                c_nonlin2 * grad_h**2 - 
                c_nonlin1 * np.abs(grad_h))
        
        # Noise term \eta(x,t) scaled by sqrt(dt) for Euler-Maruyama
        noise = noise_strength * np.random.normal(0, 1, L)
        
        # Update h
        h = h + dhdt * dt + noise * np.sqrt(dt)
        
        # Calculate roughness: standard deviation of height
        roughness[i] = np.std(h)
        
        if i >= steps - tail_steps:
            slopes_tail.append(grad_h)
            
    # Final state metrics
    final_h = h - np.mean(h) # Center for better visualization
    
    return roughness, final_h, np.concatenate(slopes_tail)


def run_ensemble(num_runs, L, T, dt, N, a, c, dx, noise_strength):
    seeds = np.random.randint(0, 1000000, num_runs)
    args_list = [(seeds[i], L, T, dt, N, a, c, dx, noise_strength) for i in range(num_runs)]
    
    with mp.Pool() as pool:
        results = pool.starmap(run_single_sim, args_list)
        
    avg_roughness = np.mean([r[0] for r in results], axis=0)
    rep_profile = results[0][1] # Take the first run as representative
    all_slopes = np.concatenate([r[2] for r in results])
    
    return avg_roughness, rep_profile, all_slopes


def main():
    L = 256
    T = 1000
    dt = 0.001
    a = 1.0
    c = 1.0
    dx = 1.0
    noise_strength = 0.1
    num_runs = 8 # Number of parallel simulations per N
    
    Ns = [1, 2, 3, 4]
    
    steps = int(T / dt)
    times = np.arange(steps) * dt
    
    results = {}
    
    for N in Ns:
        print(f"Running ensemble for N={N}...")
        results[N] = run_ensemble(num_runs, L, T, dt, N, a, c, dx, noise_strength)
    
    # Plotting
    plots_dir = 'src/continuumPuyo/plots'
    os.makedirs(plots_dir, exist_ok=True)
    params_str = f'L{L}_T{T}_dt{dt}_a{a}_c{c}_dx{dx}_ns{noise_strength}_runs{num_runs}'
    
    # 1. Roughness vs Time
    plt.figure()
    for N in Ns:
        plt.plot(times, results[N][0], label=f'N={N}')
        
    # Plot reference lines
    t_ref = times[times > 0]
    offset_idx = len(times) // 2
    ref_N = Ns[1] if len(Ns) > 1 else Ns[0] # Try to align with one of the middle N curves
    ref_val = results[ref_N][0][offset_idx]
    t_mid = times[offset_idx]
    
    if ref_val > 0:
        plt.plot(t_ref, ref_val * (t_ref/t_mid)**(1./3.), 'k--', label=r'$t^{1/3}$')
        plt.plot(t_ref, ref_val * (t_ref/t_mid)**(1./2.), 'k:', label=r'$t^{1/2}$')
        
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('Time (t)')
    plt.ylabel('Ensemble Average Roughness W(t)')
    plt.title('Average Roughness vs Time')
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(plots_dir, f'roughness_{params_str}.png'))
    plt.close()
    
    # 2. Final Profile
    plt.figure()
    for N in Ns:
        plt.plot(np.arange(L) * dx, results[N][1], alpha=0.8, label=f'N={N}')
    plt.xlabel('Position (x)')
    plt.ylabel('Height (h - <h>)')
    plt.title('Representative Final Profiles')
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(plots_dir, f'profile_{params_str}.png'))
    plt.close()
    
    # 3. Distribution of local slopes
    plt.figure()
    for N in Ns:
        plt.hist(results[N][2], bins=100, density=True, alpha=0.7, histtype='step', linewidth=2, label=f'N={N}')
    plt.yscale('log')
    plt.xlabel('Local Slope (grad h)')
    plt.ylabel('Density')
    plt.title('Distribution of Local Slopes (Tail 10% steps)')
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(plots_dir, f'slopes_{params_str}.png'))
    plt.close()


if __name__ == "__main__":
    main()
