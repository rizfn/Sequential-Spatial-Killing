import matplotlib.pyplot as plt
import numpy as np

def viz_2D():
    L = 128
    N_list = [2, 3]
    steps = 4096

    for N in N_list:
        file = f"src/edenPuyoPuyo/evaporation/outputs/survival/2D/L_{L}_N_{N}_steps_{steps}.tsv"
        sim, t_death = np.loadtxt(file, delimiter="\t", skiprows=1, unpack=True)
        t_death = t_death.astype(int)
        N_sims = len(t_death)

        times = np.arange(steps)
        # For each t, count how many sims have t_death > t or t_death == -1
        survivors = np.array([(t_death > t).sum() + (t_death == -1).sum() for t in times])
        survival_prob = survivors / N_sims

        plt.plot(times, survival_prob, label=f"N={N}")

    plt.xlabel("Time")
    plt.ylabel("Survival Probability")
    plt.xscale("log")
    plt.yscale("log")
    plt.grid()
    plt.legend()
    plt.tight_layout()
    plt.show()


def viz_3D():
    L = 16
    # N_list = [2, 3, 4, 5]
    # steps = 1024
    N_list = [2]
    steps = 1024 * 4


    for N in N_list:
        file = f"src/edenPuyoPuyo/evaporation/outputs/survival/3D/L_{L}_N_{N}_steps_{steps}.tsv"
        sim, t_death = np.loadtxt(file, delimiter="\t", skiprows=1, unpack=True)
        t_death = t_death.astype(int)
        N_sims = len(t_death)

        times = np.arange(steps)
        # For each t, count how many sims have t_death > t or t_death == -1
        survivors = np.array([(t_death > t).sum() + (t_death == -1).sum() for t in times])
        survival_prob = survivors / N_sims

        plt.plot(times, survival_prob, label=f"N={N}")

    plt.xlabel("Time")
    plt.ylabel("Survival Probability")
    plt.xscale("log")
    plt.yscale("log")
    plt.grid()
    plt.legend()
    plt.tight_layout()
    plt.show()



if __name__ == "__main__":
    # viz_2D()
    viz_3D()
