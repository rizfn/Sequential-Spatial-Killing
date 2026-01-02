# Growth from Selected Strategy - Analytical Boundary

## Derivation of Critical Boundary

Probability all $K$ balls are different colors:
$$P(\text{all different}) = \frac{N!}{(N-K)! N^K}$$

Large $N$ approximation [(birthday paradox)](https://en.wikipedia.org/wiki/Birthday_problem#Approximations):
$$P(\text{all different}) = \prod_{i=0}^{K-1} \left(1 - \frac{i}{N}\right) \approx \exp\left(-\sum_{i=0}^{K-1} \frac{i}{N}\right) = \exp\left(-\frac{K(K-1)}{2N}\right)$$

Expected growth per step:
$$\langle\Delta M_{\text{growth}}\rangle = (+1) \cdot \exp\left(-\frac{K(K-1)}{2N}\right)$$

Probability of collision (birthday paradox + taylor expansion of $e^x$ for small $x$):
$$P(\text{collision}) = 1 - \exp\left(-\frac{K(K-1)}{2N}\right) \approx \frac{K(K-1)}{2N}$$

For small collision probability, most collisions are pairs. Expected loss per collision $\approx$ 2 balls removed.

Expected loss per step:
$$\langle\Delta M_{\text{loss}}\rangle = -2 \cdot P(\text{collision}) \approx -2 \cdot \frac{K(K-1)}{2N} = -\frac{K(K-1)}{N}$$

Critical boundary (growth = loss):
$$\exp\left(-\frac{K(K-1)}{2N}\right) = \frac{K(K-1)}{N}$$

Dilute limit ($K \ll \sqrt{N}$), using $e^{-x} \approx 1 - x$:
$$1 - \frac{K(K-1)}{2N} = \frac{K(K-1)}{N}$$

$$1 = \frac{K(K-1)}{N} + \frac{K(K-1)}{2N}$$

$$1 = \frac{3K(K-1)}{2N}$$

$$N_c = \frac{3K(K-1)}{2}$$