## A single-player Ehrenfest game

Imagine you have an urn, with multiple balls inside. Balls can be of different colours (set by the parameter `N`). Then, every time step, we pull out `k` balls from the urn, randomly. If there are duplicates of the same color, we remove them, and place the remaining balls back inside. If there are no duplicates, we place a new ball inside the urn.

Essentially, we're "adding" mass if there are no duplicates, and "removing" it if there are. While the drawing itself is random, the player is allowed to use different strategies on **what colour ball** they add to the urn.

### No strategy:

When it's time to add, the agent adds a ball of a random color.

The growth rate is $P(\text{all different})$

The death rate is $P(\text{duplicates})×E[\text{balls removed∣duplicates}]$

You're right! Let me refine the analysis to account for the expected number of balls lost.

**More precise analysis:**

When drawing K balls from N colors, the expected change per step is:

$$\Delta E = P(\text{all different}) \times (+1) - P(\text{duplicates}) \times E[\text{balls removed} | \text{duplicates}]$$

For small K/N, when duplicates occur, the expected number of balls removed is approximately **2 balls** (the most common scenario is a single pair of duplicates). However, for larger K, we can get more duplicates.

A better approximation for the expected loss is $\approx 2K/\sqrt{N}$ (the expected number of colliding balls scales with the collision probability).

Setting the transition at $\Delta E = 0$:

$$P(\text{all different}) \approx 2 \cdot P(\text{duplicates})$$
$$e^{-K^2/(2N)} \approx 2(1 - e^{-K^2/(2N)})$$

Solving this gives:
$$e^{-K^2/(2N)} = \frac{2}{3}$$
$$\frac{K^2}{2N} = \ln\left(\frac{3}{2}\right) \approx 0.405$$

Therefore:

$$\boxed{K_{\text{critical}} \approx 0.90\sqrt{N}}$$



