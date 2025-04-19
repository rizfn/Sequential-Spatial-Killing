## Percolation Power-Law Populations:

$N$ species, on a system size of $L$. Now, each species doesn't have an equal population, but rather, abundances are power-law distributed.

Note that there's a key distinction in the definition of the exponent.

If you have many species, each with their own abundances, you can plot the abundance on the x-axis, and the number of species having that particular abundance on y. Then, the slope of that line (on a log-log plot) can be $-\tau$.

However, there's also the exponent for the zipf-like distribution $-\alpha$: if you rank the species in order by the number of organisms each has, and plot that, you get a different exponent.

The Zipf exponent $\alpha$ (abundances vs rank) can be linked to the previously discussed one (Number of species vs abundance) by the following relationship:

$$\alpha = \frac1\tau + 1$$

This is due to the axes first being flipped (so $\tau$ goes to $1/\tau$), and then the rank plot being a cumulative distribution function (CDF) of the original distribution (so you need to add 1).

Or in other words, if the original has slope $-\tau$, the rank plot has slope $-1/\tau - 1$.

This can be seen by this script: The cumulative effect can be seen without dividing by the bin widths.

```python
import numpy as np
import matplotlib.pyplot as plt

species = np.arange(1, 101)
tau = 1.5
weights = species ** (-1/tau - 1)
weights /= np.sum(weights)  # Normalize to sum to 1

samples = np.random.choice(species, size=10_000_000, p=weights)

species, pops = np.unique(samples, return_counts=True)

plt.plot(species, pops, linestyle='', marker='x')
plt.plot(species, pops[0]*species**(-1/tau - 1), 'r-', label='Zipf distribution fit')
plt.yscale('log')
plt.xscale('log')
plt.show()

# log_pops = np.log10(pops)
# hist, edges = np.histogram(log_pops, bins=20)
# log_hist = np.log10(hist)

# plt.plot(edges[1:], log_hist, 'o', label='Log-log histogram of populations')
# plt.plot(edges[1:], edges[1:]*(-tau) + (log_hist[0]-edges[1]*-tau), 'r-', label='Zipf distribution fit')

hist, bin_edges = np.histogram(pops, bins=np.geomspace(pops.min(), pops.max(), 20))
bin_widths = bin_edges[1:] - bin_edges[:-1]
hist = hist / bin_widths
plt.plot(bin_edges[1:], hist, 'o', label='Log-log histogram of populations')
plt.plot(bin_edges[1:], bin_edges[1:]**(-tau) * hist[0], 'r-', label='Zipf distribution fit')
plt.yscale('log')
plt.xscale('log')
plt.show()
```