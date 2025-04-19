# Sequential Spatial Killing

### Basic percolation (random initial conditions)

A simple model. We have $N$ different species on a lattice of length $L$.

We randomly fill the lattice with an equal number of each species. We then create a list of all the species, which is the order in which we kill.

We then sequentially kill the system, one species at a time. Say, we first start by killing green: then, from the edges of the graph, we feed in a "green" killer (eg: antibiotic, phage, etc) and it kills any greens or clusters of greens it can find, but is stopped by other colours. Then, we go on to the next color (red), and repeat the process.