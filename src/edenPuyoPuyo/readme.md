## Eden Growth Puyo Puyo

We have similar dynamics as the puyo-puyo model, but on a real $L\times L$ 2D grid, instead of 1+1D with a directed axis.

We have a 2D grid, and seed a single puyo in the center. Then, every time step, it "grows" as an Eden growth model: it "replicates" into one of it's neighbours. This can be done in two ways:

1. **Bond Growth**: Each "edge" between a cell and it's empty neighbour has an equal probability to be filled
2. **Site Growth**: Each empty neighbour has an equal probability to be filled

The two are not equivalent (ADD FIGURE).

The main difference from the Eden model is that we add **death** into the system. Having cells die randomly causes the system to be mapped to the Directed Percolation universality class.

Instead, we say that every time a cell replicates, it's "offspring" is chosen to be (potentially) a different species, selected from a pool of $N$ total species. If, at any point, a cluster of cells of the same species is formed, it is removed.

Removal could leave a floating cluster of cells. We need to decide on how to deal with them:
1. Attract towards the center of gravity?
2. Attract towards the cluster?
3. Attract towards the nearest filled state?
4. Only move horizontally and vertically? (may not work?)

The problem with all the above options is they can be order-dependent: our final configuration depends on the order in which floating clusters "fall".

This leaves another, simpler option:

5. Evaporate (or remove) the floating clusters

However, this restricts avalanches in the system, and we will not be able to replicate the power-law avalanche distribution of the [Puyo Puyo model](../puyopuyo/periodicCpp/).
