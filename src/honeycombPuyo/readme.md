### Honeycomb Lattice

Is the 'power law' of around 3 universal? To test it, I'm running simulations similar to those in [probabilityPuyoPuyo](../probabilityPuyoPuyo/) here, except on a different lattice.

The main difference is that now, alternate rows are shifted vertically. In other words,

```
    A
B   A   C
B   x   C
D   x   E
D   F   E
    F   
``` 

Where `x` has 6 neighbors. Each cell is placed vertically, and only falls vertically, but avalanches are still possible.

There are two slight optimizations made to the code as compared to that in 
[probabilityPuyoPuyo](../probabilityPuyoPuyo/):

1. The `movedSites` is now no longer a boolean vector, but rather modified to use a buffer. In the case of very few eliminations, such as large `N`, there is a slight speedup.
2. The random colors and columns are not drawn every timestep, but instead predrawn and referred to.