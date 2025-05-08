## Puyo-Puyo Periodic

The same dynamics as earlier. A system size of $L$, and $N$ different "colours" or species. One by one, drop a single "puyo" (cell) into a random column on the grid. If a cluster of same-coloured puyos is formed, it's removed and the puyos above fall down.

The main difference is we use **periodic boundaries** in L, meaning that `grid[0]` neighbours `grid[L-1]`. This removes the edge effects in the previous system ([`puyopuyo/cpp/`](../cpp/)): we had the two edges being more "protected" than interior columns, as they could only be attacked from one side, which caused roughness and the max local slope to grow unbounded.