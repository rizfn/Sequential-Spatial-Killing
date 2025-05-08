## PuyoPuyo Model

C++ implementation of the PuyoPuyo model. We ignore the ballistic deposition case, as this seems more interesting.

We have a system size of $L$, and $N$ different "colours" or species. One by one, drop a single "puyo" (cell) into a random column on the grid. If a cluster of same-coloured puyos is formed, it's removed and the puyos above fall down.

