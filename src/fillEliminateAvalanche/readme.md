## Sequential Filling & Elimination

A simple model!

We have a lattice, say $L\times L$ in 2D. Then, we "fill" it: we place a random species from a list of $N$ species in each site.

Next, we **eliminate**: We look through the lattice for connected clusters of the same species. If we find one, we kill all species in the cluster.

Afterwards, we simply fill all the empty sites with random species again, and repeat. Does the pattern ever coarsen and freeze? Or is it constantly stuck in the fill/eliminate cycle? It depends on $N$!


#### Models:

- `noGravity`: The simplest model, described above (periodic in both directions)
- `gravity`: Periodic only in x. We fill the lattice, then enter a loop where we eliminate, fall, eliminate, etc.