## Non-integer Number of Colors!

An extension of the Puyo-Puyo model ([`puyopuyo/`](../puyopuyo/)). In the previous model, we have a "critical number" of colors somewhere between 5 and 6 (for 2D): with 5 colors, the system experiences negative drift and growth is halted at a fixed height, while with 6 colors, there's a finite positive drift and the system grows infinitely.

Now, we try to interpolate between them, extending our model to $N \in \mathbb{R}$.

As a simplest interpolation, consider the following:

- $N$ = 5.0 can be thought of as a 6 color system, with the probability of the 6th species being different, in this case, 0.
- $N$ = 6.0 is also a 6-color system, but with the 6th species being equiprobable as the other 5.

Thus, we can potentially interpolate by:
- $N$ = 5.5 can be a 6-color system, but the 6th species is only half as likely to be chosen as the other 5.

This can be generalized to any $N\in \mathbb{R}$. The floor (integer part) gives the number of colors, and the fractional part gives the relative weight of an "additional" color.
