## Multi-layered Langmuir-like adsorbption

Imagine we have a catalyst, which has an 'attractive force' which helps lone atoms bond to it. Unlike the Langmuir model, we set up a rule that once bound, an atom cannot leave until it forms a molecule (gas), at which point it can float away.

In addition, we allow for multiple layers: There is a possibility for the catalyst to 'attract' atoms even at a distance. However, if an atom is at a distance and it's possibile for it to come closer, the attractive force will pull it closer, until it bumps into another atom. All atoms are also assumed to be of similar size.

This can be modelled in way quite similar to [Puyo Puyo](../puyopuyo). We have $N$ colors, where each color represents a different atom, such as oxygen $\mathrm{O}$ or nitrogen $\mathrm{N}$. We also define a $N\times N$ matrix $J$, which defines reactivity across atoms. For instance, say we have 3 atoms of $\mathrm{O}$, $\mathrm{N}$, $\mathrm{C}$, we can come up with the following reactions to create gases:
- $\mathrm{O} + \mathrm{O} \rightarrow \mathrm{O}_2 \nearrow $
- $\mathrm{N} + \mathrm{N} \rightarrow \mathrm{N}_2 \nearrow $
- $\mathrm{O} + \mathrm{N} \rightarrow \mathrm{NO} \nearrow $
- $\mathrm{O} + \mathrm{C} \rightarrow \mathrm{CO} \nearrow $

These can be represented in the $J$ matrix as follows:

$$
J = \begin{bmatrix}
1 & 1 & 1 \\
1 & 1 & 0 \\
1 & 0 & 0
\end{bmatrix}
$$

Where $\mathrm{O}$, $\mathrm{N}$, $\mathrm{C}$ are indexed by `0`, `1` and `2` appropriately, and a $1$ in the $J$ matrix implies the presence of a reaction.

Now, we can set up the dynamics:

- Every timestep, drop a puyo
- Check from the $J$ matrix if it can react with any of it's neighbors
- If so, react and annihilate
- Have the others fall, and attempt to react all the fallen puyos
- Continue until no further annihilations, and then drop another puyo