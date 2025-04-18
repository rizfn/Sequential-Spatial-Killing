# Puyo Puyo

[See presentation here!](https://rizfn.github.io/Sequential-Spatial-Killing/presentations/puyopuyo_apr16)

$N$ different colors (species), and a system of size $L\times L$

## Ballistic Deposition

Each time-step, choose a random color, and drop it down a random column. It falls following ballistic deposition rules (i.e. it sticks to the first thing it hits, and can stick to other particles horizontally).

What's new is **annihilation**. If two species of the same color are adjacent, they annihilate each other (and any others of the same color in the cluster).

Annihilation can cause "floating" clusters, so those clusters drift downards until they stick via BD rules.

## Simple model (`gravity`)

Same as above, except things just fall vertically and don't "stick" like in ballistic deposition. You fall until you have someone below you.

The same occurs for floating clusters, each column independently falls under gravity.