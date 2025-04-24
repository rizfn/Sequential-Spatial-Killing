---
marp: true
theme: uncover
math: mathjax
paginate: true
_paginate: skip
backgroundColor: #111115
color: #ffeeee
style: |
        .columns {
          display: grid;
          grid-template-columns: repeat(2, minmax(0, 1fr));
          gap: 0.6rem;
        }
        h1, h2, h3, h4, h5, h6, strong {
          color: #fa99a5;
        }
        section::after {
          color: #998888;
          text-shadow: 0 0 1px #ff9999;
        }

---

![bg right:39.6% fit](../puyopuyo_apr16/fig/PuyoPuyoArt.jpg)

# PuyoPuyo

$\\$

Riz Fernando Noronha

---

<video src="../puyopuyo_apr16/fig/pom_pom_party.mp4" controls width="100%"></video>

---

### Simplify!

$\\$

- Drop one puyo at a time

- 2 is enough to trigger annihilation

- Choose a random column (no player input!)

---

<iframe width="100%" height="100%" src="https://rizfn.github.io/Sequential-Spatial-Killing/visualizations/puyopuyo/gravity" style="border: 1px solid #888888">

</iframe>

---

![bg right:25% fit](fig/1D_schematic.png)

# 1D Model

$\\$

Largest Cluster = 2

No falling!

Either add or delete each time step

Random walk!

---

Master equation:

![bg right:25% fit](fig/1D_schematic.png)

$$
\begin{align*}
P(h,\, t) =\, & P(h-1,\, t-1)\cdot\left(1-\frac{1}{N}\right) \\
&+ P(h+1,\, t-1) \cdot\frac{1}{N}
\end{align*}
$$

$\\$

We can derive:

$$
\frac{\partial P}{\partial t} = \underbrace{\left(\frac12\right)}_{\color{#fa99a5}{\boldsymbol{D}}} \frac{\partial^2 P}{\partial h^2} + \underbrace{\left(\frac2N-1\right)}_{\color{#fa99a5}\boldsymbol{-v_d}} \frac{\partial P}{\partial h} \\
$$

---

### 1D Mass/Height vs Time

![width:700px](fig/1D_gravityMassVsTime.png)

---

### 2D Mass, Height vs Time

![width:1150px](fig/2D_gravityMassVsTime_L_128.png)

---

### 3D Mass, Height vs Time

![width:1150px](fig/3D_gravityMassVsTime_L_24.png)

---

### Magic Number

$\\$

2D: **6** species can grow

3D: **9** species can grow

4D: **13** species can grow

---

![bg right:30% fit](fig/number_of_neighbours.png)

### Magic Number (?)

$\\$

Average ~2 neighbours
Disregard correlations
Assume no avalanches
$P(\text{same})$, $P(\text{different})$

0 drift at 5.4!
(Need to try 3D)


---

### Finite Size Scaling 

$\\$

![width:1150px](fig/2DMassFiniteSizeScaling.png)

---

### Avalanches (2D)

1. Number of **clusters** eliminated per 'added puyo'

2. Number of **puyos** eliminated per 'added puyo'

![width:1000px](fig/2D_avalanche_L_64_loglinear.png)

---

### Avalanches (2D)

![width:1150px](fig/2D_avalanche_L_64.png)

---

![bg right:36% fit](fig/sandpile_2e18_corner.png)

### Sandpile Models

$\\$

The original SOC model

Power-law with exponent $\color{#fa99a5}\boldsymbol{\tau\,}$**=1.5**

Difficult to push to $\tau\,$>2!

---

![width:1100px](fig/bipartiteMergingAnnihilation.png)

$\\$

Possibly linked to **bipartite networks**?

Merging donor+reciever is an "emission"

Emissions follow power law, $\tau\,$=3

---

![bg right:60% fit](fig/2D_roughness_L_128.png)

### Roughness (2D)

$\\$

Family–Vicsek scaling only for *growing systems*.

KPZ $\color{#fa99a5}\boldsymbol{\beta}\,$**=1/3** (same!)

---

### Why KPZ-like?

Local slopes must be restricted!

Taller towers are more vulnerable

![width:400px](fig/vulnerability_to_attacks.png)


---

## Future work

$\\$

Length scaling with roughness

Explain magic number in $d$ dimensions

Explore universal scaling

Measure slope distributions

Start from random lattice, then eliminate

Real (undirected 2D model)

Optimum number of chains

Abelian?

