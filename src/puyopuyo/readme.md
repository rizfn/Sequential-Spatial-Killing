# Puyo Puyo

[See presentation here!](https://rizfn.github.io/Sequential-Spatial-Killing/presentations/puyopuyo_apr24)

$N$ different colors (species), and a system of size $L\times L$

## Ballistic Deposition

Each time-step, choose a random color, and drop it down a random column. It falls following ballistic deposition rules (i.e. it sticks to the first thing it hits, and can stick to other particles horizontally).

What's new is **annihilation**. If two species of the same color are adjacent, they annihilate each other (and any others of the same color in the cluster).

Annihilation can cause "floating" clusters, so those clusters drift downards until they stick via BD rules.

## Simple model (`gravity`)

Same as above, except things just fall vertically and don't "stick" like in ballistic deposition. You fall until you have someone below you.

The same occurs for floating clusters, each column independently falls under gravity.

### 1D Gravity:

$N$ species, so the probability of a given color is $1/N$

Let $P(h,t)$ represent the probability of having height $h$ at time $t$

$$
    P(h,t) = P(h-1,t-1)\left(1-\frac{1}{N}\right) + P(h+1,t-1) \frac{1}{N} 
$$

Let $1/N$ be called $\alpha$
$$
    P(h,t) - P(h, t-1) = (1-\alpha)P(h-1,t-1) + \alpha P(h+1,t-1) - P(h,t-1) \\
    P(h,t) - P(h, t-1) = (1-\alpha)P(h-1,t-1) + \alpha P(h+1,t-1) - [(1-\alpha) P(h,t-1) + \alpha P(h, t-1)] \\
    P(h,t) - P(h, t-1) = (1-\alpha)[P(h-1,t-1) - P(h, t-1)] + \alpha [P(h+1,t-1)-P(h,t-1)] 
$$
Let $Q(h,t) = P(h, t) - P(h-1, t)$
$$
    P(h,t) - P(h, t-1) = (1-\alpha)[-Q(h,t-1)] + \alpha [Q(h+1,t-1)] \\
$$
Let $\beta = \frac12 - \alpha$
$$
    P(h,t) - P(h, t-1) = \left(\beta + \frac12\right)[-Q(h,t-1)] + \left(\frac12-\beta\right) [Q(h+1,t-1)] \\
    P(h,t) - P(h, t-1) = \left[-\beta Q(h,t-1) - \frac12 Q(h, t-1)\right] + \left[\frac12Q(h+1,t-1) - \beta Q(h+1,t-1)\right] \\
    P(h,t) - P(h, t-1) = \frac12\left(Q(h+1, t-1)-Q(h,t-1)\right) - \beta(Q(h+1,t-1)+Q(h,t-1)) \\
$$
Now, looking at them as derivatives:
$$
\frac{\partial P(h,t)}{\partial t} = \frac12 \frac{\partial Q(h,t)}{\partial h} - \beta \left[P(h+1, t-1) - P(h, t-1) + P(h, t-1) - P(h-1, t-1)\right] \\
\frac{\partial P(h,t)}{\partial t} = \frac12 \frac{\partial^2 P(h,t)}{\partial h^2} - \beta \left[P(h+1, t-1) - P(h-1, t-1)\right] \\
\frac{\partial P(h,t)}{\partial t} = \frac12 \frac{\partial^2 P(h,t)}{\partial h^2} - 2\beta \frac{\partial P(h,t)}{\partial h} \\
\frac{\partial P(h,t)}{\partial t} = \frac12 \frac{\partial^2 P(h,t)}{\partial h^2} - 2\left(\frac12 - \alpha\right) \frac{\partial P(h,t)}{\partial h} \\
\frac{\partial P(h,t)}{\partial t} = \frac12 \frac{\partial^2 P(h,t)}{\partial h^2} + \left(2\alpha-1\right) \frac{\partial P(h,t)}{\partial h} \\
\frac{\partial P(h,t)}{\partial t} = \frac12 \frac{\partial^2 P(h,t)}{\partial h^2} + \left(\frac2N-1\right) \frac{\partial P(h,t)}{\partial h} \\
$$

Master equation:

$$
\frac{\partial P}{\partial t} = \frac12 \frac{\partial^2 P}{\partial h^2} + \left(\frac2N-1\right) \frac{\partial P}{\partial h} \\
$$

---

#### Try 2D:

You have, on average 2 neighbours (one below, one to the side). The probability of them being the same is $1/N$. The probability of them being different is $1-1/N$.

If they're the same, you can grow by dropping a new species that's different, which has probability $1-1/N$. If they're different, the dropped cell must b different from both neighbours, so you have probability $1 - 2/N$

$$
 P_{\text{growth}} = \frac{1}{N} \cdot \left(1 - \frac{1}{N}\right) + \left(1 - \frac{1}{N}\right) \cdot \frac{N - 2}{N}
$$
Combine the terms: 
$$ P_{\text{growth}} = \frac{N - 1}{N^2} + \frac{(N - 1)(N - 2)}{N^2} $$
Factorize: 
$$ P_{\text{growth}} = \frac{N - 1}{N^2} \left[1 + (N - 2)\right] $$
Finally,
$$ P_{\text{growth}} = \left(\frac{N - 1}{N}\right)^2 $$


---

#### 3D:

On average, you have 3 neighbours (one below, and 2 on the same layer). Not exact, as it's biased lower... but anyway.

$$ P_\text{growth} = \left[\text{all same col}\right]\left[P(\text{different from 1 col})\right] + \left[\text{2 same col}\right]\left[P(\text{different from 2 col})\right] + \left[\text{all diff col}\right]\left[P(\text{different from 3 col})\right]  $$

$$ P_\text{growth} = \left[\frac1{N^2}\right]\left[1-\frac1N\right] + \left[\frac{3(N-1)}{N^2}\right]\left[1-\frac2N\right] + \left[\frac{(N-1)(N-2)}{N^2}\right]\left[1-\frac3N\right]  $$

$$ P_\text{growth} = \frac{N-1}{N^3} \big(1+3(N-2) + (N-2)(N-3) \big)$$