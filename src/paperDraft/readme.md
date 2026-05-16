## Some math for a continuum version

### No correlations, no avalanches

Let's say we look at 3 cases: a "top" point (a column with shorter neighbours), a "bottom" point (a column lower than its neighbours) and a flat point (equal height to its neighbours). Let the probability of matching be $p=1/N$.

- Top: The height grows if it lands on the top, and doesn't annihilate (probability $1-p$). It reduces if it annihilates the top (probability $p$) or, if it drops in an adjacent cell, the neighbour in the tall column ($p+p$). Thus, $v_\text{peak}=(1-p)-p-2p = 1-4p$
- Bottom: Neighbours don't matter. It'll grow if it lands on top, and doesn't eliminate any of its three neighbours (probability $1-3p$). It'll fall if it eliminates the bottom with probability $p$. Then, $v_\text{valley} = 1-3p-p = 1-4p$
- Flat: It'll only change the column's height if it eliminates (probability $p$) or doesn't ($p-1$), which gives $v_\text{flat} = (1-p)-p = 1-2p$

In addition, there are three more cases: A flat-lower, a flat-higher, and a lower-higher boundary:
- Flat-lower: Can be eliminated with probability $2p$, can grow with probability $(1-p)$, giving $v=1-3p$
- Flat-higher: Can be eliminated with probability $p$, can grow with probability $1-2p$. $v=1-3p$
- Lower-higher: Eliminated with probability $2p$, can grow with $(1-2p)$: $v=1-4p$

There is some symmetry here, where flat-lower = flat-higher, and peak=valley. This is not true for avalanches (as they allow a peak to fall successively based on it's height), but if they're rare, we can approximate without them.

Drift is maximized at $1-2p$ on a flat interface where $m=0$, and drops down to $1-4p$. If we treat the local slope $m=\nabla h$ as a variable, the we can discuss a velocity profile $V(m)$.

Based on earlier, $V(0) = 1-2p$, and $V(m) = V(-m)$. As such, we can taylor expand it as a parabola:

$$V(m) \approx V(0) + V^\prime(0) m + \frac{V^{\prime\prime}(0)}{2}m^2$$

As the velocity decreases as the slope increases (dropping from $1-2p$ to $1-4p$), the maxima is at 0, and so $V^{\prime\prime}(0)<0$. The decrease is linked to $p$, and so we can assume that $ \frac{V^{\prime\prime}(0)}{2} = c p$. Then,

$$ V(m) \approx (1-2p) - cpm^2 $$

However, there's also noise, and since $m=\nabla h$ we can write

$$\frac{\partial h}{\partial t} = V(\nabla h) + \eta(x, t)$$
$$ \frac{\partial h}{\partial t} = (1-2p) - cp(\nabla h)^2 + \eta(x, t) $$

As $p=1/N$, we get

$$ \frac{\partial h}{\partial t} = \left(1-\frac2N\right) - \frac{c}N\left(\nabla h\right)^2 + \eta(x, t) $$


---

### Add avalanches

Now, how do avalanches modify this? Suppose we have three columns, $a$, $b$, and $c$. Both $b$ and $c$ are tall, while $a$ is low, forming a deep valley on the left. If we set the valley height to 0, the relative heights can be expressed as $h_b$ and $h_c$.

When a puyo drops into the left valley, it lands at the bottom of $b$'s exposed left wall. With probability $p$, it annihilates that bottom block.
Because of gravity, the remaining $h_b - 1$ blocks in $b$ will fall by 1 unit. As it falls, it slides against column $c$. The size of the shared boundary between them during this slide is $\min(h_1-1, h_2)$.

Each pair of adjacent blocks during the slide has a probability $p$ of matching. In addition, if $h_b>1$ (otherwise there's no avalanche), there's an additional link between the new bottom of $h_b$, and what's below it. Thus, the expected number of secondary eliminations triggered by the avalanche is $l= p \times \min(h_b, h_c)$. This creates a cascade with tertiary eliminations, etc, but for now, we only consider secondary ones.

What does this mean for the continuum limit?

Column $a$ only loses height through avalanches when a puyo falls in an adjacent column that's lower than it. Let's think about $h(x)$, and consider the case where the slope $m=\nabla h > 0$. This means that $h(x-dx) < h(x) < h(x+dx)$. Assuming a puyo is dropping in either neighbour (as there's 1 puyo dropped everywhere each timestep, on average). On the right, it does nothing: on the left, it causes an annihilation on the wall with probability $p$, and then can trigger new matches. The number of matches is approximately $h(x+dx) - h(x) \approx \lvert \nabla h \rvert$. Thus, the additional annihilation term is $-p\lvert\nabla h\rvert$. If we consider the opposite ($m<0$), we get the same result: the system is still symmetric to slopes, i.e, $V(-m)=V(m)$.

We can work through it a bit more carefully, starting from the discrete case, and only looking at a single (secondary) avalanche. Again, column $i$ only loses height when a puyo falls into an adjacent column lower than it. The expected erosion rate $v_i$ is based on the exposed walls of it: there's a probability $p$ for the initial annihilation, and then a $p$ for subsequent erosions, so that

$$ v_{\text{erosion}}(i) = - p^2\left[ \max(0, h_i - h_{i-1}) + \max(0, h_i - h_{i+1}) \right] $$

$\max$ can be written as

$$ \max(0, z) = \frac{z+\lvert z \rvert}2 $$

Thus,

$$ v_{\text{erosion}}(i) = - \frac{p^2}2\left[ (h_i - h_{i-1}) + |h_i - h_{i-1}| + (h_i - h_{i+1}) + |h_i - h_{i+1}| \right] $$

In the continuum case, $h_i - h_{i-1} \approx \nabla h dx$, and $h_{i-1} - 2h_i + h_{i+1} \approx \nabla^2 h \cdot dx^2$. If we let $dx=a$, the lattice spacing scale where $x_\text{continuous}=i\cdot a$, then

$$ v_{\text{erosion}} \approx p^2\left[ -a \lvert\nabla h\rvert + \frac{a^2}2\nabla^2 h \right] $$


That can be added into our main PDE:

$$ \frac{\partial h}{\partial t} = \left(1-\frac2N\right) + \frac{a^2p^2}2\nabla^2h - \frac{c}N\left(\nabla h\right)^2 - ap^2\lvert\nabla h \rvert + \eta(x, t) $$

Note that the fully coarse-grained model looks at $\lim_{a\to0}$. In other words, by only considering the secondary avalanche, avalanches seem to be small-scale effects that disappear under renormalization.

### Avalanche cascades

A lot can happen after the second avalanche. It can cause the eliminated sites to propagate both horizontally as well as vertically. Horizontally, column $c$ could fall from some eliminations, triggering a match with column $d$. While it's a lower-probability event, the new bottom of column $d$ may also match with the site below, triggering it to avalanche deeper, where it can even effect $a$ again. Fully mapping out the possibilities and their probabilities is challenging.

That said, we can make a much simpler approximation. If we assume that the average number of sites that slip stays the same (a bold assumption!) we can write that

$$ \langle \text{Cascade Size} \rangle \approx 1 + p|\nabla h| + (p|\nabla h|)^2 + \ldots = \frac{1}{1 - p|\nabla h|} $$

If $p|\nabla h| <1$, the cascade sizes are finite. However, at $p|\nabla h|\approx1$, it diverges, signifying long-range correlations. As $p=1/N$, that implies that the local slopes are of the order $|\nabla h|\approx N$, which may be possible?