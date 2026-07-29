# What variable linearises $v(N)$? (Not an entropy.)

**Short answer.** The right variable is the **deposit-meets-pile collision
probability** $\sum_i p_i\rho_i$ — the chance that a landing block matches a
block already in the pile. It reduces to exactly $1/N$ at integer $N$, and it
straightens the whole curve. **It is not an information entropy**: no Rényi/Hill
order works, and the reason is structural, not a matter of picking the right $q$.

**Can $\rho$ be derived rather than measured? Yes.** Removal is *deposition-driven*,
$r_i = k\,p_i\rho_i$ (measured directly, $R^2=0.998$), which with the exact
steady-state balance gives $\rho_i = p_i/(v+kp_i)$ — the measured composition to
**1–2%**, impurity enrichment included. So the pile probabilities *are* derivable
from the drop probabilities (§5).

**But that still does not give you $v$**, because $v=1-\sum_i r_i$ is a near-total
cancellation ($\sum r_i=0.94$ at $N=6$): a 2.6% error in $r$ becomes a **42%**
error in $v$. Deriving $\rho$ and measuring $\rho$ are therefore very different
things here.

## 1. The thing to explain: $v(N)$ is a staircase

$L=1024$, $5\times10^4$ steps, 16 sims, warmup $10^4$, zero ceiling hits.
$v=\mathrm{d}(M/L)/\mathrm{d}t$ ($M/L$ is exactly the mean height — gravity
leaves no holes, so mass and height are the same observable, verified to $10^{-7}$
against `avalancheDist`'s independent lattice scan).

Inside each interval $v$ jumps almost discontinuously at $f=0^+$ and then flattens:

| $N$ | 6.0 | 6.1 | 6.2 | 6.5 | 6.9 | 7.0 |
|---|---|---|---|---|---|---|
| $v$ | 0.0582 | 0.1404 | 0.1671 | 0.1991 | 0.2096 | 0.2101 |
| $\Delta v$ | — | **+0.0822** | +0.0267 | +0.0062 | +0.0004 | +0.0005 |

The first tenth of the interval carries **more than half** the total rise; the
last tenth carries 0.3%. A factor of $\sim160$ in slope across one interval.
That is why $v$ vs non-integer $N$ is not a straight line.

## 2. The benchmark, and the answer

The figure of merit is the RMS residual of a straight-line fit. The target is
**not zero** — it is how straight $v$ vs $1/N$ is over integers, which is the
best the model itself offers:

| variable | RMS residual | $\times$ noise floor |
|---|---|---|
| $1/N$, **integer $N$ only** (benchmark) | 0.00259 | 13 |
| $\sum_i p_i\rho_i$, **all 46 $N$** | **0.00235** | **12** |
| $\sum_i p_i^2$ (deposit–deposit), all $N$ | 0.01095 | 55 |
| $\sum_i \rho_i^2$ (pile–pile), all $N$ | 0.00569 | 28 |
| $1/{}^1\!D$ (Shannon), all $N$ | 0.00544 | 27 |
| $1/N$, all $N$ | 0.01650 | 82 |

Noise floor $\approx2\times10^{-4}$ (spread across sims).

**$\sum_i p_i\rho_i$ over the entire range is as linear as $1/N$ is over the
integers alone.** It also reduces to $1/N$ at integer $N$ to machine precision
(checked: $|\sum p_i\rho_i - 1/N| \le 1.4\times10^{-17}$ at $N=6..10$), so it is
the correct generalisation rather than a competing fit. Best line:

$$v = 1.082 - 6.120\sum_i p_i\rho_i$$

Note $v$ is linear in the collision probability only to $\sim1\%$ absolute — the
benchmark itself sits 13$\times$ above the noise floor and has *systematic*
curvature (integer residuals $-0.0029, +0.0028, +0.0027, +0.0004, -0.0031$). So
"$v$ is linear in $1/N$ at integer $N$" is a good approximation, **not an
identity**. Don't over-read either fit.

## 3. Is it information entropy? No.

Scanning every Rényi order $q$ (Hill number ${}^q\!D=\exp H_q$, $q=1$ Shannon,
$q=2$ collision/Simpson):

- best order is $q\approx0.70$, RMS $0.00388$ — **still 19$\times$ noise**, and
  worse than $\sum p_i\rho_i$;
- $q=0.70$ is not a principled value. It is a fudge, not a principle;
- the two *natural* orders both do badly: Shannon 27$\times$, collision 55$\times$.

This was worth testing because it looks so promising: for uniform weights
$\sum_i p_i^2 = 1/N$ exactly, so the collision entropy $H_2$ *does* reproduce
$1/N$ at integer $N$. It simply fails off the integers, and §4 says why.

## 4. Hand-waving argument

**Why $1/N$ works at integer $N$.** A deposited block adds to the pile unless it
annihilates, and it annihilates iff it matches a neighbour. Its neighbours are
*pile* blocks. So the elementary event is: a block drawn from the deposition
distribution $p$ meets a block drawn from the pile distribution $\rho$, and

$$P(\text{match}) = \sum_i p_i \rho_i$$

Mass balance (`../avalancheScaling`) gives $v = 1 - \langle s\rangle f_{\rm active}$
with $f_{\rm active}\propto P(\text{match})$, so to leading order
$v \approx 1 - c\sum_i p_i\rho_i$. At integer $N$ the $N$ species are related by
a permutation symmetry, which forces $\rho = p = $ uniform, hence

$$\sum_i p_i \rho_i = N\cdot\tfrac1N\cdot\tfrac1N = \tfrac1N \qquad\Rightarrow\qquad v = 1 - c/N$$

So $1/N$ was never fundamental. It is the collision probability *in disguise*,
and only because symmetry makes the two coincide.

**Why it breaks off the integers.** At $N=n+f$ the symmetry is gone and the rare
species is **enriched**: measured $\rho_{\rm imp}/(f/N)$ across $N\in[6,7]$ is
**4.03$\times$** at $f=0.1$, falling to 1.00$\times$ at $f=1$ (where the impurity is
an ordinary species and symmetry is restored). It dilutes the majority species,
suppresses matching, and $v$ jumps.

The *mechanism* of the enrichment is §5's closure, and it is **not** a
$\sqrt{f}$ cusp:

$$\rho_{\rm imp} = \frac{p_{\rm imp}}{v + k\,p_{\rm imp}},\qquad p_{\rm imp}=f/N$$

a smooth crossover that is linear in $f$ below $f^* = Nv/k$ and saturates at
$1/k$ above it. $f^*$ is small precisely when $v$ is small, so near $N_c$ the rise
becomes arbitrarily sharp *without any singularity*. Since $f^*$ grows with $N$
(0.156, 0.329, 0.513, 0.699 at $N=6.1,7.1,8.1,9.1$), more of each interval sits in
the linear regime as $N$ grows, so a power-law fit must return an exponent
drifting toward 1 — measured $b = 0.329, 0.467, 0.564, 0.632$ across
$[6,7],[7,8],[8,9],[9,10]$. **RETRACTED:** the earlier $f\sim\rho^2 \Rightarrow
\rho\sim\sqrt f$ pair-balance argument (also in `../avalancheScaling`) is directly
refuted — see §5. It gets the *sharpness* right for the wrong reason, and it
cannot explain the drifting $b$, which the crossover does.

**Why no entropy can work.** Every entropy — Shannon, Rényi, Hill, any order — is
a functional of the **deposition** distribution $p$. It measures the diversity of
what you *put in*. But the dynamics collides against what is *already there*, and
self-poisoning makes those two differ by a factor $1/(v+kp_i)$ that depends on
$v$ — i.e. on the *answer*. An entropy is a fixed functional of $p$ with no
knowledge of $v$, so it cannot produce a family of curves whose shape changes with
the very quantity being predicted. The failure is structural: it is not that we
picked the wrong $q$.

## 5. What *is* $\sum_i p_i\rho_i$, and can it come from $p$?

$p_i$ is the probability of **dropping** species $i$; $\rho_i$ is the fraction of
the **pile** that is species $i$. So $\sum_i p_i\rho_i$ is the probability that a
freshly dropped block and a randomly chosen pile block carry the same species —
the chance the newcomer finds something to annihilate against. At $N=6.5$:

| | $p_i$ (drop) | $\rho_i$ (pile) | ratio |
|---|---|---|---|
| each major | 0.1538 | 0.1465 | 0.95 |
| impurity | 0.0769 | 0.1211 | **1.57** |

The impurity is enriched, the majors depleted, so $\sum p_i^2 = 0.1479$
overestimates the true $\sum p_i\rho_i = 0.1445$ by 2.4%.

### There is an exact constraint linking them

In steady state the composition is fixed while the pile grows at $v$ per
deposition. Per species: added $-$ removed $=$ what the growing pile carries,

$$\boxed{\;p_i - r_i = v\,\rho_i\;}$$

with $r_i$ the removal rate of species $i$. Summing recovers
$v = 1-\sum_i r_i$, consistent with the mass balance of `../avalancheScaling`.
Exact, but not closed: $r_i$ needs the dynamics.

### Closing it: removal is deposition-driven

Nothing happens in this model until a block is **dropped**. Species $i$ is removed
when species $i$ is deposited (rate $p_i$) *and* it lands next to a match
(probability $\propto\rho_i$). So the closure is $r_i = k\,p_i\rho_i$ — **not**
$k\rho_i^2$ (two pile blocks meeting), which is what an equilibrium pair-reaction
intuition suggests and which is wrong here because pile blocks do not move on
their own.

Recording $r_i$ per species and fitting $r_i = k\,x$ through the origin over all
46 $N$ and every species (322 points) decides it directly:

| closure for $r_i$ | $k$ | $R^2$ | rel. RMS |
|---|---|---|---|
| $\rho_i^2$ (pair meeting) | 5.59 | 0.965 | 27% |
| **$p_i\rho_i$ (deposition-driven)** | **5.51** | **0.9977** | **2.6%** |
| $p_i^2$ | 5.37 | 0.980 | 13% |

With the exact balance this gives, with **one** constant,

$$\boxed{\;\rho_i = \frac{p_i}{v + k\,p_i}\;}$$

and it reproduces the pile composition from the drop probabilities:

| $N$ | $\rho_{\rm imp}$ meas | pred | err | $\rho_{\rm maj}$ meas | pred | err |
|---|---|---|---|---|---|---|
| 6.1 | 0.0660 | 0.0700 | +6.1% | 0.1557 | 0.1550 | −0.4% |
| 6.5 | 0.1211 | 0.1228 | +1.4% | 0.1464 | 0.1462 | −0.2% |
| 8.5 | 0.0834 | 0.0832 | −0.3% | 0.1146 | 0.1146 | +0.0% |

**So yes — the probabilities inside the system are derivable from the drop
probabilities.** The enrichment, the staircase and the drifting $b$ (§4) all follow.

### …and yet $v$ still does not follow, because of a cancellation

Imposing $\sum_i\rho_i=1$ closes the system and predicts $v$ from $p$ with no
simulation at all. It is poor: RMS $0.01403$ (**70$\times$ noise**), $N_c$
predicted $k=5.50$ vs measured **5.0765** (8.3% off).

This is not a failure of the closure — it is amplification. $v$ is what
*survives a near-total cancellation*:

| $N$ | $\sum_i r_i$ | $v = 1-\sum r_i$ | 2.6% error in $r$ $\to$ | as % of $v$ |
|---|---|---|---|---|
| 6.0 | 0.9418 | 0.0582 | 0.0245 | **42%** |
| 7.0 | 0.7898 | 0.2102 | 0.0205 | 10% |
| 10.0 | 0.5326 | 0.4674 | 0.0138 | 3% |

At $N=6$, 94% of deposited mass is removed again and $v$ is the 6% left over. A
closure good to 2.6% in $r_i$ therefore lands a 42% error on $v$ — and it gets
worse toward $N_c$, where $v\to0$ and the relative error diverges. To predict $v$
to 10% at $N=6$ you would need $r_i$ to **0.6%**; at $N_c$, to arbitrary precision.

**This is the honest answer to "can it be described from the drop
probabilities?"**: the composition, yes, to 1–2%. The velocity, no — not because
the physics is opaque, but because $v$ is a small difference of large numbers and
no closure of this kind will ever be accurate enough near $N_c$. Measuring $\rho$
sidesteps the cancellation entirely, which is why $\sum p_i\rho_i$ with measured
$\rho$ reaches 12$\times$ noise while the derived version reaches only 70$\times$.

## 6. Can it be rescued as a function of $N$ alone? Apparently not.

The other repair is to *predict* $\rho$ empirically and substitute. Fitting
$\rho_{\rm imp}=Af^b$ and forming $\sum_i p_i \rho_i^{\rm model}$ gives RMS
$0.00597$ (**30$\times$ noise**) — worse than the entropies it was meant to beat.

The reason: **$b$ is not universal.** Per-interval fits give

| interval | [6,7] | [7,8] | [8,9] | [9,10] |
|---|---|---|---|---|
| $b$ | 0.329 | 0.467 | 0.564 | 0.632 |

so there is no single exponent, and a one-law model cannot cover the range. (The
drift toward $\approx0.5$ at larger $N$ is consistent with the $\sqrt{f}$
argument becoming accurate when the impurity is dilute in a large pile; the small-$N$
intervals are where cascades matter most and the naive pair-balance is worst. Not
tested.)

**Conclusion: as a function of $N$ alone, no linearising transformation was
found, and the $\sqrt{f}$ cusp is a structural reason to doubt a simple one
exists.** $\sum_i p_i\rho_i$ linearises $v$, but it is a self-consistency
relation — you must measure the pile to evaluate it.

## Files

- `velocity.cpp` — dynamics **identical** to `../criticalScaling/criticalScaling.cpp`
  and `../avalancheScaling/avalancheDist.cpp` (verified: final $M/L$ matches
  `avalancheDist`'s `pile_total/L` exactly, incl. fractional $N$). Records
  $M/L$ at 200 **linearly** spaced post-warmup times (v is a straight-line fit,
  so linear spacing, unlike the log spacing used for critical scaling) plus the
  final pile composition. CLI: `L N steps sim H [warmup]`.
- `common.py` — loaders, `weights`, `hill`, `candidates`, `linearity`.
- `figures.py` — `staircase`, `collapse`, `entropy_scan`, `impurity_enrichment`.
- `run_sweep.py` — the sweep (N = 5.5…10.0 step 0.1, 16 sims).
- `closure.py` — the exact steady-state balance p_i - r_i = v.rho_i, the pair
  closure r_i = k.rho_i^2, and the test of whether it predicts v from p alone.

```sh
g++ -O3 -march=native -std=c++17 -o velocity velocity.cpp
python run_sweep.py | xargs -P 16 -L 1 ./velocity > /dev/null
python figures.py
```

## Caveats

- Range is $N\in[5.5,10]$ at one $L=1024$. $v$ is a bulk quantity far from $N_c$,
  but its $L$-dependence was **not** checked here — and near $N_c$ it is severe
  (see `../criticalScaling`). Do not extend this analysis below $N\approx5.5$
  without an $L$ check.
- $\rho$ is the composition at the final time only, not a time average.
- The mass-balance step ($f_{\rm active}\propto P(\text{match})$, $\langle s\rangle$
  roughly constant) is hand-waving, as requested — $\langle s\rangle$ does drift
  (3.27→3.37 over $N=5.075..6$, `../avalancheScaling`), which is likely part of
  the residual 1% non-linearity.
