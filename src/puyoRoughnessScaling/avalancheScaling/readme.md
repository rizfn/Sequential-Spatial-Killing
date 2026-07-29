# Avalanches, and the growth-arrest transition at $N_c$

**There is a genuine critical point at $N_c=5.075(10)$ — and the avalanches do
not care.** $\langle s^2\rangle/\langle s\rangle$ is flat in $L$ from 128 to
4096 *at $N_c$ itself*. The reason is structural, and it is the main result
here: arrest is a condition on the **first** moment of $P(s)$, which needs no
divergence to satisfy.

## The transition is real

$N$ is continuous (the last species carries the fractional weight, as in
`probabilityPuyoPuyo/onlyAvalanche2D.cpp`), so we can tune through it. The
growth velocity $v=\mathrm{d}\langle h\rangle/\mathrm{d}t$ is an order
parameter and vanishes continuously (`transition_vN`, $L=512$, 8 sims):

| $N$ | 5.04 | 5.06 | 5.07 | 5.08 | 5.10 | 5.12 | 5.20 |
|---|---|---|---|---|---|---|---|
| $v$ | 0.00000(7) | −0.00004(33) | −0.0002(9) | 0.0022(12) | 0.0074(10) | 0.0112(21) | 0.0265 |

so $N_c=5.075\pm0.010$. Below it the pile arrests at $\langle h\rangle\approx5$;
above it it grows freely. At $N_c$ the pile still grows, but as a **sublinear
power law $\langle h\rangle\sim t^{0.32}$** rather than $\sim vt$ — the scaling
signature of a genuine critical point. So the criticality is real, and it lives
in the height/time sector.

## But the avalanches are not scale-free — not even at $N_c$

At a critical point $\xi\to\infty$, so the *only* remaining cutoff is $L$ and
$\langle s^2\rangle/\langle s\rangle\sim L^{D}$ must grow without bound. It does
not (`critical_test`):

| $L$ | 128 | 256 | 512 | 1024 | 2048 | 4096 |
|---|---|---|---|---|---|---|
| $N=5.075\;(=N_c)$ | 8.466 | 8.573 | 8.646 | 8.680 | 8.643 | 8.652 |
| $N=5.080$ | 8.488 | 8.625 | 8.664 | 8.689 | 8.723 | 8.732 |

A factor of 32 in $L$, sitting on the critical point, and the avalanche scale
does not move. (Steady state is verified, not assumed — `outputs/moments/` shows
the plateau is reached by $t\approx7000$ even at $N_c$, where critical slowing
down is a real worry, and $\mathrm{WARMUP}=8192$.)

## Why: mass balance fixes the *first* moment

Each deposition adds exactly one block and removes, on average,
$f_{\rm active}\langle s\rangle$ — the probability it triggers anything times the
mean cascade mass. So, exactly,

$$\boxed{\;v = 1 - \langle s\rangle\, f_{\rm active}\;}$$

and arrest ($v=0$) is precisely $\langle s\rangle f_{\rm active}=1$. Measured
(`mass_balance`):

| $N$ | $\langle s\rangle$ | $f_{\rm active}$ | $\langle s\rangle f_{\rm active}$ | $1-\langle s\rangle f$ | $v$ measured |
|---|---|---|---|---|---|
| 5.075 | 3.269 | 0.3057 | **0.9993** | 0.0007 | 0.0000 |
| 5.200 | 3.290 | 0.2958 | 0.9734 | 0.0266 | 0.0265 |
| 6.000 | 3.373 | 0.2792 | 0.9417 | 0.0583 | 0.0612 |

The identity holds to three decimals. Now the punchline: $\langle s\rangle$ is
the **first** moment, and it is dominated by the swarm of tiny events
($P(s=2)\approx0.74$). Tuning it to exactly 1 requires no diverging cutoff
whatsoever — $\langle s\rangle$ stays $\approx3.3$ while $f_{\rm active}$ drifts.
The cutoff $\langle s^2\rangle/\langle s\rangle$ is a *different* moment ratio
and is free to remain finite. That is exactly what it does.

**Contrast the BTW sandpile**, where this is the whole ballgame: dissipation
happens only at the open boundary, so in steady state every added grain must
random-walk *out*, which forces $\langle s\rangle\sim L^{2}$ — a diverging first
moment, hence forced scale-free avalanches. Puyo dissipates **in the bulk**
(a matched cluster vanishes wherever it forms), so mass balance is satisfied
locally with $\langle s\rangle = 1/f_{\rm active} = O(1)$. Nothing forces a
divergence, and none appears. This is the structural reason the model is not SOC.

## The $N=6$ peak is a fractional-species artifact

The avalanche scale vs $N$ is a **sawtooth with maxima at integer $N$**
(`speciesSawtooth/cutoff_vs_N`, $L=1024$):

| $N$ | 5.075 | 5.5 | **6.0** | 6.25 | 6.5 | **7.0** | 7.5 | **8.0** |
|---|---|---|---|---|---|---|---|---|
| $\langle s^2\rangle/\langle s\rangle$ | 8.68 | 10.21 | **10.59** | 8.43 | 8.80 | **8.97** | 7.85 | **7.97** |

### Why: the fractional species is a self-poisoning impurity

At $N=n+f$ the last species carries weight $f$, so it is deposited with
probability $f/N$. Being rare, it seldom finds a same-species neighbour -- so it
is seldom eliminated, and therefore **accumulates**.

Balance its rates. It arrives at rate $\propto f$. It can only leave in pairs
(two impurities must meet), so it departs at rate $\propto\rho^{2}$ in its pile
density $\rho$. Steady state $f\sim\rho^{2}$ gives

$$\rho\sim\sqrt{f},\qquad \text{enrichment}=\frac{\rho}{f/N}\sim f^{-1/2}$$

Measured over one tooth ($L=512$, `speciesSawtooth/sawtooth_mechanism`):

$$\rho\sim f^{0.415},\qquad \text{enrichment}\sim f^{-0.555}$$

against the predicted $\pm1/2$. The enrichment reaches **5.6$\times$** at
$N=6.05$ (deposited at 0.83% of blocks, but 4.7% of the pile) and decays to
$1.0$ at $N=7$, where the impurity has become an ordinary species.

This explains the sawtooth's **asymmetry**. $\rho\sim\sqrt{f}$ has *infinite
slope* at $f=0$: an arbitrarily rare impurity still builds up a
disproportionate frozen-defect density, which fragments clusters. So the scale
plummets the instant you leave an integer ($10.52\to8.61$ between $N=6.00$ and
$6.05$), bottoms out near $f\approx0.1$, then recovers slowly as the impurity
becomes a real species -- reaching a fresh maximum at the next integer. Both
ends of a tooth are integers ($f=0$ and $f=1$ are the same uniform system), so
the dip must vanish at both, which is exactly what is seen.

The integer peaks themselves decay with $N$ (10.59, 8.97, 7.97) simply because
matches get rarer. So the global maximum at $N=6$ is *the largest integer peak
above $N_c$* -- the compromise between wanting low $N$ (dense matching) and
needing $N>N_c$ (a pile at all). It is not a critical point: every column above
saturates in $L$.

## Off-critical $L$ scaling ($N=6$)

For completeness, the original family: 16 sizes, $L=16\ldots4096$. Everything
converges by $L\approx512$, and the cutoff probe converges to a **finite** limit

$$\frac{\langle s^2\rangle}{\langle s\rangle} = 10.70(3) - 25.2\,L^{-0.76(3)},$$

with active fraction $0.2792$ at *every* $L$ to four digits. $P(s)$ collapses
with no rescaling at all. This is what off-critical looks like: $\xi$ finite,
$L\gg\xi$, no $L$ dependence.

$s_{\max}$ *does* grow with $L$, and that is a trap: it is a sample-size effect
(more depositions reach further into a *fixed* tail). It is visibly
non-monotonic here because `NSIMS` varies with $L$ — the tell that it tracks
sample count, not physics. $\langle s^2\rangle/\langle s\rangle$ is the honest
probe, and needs no $\tau$ fit: for $P\sim s^{-\tau}f(s/s_c)$,
$\langle s^k\rangle\sim s_c^{k-\tau+1}$, so the ratio $\sim s_c$ with $\tau$
cancelling.

## Shape of P(s) at integer N: one fixed shape, only the cutoff moves

*(This supersedes both the earlier "Weibull, and why" and the intermediate
"stiffens with N" analyses. `integerN.py`, full 1024-sim family. Fractional N is
a different system — its rare partial species is a self-poisoning impurity, the
sawtooth above — so this family is **integer N only**: $N=6,7,\dots,20$, all at
$L=1024$.)*

**The fix is the fitting method, not the functional form.** The tail is
$$P(s)\;\propto\;s^{-\tau}\,\exp\!\left[-(s/s_0)^{b}\right],\qquad s\gtrsim20,$$
the same three-parameter power-law-times-stretched-exponential the readme called
"unidentifiable." It is only unidentifiable under **least-squares on the binned
PDF**, where a free amplitude $A$ absorbs whatever $\tau$ and $b$ do and the fit
runs to $\tau=-9.57$. Fit instead by **discrete maximum likelihood on the raw
integer histogram**,
$$\ln\mathrm{Lik}=\sum_s \mathrm{count}[s]\,\ln p(s),\qquad p\ \text{normalized on }s=s_{\min}\dots\text{CAP},$$
and the amplitude is no longer free — it *is* the normalization — so the
$\tau\!\leftrightarrow\!b$ ridge is lifted and the fit is sharp. At $N=6$ it
lands at $\tau=2.02(11)$, $b=0.617(46)$, $s_0=16.7$, tracking $P(s)$ to a
worst-case factor $\approx1.25$ over $\sim8$ decades.

What the tail **is** — a sub-exponential cutoff — is clear. A pure power law is
off by $97\times$; a plain exponential cutoff ($b=1$) is rejected
($\Delta\ln\mathrm{Lik}=-297$ at $N=6$); a compressed cutoff ($b=1.76$) is
rejected ($-2.1\times10^4$, worst factor $9.7$). So the tail is fatter than
exponential.

**But — honest caveat — the power-law *prefactor* is not required by $P(s)$
itself.** Fit the tail with a *pure* stretched exponential $\exp[-(s/s_0)^b]$ (no
prefactor) and it tracks $P(s)$ to a worst-case factor $1.1$–$1.3$, against
$1.0$–$1.1$ for the power-law$\times$stretched form; the power law buys only
$\Delta\ln\mathrm{Lik}\approx30$–$270$. That is the *same* order of "improvement"
this section dismisses as a sample-size effect when it freezes the exponents
below, so consistency forbids leaning on it here: **$P(s)$ alone does not
establish the prefactor, and $b$ is not even defined without it** (pure stretched
gives $b\approx0.33$; with the prefactor $b\approx0.6$). The justification for the
power law comes not from $P(s)$ but from the extent distribution $P(w)$ — where
it is unambiguous (next). The quoted $\tau=2.03,\,b=0.62$ is the shape that form
takes once the prefactor is granted.

### Freezing both exponents costs nothing on the genuine tail

A per-$N$ free fit of $(\tau_N,\,s_{0,N},\,b_N)$ makes $\tau$ *and* $b$ appear to
rise with $N$ ($b$ even crossing 1 near $N\approx11$). That is a
$\tau\!-\!b\!-\!s_0$ **ridge artifact**: the three parameters trade off inside
the abundant bins just above $s_{\min}$, and with $\sim10^9$ tail counts the fit
exploits sub-percent curvature there and books it as moving exponents. On the
genuine tail the motion vanishes. Fit the same data four ways by the same MLE —
full per-$N$; one shared $b^*$; one shared $\tau^*$; one shared $(\tau^*,b^*)$
with only $s_0$ free — and raise $s_{\min}$ to strip the near-threshold bins:

| $s_{\min}$ | freeze $b$ | freeze $\tau$ | freeze **both** | worst-factor, freeze both |
|---|---|---|---|---|
| 30 | $-2986$ | $-3949$ | $-5498$ | 1.16 |
| 50 | $-131$ | $-146$ | $-206$ | 1.10 |
| 70 | $-20$ | $-22$ | $-37$ | 1.14 |
| 100 | $-6$ | $-6$ | $-11$ | 1.13 |

($\Delta\ln\mathrm{Lik}$ against the full per-$N$ fit.) Freezing *both* exponents
costs 5500 nats at $s_{\min}=30$ but only 11 at $s_{\min}=100$, while the tail
stays tracked to $\sim1.1\times$ throughout — the "cost" was always the
near-$s_{\min}$ bins, never the tail. So the honest description is a **single
universal shape**, $\tau=2.03,\ b=0.62$, whose only $N$-dependence is the cutoff:

| $N$ | 6 | 7 | 8 | 9 | 10 | 11 | 12 | 14 | 16 | 20 |
|---|---|---|---|---|---|---|---|---|---|---|
| $s_0(N)$ | 16.9 | 13.0 | 10.7 | 9.3 | 8.4 | 7.8 | 7.3 | 6.6 | 6.2 | 5.7 |
| $\langle s^2\rangle/\langle s\rangle$ | 10.6 | 9.0 | 8.0 | 7.3 | 6.9 | 6.5 | 6.3 | 5.9 | 5.6 | 5.3 |

$s_0$ tracks the fit-free cutoff $\langle s^2\rangle/\langle s\rangle$, worst-factor
is $1.03$–$1.10$ at every $N$, and all ten $N$ collapse onto the one master curve
$\exp[-(s/s_0)^{b}]$ after rescaling by $s^{\tau}$ and $s/s_0(N)$
(`integerN/tail_collapse_fixedshape`).

### The mechanism: $b=1/d$, the composition law resurrected

$P(s)$ is not fundamental — it is inherited from the cascade's spatial **extent**
$w$ (distinct columns eliminated), through two measured facts (`extentMechanism.py`,
`integerN/extent_mechanism`):

1. **Geometry.** Mass grows with extent as a power, $\langle s\,|\,w\rangle\sim
   w^{d}$ with $d\approx1.71$ **flat in $N$**. A width-$w$ cascade digs a depth
   $\sim w^{d-1}$, i.e. it is a compact 2D-ish patch of surface, not a thin sliver.
2. **Spreading.** The extent is a **truncated power law**,
   $P(w)\sim w^{-a}e^{-w/w_*}$ ($a\approx3.0\to2.2$, $w_*\approx4.9\to1.8$), and
   *this* is where the power law is data-forced, not assumed: by discrete MLE a
   plain exponential is off by $10^4$–$10^5$, and even a stretched exponential or
   a lognormal in $w$ is rejected (worst factor $2$–$5$ vs $1.1$ for
   $w^{-a}e^{-w/w_*}$, $\Delta\ln\mathrm{Lik}\sim10^3$–$10^4$). The power-law form
   means the lateral spread is (weakly near-critical) **branching**: an eliminated
   column can topple its neighbours, a scale-free process, cut off at a finite
   correlation width $w_*$ set by how often a neighbour actually matches. This
   is the honest source of the $P(s)$ prefactor — $P(s)$ alone
   cannot see it (above), but through $s\sim w^d$ the power law in $P(w)$ becomes
   the power law in $P(s)$.

Change variables $s=Cw^{d}$ in $P(w)$:

$$P(s)\sim s^{-\tau}\exp\!\big[-(s/s_0)^{b}\big],\qquad \boxed{b=\tfrac1d},\quad s_0=C\,w_*^{\,d},\quad \tau=\tfrac{a+d-1}{d}.$$

So the **stretch exponent is the reciprocal of the mass–extent exponent**,
$b=1/d\approx0.58$ (fitted $0.62$) — $N$-independent because the geometry $d$ is.
The only $N$-knob is $w_*$ (per-column match rate $\sim1/N$), setting
$s_0=C\,w_*^{\,d}$:

| $N$ | 6 | 8 | 10 | 12 | 16 | 20 |
|---|---|---|---|---|---|---|
| $d$ | 1.66 | 1.71 | 1.72 | 1.73 | 1.71 | 1.77 |
| $1/d$ | 0.60 | 0.58 | 0.58 | 0.58 | 0.59 | 0.57 |
| $w_*$ | 5.3 | 3.2 | 2.7 | 2.9 | 1.7 | 1.6 |

($w_*^{\,d}$ recovers $s_0$ at low $N$: $5.3^{1.66}=15.9$ vs $s_0=16.9$; the
high-$N$ $w_*$ is unreliable, the exponential barely resolved above the $w=1,2$
core.)

This **voids the "$b=1/d$ is refuted" claim of the intermediate draft**, which
rested entirely on "$b$ rises with $N$ while $d$ is flat" — but that $b(N)$ rise
was the fitting artifact above. With $b$ correctly frozen, flat $1/d=$ flat $b$
is the *generic* law, not an $N=6$ coincidence.

### The branching exponent, and the roughness is decoupled

What stays soft is $\tau$, the **mass** exponent. Take it not from $P(s)$ (where
$\tau\!\leftrightarrow\!$cutoff is degenerate) but from the extent, via
$\tau=(a+d-1)/d\approx2.0$–$2.2$. **This $\tau$, not $a$, is the branching
exponent**: the cascade of eliminations is a branching process whose natural size
is the mass, and $\tau\approx2$ sits at the upper edge of the branching range.
Mean-field Borel–Tanner ($\tau=3/2,\ b=1$) is rejected precisely because the
branching is embedded in $1{+}1$D (which pushes $\tau$ above the mean-field $3/2$)
and its cutoff is stretched ($b=1/d$), not exponential. The large *extent*
exponent $a\approx2.8$ is not a second puzzle — it is just $\tau$ mapped through
the fractal, $a=d(\tau-1)+1$.

**Cutoff scaling, and no roughness link.** The one $N$-knob
$s_0$ scales as the match probability, $s_0\sim1/N$ (from high-moment ratios;
$\langle s^2\rangle/\langle s\rangle$ reads a shallower slope only because it is
floor-contaminated by the tiny-event bulk). Decompose $s_0=C(N)\,w_*^{\,d}$, with
$C$ the prefactor of $\langle s\,|\,w\rangle=C\,w^{d}$ measured per $N$ — so *not*
collinear with the cutoff. $C$ is essentially flat ($\sim N^{-0.2}$) while the
slope-distribution width grows $\sigma\sim N^{1.34}$: $C\sim\sigma^{-0.14}$ with no
real trend. So the **roughness is decoupled from the size distribution** — all of
the $N$-dependence is the branching cutoff $w_*$, and the surface roughness is a
*parallel* consequence of $N$, not a driver. (A lateral slope *correlation length*
— not measured here — is the only roughness observable that could still enter; the
width and tail-rate do not.)

Below $s\approx18$ a strong **even/odd parity oscillation** (even masses
favoured — the elementary cluster is a pair) means no smooth form applies; all
fits start at $s_{\min}=50$ (past both the parity band and the near-threshold
ridge), and $P(s=2)\approx0.74$ is a discrete elementary event, not an avalanche
in any scaling sense.

## Definitions

One avalanche = the whole cascade triggered by a single deposition:
**mass** $s$ (blocks eliminated), **clusters** $n$ (components eliminated),
**duration** $d$ (chain generations). Depositions eliminating nothing are
counted (`drops_counted`) but not histogrammed, so $P(s)$ is normalized over
*active* avalanches.

## Files

- `avalancheDist.cpp` — the simulation. Same dynamics as
  `probabilityPuyoPuyo/onlyAvalanche2D.cpp` (including **continuous $N$**, which
  is what lets us reach $N_c$), carrying the `../slopeDistFast.cpp` optimizations
  (column heights, moved-site seeds instead of $H\times L$ scans,
  generation-stamped BFS, dirty-column gravity), so per-drop cost is
  $O(\text{active region})$ and a tall box costs only memory. `lat` is `uint8_t`.
  Histograms accumulate in RAM and dump once: output is $O(s_{\max})$, not
  $O(L\cdot\text{steps})$. ~12.8M depositions/s (~6.8M at $L=4096$, out of cache).
- `common.py` — sweep definition (single source of truth), pooling/log-binning,
  `moments_vs_t` (steady-state check), `velocity` (order parameter),
  `slope_resolved` / `composition` (mechanism), `weibull_slope`.
- `run_sweep.py` / `figures.py`.
- `run_integerN.py` / `integerN.py` — the integer-$N$ family (1024 sims/N at
  $L=1024$). `integerN.py` fits the one universal shape by discrete MLE (shared
  $\tau,b$; only $s_0(N)$ free — the amplitude *is* the normalization, which lifts
  the $\tau\!\leftrightarrow\!b$ ridge) and makes `avalanche_pdf_integerN` and
  `tail_collapse_fixedshape`. Pooling 1024 files/N is slow, so it caches the
  pooled histograms to a hidden `.integerN_cache.npz` (delete to repool). The
  discarded per-$N$ / single-exponent comparisons live in the git history.
- `extentMechanism.py` — the $b=1/d$ mechanism: measures $\langle s\,|\,w\rangle
  \sim w^{d}$ and $P(w)\sim w^{-a}e^{-w/w_*}$ from `outputs/slopeResolved/`, makes
  `integerN/extent_mechanism`; also the honesty checks `test_ps_prefactor` (is the
  $P(s)$ power law required?) and `test_pw_form` (is $P(w)$ a power law $\times$ exp?).

### outputs/
| folder | contents |
|---|---|
| `avalancheDist/` | the $s$ / $n$ / $d$ histograms, over $t\ge$ WARMUP |
| `moments/` | moments of $s$ + mean height in log-spaced time windows, over ALL $t$ — steady-state check and $v(N)$ |
| `slopeResolved/` | $\langle s\,\vert\,m\rangle$, $\langle s\,\vert\,w\rangle$, and the final species composition — the mechanism tests |

### plots/
Themed subfolders (`figures.py`'s `save()` takes the theme):
| folder | contents |
|---|---|
| `criticality/` | `transition_vN`, `critical_test`, `mass_balance` |
| `mechanism/` | `mechanism_slope_vs_extent` |
| `speciesSawtooth/` | `cutoff_vs_N`, `sawtooth_mechanism` |
| `finiteSize/` | `cutoff_vs_L`, `steady_state` |
| `sizeDistribution/` | `avalanche_pdf_vs_L`, `avalanche_pdf_vs_N`, `clusters_duration` |
| `integerN/` | `avalanche_pdf_integerN`, `tail_collapse_fixedshape`, `extent_mechanism` |

## Reproducing

```sh
g++ -O3 -march=native -std=c++17 -o avalancheDist avalancheDist.cpp
python run_sweep.py | xargs -P 12 -L 1 ./avalancheDist > /dev/null
python figures.py

# the integer-N family (separate, heavier: 1024 sims/N for the tail statistics)
python run_integerN.py | xargs -P 24 -L 1 ./avalancheDist > /dev/null
python integerN.py         # the fixed-shape fit + collapse
python extentMechanism.py  # the b = 1/d mechanism + honesty checks
```

The main sweep is 1480 sims; **3.5 min wall** on 12 cores (42 core-min). `-P` is
bounded by memory rather than cores: a job holds $5HL$ bytes, ~212 MB for the
largest, and jobs are emitted biggest-box-first. `box_H` is sized from the
measured $v(N)$ (which is ~0 at $N_c$ but 0.216 at $N=7$, so a single ratio would
either waste memory or hit the ceiling); every run reports `ceiling_hits`, 0 for
all 1480. The integer-N family adds ~9k sims (~35 min on 24 cores, 0 ceiling
hits), each holding $5HL$ bytes with $H\propto v(N)$, ~137 MB at $N=20$.

## Caveats

- $N_c=5.075(10)$ comes from $v(N)$ at $L=512$ with 8 sims; no finite-size
  scaling of $N_c$ itself was done. Since we bracket it by $\pm0.01$ and test
  *both* 5.075 and 5.080 with identical flat results, the conclusion is not
  sensitive to this.
- $\langle h\rangle\sim t^{0.32}$ at $N_c$ is a single-$N$, single-$L$ fit over
  ~2 decades. It is clearly sublinear and unbounded; the exponent itself is
  provisional (suggestively $\approx1/3$).
- The ratio argument assumes $\tau<3$ (see `cutoff_vs_L` discussion). The primary
  evidence is the raw $P(s)$ collapse, which constrains the whole distribution.
