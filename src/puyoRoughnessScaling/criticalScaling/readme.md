# The growth-arrest transition is an interface *unbinding* transition

**Verdict up front: this is real physics, it is better than the previous
description of it, and it is almost certainly *not* a PRL.** The reason is not
that the transition is uninteresting — it is that the transition already has a
name, a universality class, a founding PRL from 1997 and a 2009 review; and the
exponents that could show we are *outside* that class are, in this model, not
measurable in the window we can reach. Details below, including the numbers
that would be needed to overturn this.

## 1. What the transition actually is

Below $N_c$ the pile does **not** simply "arrest at $\langle h\rangle\approx5$".
The arrest height **diverges** as $N\to N_c$, and it stays roughly proportional
to the interface width. So the surface is *bound to the substrate* at a finite
binding length $h^*$, and at $N_c$ it *comes unbound* and starts to move. That
is an interface unbinding (wetting) transition, and $h^*$ is the diverging
length that was previously missing from the picture.

$L=4096$ (except $N=5.075$: $L=8192$), $t=2\times10^6$, no ceiling hits:

| $N$ | 5.000 | 5.020 | 5.040 | 5.050 | 5.060 | 5.065 | 5.070 | 5.072 | 5.075 |
|---|---|---|---|---|---|---|---|---|---|
| $h^*$ | 5.57 | 9.85 | 14.57 | 18.14 | 24.42 | 29.95 | 41.43 | 51.03 | **97.22** |
| $W^*$ | 3.72 | 4.83 | 6.28 | 7.50 | 9.43 | 11.33 | 15.23 | 17.80 | 31.81 |
| $h^*/W^*$ | 1.50 | 2.04 | 2.32 | 2.42 | 2.59 | 2.64 | 2.72 | 2.87 | 3.06 |

$h^*/W^*$ is *not* constant — it drifts from 1.5 and appears to approach $\approx3$
near $N_c$. The drift lives entirely in the small-$h^*$ points, which is what
§4 says it should be. Above $N_c$ the picture inverts: $\langle h\rangle$ grows
linearly while $W$ saturates, i.e. an ordinary free KPZ interface that has
detached from the wall.

## 2. $N_c$ is bracketed, not fitted — and the old value was biased

At $L=4096$ every $N\le5.075$ arrests (late log-slope $\le0.07$) and every
$N\ge5.078$ grows linearly (late log-slope $\ge0.85$). So

$$N_c = 5.0765 \pm 0.0015$$

Bracketing rather than fitting matters: it is the *only* reason $\psi$ below is
identifiable at all (§3).

**The previous $N_c=5.075(10)$ was measured at a single $L=512$ and is biased.**
At $N=5.075$ the behaviour depends entirely on $L$ (all with adequate box
height, zero ceiling hits):

| $L$ | 512 | 1024 | 2048 | 4096 | 8192 |
|---|---|---|---|---|---|
| $\langle h\rangle(2\times10^6)$ | 804 | 274 | 112 | 92.3 | 97.2 |
| late log-slope | 0.94 | 0.78 | 0.38 | **0.07** | **0.09** |
| verdict | grows | grows | grows | **arrests** | **arrests** |

$h^*$ converges to $\approx95$ for $L\ge4096$, so the arrest survives the
thermodynamic limit — but at $L=512$, the size the old $N_c$ came from, the
pile still grows. The old measurement read $v(5.075)=-0.0002\pm0.0009$ and
called it zero; the true value there is $\sim4\times10^{-4}$, i.e. *below their
own noise floor*. A velocity that vanishes only as $L\to\infty$ cannot be
located from one $L$.

This also **retracts an intermediate claim made during this analysis**: the
clean power law $v(N{=}5.075,L)\sim L^{-1.2}$ over $L=64..1024$ is *not*
critical scaling, because $N=5.075$ is in the bound phase, where $v(L)\to0$
regardless. It was measured at the wrong $N$.

## 3. The binding-length exponent

Fitting $h^*\sim(N_c-N)^{-\psi}$ over the points with $h^*>\ell^*$ (see §4),
with $N_c$ **pinned by the bracket**:

$$\psi = 0.58 \pm 0.03\ \text{(fit)} \pm 0.08\ \text{(from } N_c)$$

The fit is excellent ($R^2=0.9999$) and the pairwise local slopes show no drift
(0.62, 0.56, 0.57, 0.59). Sensitivity to the bracket: $\psi=0.50$ at
$N_c=5.0760$, $0.58$ at $5.0765$, $0.71$ at $5.0775$.

**This only works because $N_c$ is bracketed.** Left free, the 3-parameter fit
is *unidentifiable* — it runs to $N_c=5.0945$, $\psi=1.19$, contradicting the
bracket, and a logarithm $h^*=A-B\ln(N_c-N)$ fits as well ($R^2=0.997$). This is
the same failure mode as the 4-parameter stretched exponential in
`../avalancheScaling`. Never fit $N_c$ and $\psi$ together here.

## 4. Why the exponents cannot settle the universality class

Standard scaling gives $h^*\sim\xi_\perp^{\alpha}$, so $\psi=\nu_\perp\alpha$.
With KPZ $\alpha=1/2$:

| class | $\nu_\perp$ | predicted $\psi$ | vs measured 0.58(9) |
|---|---|---|---|
| bEW (equilibrium) | 2/3 | 0.333 | **excluded** |
| DP | 1.097 | 0.548 | consistent |
| measured | — | 0.58(9) | — |

The DP row is *suggestive only* — I could not extract published $\nu_\perp$ for
bKPZ$^-$, so this is my scaling relation applied to a guessed $\nu_\perp$, not a
class identification. Do not quote it.

**The hard obstruction is scale separation.** From `roughness-not-kpz`, the model
has a UV cutoff $\ell^*(N)\sim N^{1.6}$ — the Gaussianisation length below which
*no continuum, hence no KPZ, exists*. Near $N_c$, $\ell^*\approx20$. The last
column of §1 is $h^*/\ell^*$: it reaches only **4.9** at the closest $N$ we can
equilibrate, and is $<1$ for $N\le5.05$. The bound interface is a handful of UV
cutoffs thick. This is precisely the trap that produced the bogus $\alpha(N)$
flow in the roughness work, in a new observable.

Compounding it: $N_c=5.0765$ sits only **0.077 above the integer $N=5$**, where
the $\sqrt{f}$ self-poisoning-impurity cusp lives (`avalanches-not-critical`:
$\rho\sim f^{0.415}$). The bound-phase window is squeezed between two
non-analyticities of comparable distance. The far points are cusp-contaminated;
the near points are UV-limited.

So $\psi=0.58$ is an **effective exponent**, not an asymptotic one.

## 5. The novelty problem: this class already exists

Non-equilibrium wetting *is* the bounded-KPZ (bKPZ) equation,

$$\partial_t h = a - \frac{\mathrm{d}V}{\mathrm{d}h} + \sigma\nabla^2h + \lambda(\nabla h)^2 + \zeta$$

The 2009 review ([arXiv:0908.3068](https://arxiv.org/abs/0908.3068)) describes
the field as "the KPZ equation with a soft-wall potential", realised by
"microscopic models in the KPZ universality class featuring particle evaporation
and deposition near hard walls". That is a literal description of this model:
KPZ interface (established asymptotically in `roughness-not-kpz`), deposition,
annihilation-as-evaporation, hard floor. Founding paper: Hinrichsen, Livi,
Mukamel & Politi, **PRL 79, 2710 (1997)**
([cond-mat/9706031](https://arxiv.org/abs/cond-mat/9706031)).

Classes are set by the sign of $\lambda$: bKPZ$^-$ ($\lambda<0$), bKPZ$^+$
($\lambda>0$), bEW ($\lambda=0$, detailed balance). This model's own gradient
expansion (`../../paperDraft/readme.md`) gives $V(m)=(1-2p)-cpm^2$, i.e.
$\lambda=-cp<0$ — growth is *fastest on a flat surface* — placing it in
**bKPZ$^-$**. Caveat: $\lambda$'s sign here is taken from the model's derivation,
**not measured**; the direct test is to tilt the interface (helical BCs) and fit
$v(m)$.

And the review quotes bKPZ$^-$ inheriting $\alpha=1/2$, $z=3/2$, which predicts
$\langle h\rangle\sim t^{\alpha/z}=t^{1/3}$ at $N_c$. The bEW analogue checks
out ($\alpha=1/2$, $z=2\Rightarrow\gamma=1/4$, the review's quoted value). **So
the old "$\langle h\rangle\sim t^{0.32}$ at $N_c$" was the textbook bKPZ
prediction all along** — it was measured, but not recognised.

Our own $t$-dependence is too noisy to confirm $1/3$: at $N=5.075$ the local
slope bounces (0.43, 0.14, 0.24, 0.02, 0.06, 0.33) on 4 sims, and no $N$ shows a
sustained $t^{1/3}$ — because no accessible $N$ sits close enough to $N_c$ with
$L$ large enough.

## 6. So: PRL?

**No.**

- The phenomenon is a **known universality class** with a 1997 PRL and a
  review. "We found nonequilibrium wetting in a new lattice model" is a PRE.
- The one claim that *would* be a PRL — *these exponents are not bKPZ$^-$* —
  needs $\psi$ to $\pm0.02$ asymptotically. We have $\pm0.09$, effective, at
  $h^*/\ell^*\le5$, with a cusp 0.077 away. The measurement cannot currently
  distinguish "new class" from "bKPZ$^-$ plus corrections to scaling".
- The genuinely novel ingredients are microscopic, and universality is designed
  not to care about them: evaporation here is *collective* (whole matched
  clusters plus cascades) rather than single-particle, and the wall potential is
  *emergent* rather than imposed. But the cascades are short-ranged
  ($w_0\approx4.8$ columns, exponential cutoff — `avalanches-not-critical`), so
  by the standard argument they should not change the class. Verifying that a
  collective-evaporation model still lands in bKPZ$^-$ is a real but modest
  contribution.

**What would change the verdict.** Either (a) show $\lambda$'s sign by tilting and
find bKPZ$^+$ behaviour where bKPZ$^-$ is expected, or (b) reach
$h^*/\ell^*\gtrsim10$ and show $\psi$ misses bKPZ$^-$. (b) needs $h^*\gtrsim200$:
$N_c$ known to $10^{-4}$, $L\gtrsim4\times10^4$ ($\xi_\parallel\sim h^{*2}$),
$t\gtrsim10^7$ — roughly 10 core-hours/sim, so ~feasible, but only *after* the
control-parameter problem is fixed.

**The control parameter is the real blocker.** Tuning $N$ via a fractional
species puts $N_c$ next to the integer-$N$ cusp *and* couples the impurity
density to the distance from criticality ($\rho\sim\sqrt{N-5}$ varies across the
whole scaling window). A knob without that pathology — e.g. integer $N=5$ with a
tunable annihilation probability — would open the window and decouple the two
singularities. That is a *different model*, so it is flagged here rather than
done.

## Files

- `criticalScaling.cpp` — dynamics **bit-identical** to
  `../avalancheScaling/avalancheDist.cpp` (verified: final $\langle h\rangle$
  matches `pile_total/L` exactly on 3 parameter sets incl. fractional $N$); only
  the recording differs. Samples $\langle h\rangle$, $W$, active fraction and
  $\langle s\rangle$ at 400 log-spaced times. CLI: `L N steps sim H`.
- `common.py` — loaders, `local_slope`, `velocity`, `linear_frac`,
  `ceiling_hits`, collapse cost.
- `figures.py` — `growth_vs_N`, `unbinding`, `finite_size`, `psi_fit`,
  `scale_audit`.

```sh
g++ -O3 -march=native -std=c++17 -o criticalScaling criticalScaling.cpp
python figures.py
```

## Caveats

- **Always check `ceiling_hits`.** A too-short box caps the pile and *fakes an
  arrest*. The first $L=1024$, $N=5.075$ batch hit the ceiling $7\times10^8$
  times ($H=500$) and had to be discarded; `plateau()` in `figures.py` refuses
  ceiling-contaminated runs for this reason.
- The near-critical points rest on 4 sims and fluctuate heavily; $h^*(5.072)$ in
  particular is noisy. More sims would tighten $\psi$ but cannot fix the
  UV-cutoff limit, which is the binding constraint.
- $\ell^*\approx20$ near $N_c$ is *extrapolated* from $\ell_G=19.8$ at $N=6$; it
  was not measured at $N\approx5.08$. Since $\ell^*\sim N^{1.6}$ is increasing,
  20 is if anything an over-estimate near $N_c$, so $h^*/\ell^*$ is a mild
  under-estimate — this does not rescue the scale separation.
- $\lambda<0$ is inferred from the model's gradient expansion, not measured.
