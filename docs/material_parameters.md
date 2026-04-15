# SS316L solid-phase material priors (room temperature)

This document proposes **starting values** and **bounds** for the **solid-phase**
material parameters used by the LS-EVPFFT solver in this repo (`mimosa`).

This Markdown file is **generated** by:

* `tmp/generate_ss316l_solid_priors_doc.py`

Do not hand-edit `docs/material_priors_ss316l_solid.md`; instead, edit the
generator script and rerun it.

It is written to support the near-term goal:

1. Run stable creep simulations (high loads up to ~575 MPa, large time increments) with a dummy gas phase.
2. Calibrate *solid* material parameters first (AM SS316L, room temperature), before widening the calibration space to gas parameters.

Although the gas phase strongly affects numerical stability, the **first task**
is to choose solid parameters that are physically reasonable for **high-yield AM
316L** at room temperature and are good priors for Bayesian optimization.

This document defines all concepts that are used (Taylor factor, CRSS, rate sensitivity, Voce hardening, etc.) and explains why the proposed numerical ranges are sensible.

---

## 0. TL;DR: recommended “v2” priors (copy/paste)

This section is meant to be immediately actionable: it lists recommended
starting values and calibration bounds, then later sections justify them in
depth.

### 0.1 Solid elasticity (cuel.sx / SolidElastic)

The solid phase uses **cubic** elastic constants (C11, C12, C44). Start by
**fixing** elasticity (see Section 5).

Repo-default values (MPa):

| Parameter | Start (MPa) | Suggested bounds (MPa) |
|---|---:|---:|
| C11 (`c11`) | 206000 | 195000–215000 |
| C12 (`c12`) | 133000 | 125000–145000 |
| C44 (`c44`) | 119000 | 115000–130000 |

### 0.2 Solid slip kinetics + hardening (cupl2.sx / SolidPlastic)

These parameters define FCC slip kinetics + extended Voce hardening.

Recommended starting values (AM 316L, room temperature):

| Parameter | Start | Units | Meaning |
|---|---:|---|---|
| `tau0xf` | 185 | MPa | initial slip resistance (forward) |
| `rb` | 0.95 | - | ratio `tau0xb/tau0xf` |
| `tau0xb` | 175.75 | MPa | initial slip resistance (backward) |
| `tau1x` | 200 | MPa | Voce saturation increment |
| `thet0` | 5000 | MPa | initial hardening slope vs accumulated slip |
| `thet1` | 300 | MPa | long-term hardening slope vs accumulated slip |
| `nrsx` | 20 | - | stress exponent ($m=1/n=0.05$) |
| `gamd0x` | 1e-4 | 1/s | reference shear rate |
| `hselfx` | 1.0 | - | self hardening coefficient |
| `hlatex` | 1.4 | - | latent hardening coefficient |

### 0.3 Recommended calibration bounds (solid plasticity)

The following bounds are tuned for AM 316L at room temperature and for stable,
well-scaled Bayesian optimization.

| Variable (calibration space) | Low | High | Notes |
|---|---:|---:|---|
| `tau0xf` (MPa) | 130 | 250 | matches AM yield scales via Taylor-factor reasoning |
| `rb` (-) | 0.90 | 1.00 | keep forward/backward nearly symmetric |
| `tau1x` (MPa) | 80 | 300 | allows $tau0+tau1$ saturation scale up to ~550 MPa |
| `thet0` (MPa) | 2000 | 9000 | controls how fast saturation is approached |
| `thet1` (MPa) | 50 | 800 | long-term hardening slope |
| `log10_gamd0x` | -7.5 | -3.0 | $gamd0x$ in [3e-8, 1e-3] 1/s |
| `log10_m` | -1.70 | -1.00 | $m$ in [0.02, 0.10] => $n$ in [10, 50] |

Notes on log-parameterization (this repo’s convention):

* `log10_gamd0x` means $\\log_{10}(gamd0x)$, so $gamd0x = 10^{log10\\_gamd0x}$.
* `log10_m` means $\\log_{10}(m)$, where $m = 1/nrsx$. So:
  - $m = 10^{log10\\_m}$
  - $nrsx = 1/m$

Staged calibration is often helpful (see Section 7.3):
* Stage 1: fix `log10_m` (or keep it narrow); calibrate the rest.
* Stage 2: widen `log10_m` and refine.

---

## 1. Context: what we are calibrating

### 1.1 The constitutive model in the solver (high level)
The LS-EVPFFT solver implements an **elasto-viscoplastic crystal plasticity** model in an FFT (Fast Fourier Transform) spectral solver.

Key ideas:

* **Elasticity**: the material responds elastically (reversible) at small stresses.
* **Plasticity by slip**: permanent deformation occurs by **dislocation slip** on specific crystallographic slip systems.
* **Viscoplasticity**: slip rates follow a *rate-dependent* power law, which provides a time scale and allows creep.
* **Hardening**: as plastic deformation accumulates, the slip resistance increases (strain hardening).

In this repo, the solid phase is treated as an FCC (face-centered cubic) crystal with {111}<110> slip systems.

### 1.2 Where parameters live

Solid material properties are provided via two files:

* `cuel.sx`: elastic constants (stiffness tensor)
* `cupl2.sx`: plasticity (slip kinetics + hardening)

Examples in the repo:

* `tests/creep_stress/cuel.sx`
* `tests/creep_stress/cupl2.sx`

Python-side defaults and calibration bounds live in:

* `lsevpfft_utils/config.py` (`SolidElastic` and `SolidPlastic`)

### 1.3 What “solid-only calibration” means here

“Solid-only calibration” means:

* the solid phase parameters (`cuel.sx`, `cupl2.sx`) are calibrated
* gas-phase parameters (e.g. `complgas`, gas kinetics for `igas=2`) are held fixed or disabled

Even if production simulations must include a gas phase for slip-free boundaries, it is still valuable to first establish a realistic solid prior for AM 316L.

---

## 2. Definitions: stress, slip, CRSS, and Taylor factor

### 2.0 Acronyms (defined once)

* **AM**: additive manufacturing
* **LPBF / SLM**: laser powder bed fusion / selective laser melting
* **CP**: crystal plasticity
* **CRSS**: critical resolved shear stress (slip resistance scale)
* **EVP**: elasto-viscoplastic
* **FFT**: fast Fourier transform

### 2.1 Resolved shear stress and slip resistance

In crystal plasticity, each slip system s feels a **resolved shear stress** tau_s, computed by projecting the local stress tensor onto the slip direction and slip-plane normal.

Plastic slip on that system occurs when |tau_s| is sufficiently large relative to a resistance g_s.

The resistance g_s is often called **CRSS** (Critical Resolved Shear Stress) in rate-independent crystal plasticity. In rate-dependent models, g_s plays the same role but slip can occur for tau_s < g_s at reduced rate.

### 2.2 Macroscopic yield stress vs microscopic slip resistance

In an FCC polycrystal, a common approximation connects macroscopic yield stress sigma_y to a characteristic microscopic slip resistance tau0 via:

    $sigma_y  \sim  M * tau0$

where:

* M is the **Taylor factor**.

#### What is the Taylor factor M?
The Taylor factor is a dimensionless number that maps macroscopic strain to the amount of crystallographic slip needed to accommodate that strain in a polycrystal.

Intuition:

* If the crystal orientations are random and deformation requires “a lot” of slip activity, M is larger.
* If the texture is favorable for slip, M can be smaller.

For FCC polycrystals under uniaxial tension, typical M values are roughly:

* M ~ 3.0 to 3.3 for many textures (often ~3.06 for random textures in classic Taylor theory, but AM textures can shift it).

#### Using AM yield strengths as a target

Open-access AM 316L tensile datasets frequently report 0.2% proof stresses in the
~400–650 MPa range depending on process/orientation/heat treatment (Section
9A.2 shows an open-access table).

Using $sigma_y \sim 450$–650 MPa and $M \sim 2.8$–3.3 gives:

    $tau0 \sim sigma_y / M \sim 140$–230 MPa

This estimate is not exact (texture, grain interactions, hardening, and rate sensitivity all matter), but it is an excellent guide for **reasonable priors**.

---

## 3. Definitions: rate-dependent slip kinetics (creep-capable)

### 3.1 Slip rate power law

In this repo/solver, slip rate is modeled by a power law of the form:

    $ \dot{\gamma} = \dot{\gamma}_0 \left(\frac{|\tau|}{g}\right)^n \mathrm{sign}(\tau) $

Repo parameter names:
* $\dot{\gamma}_0$ is `gamd0x` (1/s)
* $n$ is `nrsx` (dimensionless)

Parameters:

* `gamd0x` (1/s): reference shear rate
* `nrsx` (dimensionless): stress exponent
* tau: resolved shear stress
* g: current slip resistance

This form is widely used in EVPFFT and related crystal plasticity literature because it is simple, robust, and naturally supports creep.

### 3.2 Rate sensitivity m

It is common to define the **rate sensitivity**:

    $m = 1 / n$

In this repo’s parameterization:

* n == `nrsx`
* m == 1 / `nrsx`

Interpretation:

* small m (large n): nearly rate-independent behavior (sharp yielding)
* large m (small n): more rate-sensitive / more viscous behavior

Typical “reasonable” room-temperature values used in many FCC crystal plasticity models are around:

* n ~ 10 to 50   (m ~ 0.1 to 0.02)

We will recommend a baseline around n=20 (m=0.05) as a stable starting point.

### 3.3 Why (gamd0x, nrsx) matter for creep

At fixed stress ratio |tau/g|, increasing `gamd0x` increases slip rate linearly.

For the exponent, consider two cases at the same ratio r = |tau/g| < 1:

* If n is large, r^n is very small => creep is very slow.
* If n is smaller, r^n is less small => creep is faster.

So creep rate is very sensitive to `nrsx` when r is not close to 1.

Practical calibration implication:

* If your model is “too stiff” in creep (too little strain), it is often better to adjust `gamd0x` and/or `nrsx` than to push tau0xf unrealistically low.
### 3.4 Physical meaning of `gamd0x` (reference shear rate)

A useful way to read the slip-rate power law is:

$$
\dot{\gamma} = \dot{\gamma}_0 \left(\frac{|\tau|}{g}\right)^n \mathrm{sign}(\tau)
$$

When $|\tau| = g$, the magnitude of the slip rate is
$|\dot{\gamma}| = \dot{\gamma}_0$. In other words:

* `gamd0x` ($\dot{\gamma}_0$, 1/s) sets the **intrinsic timescale** of
  viscoplastic flow once a slip system is “at its resistance”.

More physically, $\dot{\gamma}_0$ is a **lumped mobility prefactor**. In a
dislocation-based picture, plastic shear rate on a slip system is often written
using the Orowan relation:

$$
\dot{\gamma} = \rho_m \, b \, v
$$

where $\rho_m$ is the mobile dislocation density, $b$ is the Burgers vector,
and $v$ is the average dislocation velocity. A power-law crystal-plasticity
model does not explicitly evolve $(\rho_m, v)$, so $\dot{\gamma}_0$ effectively
controls “how fast dislocations can move” at the reference stress ratio.

Practical implications (calibration + numerics):

* **Units:** 1/s, so $\dot{\gamma}_0 \, \Delta t$ is a characteristic shear
  increment scale per macrostep.
* **Creep/relaxation timescale:** increasing `gamd0x` speeds up creep roughly
  linearly at fixed $|\tau|/g$; decreasing it slows creep without changing the
  elastic stiffness.
* **Trade-offs / identifiability:** over limited windows, similar strain–time
  curves can be produced by trading higher `tau0xf` against higher `gamd0x`, or
  by changing `nrsx` (since $(|\tau|/g)^n$ is extremely sensitive when $|\tau|/g$
  is not close to 1).
* **Large-$\Delta t$ stability:** for very large time increments, avoid regimes
  where the implied slip increment per macrostep becomes enormous. A crude
  sanity check is to keep $\dot{\gamma}_0 \, \Delta t \ll 1$ for the padding
  layer transient. Note that the *actual* $|\dot{\gamma}|$ can be much larger
  than $\dot{\gamma}_0$ when $|\tau|/g > 1$.


---

## 4. Definitions: hardening and the (extended) Voce law

### 4.1 Why hardening is needed

Without hardening, once slip begins the material can flow too easily and would under-predict the rising stress vs plastic strain seen in real metals.

Hardening represents the increasing difficulty of slip due to dislocation interactions and evolving substructure.

### 4.2 Extended Voce hardening (Tome et al., 1984)

The classic Voce hardening has an exponential saturation form. Tome et al. (1984) proposed an **extended Voce** form that adds a linear-asymptote term so hardening can continue at large strain.

In a scalar form consistent with how these parameters are used in many EVPFFT implementations:

$$g(\Gamma) = \tau_0 + (\tau_1 + \theta_1\Gamma)\left(1-\exp(-\theta_0\Gamma/\tau_1)\right)$$

Where:

* Gamma is accumulated slip (or accumulated shear)
* tau0 is the initial slip resistance
* tau1 controls the “distance” to saturation for the exponential part
* thet0 is the initial hardening rate
* thet1 is the asymptotic hardening slope at large strain

Key limiting behaviors:

* At Gamma ~ 0:

$$dg/d\Gamma \approx \theta_0$$

* At large Gamma:

    $$g \approx \tau_0 + \tau_1 + \theta_1\Gamma$$

So tau1 shifts the saturation strength, while thet1 controls whether hardening continues linearly.

#### Practical tuning intuition

The characteristic accumulated slip to approach the exponential saturation is:

    Gamma* ~ tau1/thet0

Example: tau1=200 MPa and thet0=5000 MPa => Gamma* ~ 0.04.

### 4.3 Self and latent hardening

Slip systems can harden due to slip on themselves (self hardening) and due to slip on other systems (latent hardening).

In this repo the parameters are:

* `hselfx`: self hardening coefficient
* `hlatex`: latent hardening coefficient

For FCC, it is common to use:

* hself = 1.0
* hlate ~ 1.4

This is also the legacy/default choice in this repo and is consistent with many FCC CP parameterizations.

---

## 5. Elastic constants priors (solid): recommended values

### 5.1 Current repo values
In `tests/creep_stress/cuel.sx`, the solid elastic constants are:

* C11 = 206000 MPa
* C12 = 133000 MPa
* C44 = 119000 MPa

This is a standard cubic stiffness matrix representation.

### 5.2 Literature support

These exact constants are cited in the AM 316L anisotropy paper included in this repo:

* Charmi et al. (2021), Table 2: C11=206 GPa, C12=133 GPa, C44=119 GPa.

They are also consistent with the open-access NIST IR 83-1690 report by Ledbetter et al.

Note: some open-access CP papers use a very similar constant set (e.g., C11~204.6,
C12~137.7, C44~126.2 GPa). These differences are small (few percent) compared to
uncertainties in plasticity and creep parameters.

### 5.3 Recommendation

**Recommendation:** treat these elastic constants as fixed for initial calibrations.

Rationale:

* You stated they are verified in literature.
* For calibrating creep/hardening and rate effects, elastic constants are generally not the dominant uncertainty.
* Reducing calibration dimensionality improves Bayesian optimization robustness.

---

## 6. Solid plasticity priors: recommended starting values (AM 316L, sigma_y ~ 500 MPa)

This section proposes a “v2” starting point for calibration.

### 6.1 Repo parameter names (as in `cupl2.sx`)
From `tests/creep_stress/cupl2.sx`, the line:

    300.0  285.0   80.0   4000.0  200.0  0  tau0xf,tau0xb,tau1x,thet0,thet1,hpfac

maps to:

* tau0xf: initial slip resistance (forward)
* tau0xb: initial slip resistance (backward)
* tau1x: Voce saturation parameter
* thet0: initial hardening rate
* thet1: asymptotic hardening rate

And the kinetics line:

    1   12   20.00000000   1.00000000e-04   0.0  0

maps to:

* nrsx = 20
* gamd0x = 1e-4 1/s

### 6.2 Key experimental context: AM high yield
You confirmed the material is **high yield** with sigma_y ~ 500 MPa at room temperature.

Using the Taylor-factor estimate discussed earlier (sigma_y ~ M*tau0), that suggests tau0 on the order of ~150–200 MPa.

### 6.3 Recommended starting values

Recommended starting values for the solid phase:

| Parameter | Start | Units | Meaning |
|---|---:|---|---|
| tau0xf | 185 | MPa | initial slip resistance (forward) |
| rb | 0.95 | - | ratio tau0xb/tau0xf |
| tau0xb | 175.75 | MPa | initial slip resistance (backward) |
| tau1x | 200 | MPa | Voce saturation increment |
| thet0 | 5000 | MPa | initial hardening rate |
| thet1 | 300 | MPa | asymptotic hardening slope |
| gamd0x | 1e-4 | 1/s | reference shear rate |
| nrsx | 20 | - | stress exponent (m=0.05) |
| hselfx | 1.0 | - | self hardening |
| hlatex | 1.4 | - | latent hardening |

Notes:

* These values are intentionally conservative for AM 316L (sigma_y ~500 MPa) and aim to keep the simulation in a realistic strength range without being excessively strong.
* They are also close to parameter magnitudes used in AM 316L crystal plasticity studies in the repo (e.g., Somlo et al. 2022 uses slip resistance magnitudes ~O(200 MPa) for FCC 316L).

### 6.4 Why start with n=20 (m=0.05) if we’re unsure?
You said you “don’t know” what rate sensitivity to target. That is normal: rate sensitivity and reference shear rate are often hard to uniquely identify from limited tests.

I recommend:

* Start with a moderate, common baseline: m=0.05 (n=20)
* Allow calibration bounds to explore reasonable deviations

This strikes a balance:

* If we force n too large (m too small), creep can become unrealistically slow unless gamd0x is pushed large.
* If we force n too small (m too large), the model can become too viscous and “smear” the yield behavior.

---

## 7. Solid plasticity priors: recommended bounds for calibration

These bounds are designed for:

* room temperature
* AM 316L with sigma_y ~500 MPa
* robust Bayesian optimization (avoid exploring wildly unrealistic regions)

### 7.1 Proposed bounds

| Variable (calibration space) | Low | High | Notes |
|---|---:|---:|---|
| tau0xf (MPa) | 130 | 250 | corresponds roughly to sigma_y ~ 360–825 MPa for M~2.8–3.3 |
| rb (-) | 0.90 | 1.00 | keep forward/backward nearly symmetric |
| tau1x (MPa) | 80 | 300 | reasonable Voce saturation increment |
| thet0 (MPa) | 2000 | 9000 | covers a wide range of initial hardening |
| thet1 (MPa) | 50 | 800 | keep thet1 <= thet0 (enforced by code) |
| log10_gamd0x | -7.5 | -3.0 | gamd0x in [3e-8, 1e-3] 1/s |
| log10_m | -1.70 | -1.00 | m in [0.02, 0.10]  => n in [10, 50] |

Notes:

* Calibration is done in log-space for `gamd0x` and `m` to keep the optimizer
  well-scaled. Convert back via:
  - $gamd0x = 10^{log10\\_gamd0x}$
  - $m = 10^{log10\\_m}$ and $nrsx = 1/m$

### 7.1A Bounds justified by open-access AM 316L CP parameter sets

Two open-access CP papers provide useful numerical anchors for AM 316L at/near
room temperature:

| Source | AM yield (MPa) | tau0 (MPa) | tau_sat (MPa) | tau_sat - tau0 (MPa) | gamma0 (1/s) | m | n=1/m | q |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| Grilli et al. (2022) | (not reported) | 213 | 302 | 89 | 1e-7 | 0.1 | 10 | 1.4 |
| Wang et al. (2024) | 516 | 186 | 420 | 234 | 1e-3 | 0.05 | 20 | 1.4 |

How this supports our bounds (qualitatively):

* **tau0xf:** both papers use tau0 in ~186–213 MPa. Combined with the Taylor
  estimate $tau0 \\sim sigma_y/M$ (Section 2.2), a bound of 130–250 MPa is wide
  enough to cover texture/process variation without letting the optimizer spend
  too much time at unrealistically low/high strengths.
* **tau1x:** if we interpret tau1 as a saturation increment scale, then
  tau_sat - tau0 ranges from ~90 to ~235 MPa in these two fits. Hence allowing
  tau1 up to ~300 MPa is defensible for strong AM hardening.
* **gamd0x:** gamma0 differs by **4 orders of magnitude** across the two fits
  (1e-7 vs 1e-3), which strongly motivates calibrating this parameter in
  log-space over a wide range.
* **m/nrsx:** both fits sit in a compact band (m=0.05–0.1, i.e. n=10–20). We
  still allow modest expansion (m down to 0.02, n up to 50) for flexibility in
  creep calibration.

Important caveat:

* Both papers use a *different* saturation hardening law (with h0 and exponent
  a) than this solver’s extended Voce form. So you should not expect a
  one-to-one mapping for (thet0, thet1). We therefore bound (thet0, thet1) based
  on solver behavior and typical CP magnitudes (thousands of MPa), rather than
  directly matching h0.

### 7.2 Rationale for narrowing tau0xf vs legacy bounds
Earlier bounds in the repo allowed tau0xf up to 300 MPa, which implies sigma_y near or above 900 MPa for M~3, which is far above your confirmed ~500 MPa yield.

Allowing overly high tau0xf can:

* make early calibration runs systematically too stiff
* force optimizer to compensate by pushing `gamd0x`/`m` to extreme values
* increase the chance of mis-calibrating for numerical artifacts rather than physics

So for initial calibration windows, narrowing tau0xf is usually beneficial.

### 7.3 Identifiability note: (tau0xf, gamd0x, m)
Be aware that tau0xf, gamd0x, and m can trade off to produce similar creep strain over a limited time window.

Practical tip (staged calibration):

1. First stage: fix m (or keep it in a narrower range), calibrate tau0xf, Voce params, and gamd0x.
2. Second stage: widen m bounds and refine.

This is not required, but it can reduce Bayesian optimization sensitivity.

---

## 8. Recommended “v2” priors in repo terms (what to set in code)

In `lsevpfft_utils/config.py`, the following default values are recommended for the solid phase:

* DEFAULT_tau0xf = 185.0
* DEFAULT_rb = 0.95
* DEFAULT_tau1x = 200.0
* DEFAULT_thet0 = 5000.0
* DEFAULT_thet1 = 300.0
* DEFAULT_log10_gamd0x = log10(1e-4)
* DEFAULT_log10_m = log10(0.05)

And the `SolidPlastic.bounds()` method should match the bounds listed in Section 7.

---

## 9. References (from this repo)

1. **Charmi et al. (2021)**, “Mechanical anisotropy of additively manufactured stainless steel 316L: An experimental and numerical study”
   * File: `docs/papers/AM/anisotropy of AM 316L.pdf`
   * Used for: single-crystal elastic constants and parameter magnitude context.

2. **Tome et al. (1984)**, “The relation between macroscopic and microscopic strain hardening in f.c.c. polycrystals”
   * File: `docs/papers/hardening/hardening-fcc-polycrystals extended voce.pdf`
   * Used for: extended Voce form and interpretation.

3. **Somlo et al. (2022)**, “Anisotropic yield surfaces of additively manufactured metals simulated with crystal plasticity”
   * File: `docs/papers/AM/crystal plasticity sims anisotropic yield surfaces 316L.pdf`
   * Used for: order-of-magnitude slip resistance and hardening parameter magnitudes for AM 316L.

4. **Repo defaults and file formats**
   * `tests/creep_stress/cuel.sx`
   * `tests/creep_stress/cupl2.sx`
   * `lsevpfft_utils/config.py`

---

## 9A. Internet additions (open-access)

This section collects open-access priors found online that are directly useful
for AM 316L parameter magnitudes and bounds.

### 9A.1 Grilli et al. (2022) CPFEM for SLM 316L residual stress

**Source:** N. Grilli et al., *Computational Mechanics* (2022),
"Crystal plasticity model of residual stress in additive manufacturing using the
element elimination and reactivation method".

* DOI: https://doi.org/10.1007/s00466-021-02116-z
* Open-access PDF: https://link.springer.com/content/pdf/10.1007/s00466-021-02116-z.pdf

They use a power-law slip kinetics form equivalent to:

    gamma_dot = gamma0 * (|tau/tau_c|)^(1/m) * sign(tau)

So their exponent corresponds to:

    n = 1/m

Reported (Table 2, around room temperature T0~303 K):

* tau_c(T0) = 213 MPa
* tau_sat = 302 MPa
* gamma0 = 1e-7 1/s
* m = 0.1 (=> n=10)
* h0 = 3839 MPa, a = 2.5, q_r = 1.4
* Elastic constants: C11=204.6 GPa, C12=137.7 GPa, C44=126.2 GPa

Mapping notes to our priors:
* tau_c supports `tau0xf ~ 150–220 MPa`.
* tau_sat supports allowing `(tau0xf + tau1x)` to be at least ~300 MPa.
* gamma0 supports extending `gamd0x` bounds down to ~1e-7.
* m=0.1 supports `nrsx ~ 10` being plausible.

### 9A.2 Wang et al. (2024, MDPI Materials) AM yield + CP parameter magnitudes

**Source:** H. Wang, H.-W. Lee, M. T. Tran, D.-K. Kim, *Materials* (2024),
"Microstructure-Based Modeling of Deformation and Damage Behavior of Extruded and
Additively Manufactured 316L Stainless Steels".

* DOI: https://doi.org/10.3390/ma17102360
* Open-access (MDPI): https://www.mdpi.com/1996-1944/17/10/2360

Two especially useful tables:

* Table 2 reports their LPBF sample yield strength (0.2% proof): 516 MPa, and
  cites several other LPBF yield strengths: 590, 648, 412, 596, 583, 573 MPa.
* Table 3 reports CP hardening magnitudes for their AM sample (LA):
  - tau0 ~ 186 MPa
  - tau_sat ~ 420 MPa
  - h0 ~ 1300 MPa, a=2.5, q=1.4
  - slip kinetics: gamma0=0.001 1/s, m=0.05 (=> n=20)

Mapping notes to our priors:
* tau0 ~ 186 MPa strongly supports `tau0xf ~ 170–220 MPa`.
* tau_sat ~ 420 MPa motivates allowing `tau1x` high enough that
  `tau0xf+tau1x` can reach ~400–450 MPa for strong AM hardening.
* m=0.05 is exactly our default choice (`nrsx=20`), and gamma0=1e-3 supports the
  high end of our `gamd0x` bounds.

### 9A.3 NIST IR 83-1690 (Ledbetter et al.) monocrystal elastic constants for 316

**Source:** H. Ledbetter et al., NBS (NIST) Interagency Report 83-1690,
"Predicted single-crystal elastic constants and anisotropy ratios for
polycrystalline stainless steel 316" (open PDF).

* PDF: https://nvlpubs.nist.gov/nistpubs/Legacy/IR/nbsir83-1690.pdf

They report predicted monocrystalline constants:
* C11 = 206 GPa, C12 = 133 GPa, C44 = 119 GPa

This supports using the repo’s current elastic constants as the default.

---
## 9B. Gas/padding layer (“dummy gas” / “damper”) parameters and stability diagnostics

This repo uses a thin padding layer (“gas voxels”) on two free surfaces to
approximate slip-free / traction-free boundary conditions while retaining the
FFT solver’s periodic machinery.

Important: this “gas” is **not** physical air. It is a **numerical
regularization layer**, and its parameters are best treated as stability knobs.

### 9B.1 Two solver modes: `igas=1` vs `igas=2`

The solver supports two distinct padding-layer behaviors:

* `igas=1` (**dummy compliance operator**): the padding layer skips elastic and
  plastic state updates. A linear operator is applied that is controlled mainly
  by `complgas`.
* `igas=2` (**damper phase**): the padding layer uses **isotropic elasticity**
  plus the **same power-law slip kinetics** used for the solid. However,
  **hardening evolution is disabled** in the padding layer (it can dissipate,
  but it does not “lock up” by accumulating hardening).

### 9B.2 `igas=1`: what `complgas` does

In repo terms, `complgas` is a scalar compliance factor (a fraction of solid
compliance) used in the dummy-gas operator.

Qualitatively:

* smaller `complgas` -> a **softer / more compliant** padding layer (closer to
  traction-free), but typically **more ill-conditioned** numerically
* larger `complgas` -> a **stiffer** padding layer (less slip-free), usually
  more stable

Rule of thumb:

* reducing `complgas` from 0.01 to 0.0001 makes the padding layer about 100x
  “softer”, which can increase deformation-gradient updates in gas voxels and
  accelerate `det(F)` drift toward zero.

### 9B.3 `igas=2`: gas elasticity parameters (`E`, `nu`) and derived moduli (`G`, `K`)

For `igas=2`, the padding layer uses isotropic elasticity specified by:

* $E$ (Young’s modulus, MPa)
* $\nu$ (Poisson ratio, dimensionless)

The solver converts these to Lamé parameters:

$$
G = \mu = \frac{E}{2(1+\nu)}
$$

$$
\lambda = \frac{2G\nu}{1-2\nu}
$$

and the bulk modulus is:

$$
K = \lambda + \frac{2}{3}G = \frac{E}{3(1-2\nu)}
$$

Interpretation:

* $G$ (shear modulus) controls resistance to **shape change** / shear.
* $K$ (bulk modulus) controls resistance to **volume change**, i.e. how strongly
  the layer resists `det(F)` drifting away from 1.

Current repo defaults (intended as a stable starting point, not a physical gas):

* $E_\mathrm{gas} = 2000$ MPa
* $\nu_\mathrm{gas} = 0.49$

This corresponds to approximately:

* $G_\mathrm{gas} \approx 6.7 \times 10^2$ MPa
* $K_\mathrm{gas} \approx 3.3 \times 10^4$ MPa

Practical guidance:

* Choosing $E_\mathrm{gas}$ as roughly **0.1%–2%** of the metal’s $E$ keeps the
  padding layer much softer than the solid while remaining numerically
  meaningful.
* Choosing $\nu_\mathrm{gas}$ closer to 0.5 increases $K/G$ (more volume
  resistance), which can help prevent `det(F)` collapse, but excessive
  near-incompressibility can worsen conditioning. A sweep over
  $\nu_\mathrm{gas} \in [0.45, 0.49]$ is often more informative than fixing
  exactly 0.49.
* If you prefer to reason in terms of a target shear modulus, you can set:
  $E = 2G(1+\nu)$.

### 9B.4 `igas=2`: gas slip kinetics parameters (`tau0xf`, `rb`, `gamd0x`, `nrsx`)

The padding layer uses the same power-law kinetics:

$$
\dot{\gamma} = \dot{\gamma}_0 \left(\frac{|\tau|}{g}\right)^n \mathrm{sign}(\tau)
$$

For `igas=2`, hardening evolution is disabled in the padding layer, so $g$ is
effectively set by the *initial* slip resistance parameters (and any
initialization). The “hardening” parameters (`tau1x`, `thet0`, `thet1`, latent
and self hardening) are written/read but are not used dynamically in the
padding layer in this repo and should not be calibrated for the gas.

Parameter roles:

* `tau0xf` (MPa): forward initial slip resistance. Larger -> stiffer padding
  (less slip), smaller -> more compliant (more slip).
* `rb` (dimensionless): backward/forward resistance ratio,
  `tau0xb = rb * tau0xf`. Values < 1 make reverse slip easier than forward. For
  a purely numerical damper, `rb = 1` is usually the safest default unless
  cyclic loading suggests otherwise.
* `nrsx` (or $m = 1/n$): rate sensitivity. Smaller $n$ (larger $m$) gives a more
  viscous, smoother response and can reduce abrupt kinematic bursts.
* `gamd0x` ($\dot{\gamma}_0$, 1/s): reference shear rate / timescale (see
  Section 3.4).

Current repo defaults (numerical starting point):

* `gas tau0xf = 10` MPa
* `gas rb = 1.0`
* `gas nrsx = 5` (so $m = 0.2$)
* `gas gamd0x = 2e-5` 1/s

Suggested stability sweep ranges (non-physical numerical tuning):

* `tau0xf`: 5–50 MPa
* `nrsx`: 3–10 (so $m \approx 0.33$–0.1)
* `gamd0x`: $10^{-6}$–$10^{-4}$ 1/s
* `rb`: fix at 1.0 initially

Heuristic for large-$\Delta t$ runs:

* Choose `gamd0x` so that $\dot{\gamma}_0 \Delta t$ is not $O(1)$ (e.g.
  $\Delta t = 2000$ s and $\dot{\gamma}_0 = 10^{-5}$ 1/s gives
  $\dot{\gamma}_0 \Delta t = 0.02$). Remember the actual slip rate scales with
  $(|\tau|/g)^n$, so if $|\tau|$ exceeds $g$ the slip rate can still become
  large.

### 9B.5 Interpreting gas-layer stability diagnostics in `output.txt`

The padding-layer instabilities we see (drifting `det(F)`, NaNs, LAPACK
failures) are closely tied to a small set of quantities printed each macrostep:

* **`detFmin gas`**: minimum $\det(F)$ among gas voxels. $\det(F)$ is the local
  volume ratio. $\det(F) \le 0$ implies inversion (non-physical) and often
  triggers linear-algebra failures. A monotone decay of `detFmin gas` toward
  zero is a strong early-warning signal.
* **`max|c066mod| gas`**: maximum magnitude of the *modified* stiffness used by
  the constitutive update. In the solver, the reference stiffness is pushed
  forward by the current deformation gradient:

  $$
  C^\mathrm{mod}_{ijkl} =
  \frac{1}{\det F}\, C_{imkn}\, F_{jm}\, F_{ln}.
  $$

  The explicit $1/\det(F)$ factor means shrinking $\det(F)$ can blow up entries
  in $C^\mathrm{mod}$ even if the underlying material is “soft”. So very large
  `max|c066mod| gas` is primarily a *conditioning / kinematics* indicator, not a
  statement that “the gas got stiffer”.
* **`tr(L)*dt gas`**: min/max of the trace of the velocity-gradient increment.
  For an exponential update of the deformation gradient,

  $$
  \det(\exp(L\,\Delta t)) = \exp(\mathrm{tr}(L\,\Delta t)),
  $$

  so `tr(L)*dt` is a *log-volume increment*. Large negative values predict
  strong volumetric compression in one macrostep (and hence decreasing
  $\det(F)$); large positive values predict strong expansion.
* **`max||W||*dt gas`**: with $W=\frac{1}{2}(L-L^T)$, this is a bound on the
  rotation increment per macrostep. Very large rotations often coincide with
  poor Newton convergence and large deformation-gradient updates.

Practical use:

* If `detFmin gas` is decaying and `tr(L)*dt gas` has large negative excursions,
  increasing $K_\mathrm{gas}$ (via higher $\nu_\mathrm{gas}$ and/or $E_\mathrm{gas}$),
  reducing $\Delta t$, or making the damper kinetics more viscous (smaller $n$,
  smaller `gamd0x`, larger `tau0xf`) can help.
* If `max|c066mod| gas` is exploding while `detFmin gas` is falling, the issue
  is kinematic/volumetric collapse rather than “too much shear traction”; focus
  on volume control ($K_\mathrm{gas}$) and time-increment control.

---

## 9C. A more interpretable hardening reparameterization

The triplet `tau1x`, `thet0`, `thet1` is the correct solver input, but it is
not the most intuitive way to think about hardening. The reason is that these
three numbers do not each control one clean visual aspect of the constitutive
response. In practice:

* `tau1x` and `thet0` interact strongly,
* `thet1` is easiest to interpret as a *fraction* of the early hardening rate,
  not as an isolated MPa-valued slope, and
* the experimentally visible "knee" in a stress-strain or creep curve is more
  naturally described by a *slip-accumulation scale* than by the raw pair
  (`tau1x`, `thet0`).

For interpretation, it is therefore often better to think in terms of the
following combinations:

* a **saturation increment** or **saturation level**,
* a **slip-accumulation scale** for when the curve bends over,
* and a **tail ratio** that says how much hardening remains after the
  saturating part has mostly played out.

### 9C.1 Start from the scalar hardening form

Section 4.2 introduced the scalar law

$$
g(\Gamma)=\tau_0 + (\tau_1 + \theta_1 \Gamma)
\left(1-\exp\left(-\theta_0 \Gamma/\tau_1\right)\right)
$$

where:

* $\tau_0$ is the initial slip resistance,
* $\Gamma$ is **accumulated slip**,
* $\tau_1$ is the saturating hardening increment,
* $\theta_0$ is the initial hardening slope,
* $\theta_1$ is the late hardening slope.

That form is already enough to show the key structure:

* `tau1x` sets a **vertical scale** for the saturating part,
* `thet0` and `tau1x` together set a **horizontal scale** for how fast the knee
  is reached,
* `thet1` sets the **late-stage slope** after the initial saturation-like part.

The exact source code uses a discrete update in
`various_functions.f90::harden`, not this closed-form scalar expression, but
the same combinations control the response there too. The source first
accumulates slip through

$$
\Delta \Gamma = \sum_\alpha |\dot{\gamma}^\alpha| \Delta t,
\qquad
\Gamma_{n+1}=\Gamma_n+\Delta\Gamma,
$$

and stores the accumulated quantity in `gacumgr`. It then builds a mode-level
hardening factor containing the combinations

$$
\frac{\theta_0}{\tau_1},
\qquad
\tau_1,
\qquad
\theta_1,
$$

before distributing that hardening to each slip system through the latent
hardening matrix. So the interpretation below is not merely a documentation
convenience; it matches the actual parameter combinations used by the source.

### 9C.2 The more interpretable combinations

The simplest useful reparameterization is:

$$
\Delta \tau_{\mathrm{sat}} = \tau_1
$$

$$
\tau_{\mathrm{sat}} = \tau_0 + \tau_1
$$

$$
\Gamma_{\mathrm{sat}} = \frac{\tau_1}{\theta_0}
$$

$$
r_{\mathrm{tail}} = \frac{\theta_1}{\theta_0}
$$

with the inverse map

$$
\tau_1 = \Delta \tau_{\mathrm{sat}} = \tau_{\mathrm{sat}} - \tau_0
$$

$$
\theta_0 = \frac{\Delta \tau_{\mathrm{sat}}}{\Gamma_{\mathrm{sat}}}
$$

$$
\theta_1 = r_{\mathrm{tail}} \, \theta_0
$$

This is useful because each new quantity has a clear visual meaning:

* $\Delta \tau_{\mathrm{sat}}$ or $\tau_{\mathrm{sat}}$: the **vertical**
  hardening scale of the saturating part,
* $\Gamma_{\mathrm{sat}}$: the **horizontal** scale in accumulated slip for the
  knee/bend-over,
* $r_{\mathrm{tail}}$: the ratio of late hardening slope to early hardening
  slope.

The word "saturation" needs one careful qualification. In the present model,
`thet1` may be nonzero, so the hardening does not necessarily flatten to a true
constant plateau. Rather:

* if `thet1 = 0`, the Voce part approaches a true plateau,
* if `thet1 > 0`, the Voce part bends toward a plateau-like level and then the
  curve keeps rising with a residual tail slope.

So the most precise language is:

* `tau1x` is the **saturation increment of the Voce-like part**, and
* `tau_sat = tau0 + tau1` is the **plateau level that would be approached if the
  tail term were absent**.

### 9C.3 What "extra hardening" means here

It is better to say "extra resistance above the initial CRSS" than "extra
hardening above the base hardening". The base state is the initial slip
resistance $\tau_0$, not some other hardening slope.

So when we say that `tau1x` controls how much extra hardening is obtained before
the curve bends over, what we mean is:

* the slip resistance starts near $\tau_0$,
* the Voce-like saturating contribution raises it by roughly `tau1x`,
* and therefore the saturation-like level of the exponential part is roughly
  $\tau_0 + \tau_1$.

In the special case `thet1 = 0`, the large-$\Gamma$ limit is simply

$$
g(\Gamma \to \infty) \to \tau_0 + \tau_1.
$$

That is the cleanest place to see what `tau1x` does: it is the **difference
between the initial CRSS and the plateau CRSS** of the saturating part.

With `thet1 > 0`, the late response is instead

$$
g(\Gamma) \approx \tau_0 + \tau_1 + \theta_1 \Gamma
\quad \text{for large } \Gamma.
$$

So `tau1x` still sets the vertical offset of the saturation-like part, but the
total resistance can continue to rise beyond that because the tail term is still
adding hardening.

This is why `tau_sat` and `r_tail` are useful to think about together:

* `tau_sat` tells you where the exponential part wants to level off,
* `r_tail` tells you how strongly the law refuses to truly plateau.

### 9C.4 Is $\Gamma_{\mathrm{sat}}$ slip or strain?

It is **accumulated slip**, not macroscopic strain.

That distinction matters. In the solver, the quantity that evolves in the
hardening law is

$$
\Gamma = \int \sum_\alpha |\dot{\gamma}^\alpha| \, dt,
$$

which is a scalar sum of the absolute shear accumulated on all active slip
systems. This is the quantity stored in `gacumgr` and updated in
`harden(...)`. It is not the same as:

* axial engineering strain,
* equivalent von Mises plastic strain,
* or any single component of the strain tensor.

So the clean interpretation is:

* $\Gamma_{\mathrm{sat}} = \tau_1/\theta_0$ is the **accumulated-slip scale**
  over which the saturating part of the hardening law is approached.

Why, then, does it often behave like a "strain scale" in a macroscopic stress-
strain plot?

Because along a fixed loading path in a polycrystal, accumulated slip and
macroscopic plastic strain are often approximately proportional up to an order-
one geometric factor. In the simplest Taylor-style picture for uniaxial loading,
one can think very roughly in terms of

$$
\varepsilon_p \sim \frac{\Gamma}{M_\mathrm{eff}},
$$

with an effective Taylor factor $M_\mathrm{eff}$ of order 3 for FCC
polycrystals. This is only an approximation, but it explains the common
observation:

* increasing $\Gamma_{\mathrm{sat}}$ tends to push the macroscopic knee to
  larger strains,
* decreasing $\Gamma_{\mathrm{sat}}$ makes the curve bend over earlier.

So:

* in the **solver**, it is a slip-accumulation scale;
* in the **experimentally visible macro curve**, it often appears as an
  effective strain scale.

That is why the name `gamma_sat` or `Gamma_sat` is safer and more accurate than
calling it a strain directly.

### 9C.5 How the tail ratio implements late hardening

The most transparent way to understand the tail is through the large-$\Gamma$
limit:

$$
g(\Gamma) \approx \tau_0 + \tau_1 + \theta_1 \Gamma.
$$

So the late-stage slope is simply

$$
\frac{dg}{d\Gamma} \to \theta_1.
$$

At the start of hardening, the slope is instead

$$
\left.\frac{dg}{d\Gamma}\right|_{\Gamma=0} \approx \theta_0.
$$

That makes

$$
r_{\mathrm{tail}} = \frac{\theta_1}{\theta_0}
$$

an especially useful dimensionless summary:

* `r_tail = 0`: true plateau of the Voce-like part
* small `r_tail`: strong early hardening, weak late tail
* larger `r_tail`: less contrast between the early slope and the late slope

The source code enforces `thet1 <= thet0`, which means

$$
0 \le r_{\mathrm{tail}} \le 1
$$

for physically admissible input in this workflow.

In implementation terms:

* `thet1` appears directly in the discrete mode-level hardening factor `voce`,
* one term is proportional to $\Delta\Gamma \, \theta_1$,
* and additional terms proportional to `thet1` remain even after the
  exponential part has mostly decayed.

That is exactly how the code gives the model a persistent late-stage slope.

So `r_tail` is not merely a cosmetic ratio. It tells you, in one dimensionless
number, how much of the initial hardening rate survives into the late tail.

### 9C.6 Why this basis is easier to calibrate

The raw triplet (`tau1x`, `thet0`, `thet1`) mixes three effects:

* a vertical scale,
* a horizontal slip scale,
* a late-tail slope.

The reparameterized set separates those roles more cleanly:

* $\Delta \tau_{\mathrm{sat}}$ or $\tau_{\mathrm{sat}}$: "how high does the
  saturation-like part go?"
* $\Gamma_{\mathrm{sat}}$: "how much accumulated slip does it take to get
  there?"
* $r_{\mathrm{tail}}$: "how much hardening remains after the knee?"

This is especially helpful in optimization because the macroscopic curves often
identify these combinations more directly than the raw parameters. For example:

* two parameter sets with different (`tau1x`, `thet0`) can produce a similar
  knee location if `tau1x/thet0` is similar;
* two parameter sets with different `thet0` and `thet1` can produce a similar
  late-to-early slope contrast if `thet1/thet0` is similar.

So the reparameterized basis is often closer to what the data actually "sees".
It does not change the solver physics at all; it only changes how we think about
and potentially optimize the same hardening law.

---


## 10. Change log

* 2026-02-17: initial version written (solid priors for AM SS316L, room temperature).
* 2026-02-17: expanded doc with v2 priors + open-access sources (Grilli 2022, Wang 2024, NIST IR 83-1690); widened kinetics bounds and updated tau1 recommendation.
* 2026-02-26: documented gas/padding-layer parameters (igas=1 vs igas=2) and expanded the physical interpretation of gamd0x; added notes on stability diagnostics.
* 2026-03-20: added a solver-consistent interpretation of the hardening triplet (`tau1x`, `thet0`, `thet1`) in terms of saturation increment/level, accumulated-slip scale, and tail ratio.
