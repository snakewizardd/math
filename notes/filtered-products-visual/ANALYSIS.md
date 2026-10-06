# Filtered products: a read-through of the paper

Paper: `preprints/Filtered-products-and-boundary-preserving-compression-in-complex-cobordism-September-25-2026/`
(source in `build/sections/`). Companion:
`preprints/Radius-of-comparison-equals-half-the-mean-dimension-September-25-2026/`.

This is a careful reading, not a referee report or a formal check. Every
inequality in Sections 2 to 5 that can be checked by hand was checked.
The cited external results (Mathew, Burklund, Basterra-Mandell, Lurie,
Freudenthal) were taken as given. Lemma "Thickened complements" and
Lemma "Compatible local duality" in Section 3 were read for their
statements and how they are used, not checked line by line.

## 1. What the paper is for

The companion proves `rc(C(X) rtimes_h Z) = (1/2) mdim(X, h)` for minimal
`h`. Niu had already proved `<=` (companion intro, eq.
`known-upper-bound`). The open direction is `>=`. This paper supplies
the one geometric tool it needs:

> If a cube-valued map has a sparse bundle embedding of complex rank
> strictly below `2pn`, then some large power of it can be moved so that
> a fixed fraction `a/8` of the slots lie on a cube face. Every slot
> already on a face stays where it is.

The companion turns "slots on a face" into "lower covering order". It
uses a comparison bound `r` to build the bundle embedding, and it gets
a contradiction whenever a cover rate `d > 2r`.

## 2. Pipeline

```
bundle embedding, rank K < 2pn          (Section 5, Lemma torus-mu)
   |  Chern-class computation in MU, square-zero classes on tori
   v
MU point test vanishes on Z_{an}         omega_MU(f, Z_an) = 0
   |  Theorem 1.1 (filtered products) + Lemma batch-vanishing
   v
sphere point test of f^{x ell} vanishes on Z_H,  H = floor(a n ell / 4)
   |  padding (Theorem boundary-compression): old output + one new batch
   v
stable avoidance (Section 3): a real map within 1/2 that avoids Z_{jH-1}
   |  boundary restoration
   v
>= jH face slots, old face slots exact
   |  finite-parameter reduction (Section 5): compact Y -> polyhedron P
   v
Theorem 1.3 on any compact metrizable Y
```

## 3. Section by section

### Section 2: filtered products

- **Compatible Adams lift.** The difficulty is to get one lift through
  `I_R^{^n}` on the whole union `U_{ell,n}` that agrees on every overlap.
  The paper builds it in `Fun(Delta, Sp)` using the cobar object
  `R^{^(j+1)}`. The step that makes this work is the simplex-slicing
  lemma: one simplex coordinate is shared out among the `ell` factors so
  that at most `j` of them are away from a vertex. In degree `j < n`, at
  least one of the `>= n` factors in `G` is therefore at a vertex and
  collapses. The restriction to the punctured cube uses
  `S -> [|S|-1]`, which is a functor into `Delta_{<= n-1}`. Only
  restriction of a null cone is used, so no cofinality argument is
  needed. I agree with this.
- **Slicing lemma, checked numerically.** `slicing_check.py` tests the
  three claims on 5,400 points, including points with breaks exactly on
  slice endpoints, and on 221,400 (point, monotone map) pairs:
  - the values are simplex points (error `9e-16`);
  - there are at most `j` non-vertices, and the sample reaches exactly `j`;
  - naturality holds (error `4e-15`).
  See `slicing_check_log.txt`.
- **Getting from Adams filtration to vanishing.** At a fixed prime,
  Mathew's sublinear `f_R(k) = o(k)` beats `n_ell = ceil(ell/2)`. The
  cell bound is `D_ell = (e-d) ell`. At large primes, `I_BP` is
  `(2p-3)`-connective, so connectivity alone kills the class once
  `(2p-3)/2 > e-d`. This leaves finitely many primes, each with an
  eventual threshold. The local-to-global step at the end is standard
  for finitely generated groups.
- **Small notes.**
  - The large-prime lemma uses `BP` as a unit-preserving retract of
    `MU_(p)`. Quillen's idempotent already gives this, so the E_4
    splitting cited is more than needed (not wrong).
  - `H_*(BP; F_p)` starts at `xi_1` in degree `2p-2`, so the stated
    connectivity is sharp.
  - In the application (`E = Q_c(f, T(t,n))`, `d = tn`), `e - d` is
    exactly `dim P`. The prime cutoff is therefore `(2p-3)/2 > dim P`,
    and the nilpotence input is evaluated at `k_ell = dim P * ell + 1`.
    The paper does not say this; it is a useful way to read the
    dependence.

### Section 3: point tests, stable avoidance, restoration

- **Point test.** This is the relative class of "`f(x) - z` near 0" on
  `P x Z` modulo the far set. Products of tests are smash products of
  tests (Lemma restrictions-and-products). That is what lets Section 2
  apply.
- **Stable avoidance.** Simplices of `P` get nested boxes of radius
  `2c` to `3c`, with `c = 1/32`. Local duality makes the obstruction to
  lifting into `box \ Z` a pullback of the one global point test. One
  global nullhomotopy gives a compatible stable lift. Freudenthal plus
  `(q-2)`-connected complements then gives actual maps when
  `dim P <= 2q-5`. The displacement is below `4c = 1/8`, which is well
  under the `1/2` claimed. I checked the box arithmetic.
- **Boundary restoration.** This is elementary, and simulated in
  `restoration_padding_sim.py`.
  - The cutoff at `3/4` plus a displacement `< 1/2` means the
    interpolation cannot create a new zero.
  - Dividing by `max(delta, |.|)` is continuous at 0.
  - An old face slot has norm 1, which is at least `delta`, so it is
    fixed exactly.

### Section 4: batch vanishing and padding

- **Batch vanishing.** A group outside `Z_{an}` costs `> an` non-zero
  slots, so at most `H/(an) <= ell/4` groups can be outside. The paper
  needs `ceil(ell/2)` inside and in fact gets `> 3 ell/4`. This is
  slack, not an error.
- **Padding.** This is the idea in the paper that is not routine. A
  one-shot argument on `P^u` fails because `dim P^u = u dim P` grows
  faster than the stable range `~ 2tun` once `dim P >= 2tn`. The
  `restoration_padding_sim.py` ledger shows the crossover exactly at
  `2tn`. The fix is to replace the old parameters by the old **output**
  `Q_s(t, N_{j-1})`, whose dimension is at most `tN_{j-1}`. Only one
  fresh batch `P^ell` is added. The old factor is glued in by the
  observation `||y - z_o|| >= 1` on `A_j`, which is right because a
  boundary slot of `y` is paired with a zero parameter.
- **Arithmetic.**
  - `t(jH-1) < (a/4) D_j` and `dim W_j < (3/2) D_j` both hold.
  - `2q_j - 5 > (3/2) D_j` needs `D_j >= 20`, which is the reason for
    the `20` in the choice of `B`.
  - `H/b >= 3a/16` needs `ab >= 16`, which is the reason for that
    condition.
  - `j >= 2B` gives the factor `2/3`.

  The sim checks the ledger for `dim P` from 10 to 10,000.
- **Slack in `a`.** Section 4 only needs `a < 1/2`. The `1/16` is either
  inherited from the companion's budget or chosen conservatively. In
  the sim, the fraction `a/8` is reached around `j ~ B`, not `2B`.

### Section 5: torus criterion and compact inputs

- **Framed torus.** `(T^2)^p` sits in `R^{2p+1}`, built by iterated
  tube boundaries. The line bundles are pulled back from `T^2 -> S^2`,
  and `<prod c_1, [L]> = 1`. This is standard.
- **Cobordism criterion.**
  - Each `x_i = c_1^{MU}(lambda_i)` comes from a 2-torus, and
    `MU^4(T^2) = 0`, so `x_i^2 = 0`.
  - If `E` embeds in `C^K` with `K < 2r`, the complement `F` has rank
    `K - r < r`.
  - Then `c(F) = prod (1 - x_i v)`, and the `v^r` coefficient gives
    `prod x_i = 0` on the tube.
  - Multiplying by the framed Thom class of `L^n` gives the `MU` point
    class, checked against `H Z`.
  - The degrees match: `2pn + n = (2p+1) n = M`.

  The Tietze extension plus compactness argument, which extends the
  embedding from the locus to a thin tube, is correct as written. The
  reason to work in `MU` rather than `H Z` is that Section 2 needs an
  `MU`-null input. This is the stated improvement over the
  Hirshberg-Phillips Chern obstruction.
- **Finite-parameter reduction.**
  - Stone-Weierstrass is applied on `Y x Z`, with the coordinates of `f`
    kept as the first coordinates of `j_0`, so that `f = f_P o j`
    exactly.
  - The "bad set" is compact and its projection misses `C = j_0(Y)`,
    because the bundle depends only on `f_0(u) - z`.
  - `P` is the union of the small cubes meeting `C`.

  I agree with this argument, and it is the cleanest part of the paper.
- **Sharpness remark.** For the identity on `C(2p+1, u)` with `K = 2p`,
  the argument goes as follows:
  - Any `G` with at least one face slot maps into `\partial B`.
  - On `\partial B`, the straight-line homotopy keeps a fixed face slot
    on its face, so `G|_{\partial B}` has degree 1.
  - But `G` extends that restriction over the ball, which forces
    degree 0.

  This is correct. It shows `K <= 2pn` cannot replace `K < 2pn` in
  general. It says nothing about the best fraction under the strict
  inequality, and the paper says so.

## 4. Where the 1/2 comes from

A slot has real dimension `t = 2p+1`. The torus `L` carries complex rank
`p`. The criterion allows a target rank `< 2pn`, which is twice the
source rank. In the companion, a cover rate `d` gives about
`g ~ d k/(t+1)` slots per block. The comparison bound pays
`R ~ (2 - 10a) p g`, so the excess trace is `~ p d /(t+1)`. Choosing
`p` large makes `p/(2p+2)` as close to `1/2` as needed, and that
contradicts `d > 2r`.

So the `1/2` is "rank `p` inside real dimension `2p+1`, with target
budget twice the source rank". It is not "pin half of the shift
parameters". The first-pass sim's `d/2` ledger should be read as a
picture only.

## 5. Assessment

- **Correctness.** I found no error. Every quantitative inequality I
  checked holds, often with room to spare: `a < 1/2` would do in
  Section 4, the batch step has `3 ell/4` against the `ell/2` needed,
  and the displacement is `1/8` against `1/2`.
- **Load-bearing external input.** The sublinear quantitative
  nilpotence of Mathew and Burklund. If that input failed, Theorem 1.1
  would fail with it, and so would the paper's main theorem.
- **Effectiveness.** The `ell` from Theorem 1.1 is not explicit, so `u`
  is not explicit. The companion only needs existence, so this costs
  nothing there. A clean concrete project would be an explicit `ell`
  for one small case, for example `P` a point, `p = 1`, `n = 1`.
- **Novel parts.** The slicing lemma, the padding trick, and the
  `MU` square-zero computation on a tube. The rest is careful assembly
  of standard tools.
- **Points to watch in the companion.** The companion has to produce
  the bundle embedding *continuously in the sparse parameter `z`*,
  including as its support changes. That is a strong hypothesis
  (companion Prop. `comparison-bundle`), and it is where the
  operator-algebra content sits. I did not audit it here.

## 6. Files in this folder

| file | what it does |
|---|---|
| `cube_sim.py`, `sim_log.txt` | first pass (Grok): face push picture, Remark 5.4 picture |
| `restoration_padding_sim.py`, `restoration_padding_log.txt` | boundary restoration lemma as written, padding ledger, the 1/2 table |
| `slicing_check.py`, `slicing_check_log.txt` | numeric check of the simplex-slicing lemma |
| `README.md` | the note, with corrections to the first pass |
| `ANALYSIS.md` | this file |
