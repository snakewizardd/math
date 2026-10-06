# Filtered products, first visual

Notes on OpenAI Math Release preprint
`Filtered products and boundary-preserving compression in complex cobordism`
(25 September 2026), plus a toy simulation of the geometric move.

This is not a formalization, not a referee report, and not a proof of the
companion radius-of-comparison claim. It is the first pass at showing the
move the lemma is for.

## What the paper is

An informal proof of a specialist compression lemma. Author line is OpenAI.
The folder upstream has the PDF and a build directory, not a Lean certificate
next to this manuscript.

Theorem 1.1: if a map is MU-null on a finite pointed subcomplex, then for all
large smash powers it is stably null on the union of factors with a fixed
positive proportion of coordinates in that subcomplex. One nullhomotopy on
the union, intersections included.

Theorem 1.3: under the bundle-rank window `pn <= K < 2pn`, a map from any
compact metrizable space into a product of cubes can be nudged, in
arbitrarily large Cartesian powers, so a positive fraction of slots land on
a face, and every slot already on a face stays fixed. The extracted fraction
is `a/8` with `0 < a < 1/16`. Remark 5.4 shows the inequality is sharp: the
identity map of a cube has no such continuous nudge.

The companion preprint uses Theorem 1.3 toward
`rc(C(X) rtimes_h Z) = (1/2) mdim(X, h)`. This note does not prove that.

## What the sim is

`cube_sim.py` does not run the cobordism detector. It runs the cube move.

- Filled square: boundary degree is 1. Radial face push agrees with the
  identity on the boundary and is undefined at the center. A continuous
  attempt that stays defined there freezes the boundary and leaves the
  interior inside. That is the Remark 5.4 obstruction.
- Thin loop: image already misses the center, so the same push is
  continuous. Face rate went from about 0.17 to 1.00, old face slots
  frozen, and the rate stayed 1 across 1 to 16 copies. This is the easy
  case. The paper only claims a fraction `a/8` after enough copies.
- Half-parameter chart and shift window: a ledger, not a measurement.
  A shift that reveals `d` new reals per iterate has mean dimension `d`.
  The rank window needs a complementary copy of the source rank before the
  Chern product is forced to vanish, so the pin budget drawn here is `d/2`.
  A full cube pins 0. The green row is that ledger. It is not a mechanism
  that keeps a full cube thin.

## How to read the figures

`compression_move.png` is the obstruction and the easy case.
`half_parameter.png` plots the ledger. It does not prove the factor 1/2.
`shift_ledger.png` marks the half the lemma is allowed to pin. Pinning is
the output. Thinness is the hypothesis.

## Run

```
python3 cube_sim.py
```

Writes the three figures and `sim_log.txt` next to the script. Needs
numpy and matplotlib.

Run log from the container pass is in `sim_log.txt`.

## Second pass: restoration and padding

`restoration_padding_sim.py` runs two steps of the paper as written. The
first pass drew the face push. Run log is in `restoration_padding_log.txt`.

- Boundary restoration (Section 3). The paper does not push `v/|v|`. It
  pushes `v / max(delta, |v|)`, with `delta` below the `(s+1)`-st largest
  slot norm. That is continuous at 0, so slots near 0 stay inside.
  Toy: a 2-torus family, 4 slots in `[-1,1]^2`. Every slot has stable
  zeros, but no two zeros meet, so a 0.05 move avoids `Z_2`. Result: at
  least 3 boundary slots everywhere, and every old boundary slot kept
  bit for bit. Largest neighbour jump halves as the grid doubles. The
  plain radial push jumps by 2 at every grid size.
- Padding (Section 4). This is an integer ledger of the proof, not a
  measurement. It checks `2q_j - 5 > dim W_j` at every step and tracks
  `jH / (b(B+j))`. Every `dim P` from 10 to 10000 passes. A larger
  `dim P` only increases `B`, so `a/8` takes more copies. A single pass
  on `P^u` leaves the stable range once `dim P >= 2tn`. Taking the old
  output instead of the old parameters is what removes the dimension
  bound. The `ell` used is the smallest one with `ab >= 16`. The
  nilpotence cutoff in the batch-vanishing lemma is not computed.
- The 1/2. A slot of real dimension `t = 2p+1` carries complex rank `p`.
  The target has to stay below `2pn`. The companion compares a cover
  rate `d` with `p d / (t+1)`, and `p/(2p+2) -> 1/2`. That is where the
  half comes from. It does not come from pinning half of the shift
  parameters.

Corrections to the first pass:

- `thin_map` slot 2 is discontinuous at `cos(theta) = 0`. It jumps from
  `(1, 0.3, -0.2)` to about `(0, 0.35, 0.2)`. Slots 0 and 1 stay away
  from 0, so their face push was never at risk. So face rate 1.00 does
  not test the lemma.
- `power_rate` takes independent random phases. Face rate 1 for every
  `u` follows from the slot shape alone. It is not a power effect.
- The `d/2` pin ledger and `shift_ledger.png` are a picture, not
  something the paper states. See the 1/2 bullet above.
- The Remark 5.4 obstruction is a degree argument on `C(t,u)` for every
  `u`. One failed continuous attempt in `[-1,1]^2` illustrates it but
  does not prove it.

Run `python3 restoration_padding_sim.py` to write `restoration.png`,
`padding_ledger.png` and the log.
