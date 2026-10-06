#!/usr/bin/env python3
"""Second-pass sim for the filtered-products compression paper.

cube_sim.py draws the face push. This script runs two pieces of the paper
as they are written, and prints where the 1/2 actually comes from.

  1. Boundary restoration (Lemma "Boundary restoration", Section 3).
     Input: f, and a nearby g that avoids the sparse arrangement Z_s.
     Output: G with at least s+1 boundary slots, equal to f on every
     original boundary slot, continuous even where a slot of g is small.
     Toy: K = 2-torus, slots in [-1,1]^2. Each slot of f has stable zeros,
     so no small move clears every slot, but a generic move clears all
     but one at each point: g avoids Z_{n-2}, whose codimension is
     t*(n-(n-2)) = 4 > dim K = 2.

  2. Padding (Theorem "Boundary compression", Section 4).
     Exact integer ledger: b, H, B, D_j, q_j, dim W_j, the stable range
     2q_j - 5 > dim W_j, and the guaranteed boundary fraction
     jH / (b(B+j)). Compared with a one-shot attempt on P^u, which fails
     once dim P is large compared with t*n.

  3. The 1/2. A slot of real dimension t = 2p+1 carries complex rank p,
     the target rank must stay below 2pn, and the companion compares a
     cover rate d against p*d/(t+1). p/(2p+2) -> 1/2.

The cobordism and duality steps (Sections 2-3) are not simulated. The
cutoff for ell in Lemma "Batch vanishing" is not effective; the ledger
uses the smallest ell allowed by the arithmetic (a*b >= 16) only.
"""

from pathlib import Path
from math import floor
import numpy as np
import matplotlib.pyplot as plt

OUT = Path(__file__).resolve().parent
lines = []


def log(s=""):
    lines.append(s)
    print(s)


# ---------------------------------------------------------------- part 1
T_DIM, N_SLOTS = 2, 4
S = N_SLOTS - 2
ALPHA = (0.3, 1.4, 2.5)
BETA = (0.9, 2.0, 3.1)
EPS = 0.05


def toy_f(th1, th2):
    out = np.empty(th1.shape + (N_SLOTS, T_DIM))
    for j in range(3):
        out[..., j, 0] = 0.6 * np.sin(th1 - ALPHA[j])
        out[..., j, 1] = 0.6 * np.sin(th2 - BETA[j])
    out[..., 3, 0] = np.clip(1.5 * np.cos(th1), -1, 1)
    out[..., 3, 1] = 0.4 * np.sin(th2)
    return out


def toy_g(th1, th2, f):
    h = np.empty_like(f)
    for j in range(N_SLOTS):
        h[..., j, 0] = np.cos(th1 + 2 * th2 + j)
        h[..., j, 1] = np.sin(2 * th1 - th2 + 0.5 * j)
    return f + EPS * h


def slot_norm(v):
    return np.max(np.abs(v), axis=-1)


def restoration(f, g, s):
    """The proof of the Boundary restoration lemma, line by line."""
    nf = slot_norm(f)
    chi = np.clip((nf - 0.75) / 0.25, 0.0, 1.0)[..., None]
    gp = (1 - chi) * g + chi * f
    ng = slot_norm(gp)
    kth = np.sort(ng, axis=-1)[..., ::-1][..., s]
    delta = 0.5 * float(kth.min())
    G = gp / np.maximum(delta, ng)[..., None]
    return G, delta


def radial(g):
    """Face push with no cutoff, as in cube_sim.radial_to_face."""
    return g / np.maximum(slot_norm(g), 1e-300)[..., None]


def max_jump(G):
    j1 = slot_norm(G - np.roll(G, 1, axis=0)).max()
    j2 = slot_norm(G - np.roll(G, 1, axis=1)).max()
    return float(max(j1, j2))


def grid(m):
    th = np.linspace(0, 2 * np.pi, m, endpoint=False)
    return np.meshgrid(th, th, indexing="ij")


log("part 1: boundary restoration on the 2-torus, t=2, n=4, s=n-2")
TH1, TH2 = grid(400)
f = toy_f(TH1, TH2)
g = toy_g(TH1, TH2, f)
G, delta = restoration(f, g, S)
tol = 1e-9
on_f = slot_norm(f) >= 1 - tol
on_G = slot_norm(G) >= 1 - tol
sorted_g = np.sort(slot_norm(g), axis=-1)
log(f"  |g - f| max = {slot_norm(g - f).max():.3f} (lemma needs < 1/2)")
log(f"  min over K of smallest slot norm of g = {sorted_g[..., 0].min():.4f}"
    " (stable zeros, grid-limited)")
log(f"  min over K of 2nd smallest slot norm of g = {sorted_g[..., 1].min():.4f}"
    f" (> 0, so g avoids Z_{S})")
log(f"  delta = {delta:.4f}")
log(f"  min #boundary slots of G = {on_G.sum(axis=-1).min()} (claim >= {S + 1})")
log(f"  original boundary slots preserved exactly = "
    f"{bool(np.all(np.abs(G - f)[on_f] == 0))}")
log(f"  fraction of K with all {N_SLOTS} slots on boundary = "
    f"{float(np.mean(on_G.sum(axis=-1) == N_SLOTS)):.3f}")
jumps = []
for m in (100, 200, 400, 800):
    a1, a2 = grid(m)
    fm = toy_f(a1, a2)
    gm = toy_g(a1, a2, fm)
    Gm, _ = restoration(fm, gm, S)
    jumps.append((m, max_jump(Gm), max_jump(radial(gm))))
log("  max neighbour jump (grid m: restoration, radial push without cutoff)")
for m, jr, jp in jumps:
    log(f"    m={m:4d}: {jr:.4f}  {jp:.4f}")

fig, axes = plt.subplots(1, 3, figsize=(13.2, 4.2))
ext = [0, 2 * np.pi, 0, 2 * np.pi]
ax = axes[0]
im = ax.imshow(on_G.sum(axis=-1).T, origin="lower", extent=ext,
               cmap="Greens", vmin=0, vmax=N_SLOTS)
for j in range(N_SLOTS):
    small = slot_norm(g[..., j, :]) < 0.02
    ax.scatter(TH1[small], TH2[small], s=2, label=f"g slot {j} near 0")
ax.set_title(f"boundary slots of G (min {on_G.sum(axis=-1).min()})")
ax.set_xlabel("theta1")
ax.set_ylabel("theta2")
ax.legend(fontsize=6, loc="upper right", markerscale=4)
fig.colorbar(im, ax=ax, fraction=0.046)
ax = axes[1]
im = ax.imshow(slot_norm(G[..., 0, :]).T, origin="lower", extent=ext,
               cmap="viridis", vmin=0, vmax=1)
ax.set_title("|G slot 0|: interior only near its own zeros")
ax.set_xlabel(f"delta = {delta:.3f}")
fig.colorbar(im, ax=ax, fraction=0.046)
ax = axes[2]
ms = [r[0] for r in jumps]
ax.loglog(ms, [r[1] for r in jumps], marker="o", label="restoration (max{delta,|.|})")
ax.loglog(ms, [r[2] for r in jumps], marker="s", label="radial push, no cutoff")
ax.set_xlabel("grid points per side")
ax.set_ylabel("largest jump between neighbours")
ax.set_title("continuity proxy")
ax.legend(fontsize=7)
fig.suptitle("Boundary restoration: s+1 boundary slots, old boundary kept, no blow-up at 0")
fig.tight_layout()
fig.savefig(OUT / "restoration.png", dpi=140)
plt.close()


# ---------------------------------------------------------------- part 2
def padding_ledger(t, n, a, dP, J):
    ell = 1
    while a * n * ell < 16:
        ell += 1
    b = n * ell
    H = floor(a * b / 4)
    B = 1
    while t * b * (B + 1) < max(20, 2 * ell * dP):
        B += 1
    rows = []
    for j in range(1, J + 1):
        D = t * b * (B + j)
        q = D - t * (j * H - 1)
        dimW = t * b * (B + j - 1) + ell * dP
        rows.append(dict(j=j, u=ell * (B + j), D=D, q=q, dimW=dimW,
                         ok=2 * q - 5 > dimW, frac=j * H / (b * (B + j))))
    return dict(ell=ell, b=b, H=H, B=B, rows=rows)


def one_shot(t, n, a, dP, u):
    N = n * u
    Hn = floor(a * N / 4)
    D = t * N
    q = D - t * (Hn - 1)
    return 2 * q - 5 - u * dP, D


t, n, a = 3, 4, 0.06
log()
log(f"part 2: padding ledger, t={t} (p=1), n={n}, a={a}")
dPs = [10, 100, 1000, 10000]
ledgers = {}
for dP in dPs:
    L = padding_ledger(t, n, a, dP, J=4000)
    ledgers[dP] = L
    rows = L["rows"]
    first = next(r for r in rows if r["frac"] >= a / 8)
    log(f"  dim P={dP:5d}: ell={L['ell']} b={L['b']} H={L['H']} B={L['B']}"
        f"  stable range at every j={all(r['ok'] for r in rows)}"
        f"  frac>=a/8 from j={first['j']} (u={first['u']}), 2B={2 * L['B']}"
        f"  frac at j=4000: {rows[-1]['frac']:.5f}")
log(f"  a/8={a / 8:.5f}  a/4={a / 4:.5f}  H/b={ledgers[10]['H'] / ledgers[10]['b']:.5f}")
log("  one-shot on P^u (u=10^4): stable-range slack 2q-5-dim, per ambient dim")
sweep = np.unique(np.logspace(0, 4.3, 60).astype(int))
naive = []
for dP in sweep:
    slack, D = one_shot(t, n, a, int(dP), 10_000)
    naive.append(slack / D)
naive = np.array(naive)
cut = sweep[np.argmax(naive <= 0)]
log(f"    first failing dim P = {cut}  (2tn = {2 * t * n})")
for dP in dPs:
    slack, D = one_shot(t, n, a, dP, 10_000)
    log(f"    dim P={dP:5d}: {slack / D:+.3f}")
padded = []
for dP in sweep:
    rows = padding_ledger(t, n, a, int(dP), J=50)["rows"]
    padded.append(min((2 * r["q"] - 5 - r["dimW"]) / r["D"] for r in rows))

fig, axes = plt.subplots(1, 2, figsize=(11.6, 4.2))
ax = axes[0]
for dP in dPs:
    rows = ledgers[dP]["rows"]
    ax.semilogx([r["u"] for r in rows], [r["frac"] for r in rows], label=f"dim P = {dP}")
ax.axhline(a / 8, color="k", ls="--", lw=1, label="a/8 (theorem)")
ax.axhline(a / 4, color="gray", ls=":", lw=1, label="a/4")
ax.set_xlabel("copies u")
ax.set_ylabel("guaranteed boundary fraction jH / (b(B+j))")
ax.set_title("padding: dimension costs delay, not failure")
ax.legend(fontsize=7)
ax = axes[1]
ax.semilogx(sweep, naive, marker=".", label="one shot on P^u")
ax.semilogx(sweep, padded, marker=".", label="padded W_j, worst j")
ax.axhline(0, color="k", lw=1)
ax.axvline(2 * t * n, color="gray", ls=":", lw=1)
ax.set_xlabel("dim P")
ax.set_ylabel("(2q - 5 - dim) / ambient dim")
ax.set_title("stable avoidance range")
ax.legend(fontsize=7)
fig.tight_layout()
fig.savefig(OUT / "padding_ledger.png", dpi=140)
plt.close()


# ---------------------------------------------------------------- part 3
log()
log("part 3: where the 1/2 comes from")
log("  p  t=2p+1  rank/slot  target bound  p/(t+1)")
for p in (1, 2, 4, 8, 16, 64):
    tt = 2 * p + 1
    log(f"  {p:<3d}{tt:<8d}{p:<11d}<2p per slot  {p / (tt + 1):.4f}")
log("  companion: pick p with p*d/(t+1) > r, so r < d/2 is contradicted")

(OUT / "restoration_padding_log.txt").write_text("\n".join(lines) + "\n")
