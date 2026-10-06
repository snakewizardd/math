#!/usr/bin/env python3
"""Toy sim of the paper's geometric payload.

The cobordism detector is not simulated. What is simulated is the move the
detector is for, and the one case where the move is impossible:

  slots in cubes, boundary = touching a face, preservation = do not edit a
  slot that already touches a face.

Theorem shape being illustrated:
  if the map is thinner than twice its parameter rank, large powers can be
  nudged so a positive fraction of slots hit a face, old face-slots frozen.
  at equality (filled cube, identity) no continuous nudge does this.
"""

from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

OUT = Path(__file__).resolve().parent
OUT.mkdir(parents=True, exist_ok=True)
rng = np.random.default_rng(7)

def on_boundary(v, tol=1e-9):
    return np.max(np.abs(v)) >= 1 - tol

def radial_to_face(v):
    """Pointwise face-projection. Undefined / discontinuous at 0."""
    m = np.max(np.abs(v))
    if m < 1e-12:
        return None
    return np.clip(v / m, -1, 1)

n_side = 81
xs = np.linspace(-1, 1, n_side)
XX, YY = np.meshgrid(xs, xs)
pts = np.stack([XX.ravel(), YY.ravel()], axis=1)
boundary_mask = (np.abs(pts[:, 0]) > 1 - 1e-9) | (np.abs(pts[:, 1]) > 1 - 1e-9)

radial = []
died_at_center = False
for p in pts:
    q = radial_to_face(p)
    if q is None:
        died_at_center = True
        radial.append(np.array([np.nan, np.nan]))
    else:
        radial.append(q)
radial = np.array(radial)
b_ok = np.allclose(radial[boundary_mask], pts[boundary_mask])

edge = np.linspace(-1, 1, 400)
bd = np.concatenate([
    np.stack([edge, np.full_like(edge, -1)], 1),
    np.stack([np.full_like(edge, 1), edge], 1),
    np.stack([edge[::-1], np.full_like(edge, 1)], 1),
    np.stack([np.full_like(edge, -1), edge[::-1]], 1),
])
ang = np.unwrap(np.arctan2(bd[:, 1], bd[:, 0]))
degree = (ang[-1] - ang[0]) / (2 * np.pi)

def continuous_attempt(p, hole=0.35):
    m = np.max(np.abs(p))
    if m <= hole:
        return p
    q = p / m
    w = (m - hole) / (1 - hole)
    return (1 - w) * p + w * q

attempt = np.array([continuous_attempt(p) for p in pts])
still_interior = np.mean([not on_boundary(p) for p in attempt[~boundary_mask]])
boundary_frozen = np.allclose(attempt[boundary_mask], pts[boundary_mask])

theta = np.linspace(0, 2 * np.pi, 240, endpoint=False)

def thin_map(th):
    slots = np.zeros((3, 3))
    slots[0] = np.array([0.45 * np.cos(th), 0.45 * np.sin(th), 0.15])
    slots[1] = np.array([0.2, 0.55 * np.cos(2 * th), 0.25 * np.sin(th)])
    if np.cos(th) > 0:
        slots[2] = np.array([1.0, 0.3 * np.sin(th), -0.2])
    else:
        slots[2] = np.array([0.4 * np.cos(th), 0.35, 0.2 * np.sin(th)])
    return slots

before = np.array([thin_map(th) for th in theta])

def compress_slots(slots):
    out = slots.copy()
    for i, v in enumerate(slots):
        if on_boundary(v):
            continue
        q = radial_to_face(v)
        out[i] = q
    return out

after = np.array([compress_slots(s) for s in before])

def slot_hit_rate(arr):
    hits = [on_boundary(v) for sample in arr for v in sample]
    return float(np.mean(hits))

rate_before = slot_hit_rate(before)
rate_after = slot_hit_rate(after)
preserved = True
for b, a in zip(before, after):
    for vb, va in zip(b, a):
        if on_boundary(vb) and not np.allclose(vb, va):
            preserved = False

def power_rate(u):
    ph = rng.uniform(0, 2 * np.pi, size=u)
    slots = np.concatenate([thin_map(p) for p in ph], axis=0)
    comp = compress_slots(slots)
    return float(np.mean([on_boundary(v) for v in comp]))

powers = [1, 2, 4, 8, 16]
power_rates = [power_rate(u) for u in powers]

d_values = np.array([1, 2, 4, 8])
steps = 32
ledger = []
for d in d_values:
    revealed = d * steps
    pinned = 0.5 * revealed
    ledger.append((int(d), revealed, pinned, 0.0))

d = 2
window = rng.random((steps, d))

fig, axes = plt.subplots(1, 3, figsize=(12.4, 4.2))
ax = axes[0]
ax.set_title("filled square, face push dies at 0")
ax.scatter(pts[~boundary_mask, 0], pts[~boundary_mask, 1], s=6, c="#9bb7d4", linewidths=0)
ax.plot(bd[:, 0], bd[:, 1], color="#1f4e79", lw=2)
ax.scatter([0], [0], s=40, c="#b00020", zorder=3)
ax.annotate("no continuous\nvalue here", (0.05, 0.08), fontsize=8, color="#b00020")
ax.set_aspect("equal")
ax.set_xlim(-1.15, 1.15)
ax.set_ylim(-1.15, 1.15)
ax.set_xlabel(f"boundary degree {degree:.0f}. radial map undefined at center.")

ax = axes[1]
ax.set_title("continuous attempt leaves a hole")
interior = attempt[~boundary_mask]
hit = np.array([on_boundary(p) for p in interior])
ax.scatter(interior[~hit, 0], interior[~hit, 1], s=7, c="#b00020", label="still inside")
ax.scatter(interior[hit, 0], interior[hit, 1], s=7, c="#2e7d32", label="reached a face")
ax.plot(bd[:, 0], bd[:, 1], color="#1f4e79", lw=2)
ax.set_aspect("equal")
ax.set_xlim(-1.15, 1.15)
ax.set_ylim(-1.15, 1.15)
ax.legend(fontsize=7, loc="upper right")
ax.set_xlabel(f"boundary frozen={boundary_frozen}. interior still inside={still_interior:.2f}")

ax = axes[2]
ax.set_title("thin loop: face push is continuous")
ax.plot(before[:, 0, 0], before[:, 0, 1], color="#9bb7d4", lw=2, label="before")
ax.plot(after[:, 0, 0], after[:, 0, 1], color="#2e7d32", lw=2, label="after")
ax.add_patch(Rectangle((-1, -1), 2, 2, fill=False, ec="#1f4e79", lw=1))
ax.scatter([0], [0], s=20, c="#b00020")
ax.set_aspect("equal")
ax.set_xlim(-1.2, 1.2)
ax.set_ylim(-1.2, 1.2)
ax.legend(fontsize=7)
ax.set_xlabel(f"face rate {rate_before:.2f} -> {rate_after:.2f}. frozen={preserved}")
fig.suptitle("What the lemma does: pin slots to a face, unless the cube is full", fontsize=12)
fig.tight_layout()
fig.savefig(OUT / "compression_move.png", dpi=140)
plt.close()

fig, axes = plt.subplots(1, 2, figsize=(10.6, 4.2))
ax = axes[0]
ax.plot(powers, power_rates, marker="o", color="#1f4e79")
ax.set_ylim(0, 1.05)
ax.set_xlabel("copies u")
ax.set_ylabel("fraction of slots on a face")
ax.set_title("thin map, after the nudge")
ax.axhline(rate_before, color="#9bb7d4", ls="--", label="before")
ax.legend(fontsize=8)
ax = axes[1]
ds = [row[0] for row in ledger]
ax.plot(ds, [row[1] / steps for row in ledger], marker="o", label="new parameters / step  (mdim)")
ax.plot(ds, [row[2] / steps for row in ledger], marker="s", label="face pins / step  (half)")
ax.plot(ds, [row[3] / steps for row in ledger], marker="x", label="pins if the cube is full")
ax.set_xlabel("parameters revealed per iterate")
ax.set_ylabel("count per iterate")
ax.set_title("why the 1/2 is the useful part")
ax.legend(fontsize=7)
fig.tight_layout()
fig.savefig(OUT / "half_parameter.png", dpi=140)
plt.close()

fig, ax = plt.subplots(figsize=(10.2, 3.2))
ax.imshow(window.T, aspect="auto", cmap="Blues", vmin=0, vmax=1)
for t in range(steps):
    ax.add_patch(Rectangle((t - 0.5, 0.5), 1, 1, fill=False, ec="#2e7d32", lw=1.2))
ax.set_yticks([0, 1])
ax.set_yticklabels(["free parameter", "pinned parameter"])
ax.set_xlabel("iterate of the shift")
ax.set_title("d=2 shift window. green row is the half the lemma is allowed to pin")
fig.tight_layout()
fig.savefig(OUT / "shift_ledger.png", dpi=140)
plt.close()

lines = [
    f"identity boundary degree = {degree:.4f} (must be 1)",
    f"radial face push undefined at center = {died_at_center}",
    f"radial agrees with identity on boundary samples = {b_ok}",
    f"continuous attempt freezes boundary = {boundary_frozen}",
    f"continuous attempt, fraction of interior still inside = {still_interior:.4f}",
    f"thin map face rate before = {rate_before:.4f}",
    f"thin map face rate after = {rate_after:.4f}",
    f"thin map froze every old face slot = {preserved}",
    "power face rates after nudge: " + ", ".join(f"u={u}:{r:.3f}" for u, r in zip(powers, power_rates)),
    "ledger (d, params in 32 steps, pins at half, pins if full):",
]
for row in ledger:
    lines.append(f"  {row}")
(OUT / "sim_log.txt").write_text("\n".join(lines) + "\n")
print("\n".join(lines))
