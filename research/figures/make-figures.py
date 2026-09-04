#!/usr/bin/env python3
"""Figures for the dev.to post. Every number plotted is computed here.

    python3 post-figures.py

Writes into images/post/.
"""
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

OUT = "images/post"
os.makedirs(OUT, exist_ok=True)

INK, MUTED, LINE = "#14181f", "#5b6472", "#dde2e9"
BAD, OK, PANEL = "#b23b32", "#2f7a54", "#f6f7f9"

plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 11,
    "axes.edgecolor": LINE, "axes.labelcolor": MUTED,
    "xtick.color": MUTED, "ytick.color": MUTED,
    "figure.facecolor": "white", "axes.facecolor": "white",
})

H = np.array([[1, 1], [1, -1]], dtype=complex) / np.sqrt(2)
I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Z = np.diag([1, -1]).astype(complex)


def rz(t, bad=False):
    s = -1 if bad else 1
    return np.diag([np.exp(-1j * s * t / 2), np.exp(1j * s * t / 2)])


def rx(t):
    c, s = np.cos(t / 2), np.sin(t / 2)
    return np.array([[c, -1j * s], [-1j * s, c]])


def save(fig, name):
    path = f"{OUT}/{name}"
    fig.savefig(path, dpi=150, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"  {path}")


# ---------------------------------------------------------------- figure 1
def fig_witness():
    """H then Rz(theta): state fidelity is cos^2, probabilities never move."""
    theta = np.linspace(0, np.pi, 400)
    fidelity, divergence = [], []

    for t in theta:
        psi = rz(t) @ H @ np.array([1, 0], dtype=complex)
        phi = rz(t, bad=True) @ H @ np.array([1, 0], dtype=complex)
        fidelity.append(abs(np.vdot(psi, phi)) ** 2)
        divergence.append(np.abs(np.abs(psi) ** 2 - np.abs(phi) ** 2).max())

    fidelity = np.array(fidelity)
    assert np.allclose(fidelity, np.cos(theta) ** 2, atol=1e-12)
    assert max(divergence) < 1e-15, max(divergence)

    fig, ax = plt.subplots(figsize=(9, 4.6))
    ax.plot(theta, fidelity, color=BAD, lw=2.6,
            label=r"state fidelity  $|\langle\psi|\phi\rangle|^2$")
    ax.plot(theta, divergence, color=OK, lw=2.6,
            label="largest difference in measured probability")

    ax.axvline(np.pi / 2, color=MUTED, ls=":", lw=1.2)
    ax.annotate("orthogonal states,\nidentical measurements",
                xy=(np.pi / 2, 0), xytext=(np.pi / 2 + 0.18, 0.30),
                fontsize=10.5, color=INK,
                arrowprops=dict(arrowstyle="->", color=MUTED, lw=1.2))

    ax.set_xticks([0, np.pi / 4, np.pi / 2, 3 * np.pi / 4, np.pi])
    ax.set_xticklabels(["0", "π/4", "π/2", "3π/4", "π"])
    ax.set_xlabel("θ")
    ax.set_ylim(-0.04, 1.06)
    ax.set_xlim(0, np.pi)
    ax.legend(frameon=False, loc="upper center", fontsize=10.5)
    ax.set_title("Two gates: H, then Rz(θ), with the sign of Rz inverted",
                 fontsize=13, fontweight="bold", color=INK, pad=12)
    ax.grid(alpha=0.25, lw=0.7)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    save(fig, "fig-witness.png")
    return fidelity


# ---------------------------------------------------------------- figure 2
def maxcut_cost(gamma, beta, bad, n=4):
    """One QAOA layer for MaxCut on a 4-cycle. Returns the expected cut value."""
    dim = 2 ** n
    edges = [(0, 1), (1, 2), (2, 3), (3, 0)]

    def op1(m, q):
        out = np.array([[1]], dtype=complex)
        for i in range(n):
            out = np.kron(out, m if i == q else I2)
        return out

    state = np.ones(dim, dtype=complex) / np.sqrt(dim)   # |+>^n

    for a, b in edges:                                    # exp(-i gamma Z_a Z_b)
        za = np.array([1 if not (i >> (n - 1 - a)) & 1 else -1 for i in range(dim)])
        zb = np.array([1 if not (i >> (n - 1 - b)) & 1 else -1 for i in range(dim)])
        sign = -1 if bad else 1
        state = state * np.exp(-1j * sign * gamma * za * zb)

    for q in range(n):                                    # mixer
        state = op1(rx(2 * beta), q) @ state

    probs = np.abs(state) ** 2
    cut = np.zeros(dim)
    for a, b in edges:
        for i in range(dim):
            if ((i >> (n - 1 - a)) & 1) != ((i >> (n - 1 - b)) & 1):
                cut[i] += 1
    return float(probs @ cut)


def fig_qaoa_grid():
    N = 80
    gammas = np.linspace(np.pi / N, np.pi, N)     # (0, pi], as in the book
    betas = np.linspace(np.pi / N, np.pi, N)
    diff = np.zeros((N, N))

    for i, b in enumerate(betas):
        for j, g in enumerate(gammas):
            diff[i, j] = abs(maxcut_cost(g, b, False) - maxcut_cost(g, b, True))

    tol = 1e-9
    frac = float((diff > tol).mean())
    print(f"  backends differ at {frac*100:.1f}% of the {N}x{N} grid")

    gt, bt = np.pi / 4, np.pi / 8
    at_tidy = abs(maxcut_cost(gt, bt, False) - maxcut_cost(gt, bt, True))
    print(f"  at gamma=pi/4, beta=pi/8 the difference is {at_tidy:.2e}")

    fig, ax = plt.subplots(figsize=(8.2, 5.4))
    im = ax.pcolormesh(gammas, betas, diff, cmap="RdPu", shading="auto")
    fig.colorbar(im, ax=ax, label="|difference in expected cut|")

    ax.plot([gt], [bt], marker="o", ms=13, mfc="none", mec=INK, mew=2.4)
    ax.annotate("γ = π/4, β = π/8\nthe angles I chose",
                xy=(gt, bt), xytext=(gt + 0.42, bt + 0.30),
                fontsize=11, color=INK, fontweight="bold",
                arrowprops=dict(arrowstyle="->", color=INK, lw=1.6))

    ax.set_xlabel("γ")
    ax.set_ylabel("β")
    ax.set_xticks([0, np.pi / 4, np.pi / 2, 3 * np.pi / 4, np.pi])
    ax.set_xticklabels(["0", "π/4", "π/2", "3π/4", "π"])
    ax.set_yticks([0, np.pi / 4, np.pi / 2, 3 * np.pi / 4, np.pi])
    ax.set_yticklabels(["0", "π/4", "π/2", "3π/4", "π"])
    ax.set_title(f"The two backends disagree at {frac*100:.1f}% of these angles",
                 fontsize=13, fontweight="bold", color=INK, pad=12)
    save(fig, "fig-qaoa-grid.png")
    return frac, at_tidy


if __name__ == "__main__":
    print("figures:")
    fig_witness()
    frac, tidy = fig_qaoa_grid()
    print(f"\nfor the post text: {frac*100:.1f}% of the grid diverges; "
          f"the chosen angles differ by {tidy:.1e}")
