#!/usr/bin/env python3
"""Compare uniform and adaptive frequency sweep S-parameters.

Usage: python3 plot_comparison.py

Reads Palace postpro CSV files from postpro/uniform/ and postpro/adaptive/.
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import re
import os


def read_palace_csv(path):
    """Read Palace port-floquet-S.csv. Returns dict of column_name -> array."""
    with open(path) as f:
        header = f.readline().strip()
    cols = [c.strip() for c in header.split(",")]
    data = np.genfromtxt(path, delimiter=",", skip_header=1)
    if data.ndim == 1:
        data = data.reshape(1, -1)
    return {name: data[:, i] for i, name in enumerate(cols)}


def extract_palace_modes(data):
    """Extract non-null |S| columns grouped by (port, m, n, pol)."""
    freq = data["f (GHz)"]
    modes = {}
    for key, vals in data.items():
        m = re.match(
            r"\|S\[P(\d+)\((-?\d+);(-?\d+)\)(\w+)\]\[(\d+)\]\| \(dB\)", key.strip()
        )
        if m and not np.all(np.isnan(vals)):
            port, mm, nn, pol = int(m.group(1)), int(m.group(2)), int(m.group(3)), m.group(4)
            label = f"P{port}({mm},{nn}){pol}"
            modes[label] = vals
    return freq, modes


def main():
    uniform_path = "postpro/uniform/port-floquet-S.csv"
    adaptive_path = "postpro/adaptive/port-floquet-S.csv"

    have_uniform = os.path.exists(uniform_path)
    have_adaptive = os.path.exists(adaptive_path)

    if not have_uniform and not have_adaptive:
        print("No Palace data found. Run the simulations first:")
        print("  palace periodic_box_uniform.json")
        print("  palace periodic_box_adaptive.json")
        return

    freq_u, modes_u = extract_palace_modes(read_palace_csv(uniform_path)) if have_uniform else (None, {})
    freq_a, modes_a = extract_palace_modes(read_palace_csv(adaptive_path)) if have_adaptive else (None, {})

    all_labels = sorted(set(modes_u.keys()) | set(modes_a.keys()))
    p1_labels = [l for l in all_labels if l.startswith("P1")]
    p2_labels = [l for l in all_labels if l.startswith("P2")]

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 10), sharex=True)
    colors = plt.cm.tab10(np.linspace(0, 1, 10))
    y_min = -30

    for ax, labels, title in [
        (ax1, p1_labels, "Reflection (Port 1)"),
        (ax2, p2_labels, "Transmission (Port 2)"),
    ]:
        plotted_vals = []
        for idx, label in enumerate(labels):
            color = colors[idx % len(colors)]

            # Skip degenerate symmetric modes (identical S-parameter values).
            a_vals = modes_a.get(label) if have_adaptive else modes_u.get(label)
            if a_vals is not None:
                is_dup = any(
                    len(prev) == len(a_vals) and np.allclose(
                        a_vals[~np.isnan(a_vals) & ~np.isnan(prev)],
                        prev[~np.isnan(a_vals) & ~np.isnan(prev)], rtol=1e-3
                    ) if np.any(~np.isnan(a_vals) & ~np.isnan(prev)) else False
                    for prev in plotted_vals
                )
                if is_dup:
                    continue
                plotted_vals.append(a_vals)

            # Adaptive sweep (solid line).
            if label in modes_a:
                vals = modes_a[label]
                mask = ~np.isnan(vals) & (vals > y_min)
                if np.any(mask):
                    ax.plot(freq_a[mask], vals[mask], "-", color=color, linewidth=1.5,
                            label=f"{label} (adaptive)")

            # Uniform sweep (markers).
            if label in modes_u:
                vals = modes_u[label]
                mask = ~np.isnan(vals) & (vals > y_min)
                if np.any(mask):
                    ax.plot(freq_u[mask], vals[mask], "o", color=color, markersize=6,
                            markerfacecolor="none", markeredgewidth=1.5,
                            label=f"{label} (uniform)")

        ax.set_ylabel("|S| (dB)")
        ax.set_title(title)
        ax.legend(fontsize=7, ncol=3, loc="best")
        ax.grid(True, alpha=0.3)
        ax.set_ylim(bottom=y_min)

    ax2.set_xlabel("Frequency (GHz)")
    fig.suptitle(
        "Dielectric Grating Floquet Port S-Parameters\n"
        "30° oblique incidence (TE), 4×1×8 cm domain, ε_slab=7",
        fontsize=11,
    )
    plt.tight_layout()
    plt.savefig("freq_sweep_comparison.png", dpi=150, bbox_inches="tight")
    print("Saved freq_sweep_comparison.png")


if __name__ == "__main__":
    main()
