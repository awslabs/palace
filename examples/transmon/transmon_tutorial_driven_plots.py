#!/usr/bin/env -S uv run --script
# coding: utf-8
# /// script
# requires-python = ">=3.11"
# dependencies = [
#   "matplotlib",
#   "numpy",
#   "pandas",
# ]
# ///

# Transmon Driven Solver — Uniform vs Adaptive

# Generates plots for the "Driven Solver: Uniform vs Adaptive" tutorial using the
# transmon model.
#
# Run `transmon_tutorial_driven.jl` first to generate the simulation data. Output plots
# are written to the `docs/src/assets/examples` folder.
#
# Run stand-alone from the palace root:
# ```
# uv run --script examples/transmon/transmon_tutorial_driven_plots.py
# ```


## Imports

import pathlib
import re

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

mpl.rcParams.update(
    {
        "figure.dpi": 150,
        "font.family": "serif",
        "font.serif": ["Times"],
        "mathtext.fontset": "cm",
        "axes.formatter.use_mathtext": True,
        "savefig.bbox": None,
    }
)

## Configuration

TRANSMON_DIR = pathlib.Path("./examples/transmon")
DOCS_ASSET_DIR = pathlib.Path("./docs/src/assets/examples")

UNIFORM_DIR = (
    TRANSMON_DIR / "postpro/transmon_tutorial_driven_rom/driven_uniform_reference"
)

ADAPTIVE_TOLS = [1e-1, 1e-2, 1e-3, 1e-4, 1e-5]
ADAPTIVE_DIRS = {
    tol: TRANSMON_DIR
    / f"postpro/transmon_tutorial_driven_rom/driven_adaptive_1e{int(round(np.log10(tol)))}"
    for tol in ADAPTIVE_TOLS
}

# Two feedline ports (excitations 1 and 2), showing the full 2-port S-matrix
PORT_PAIRS = [(1, 1), (2, 1), (1, 2), (2, 2)]


## Load Data

_RE_S_DB = re.compile(r"\|S\[(\d+)\]\[(\d+)\]\| \(dB\)")
_RE_E_COL = re.compile(r"E_(\w+)\[(\d+)\] \(J\)")
_RE_SAMPLED_FREQS = re.compile(
    r" Sampled frequencies \(GHz\): ((?s:.)*?)\n Sample error", re.MULTILINE
)
_RE_SAMPLED_ERRORS = re.compile(r"^ Sample errors: (.*)$", re.MULTILINE)


def load_port_s(postpro_dir: pathlib.Path) -> tuple[np.ndarray, dict]:
    """Load port-S.csv; return (freq_GHz, {(i, j): complex S_ij array}).

    S-parameter columns containing non-numeric values (e.g. 'inf' for passive
    reactive ports) are silently skipped.
    """
    df = pd.read_csv(postpro_dir / "port-S.csv")
    df.columns = df.columns.str.strip()
    freq = df["f (GHz)"].values
    s = {}
    for col in df.columns:
        m = _RE_S_DB.match(col)
        if m:
            i, j = int(m.group(1)), int(m.group(2))
            db = pd.to_numeric(df[f"|S[{i}][{j}]| (dB)"], errors="coerce").values
            ang = pd.to_numeric(df[f"arg(S[{i}][{j}]) (deg.)"], errors="coerce").values
            if np.any(np.isnan(db)) or np.any(np.isnan(ang)):
                continue  # skip passive/reactive ports with undefined S-parameters
            s[(i, j)] = 10 ** (db / 20.0) * np.exp(1j * np.deg2rad(ang))
    return freq, s


def load_domain_e(postpro_dir: pathlib.Path) -> tuple[np.ndarray, dict]:
    """Load domain-E.csv; return (freq_GHz, {(quantity, excitation): array}).

    Keys are tuples e.g. ("elec", 1), ("elec", 2) for two excitations.
    """
    df = pd.read_csv(postpro_dir / "domain-E.csv")
    df.columns = df.columns.str.strip()
    freq = df["f (GHz)"].values
    e = {}
    for col in df.columns:
        m = _RE_E_COL.match(col)
        if m:
            e[(m.group(1), int(m.group(2)))] = df[col].values
    return freq, e


def parse_log_pivots(log_path: pathlib.Path) -> list[np.ndarray]:
    """Parse sampled frequencies from palace.log.
    Returns list of arrays, one per excitation, in GHz."""
    if not log_path.exists():
        return []
    with open(log_path) as fp:
        log_content = fp.read()
    sections = _RE_SAMPLED_FREQS.findall(log_content)
    return [
        np.array(
            [float(f.strip()) for f in match.replace("\n", "").split(",") if f.strip()],
            dtype=np.float64,
        )
        for match in sections
    ]


def parse_log_errors(log_path: pathlib.Path) -> list[np.ndarray]:
    """Parse MRI convergence error values from palace.log.
    Returns list of arrays, one per excitation. First two entries are inf."""
    if not log_path.exists():
        return []
    with open(log_path) as fp:
        log_content = fp.read()
    sections = _RE_SAMPLED_ERRORS.findall(log_content)
    return [
        np.array(
            [float(v.strip()) for v in match.split(",") if v.strip()], dtype=np.float64
        )
        for match in sections
    ]


# Load uniform reference
freq_uniform, s_uniform = load_port_s(UNIFORM_DIR)
_, e_uniform = load_domain_e(UNIFORM_DIR)
print(f"Uniform: {len(freq_uniform)} frequency points")

# Load adaptive results
freq_adaptive: dict[float, np.ndarray] = {}
s_adaptive: dict[float, dict] = {}
e_adaptive: dict[float, dict] = {}
sampled_freqs: dict[float, list[np.ndarray]] = {}
sampled_errors: dict[float, list[np.ndarray]] = {}

for tol in ADAPTIVE_TOLS:
    d = ADAPTIVE_DIRS[tol]
    if not (d / "port-S.csv").exists():
        print(f"WARNING: {d} not found — run transmon_tutorial_driven.jl first")
        continue
    freq_adaptive[tol], s_adaptive[tol] = load_port_s(d)
    _, e_adaptive[tol] = load_domain_e(d)
    log_path = d / "palace.log"
    sampled_freqs[tol] = parse_log_pivots(log_path)
    sampled_errors[tol] = parse_log_errors(log_path)
    n_samp = sum(len(a) for a in sampled_freqs[tol])
    print(
        f"Adaptive tol={tol:.0e}: {len(freq_adaptive[tol])} frequency points "
        f"({n_samp} adaptive samples)"
    )


## Common Plotting Options

tol_palette = plt.cm.plasma(np.linspace(0.1, 0.9, len(ADAPTIVE_TOLS)))

GOLDEN_RATIO = (1 + np.sqrt(5)) / 2

_FREQ_LIM = (3.3, 6.7)
_FREQ_TICKS = np.linspace(3.5, 6.5, 7)

# Known eigenmode frequencies: see test/examples/ref/transmon/transmon_coarse
EIGENMODE_FREQS_GHZ = [4.099115457610e00, 5.603265962190e00]


def add_eigenmode_lines(ax) -> None:
    for f in EIGENMODE_FREQS_GHZ:
        ax.axvline(
            f, color="black", linewidth=0.6, linestyle=":", alpha=0.25, zorder=-4
        )


def adjust_golden_layout(fig, ax, pad_in: float = 0.1) -> None:
    dpi = fig.dpi
    fig_width = fig.get_figwidth()
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    axes_bb = ax.get_window_extent(renderer)
    tight_bb = ax.get_tightbbox(renderer)
    left_extra_in = (axes_bb.x0 - tight_bb.x0) / dpi
    bottom_extra_in = (axes_bb.y0 - tight_bb.y0) / dpi
    top_extra_in = (tight_bb.y1 - axes_bb.y1) / dpi
    margin_in = left_extra_in + pad_in
    bottom_space_in = bottom_extra_in + pad_in
    top_space_in = top_extra_in + pad_in
    axes_width_in = fig_width - 2 * margin_in
    axes_height_in = axes_width_in / GOLDEN_RATIO
    fig_height = axes_height_in + bottom_space_in + top_space_in
    fig.set_figheight(fig_height)
    fig.subplots_adjust(
        left=margin_in / fig_width,
        right=1 - margin_in / fig_width,
        bottom=bottom_space_in / fig_height,
        top=1 - top_space_in / fig_height,
    )


def adjust_golden_layout_2x2(
    fig, axes, pad_in: float = 0.1, hspace_in: float = 0.35, vspace_in: float = 0.45
) -> None:
    dpi = fig.dpi
    fig_width = fig.get_figwidth()
    nrows, ncols = axes.shape
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    ax_tl = axes[0, 0]
    axes_bb = ax_tl.get_window_extent(renderer)
    tight_bb = ax_tl.get_tightbbox(renderer)
    left_extra_in = (axes_bb.x0 - tight_bb.x0) / dpi
    ax_bl = axes[nrows - 1, 0]
    axes_bb_bl = ax_bl.get_window_extent(renderer)
    tight_bb_bl = ax_bl.get_tightbbox(renderer)
    bottom_extra_in = (axes_bb_bl.y0 - tight_bb_bl.y0) / dpi
    top_extra_in = max(
        (ax.get_tightbbox(renderer).y1 - ax.get_window_extent(renderer).y1) / dpi
        for ax in axes[0]
    )
    suptitle_h_in = 0.0
    if fig._suptitle is not None:
        suptitle_h_in = (
            fig._suptitle.get_window_extent(renderer).height / dpi + pad_in * 0.5
        )
    margin_lr_in = left_extra_in + pad_in
    bottom_space_in = bottom_extra_in + pad_in
    top_space_in = top_extra_in + pad_in + suptitle_h_in
    axes_width_in = (fig_width - 2 * margin_lr_in - hspace_in * (ncols - 1)) / ncols
    axes_height_in = axes_width_in / GOLDEN_RATIO
    fig_height = (
        top_space_in
        + nrows * axes_height_in
        + vspace_in * (nrows - 1)
        + bottom_space_in
    )
    fig.set_figheight(fig_height)
    fig.subplots_adjust(
        left=margin_lr_in / fig_width,
        right=1 - margin_lr_in / fig_width,
        bottom=bottom_space_in / fig_height,
        top=1 - top_space_in / fig_height,
        wspace=hspace_in / axes_width_in,
        hspace=vspace_in / axes_height_in,
    )


def adjust_golden_layout_1x2(
    fig, axes, pad_in: float = 0.1, hspace_in: float = 0.35
) -> None:
    """Adjust a 1×2 subplot figure: equal left/right margins, golden-ratio axes.

    axes must be a 1-D array of two Axes (from plt.subplots(1, 2)).
    """
    dpi = fig.dpi
    fig_width = fig.get_figwidth()

    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()

    # Left margin from left axes (has y-labels)
    axes_bb = axes[0].get_window_extent(renderer)
    tight_bb = axes[0].get_tightbbox(renderer)
    left_extra_in = (axes_bb.x0 - tight_bb.x0) / dpi

    # Bottom decoration
    axes_bb_b = axes[0].get_window_extent(renderer)
    tight_bb_b = axes[0].get_tightbbox(renderer)
    bottom_extra_in = (axes_bb_b.y0 - tight_bb_b.y0) / dpi

    # Top decoration: max across both axes
    top_extra_in = max(
        (ax.get_tightbbox(renderer).y1 - ax.get_window_extent(renderer).y1) / dpi
        for ax in axes
    )

    # Suptitle height
    suptitle_h_in = 0.0
    if fig._suptitle is not None:
        suptitle_h_in = (
            fig._suptitle.get_window_extent(renderer).height / dpi + pad_in * 0.5
        )

    margin_lr_in = left_extra_in + pad_in
    bottom_space_in = bottom_extra_in + pad_in
    top_space_in = top_extra_in + pad_in + suptitle_h_in

    axes_width_in = (fig_width - 2 * margin_lr_in - hspace_in) / 2
    axes_height_in = axes_width_in / GOLDEN_RATIO

    fig_height = top_space_in + axes_height_in + bottom_space_in
    fig.set_figheight(fig_height)

    fig.subplots_adjust(
        left=margin_lr_in / fig_width,
        right=1 - margin_lr_in / fig_width,
        bottom=bottom_space_in / fig_height,
        top=1 - top_space_in / fig_height,
        wspace=hspace_in / axes_width_in,
    )


## S-Parameters: Uniform Reference


def plot_smat_uniform(
    freq_uniform, s_uniform, PORT_PAIRS, DOCS_ASSET_DIR, adjust_golden_layout_2x2
):
    fig, axes = plt.subplots(2, 2, figsize=(6, 4), sharex=True, sharey=True)
    for idx, (i, j) in enumerate(PORT_PAIRS):
        row, col = divmod(idx, 2)
        ax = axes[row, col]
        if (i, j) in s_uniform:
            ax.plot(
                freq_uniform,
                20 * np.log10(np.abs(s_uniform[(i, j)])),
                "k-",
                linewidth=1,
                label="Uniform",
                zorder=10,
                marker="o",
                markersize=2,
            )
        add_eigenmode_lines(ax)
        ax.text(
            0.97,
            0.93,
            rf"$S_{{{i}{j}}}$",
            transform=ax.transAxes,
            ha="right",
            va="top",
            fontsize=12,
        )
        if row == 1:
            ax.set_xlabel("f (GHz)")
        ax.set_xlim(*_FREQ_LIM)
        ax.set_xticks(_FREQ_TICKS)
        ax.set_ylim(-25, 1)
    axes[0, 1].legend(loc="lower right", fontsize=8, frameon=False)
    fig.suptitle("S-Parameters Magnitude $\\vert S_{{ij}}\\vert$ (dB)")
    adjust_golden_layout_2x2(fig, axes, hspace_in=0.1, vspace_in=0.15)
    plt.savefig(DOCS_ASSET_DIR / "driven_ua_transmon_sparam_uniform.svg")
    plt.show()


plot_smat_uniform(
    freq_uniform, s_uniform, PORT_PAIRS, DOCS_ASSET_DIR, adjust_golden_layout_2x2
)


## S-Parameters: RMS-Normalized Absolute Error


def plot_smat_adaptive_rms(
    freq_uniform,
    s_uniform,
    s_adaptive,
    freq_adaptive,
    PORT_PAIRS,
    ADAPTIVE_TOLS,
    sampled_freqs,
    tol_palette,
    DOCS_ASSET_DIR,
    adjust_golden_layout_2x2,
):
    fig, axes = plt.subplots(2, 2, figsize=(6, 4), sharex=True, sharey=True)
    for idx, (i, j) in enumerate(PORT_PAIRS):
        row, col = divmod(idx, 2)
        ax = axes[row, col]
        if (i, j) not in s_uniform:
            continue
        ref = s_uniform[(i, j)]
        freq = freq_uniform
        scale = np.sqrt(np.mean(np.square(np.abs(ref))))  # RMS of |S|
        for k, tol in enumerate(ADAPTIVE_TOLS):
            if tol not in s_adaptive or (i, j) not in s_adaptive[tol]:
                continue
            assert np.all(freq_adaptive[tol] == freq_uniform)
            err = np.abs(s_adaptive[tol][(i, j)] - ref) / scale
            ax.plot(
                freq,
                err,
                color=tol_palette[k],
                linewidth=0.2,
                label=f"{tol:.0e}",
                marker="o",
                markersize=2,
            )
            ax.axhline(
                tol,
                0,
                0.80,
                color=tol_palette[k],
                linewidth=1,
                linestyle="--",
                zorder=-5,
            )
            pivots_list = sampled_freqs.get(tol, [])
            all_pivots = np.concatenate(pivots_list) if pivots_list else np.array([])
            if all_pivots.size:
                ax.plot(
                    all_pivots,
                    np.ones_like(all_pivots) * 1.3e-13 * (1.7**k),
                    color=tol_palette[k],
                    linewidth=0,
                    marker="D",
                    markersize=4,
                    alpha=0.75,
                )
        ax.axhline(1e-12, color="black", linewidth=1, linestyle=":", zorder=-5)
        add_eigenmode_lines(ax)
        ax.text(
            0.97,
            0.97,
            rf"$S_{{{i}{j}}}$",
            transform=ax.transAxes,
            ha="right",
            va="top",
            fontsize=12,
        )
        if row == 1:
            ax.set_xlabel("f (GHz)")
        ax.set_xlim(*_FREQ_LIM)
        ax.set_xticks(_FREQ_TICKS)
        ax.set_ylim(5e-14, 0.9)
        ax.set_yscale("log")
    axes[1, 1].legend(
        loc="upper right",
        bbox_to_anchor=(1.0, 0.88),
        fontsize=8,
        frameon=False,
        ncol=1,
        handlelength=1.0,
        handletextpad=0.4,
        columnspacing=0.3,
    )
    fig.suptitle(
        "Normalized Absolute Error in S-Parameters\n"
        "$\\vert S_{{\\mathrm{{adaptive}}}} - S_{{\\mathrm{{uniform}}}}\\vert"
        " /  \\vert\\vert S_{{\\mathrm{{uniform}}}}\\vert\\vert_{{\\mathrm{{RMS}}}}$"
    )
    adjust_golden_layout_2x2(fig, axes, hspace_in=0.1, vspace_in=0.15)
    plt.savefig(DOCS_ASSET_DIR / "driven_ua_transmon_sparam_adaptive_rms.svg")
    plt.show()


plot_smat_adaptive_rms(
    freq_uniform,
    s_uniform,
    s_adaptive,
    freq_adaptive,
    PORT_PAIRS,
    ADAPTIVE_TOLS,
    sampled_freqs,
    tol_palette,
    DOCS_ASSET_DIR,
    adjust_golden_layout_2x2,
)


## Electric Energy: Uniform Reference


def plot_energy_uniform(
    freq_uniform, e_uniform, DOCS_ASSET_DIR, adjust_golden_layout_1x2
):
    fig, axes = plt.subplots(1, 2, figsize=(6, 4), sharex=True, sharey=True)
    for col, excitation in enumerate([1, 2]):
        ax = axes[col]
        key = ("elec", excitation)
        if key in e_uniform:
            ax.plot(
                freq_uniform,
                e_uniform[key] * 1e9,
                "k-",
                linewidth=1,
                label="Uniform",
                marker="o",
                markersize=2,
            )
        add_eigenmode_lines(ax)
        ax.text(
            0.97,
            0.97,
            f"Excitation {excitation}",
            transform=ax.transAxes,
            ha="right",
            va="top",
            fontsize=9,
        )
        ax.set_xlabel("f (GHz)")
        ax.set_xlim(*_FREQ_LIM)
        ax.set_xticks(_FREQ_TICKS)
    axes[1].legend(loc="upper left", fontsize=8, frameon=False)
    fig.suptitle("Domain Energy $E_\\mathrm{elec}$ (nJ)")
    adjust_golden_layout_1x2(fig, axes, hspace_in=0.1)
    plt.savefig(DOCS_ASSET_DIR / "driven_ua_transmon_energy_uniform.svg")
    plt.show()


plot_energy_uniform(freq_uniform, e_uniform, DOCS_ASSET_DIR, adjust_golden_layout_1x2)


## Electric Energy: Adaptive Sweep


def plot_energy_adaptive_sweep(
    freq_uniform,
    e_uniform,
    e_adaptive,
    freq_adaptive,
    ADAPTIVE_TOLS,
    sampled_freqs,
    tol_palette,
    DOCS_ASSET_DIR,
    adjust_golden_layout_1x2,
):
    fig, axes = plt.subplots(1, 2, figsize=(6, 4), sharex=True, sharey=True)
    freq = freq_uniform
    for col, excitation in enumerate([1, 2]):
        ax = axes[col]
        ref_key = ("elec", excitation)
        if ref_key not in e_uniform:
            continue
        ref = e_uniform[ref_key]
        scale = np.sqrt(np.mean(ref**2))  # RMS of reference, sample-number independent
        for k, tol in enumerate(ADAPTIVE_TOLS):
            if tol not in freq_adaptive or ref_key not in e_adaptive[tol]:
                continue
            assert np.all(freq_adaptive[tol] == freq_uniform)
            err = np.abs(e_adaptive[tol][ref_key] - ref) / scale
            ax.plot(
                freq,
                err,
                color=tol_palette[k],
                linewidth=0.2,
                label=f"{tol:.0e}",
                marker="o",
                markersize=2,
            )
            ax.axhline(
                tol,
                0,
                0.8,
                color=tol_palette[k],
                linewidth=1,
                linestyle="--",
                zorder=-5,
            )
            pivots_list = sampled_freqs.get(tol, [])
            all_pivots = np.concatenate(pivots_list) if pivots_list else np.array([])
            ax.plot(
                all_pivots,
                np.ones_like(all_pivots) * 2e-13 * tol**0.25,
                color=tol_palette[k],
                linewidth=0,
                alpha=0.75,
                marker="D",
                markersize=4,
            )
        ax.axhline(
            1e-12,
            color="black",
            linewidth=1,
            alpha=1,
            marker="",
            zorder=-5,
            linestyle=":",
        )
        add_eigenmode_lines(ax)
        ax.text(
            0.97,
            0.97,
            f"Excitation {excitation}",
            transform=ax.transAxes,
            ha="right",
            va="top",
            fontsize=9,
        )
        ax.set_xlim(*_FREQ_LIM)
        ax.set_xticks(_FREQ_TICKS)
        ax.set_ylim(5e-15, 10)
        ax.set_yscale("log")
        ax.set_xlabel("f (GHz)")
    # axes[0].legend(
    #     loc="upper center", fontsize=8, frameon=False, ncol=len(ADAPTIVE_TOLS)
    # )
    axes[0].legend(
        loc="upper right",
        bbox_to_anchor=(1.0, 0.88),
        fontsize=8,
        frameon=False,
        ncol=1,
        handlelength=1.0,
        handletextpad=0.4,
        columnspacing=0.3,
    )
    fig.suptitle(
        "Error in Domain Energy\n"
        "$\\vert E_\\mathrm{elec,adaptive} - E_\\mathrm{elec,uniform}\\vert"
        " /  \\vert\\vert E_\\mathrm{elec,uniform}\\vert\\vert_\\mathrm{RMS}$"
    )
    adjust_golden_layout_1x2(fig, axes, hspace_in=0.1)
    plt.savefig(DOCS_ASSET_DIR / "driven_ua_transmon_energy_adaptive_sweep.svg")
    plt.show()


plot_energy_adaptive_sweep(
    freq_uniform,
    e_uniform,
    e_adaptive,
    freq_adaptive,
    ADAPTIVE_TOLS,
    sampled_freqs,
    tol_palette,
    DOCS_ASSET_DIR,
    adjust_golden_layout_1x2,
)
