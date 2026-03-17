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

# CPW Lumped-Port Driven Solver

# Generates some of the plots in the "Driven Solver: Uniform vs Adaptive" tutorial.
#
# Run `cpw_tutorial_lumped_driven.jl` first to generate the simulation data. This script
# then parses the data. Output plots are written to the `docs/src/assets` folder.
#
# Run this script stand-alone from the base palace folder using `uv`. This will
# automatically create a temporary environment and download dependencies
#
# ```
# uv run --script examples/cpw/cpw_tutorial_lumped_driven_plots.py
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

## Configuration: Data created by cpw_tutorial_lumped_driven.jl


# Root directory of the CPW and this script
CPW_DIR = pathlib.Path("./examples/cpw")
DOCS_ASSET_DIR = pathlib.Path("./docs/src/assets/examples")

# Uniform sweep output directory
UNIFORM_DIR = CPW_DIR / "postpro/tutorial_driven_rom/driven_uniform_reference"

# Adaptive tolerances and their output directories
ADAPTIVE_TOLS = [1e-1, 1e-2, 1e-3, 1e-4, 1e-5]
ADAPTIVE_DIRS = {
    tol: CPW_DIR
    / f"postpro/tutorial_driven_rom/driven_adaptive_1e{int(round(np.log10(tol)))}"
    for tol in ADAPTIVE_TOLS
}

# S-parameter port pairs
PORT_PAIRS = [(1, 1), (2, 1), (3, 1), (4, 1)]


## Load Data


_RE_S_DB = re.compile(r"\|S\[(\d+)\]\[(\d+)\]\| \(dB\)")
_RE_E_COL = re.compile(r"E_(\w+) \(J\)")
_RE_SAMPLED_FREQS = re.compile(
    r" Sampled frequencies \(GHz\): ((?s:.)*?)\n Sample error", re.MULTILINE
)
_RE_SAMPLED_ERRORS = re.compile(
    r"^ Sample errors: ((?s:.)*?)\n Total offline phase", re.MULTILINE
)


def load_port_s(postpro_dir: pathlib.Path) -> tuple[np.ndarray, dict]:
    """Load port-S.csv; return (freq_GHz, {(i, j): complex S_ij array})."""
    df = pd.read_csv(postpro_dir / "port-S.csv")
    df.columns = df.columns.str.strip()
    freq = df["f (GHz)"].values
    s = {}
    for col in df.columns:
        m = _RE_S_DB.match(col)
        if m:
            i, j = int(m.group(1)), int(m.group(2))
            # Convert magnitude and argument into complex number
            db = df[f"|S[{i}][{j}]| (dB)"].values
            ang = df[f"arg(S[{i}][{j}]) (deg.)"].values
            s[(i, j)] = 10 ** (db / 20.0) * np.exp(1j * np.deg2rad(ang))
    return freq, s


def load_domain_e(postpro_dir: pathlib.Path) -> tuple[np.ndarray, dict]:
    """Load domain-E.csv; return (freq_GHz, {(quantity, port): array})."""
    df = pd.read_csv(postpro_dir / "domain-E.csv")
    df.columns = df.columns.str.strip()
    freq = df["f (GHz)"].values
    e = {}
    for col in df.columns:
        m = _RE_E_COL.match(col)
        if m:
            e[m.group(1)] = df[col].values
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
    Returns list of arrays, one per excitation. First two entries are inf
    (the fmin/fmax initialization points have no prior ROM to compare against).
    """
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
        print(f"WARNING: {d} not found — run cpw_tutorial_lumped_driven.jl first")
        continue
    freq_adaptive[tol], s_adaptive[tol] = load_port_s(d)
    _, e_adaptive[tol] = load_domain_e(d)
    log_path = d / "palace.log"
    sampled_freqs[tol] = parse_log_pivots(log_path)
    sampled_errors[tol] = parse_log_errors(log_path)
    n_samp = sum(len(a) for a in sampled_freqs[tol])
    print(
        f"Adaptive tol={tol:.0e}: {len(freq_adaptive[tol])} frequency points ({n_samp} adaptive samples)"
    )


# Common Plotting Options

tol_palette = plt.cm.plasma(np.linspace(0.1, 0.9, len(ADAPTIVE_TOLS)))

GOLDEN_RATIO = (1 + np.sqrt(5)) / 2


def adjust_golden_layout(fig, ax, pad_in: float = 0.1) -> None:
    """Adjust a single-axes figure: equal left/right margins (driven by y-axis decorations),
    golden-ratio axes, and vertical spacing from actual text extents."""
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
    """Adjust an nrows×ncols subplot figure so that:
      - Left margin = y-axis decoration width + pad_in; right margin mirrors it.
      - Each axes has width/height = golden ratio.
      - Vertical/horizontal inter-axes spacing is hspace_in / vspace_in (inches).
      - A suptitle (if set) is accommodated in the top margin.

    Call after all decorations (titles, labels, suptitle) are set.
    axes must be a 2-D numpy array of Axes objects (from plt.subplots).
    """
    dpi = fig.dpi
    fig_width = fig.get_figwidth()
    nrows, ncols = axes.shape

    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()

    # Left decoration from top-left axes (has y-labels)
    ax_tl = axes[0, 0]
    axes_bb = ax_tl.get_window_extent(renderer)
    tight_bb = ax_tl.get_tightbbox(renderer)
    left_extra_in = (axes_bb.x0 - tight_bb.x0) / dpi

    # Bottom decoration from bottom-left axes (has x-labels)
    ax_bl = axes[nrows - 1, 0]
    axes_bb_bl = ax_bl.get_window_extent(renderer)
    tight_bb_bl = ax_bl.get_tightbbox(renderer)
    bottom_extra_in = (axes_bb_bl.y0 - tight_bb_bl.y0) / dpi

    # Top decoration: max across top-row axes (titles)
    top_extra_in = max(
        (ax.get_tightbbox(renderer).y1 - ax.get_window_extent(renderer).y1) / dpi
        for ax in axes[0]
    )

    # Suptitle height (if present)
    suptitle_h_in = 0.0
    if fig._suptitle is not None:
        suptitle_h_in = (
            fig._suptitle.get_window_extent(renderer).height / dpi + pad_in * 0.5
        )

    # Symmetric left/right margins
    margin_lr_in = left_extra_in + pad_in
    bottom_space_in = bottom_extra_in + pad_in
    top_space_in = top_extra_in + pad_in + suptitle_h_in

    # Golden-ratio axes size from available width
    axes_width_in = (fig_width - 2 * margin_lr_in - hspace_in * (ncols - 1)) / ncols
    axes_height_in = axes_width_in / GOLDEN_RATIO

    # Total figure height
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


# S-Parameters

### Uniform Sweep Reference Data


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
                label="Uniform Driven Solver",
                zorder=10,
                marker="o",
                markersize=2,
            )
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
        ax.set_xticks(np.linspace(2, 32, 6))
        ax.set_ylim(-105, 5)
    axes[0, 1].legend(loc="lower right", fontsize=8, frameon=False)
    fig.suptitle("S-Parameters Magnitude $\\vert S_{{ij}}\\vert$ (dB)")
    adjust_golden_layout_2x2(fig, axes, hspace_in=0.1, vspace_in=0.15)
    plt.savefig(DOCS_ASSET_DIR / "driven_ua_cpw_domain_sparam_uniform.svg")
    plt.show()


plot_smat_uniform(
    freq_uniform, s_uniform, PORT_PAIRS, DOCS_ASSET_DIR, adjust_golden_layout_2x2
)


## S-Parameters: Convergence


def plot_smat_adaptive_pointwise(
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
        for k, tol in enumerate(ADAPTIVE_TOLS):
            if tol not in s_adaptive or (i, j) not in s_adaptive[tol]:
                continue
            assert np.all(freq_adaptive[tol] == freq_uniform)
            err = np.abs(s_adaptive[tol][(i, j)] - ref) / np.abs(ref)
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
                tol, color=tol_palette[k], linewidth=1, linestyle="--", zorder=-5
            )
        ax.axhline(1e-12, color="black", linewidth=1, linestyle=":", zorder=-5)
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
        ax.set_xlim(1, 33)
        ax.set_xticks(np.linspace(2, 32, 6))
        ax.set_ylim(1.1e-13, 900)
        ax.set_yscale("log")
    axes[1, 1].legend(
        loc="center",
        bbox_to_anchor=(0.5, 0.25),
        fontsize=8,
        frameon=False,
        ncol=3,
        columnspacing=0.5,
    )
    fig.suptitle(
        "Magnitude of Pointwise Relative Error in S-Parameters\n"
        "$\\vert S_{{\\mathrm{{adaptive}}}} - S_{{\\mathrm{{uniform}}}}\\vert"
        " /  \\vert S_{{\\mathrm{{uniform}}}}\\vert$"
    )
    adjust_golden_layout_2x2(fig, axes, hspace_in=0.1, vspace_in=0.15)
    plt.savefig(DOCS_ASSET_DIR / "driven_ua_cpw_domain_sparam_adaptive_pointwise.svg")
    plt.show()


plot_smat_adaptive_pointwise(
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
        scale = np.sqrt(
            np.mean(np.square(np.abs(ref)))
        )  # RMS of |S_ref|, sample-number independent
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
                tol, color=tol_palette[k], linewidth=1, linestyle="--", zorder=-5
            )
            pivots_list = sampled_freqs.get(tol, [])
            all_pivots = np.concatenate(pivots_list) if pivots_list else np.array([])
            if all_pivots.size:
                ax.plot(
                    all_pivots,
                    np.ones_like(all_pivots) * 2e-14 * tol**0.25,
                    color=tol_palette[k],
                    linewidth=0,
                    marker="D",
                    markersize=4,
                    alpha=0.75,
                )
        ax.axhline(1e-12, color="black", linewidth=1, linestyle=":", zorder=-5)
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
        ax.set_xlim(1, 33)
        ax.set_xticks(np.linspace(2, 32, 6))
        ax.set_ylim(1.1e-16, 50)
        ax.set_yscale("log")
    axes[1, 1].legend(
        loc="center",
        bbox_to_anchor=(0.5, 0.35),
        fontsize=8,
        frameon=False,
        ncol=3,
        columnspacing=0.5,
    )
    fig.suptitle(
        "Normalized Absolute Error in S-Parameters\n"
        "$\\vert S_{{\\mathrm{{adaptive}}}} - S_{{\\mathrm{{uniform}}}}\\vert"
        " /  \\vert\\vert S_{{\\mathrm{{uniform}}}}\\vert\\vert_{{\\mathrm{{RMS}}}}$"
    )
    adjust_golden_layout_2x2(fig, axes, hspace_in=0.1, vspace_in=0.15)
    plt.savefig(DOCS_ASSET_DIR / "driven_ua_cpw_domain_sparam_adaptive_rms.svg")
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


## Electric Energy

### Uniform Sweep Reference Data


def plot_energy_uniform(freq_uniform, e_uniform, DOCS_ASSET_DIR, adjust_golden_layout):
    fig, ax = plt.subplots(1, 1, figsize=(6, 4))
    ax.plot(
        freq_uniform,
        e_uniform["elec"] * 1e12,
        "k-",
        linewidth=1,
        label="Uniform Driven Solver",
        marker="o",
        markersize=2,
    )
    ax.set_title(f"Domain Energy $E_{{\\mathrm{{elec}}}}$ ($\\mathrm{{pJ}}$)")
    ax.set_xlabel("f (GHz)")
    ax.set_xlim(1, 33)
    ax.set_xticks(np.linspace(2, 32, 6))
    ax.legend(loc="upper right", fontsize=8, frameon=False)
    adjust_golden_layout(fig, ax)
    plt.savefig(DOCS_ASSET_DIR / "driven_ua_cpw_domain_energy_uniform.svg")
    plt.show()


plot_energy_uniform(freq_uniform, e_uniform, DOCS_ASSET_DIR, adjust_golden_layout)


## Electric Energy: Convergence

### Plot Single Representative Energy Converge Plot


def plot_energy_adaptive_single(
    freq_uniform,
    e_uniform,
    e_adaptive,
    freq_adaptive,
    ADAPTIVE_TOLS,
    sampled_freqs,
    tol_palette,
    DOCS_ASSET_DIR,
    adjust_golden_layout,
):
    fig, ax = plt.subplots(1, 1, figsize=(6, 4))
    freq = freq_uniform
    ref = e_uniform["elec"]
    for k, tol in enumerate(ADAPTIVE_TOLS[:1]):
        assert np.all(freq_adaptive[tol] == freq_uniform)
        err = np.abs(e_adaptive[tol]["elec"] - ref) / np.abs(ref)
        ax.plot(
            freq,
            err,
            color=tol_palette[k],
            linewidth=0.2,
            label=f"tol={tol:.0e}",
            marker="o",
            markersize=2,
        )
        ax.axhline(tol, color=tol_palette[k], linewidth=1, linestyle="--", zorder=-5)
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
        1e-12, color="black", linewidth=1, alpha=1, marker="", zorder=-5, linestyle=":"
    )
    ax.set_xlim(1, 33)
    ax.set_xticks(np.linspace(2, 32, 6))
    ax.set_ylim(2e-14, 3)
    ax.set_yscale("log")
    ax.legend(loc="upper right", fontsize=8, frameon=False)
    ax.set_title(
        "Error in Domain Energy\n"
        "$\\vert E_{{\\mathrm{{elec,adaptive}}}} - E_{{\\mathrm{{elec,uniform}}}}\\vert"
        " /  \\vert E_{{\\mathrm{{elec,uniform}}}}\\vert$"
    )
    ax.set_xlabel("f (GHz)")
    adjust_golden_layout(fig, ax)
    plt.savefig(DOCS_ASSET_DIR / "driven_ua_cpw_domain_energy_adaptive_single.svg")
    plt.show()


plot_energy_adaptive_single(
    freq_uniform,
    e_uniform,
    e_adaptive,
    freq_adaptive,
    ADAPTIVE_TOLS,
    sampled_freqs,
    tol_palette,
    DOCS_ASSET_DIR,
    adjust_golden_layout,
)


def plot_energy_adaptive_sweep(
    freq_uniform,
    e_uniform,
    e_adaptive,
    freq_adaptive,
    ADAPTIVE_TOLS,
    sampled_freqs,
    tol_palette,
    DOCS_ASSET_DIR,
    adjust_golden_layout,
):
    fig, ax = plt.subplots(1, 1, figsize=(6, 4))
    freq = freq_uniform
    ref = e_uniform["elec"]
    for k, tol in enumerate(ADAPTIVE_TOLS):
        assert np.all(freq_adaptive[tol] == freq_uniform)
        err = np.abs(e_adaptive[tol]["elec"] - ref) / np.abs(ref)
        ax.plot(
            freq,
            err,
            color=tol_palette[k],
            linewidth=0.2,
            label=f"tol={tol:.01e}",
            marker="o",
            markersize=2,
        )
        ax.axhline(tol, color=tol_palette[k], linewidth=1, linestyle="--", zorder=-5)
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
        1e-12, color="black", linewidth=1, alpha=1, marker="", zorder=-5, linestyle=":"
    )
    ax.set_xlim(1, 33)
    ax.set_xticks(np.linspace(2, 32, 6))
    ax.set_ylim(5e-15, 3)
    ax.set_yscale("log")
    ax.legend(loc="upper center", fontsize=8, frameon=False, ncol=len(ADAPTIVE_TOLS))
    ax.set_title(
        "Error in Domain Energy\n"
        "$\\vert E_{{\\mathrm{{elec,adaptive}}}} - E_{{\\mathrm{{elec,uniform}}}}\\vert"
        " /  \\vert E_{{\\mathrm{{elec,uniform}}}}\\vert$"
    )
    ax.set_xlabel("f (GHz)")
    adjust_golden_layout(fig, ax)
    plt.savefig(DOCS_ASSET_DIR / "driven_ua_cpw_domain_energy_adaptive_sweep.svg")
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
    adjust_golden_layout,
)

## Adaptive Solver Convergence Curve
#
# MRI error indicator at each greedy sample point vs sample number, for each `AdaptiveTol`.
# The first two samples (at $f_\mathrm{min}$ and $f_\mathrm{max}$) are initialization points with
# no prior ROM, so their error is $+\infty$.


def plot_adaptive_convergence_curve(
    ADAPTIVE_TOLS, sampled_errors, tol_palette, DOCS_ASSET_DIR, adjust_golden_layout
):
    fig, ax = plt.subplots(1, 1, figsize=(6, 4))
    errs_tightest = None
    for k, tol in enumerate(ADAPTIVE_TOLS):
        errs_list = sampled_errors.get(tol, [])
        if not errs_list:
            continue
        errs = errs_list[0]  # single excitation
        finite = np.isfinite(errs)
        if not np.any(finite):
            continue
        n = np.arange(1, len(errs) + 1)
        ax.plot(
            n[finite],
            errs[finite],
            color=tol_palette[k],
            linewidth=1,
            label=f"tol={tol:.0e}",
            marker="o",
            markersize=4,
            zorder=np.log(tol),
        )
        ax.axhline(tol, color=tol_palette[k], linewidth=1, linestyle="--", zorder=-5)
        errs_tightest = errs  # retains value from last processed tolerance

    if errs_tightest is not None and len(errs_tightest) > 2:
        ax.plot(
            [1, 2, 3],
            [1, 1, errs_tightest[2]],
            marker="*",
            color="lightgrey",
            zorder=-10,
            linestyle="--",
        )

    ax.set_yscale("log")
    ax.set_xlabel("Sample number")
    ax.set_xticks(np.linspace(1, 13, 13))
    ax.legend(loc="upper right", fontsize=8, frameon=False, ncol=2)
    ax.set_title("Adaptive Solver Convergence: Error Indicator")
    adjust_golden_layout(fig, ax)
    plt.savefig(DOCS_ASSET_DIR / "driven_ua_cpw_adaptive_convergence_curve.svg")
    plt.show()


plot_adaptive_convergence_curve(
    ADAPTIVE_TOLS, sampled_errors, tol_palette, DOCS_ASSET_DIR, adjust_golden_layout
)
