#!/usr/bin/env python3
# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""grade_cells.py -- the scorer of record for a corrected thin-mesh window against its thick reference (fix round R-7).

One window = one thin CORRECTED device run graded against ONE reference per Type (SA / MS / MA_sharp / C) with the D4 bars
(SA / MS / MA_sharp 5 %, C 0.5 %; x = 100 (device / reference - 1)), PER SOURCE (electrostatic, every terminal i; decision 1210)
or PER FREQUENCY (Maxwell driven; decision 1213), and the cell letter FAIL-CLOSED over them: P only if every source /
frequency is P, F if any F, else I; the worst-|x| member is quoted with the range and the counts. The single-source /
single-frequency window is the same rule with one member.

  * ES device: the FINAL AMR cycle of the run (postpro/iterationNN whose error-indicators.csv equals the top-level one, else
    the top-level postpro: the cycle_dirs rule of record); per source i the Type value = the sum over the Type's interfaces of
    `p_surf postprocessed fixed-trace[k]` (ft = the result of record); C = terminal-C-corrected.csv `C fixed-trace[i][i]`;
    the off-diagonal entries are untargeted columns (decision 910), printed with an absolute bar vs C[i][i].
  * Maxwell device: the top-level postpro (the final solve); per frequency f the Type value = the sum over the Type's
    interfaces of `E_surf corrected (J)` / the p denominator (--p-denominator: E_norm = `E_norm corrected (J)` as the
    binary prints p_surf; E_norm-minus-E_cap = the junction-driven renormalisation of the card tools, A1.1.1); C_eff = 2
    E_norm_ft / |V_port|^2 (847 (B)) on ported lines.
  * MA_sharp of the DEVICE by convention (B) (USER 422 / AMENDMENT 16): (E_surf_MA_ft - sum_m F_m) + sum_m F_m (1 + t_m),
    F_m = the model's fabricated surface energy[k over the MA interfaces] at evaluation 0 of the source / frequency,
    t_m = the class tail of the model's name (TAILS below; a spatial cluster's own library record `MA.Summary.Deficit.
    Consistent.Median` with the bar [0, 2 x record]); the bracket from the tail bars gives INCONCLUSIVE when ft is inside the
    bar and the bracket is outside (the verdict of record).
  * the TAIL FORM (--tail-form; USER decision 1308 F-C5 (b), supervisor decisions 1310 / 1312 / 1314): `422` (the default until the
    P4 re-freeze) = the class tails of record above; `fixed-shell` = ONE estimator on both sides, both estimating the converged
    fixed-shell-R10 sharp energy S* (rings {10, 30, 50} nm = {R/190, 3R/190, R/38}, R = 1.9 um; 50 nm = t/2): the device class tails
    FIXED_SHELL_TAILS = S*_class / raw_class(library rung) - 1 (the factor mapping the library-rung raw the device applies onto S*),
    and the reference's shells REQUIRED to carry that ring set (FIXED_SHELL_RINGS_UM: every stage-E r10 reference of record does); a
    reference recorded with another ring set is REFUSED (exit 3) unless --reference-rings-information, which grades it and labels
    the MA_sharp cell INFORMATION (reference ring set != of record). The spatial clusters keep their library records under both
    forms (not re-measured under R10). Under this form a reference coarser than r5 p5 is not converged in S* (P3-FC5 STEP 10: r10 p4
    -> r5 p5 +2.6..+4.5 % on CD3), so its MA_sharp cell is labelled INFORMATION (not counted) until graded against an r5 p5-or-finer
    reference (decision 1314; reference_depth_information on a run reference; a values reference's MA_sharp_reference_depth, absent =
    "depth not recorded" = INFORMATION too, decision 1318).
  * MA_sharp of a thick REFERENCE run from its recorded MA shells (the shells catalog beside the reference config):
    `radial_ma_profile` of the library (examples/cpw3d_surface/spatial_coupon/qualify), the HEADLINE estimator -- the
    collect17 form of record, ported. A reference given as a VALUES file (p_SA / p_MS / p_MA_sharp / C per source or
    frequency; the Richardson references of refresh.json) is read as is.
  * the C FORM on a Maxwell ported line (--c-form; USER decision 1230 Q1): `static-limit` (the graded cell) = the fitted
    static limit C_s of C_eff(f) = C_s (1 + a f^2) over the frequencies, unweighted least squares on f^2 per side, with
    C_eff(f) per frequency printed beside; `direct` = C_eff(f) per frequency at the bar (decision 1214); `information` =
    the C cell is not counted (F-40). The ES C is always the diagonal capacitance per source.
  * the BAND check (PLAN-v3 4.2): x per Type inside block_C of the filed coverage-report.json (--band); every member in band.

The letters, tails, class rule and verdict are the 994 / 995 tabulators' (tabulate_rerun2.py: model_class / tail_of / TAILS /
verdict; collect17.py: ma_sharp) ported verbatim so that the 47 stage-E cells and the 995 Maxwell card are reproduced from the
stored CSVs (test_grade_cells.py). Pure Python; exit 0 when the grade was written.
"""

import argparse
import csv
import json
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
QUALIFY_DIR = os.path.normpath(os.path.join(HERE, "..", "cpw3d_surface", "spatial_coupon", "qualify"))
if QUALIFY_DIR not in sys.path:
    sys.path.insert(0, QUALIFY_DIR)
import radial_ma_profile as profile  # noqa: E402

TYPES = ("SA", "MS", "MA_sharp", "C")
D4 = {"SA": 5.0, "MS": 5.0, "MA_sharp": 5.0, "C": 0.5}
# Convention (B) class tails (fraction of the class MA) and their bars [lo, hi]: AMENDMENT 16 A16.4 (tabulate_rerun2.TAILS).
TAILS = {
    "isolated": (0.0678, 0.0600, 0.0726),
    "curved": (0.0678, 0.0600, 0.0726),
    "strip": (0.0678, 0.0600, 0.0726),
    "corner": (0.1891, 0.1369, 0.1891),
}
SPATIAL_DEFAULT_TAIL = (0.0302, 0.0200, 0.0340)
TAIL_FORMS = ("422", "fixed-shell")
# FIXED-SHELL-R10 (decision 1310 Q-A): the ring boundaries (um) both sides' Consistent estimator is read on; the last = t/2.
FIXED_SHELL_RINGS_UM = (0.01, 0.03, 0.05)
# The fixed-shell class tails (supervisor decisions 1312 / 1313 / 1314 implementing USER 1308 F-C5 (b); the P3-FC5 lane record):
# t_class = S*_class / raw_class(library rung) - 1, where S*_class = the CONVERGED fixed-shell-R10 sharp MA energy of the class
# coupon (P3-FC5 fc5/results/r10-tails.json) and raw_class(library rung) = the raw MA of the coupon rung that produced the library
# model the device applies (the tail is the factor mapping THAT raw onto S*). The bar [lo, hi] = the S* spread over the ladder's
# (h, p) points widened on both sides by the converged-reference R10 h-step FIXED_SHELL_REFERENCE_STEP; the bar applies to cells
# graded against a CONVERGED reference (r5 p5 or finer) -- an r10 p4-referenced cell carries no letter under this form (1314).
#   isolated: S* = 6.214692e-19 J (the 2D edge coupon, lc 2 / lc 1 nm x p4..p8: 10 points, min 6.211713e-19, max 6.217804e-19,
#             spread 0.098 %); raw(edge lc 2 nm, p4) = 5.652615e-19 J -> t = +9.94 % [9.70, 10.19].
#   corner:   S* = 2.646967e-21 J (the 3D 90-deg convex coupon, 10 nm x p3..p7: 5 points, min 2.634013e-21, max 2.653846e-21,
#             spread 0.749 %; the 20-nm rung's sharp is ring-1 sub-mesh noisy and not used); raw(corner 20 nm, p3) = 2.322675e-21 J
#             -> t = +13.96 % [13.21, 14.45].
#   curved / strip: the isolated tail by assumption (as the 422 form); the spatial clusters keep their own library record tails
#   under both forms (1314 Q2: "not re-measured under R10").
FIXED_SHELL_S_STAR_J = {"isolated": (6.214692e-19, 6.211713e-19, 6.217804e-19), "corner": (2.646967e-21, 2.634013e-21, 2.653846e-21)}
FIXED_SHELL_RAW_LIBRARY_RUNG_J = {"isolated": ("edge lc 2 nm p4", 5.652615e-19), "corner": ("corner 20 nm p3", 2.322675e-21)}
# The converged-reference R10 h-step (fraction of MA_sharp): the CD3 thick reference read with the R10 rings at r5 p5 -> r2.5 p5
# moves -0.16..-0.19 % on its four sources (P3-FC5 STEP 10; the r10 p4 -> r5 p5 step +2.6..+4.5 % is the unconverged rung's).
# CONVENTION (decision 1314 Q1, the numbers as ruled): the step is ADDED to the tail bar as is (a fraction of MA_sharp added to a
# fraction of raw, without the S*/raw factor 1.0994 / 1.1396 that would widen it to 0.0021 / 0.0022): lo / hi = S*_min,max / raw - 1 -/+ step.
FIXED_SHELL_REFERENCE_STEP = 0.0019


def fixed_shell_tails(reference_step=FIXED_SHELL_REFERENCE_STEP):
    """{class: (tail, lo, hi)} of the fixed-shell form from the constants above: tail = S* / raw - 1, lo = S*_min / raw - 1 - step,
    hi = S*_max / raw - 1 + step."""
    out = {}
    for cls, (s_mean, s_min, s_max) in FIXED_SHELL_S_STAR_J.items():
        raw = FIXED_SHELL_RAW_LIBRARY_RUNG_J[cls][1]
        out[cls] = (s_mean / raw - 1.0, s_min / raw - 1.0 - reference_step, s_max / raw - 1.0 + reference_step)
    out["curved"] = out["strip"] = out["isolated"]
    return out


FIXED_SHELL_TAILS = fixed_shell_tails()
# The MA shells catalog's attribute numbering (collect17.B): plane j's top face B + 100 + 20 j + 1, its foot rings
# B + 200 + 20 j + k, its top-side rings B + 300 + 20 j + k.
SHELLS_BASE = 500
C_FORMS = ("direct", "static-limit", "information")
P_DENOMINATORS = ("E_norm", "E_norm-minus-E_cap")
# The two implementation tolerances of the scorer of record (R-7; REFREEZE-ROWS-P1.md; review d1230 MINOR-3). Neither affects a
# record: the frequencies of record are exact decimals and no cell sits within 1e-12 pt of a band edge.
# FREQ_TOL_GHZ: a CSV row belongs to frequency F when |f - F| < 1e-6 GHz (the frequency-row match of the Maxwell readers).
FREQ_TOL_GHZ = 1e-6
# BAND_TOL_PT: a cell is in band when low - 1e-12 <= x <= high + 1e-12 pt (the band-edge comparison of band_check).
BAND_TOL_PT = 1e-12


class MissingInput(Exception):
    """A required input of the grade is absent (exit code 2; the file is named)."""


class ReferenceRingsNotOfRecord(Exception):
    """--tail-form fixed-shell on a reference whose recorded MA shell rings are not FIXED_SHELL_RINGS_UM (exit code 3; decision 1310 Q-D)."""


INFORMATION_RINGS = "INFORMATION (reference ring set != of record)"


def read_csv(path):
    """Palace CSV -> list of dict rows with stripped headers and string values ('#' comment lines skipped)."""
    if not os.path.isfile(path):
        raise MissingInput(path)
    with open(path, newline="") as handle:
        rows = [r for r in csv.reader(handle) if r and not r[0].lstrip().startswith("#")]
    header = [h.strip() for h in rows[0]]
    return [dict(zip(header, [v.strip() for v in r])) for r in rows[1:]]


def fnum(value):
    return float(value)


def pct(a, b):
    return 100.0 * (a / b - 1.0)


def verdict(value, limit, bracket=None):
    """tabulate_rerun2.verdict: |x| > bar FAIL; ft inside with the bracket outside INCONCLUSIVE; else PASS."""
    if value != value:
        return "n/a"
    if abs(value) > limit:
        return "FAIL"
    if bracket and (abs(bracket[0]) > limit or abs(bracket[1]) > limit):
        return "INCONCLUSIVE (ft inside, bracket outside)"
    return "PASS"


def model_class(name):
    """tabulate_rerun2.model_class: the convention-(B) class of a catalog model name."""
    n = name.lower()
    if n.startswith("spatialedgecluster"):
        m = re.search(r"edgecount-(\d+)", n)
        return f"spatial cluster {m.group(1)}-edge" if m else "spatial cluster"
    if n.startswith("parallel-edge-cluster") or n.startswith("paralleledgecluster"):
        return "stack (ParallelEdgeCluster)"
    if n.startswith("isolated-edge"):
        return "isolated edge"
    if n.startswith("convex-corner") or n.startswith("concave-corner"):
        return "corner"
    if n.startswith("curved-same-conductor-strip"):
        return "curved strip"
    if n.startswith("curved-edge"):
        return "curved edge"
    if n.startswith("same-conductor-strip"):
        return "strip pair"
    if n.startswith("same-conductor-gap") or n.startswith("different-conductor-gap"):
        return "gap pair"
    return "other"


def tail_of(name, spatial_tails, tail_form="422"):
    """tabulate_rerun2.tail_of -> ((tail, lo, hi), source text); the class table by `tail_form` (TAILS = 422; FIXED_SHELL_TAILS =
    fixed-shell), the class map and the spatial clusters' library records the same under both forms."""
    table = TAILS if tail_form == "422" else FIXED_SHELL_TAILS
    cls = model_class(name)
    if cls.startswith("spatial cluster"):
        key = name.split("@")[0]
        for k, v in spatial_tails.items():
            if key.startswith(k) or k.startswith(key):
                return v, "library record (per coupon)"
        return SPATIAL_DEFAULT_TAIL, "spatial median 3.02 (no per-coupon record found)"
    src = " (FIXED-SHELL-R10: S* / raw(library rung) - 1, decisions 1312 / 1314)" if tail_form == "fixed-shell" else ""
    if cls == "corner":
        return table["corner"], ("corner set C p3 / bias bar" if tail_form == "422" else "S*(corner 10 nm ladder) / raw(corner 20 nm p3)") + src
    if cls in ("curved edge", "curved strip"):
        return table["curved"], "isolated tail" + src
    if cls in ("strip pair", "gap pair", "stack (ParallelEdgeCluster)"):
        return table["strip"], "isolated tail by assumption" + src
    return table["isolated"], ("2D straight coupon p4" if tail_form == "422" else "S*(edge lc 2 / lc 1 ladder) / raw(edge lc 2 nm p4)") + src


def spatial_tails(library_path):
    """{model Name: (median, q1, q3)} of MA.Summary.Deficit.Consistent for every library model carrying one."""
    if library_path is None:
        return {}
    with open(library_path) as handle:
        models = json.load(handle)["Models"]
    out = {}
    for m in models:
        if "MA" in m and "Summary" in m["MA"]:
            d = m["MA"]["Summary"]["Deficit"]["Consistent"]
            out[m["Name"]] = (d["Median"], d["Quartiles"][0], d["Quartiles"][1])
    return out


def model_catalog(postpro, catalog_path=None):
    """{model index: Name} from postpro/palace.json's SurfaceResponse.ModelCatalog, else from --model-catalog (a JSON list
    of {Index, Name}: the extracted catalog of a run whose palace.json was not kept)."""
    pj = os.path.join(postpro, "palace.json")
    if os.path.isfile(pj):
        with open(pj) as handle:
            cat = json.load(handle)["SurfaceResponse"]["ModelCatalog"]
        return {int(m["Index"]): m["Name"] for m in cat}, pj
    if catalog_path is None:
        raise MissingInput(pj)
    with open(catalog_path) as handle:
        cat = json.load(handle)
    return {int(m["Index"]): m["Name"] for m in cat}, catalog_path


def interfaces_by_type(config_path):
    """{'SA': [indices], 'MS': [...], 'MA': [...]} of the run's Postprocessing.Dielectric entries (the decision-344 rule:
    the window value of a Type = the sum over EVERY interface of that Type; diagnostic indices >= 3000 excluded)."""
    with open(config_path) as handle:
        cfg = json.load(handle)
    out = {"SA": [], "MS": [], "MA": []}
    for d in cfg["Boundaries"]["Postprocessing"]["Dielectric"]:
        if int(d["Index"]) >= 3000:
            continue
        out[d["Type"]].append(int(d["Index"]))
    return out, cfg


def ma_sharp_convention_b(e_surf_ma, model_energies, catalog, sp_tails, tail_form="422"):
    """The device MA_sharp (B): -> {E_out_plus_uncovered_J, models_MA_E_J, n_models, device_sharp_E_J, bracket_E_J, per_class}; the class
    tails t_m by `tail_form` (tail_of)."""
    base = e_surf_ma - sum(model_energies.values())
    cen = lo = hi = base
    per_class = {}
    for m, fm in model_energies.items():
        name = catalog.get(m, f"model{m}")
        (t, tlo, thi), srcname = tail_of(name, sp_tails, tail_form)
        cls = model_class(name)
        if cls.startswith("spatial cluster"):
            tlo, thi = 0.0, 2.0 * t
        cen += fm * (1.0 + t)
        lo += fm * (1.0 + tlo)
        hi += fm * (1.0 + thi)
        pc = per_class.setdefault(cls, {"E_J": 0.0, "tail_E_J": 0.0, "tail_src": srcname, "models": []})
        pc["E_J"] += fm
        pc["tail_E_J"] += fm * t
        pc["models"].append(name)
    return {"E_out_plus_uncovered_J": base, "models_MA_E_J": sum(model_energies.values()), "n_models": len(model_energies),
            "device_sharp_E_J": cen, "bracket_E_J": sorted([lo, hi]), "tail_form": tail_form,
            "device_tail_pct_of_raw_MA": (100.0 * (cen - e_surf_ma) / e_surf_ma) if e_surf_ma else None, "per_class": per_class}


def rings_of_record(rings_um):
    """True when a reference's MA shell ring boundaries (um, outer radii) are FIXED_SHELL_RINGS_UM (to 1e-9 um)."""
    return len(rings_um) == len(FIXED_SHELL_RINGS_UM) and all(abs(a - b) < 1e-9 for a, b in zip(rings_um, FIXED_SHELL_RINGS_UM))


def check_reference_rings(rings_um, where, tail_form, rings_information):
    """The fixed-shell form's reference-side rule (decision 1310 Q-D): the rings must be of record; otherwise refuse, or -- with
    --reference-rings-information -- grade and return the INFORMATION label for the MA_sharp cell."""
    if tail_form != "fixed-shell" or rings_of_record(rings_um):
        return None
    if not rings_information:
        raise ReferenceRingsNotOfRecord(f"{where}: recorded MA shell rings {list(rings_um)} um != FIXED_SHELL_RINGS_UM {list(FIXED_SHELL_RINGS_UM)} "
                                        "(decision 1310 Q-D: refused; --reference-rings-information grades it as INFORMATION)")
    return INFORMATION_RINGS


# Under the fixed-shell form a reference is converged in S* at r5 p5 or finer (decision 1314: the CD3 ladder read with the R10 rings
# moves +2.6..+4.5 % from r10 p4 to r5 p5 and <= 0.19 % from r5 p5 to r2.5 p5); its depth = (its first graded row height r, its order).
FIXED_SHELL_CONVERGED_R_UM = 0.005
FIXED_SHELL_CONVERGED_ORDER = 5
INFORMATION_DEPTH = "INFORMATION (reference not converged in S*: coarser than r5 p5)"
# A VALUES-file reference whose depth is not recorded (no MA_sharp_reference_depth) is not shown converged: fail-closed by label
# (decision 1318; the rings rule keeps its exit-3 / override form).
INFORMATION_DEPTH_UNRECORDED = "INFORMATION (reference depth not recorded)"


def depth_information(depth, tail_form):
    """The fixed-shell form's reference-depth rule on a depth {"r_um", "order"} (or None = not recorded): -> the INFORMATION label or
    None (converged: r_um <= FIXED_SHELL_CONVERGED_R_UM and order >= FIXED_SHELL_CONVERGED_ORDER). Always None under the 422 form."""
    if tail_form != "fixed-shell":
        return None
    if not depth:
        return INFORMATION_DEPTH_UNRECORDED
    if float(depth["r_um"]) <= FIXED_SHELL_CONVERGED_R_UM + 1e-12 and int(depth["order"]) >= FIXED_SHELL_CONVERGED_ORDER:
        return None
    return INFORMATION_DEPTH


def reference_depth_information(shells_catalog, reference_run, tail_form):
    """The fixed-shell form's reference-depth rule (decision 1314) on a RUN reference: -> ({"r_um", "order"}, the INFORMATION label or
    None); the row height from the shells catalog's RowHeightsUm[0], the order from the reference's postpro/config_resolved.json
    Solver.Order. Under the 422 form nothing is read: -> (None, None)."""
    if tail_form != "fixed-shell":
        return None, None
    with open(os.path.join(reference_run, "postpro", "config_resolved.json")) as handle:
        order = int(json.load(handle)["Solver"]["Order"])
    depth = {"r_um": shells_catalog["RowHeightsUm"][0], "order": order}
    return depth, depth_information(depth, tail_form)


def ma_sharp_from_shells(shells_catalog, p_all, e_elec, edge_out, ma_interfaces=(3,), tail_form="422", rings_information=False):
    """The thick reference's MA_sharp from its recorded MA shells (collect17.ma_sharp, A15.2 (4): rings summed over the
    planes per edge kind; the HEADLINE estimator of radial_ma_profile). p_all = {interface index: p_surf} of the reference
    row, e_elec its E_elec, edge_out = {interface: {R_um (rounded 9): E_out (J)}} of the same source / frequency,
    ma_interfaces the window's MA interface indices (the raw MA total = their sum; the record form's single interface 3).
    Under tail_form fixed-shell the recorded rings must be FIXED_SHELL_RINGS_UM (check_reference_rings)."""
    E = {k: v * e_elec for k, v in p_all.items()}
    rings = shells_catalog["Rings"]
    information = check_reference_rings([r["OuterRadius"] for r in rings], "shells catalog", tail_form, rings_information)
    nk = len(rings)
    top = {k["Ring"]: 0.0 for k in rings}
    bottom = {k["Ring"]: 0.0 for k in rings}
    top_face = {k["Ring"]: 0.0 for k in rings}
    planes = []
    for j, pl in enumerate(shells_catalog["Planes"]):
        i_top = SHELLS_BASE + 100 + 20 * j + 1
        shells = [e for e in shells_catalog["Entries"] if e["Index"] == i_top][0]["Shells"]
        prev = E[i_top]
        face = {}
        for k, d in enumerate(shells, start=1):
            e_out = edge_out[i_top][round(d, 9)]
            face[k] = prev - e_out
            prev = e_out
        for k in range(1, nk + 1):
            top_face[k] += face[k]
            top[k] += face[k] + E[SHELLS_BASE + 300 + 20 * j + k]
            bottom[k] += E[SHELLS_BASE + 200 + 20 * j + k]
        planes.append({"Plane": pl["name"], "TopFaceRings_J": face,
                       "SideFromTop_J": {k: E[SHELLS_BASE + 300 + 20 * j + k] for k in range(1, nk + 1)},
                       "SideFromFoot_J": {k: E[SHELLS_BASE + 200 + 20 * j + k] for k in range(1, nk + 1)},
                       "TopFaceTotal_J": E[i_top], "SideTotal_J": E[SHELLS_BASE + 100 + 20 * j + 2]})
    ma_total = sum(E[k] for k in ma_interfaces)

    def by_ring(d):
        return {k: (0.0 if k == 1 else rings[k - 2]["OuterRadius"], rings[k - 1]["OuterRadius"], d[k]) for k in range(1, nk + 1)}

    prof = {"top": profile.profile_of(by_ring(top)), "bottom": profile.profile_of(by_ring(bottom))}
    remainders = {}
    for name in profile.remainder_names():
        total = 0.0
        ok = True
        for kind in ("top", "bottom"):
            est = prof[kind]["Estimates"].get(profile.kind_estimator(name, kind))
            if est is None:
                ok = False
                break
            total += est["Remainder"]
        remainders[name] = total if ok else None
    tail = remainders.get(profile.HEADLINE_ESTIMATOR)
    far = ma_total - sum(top.values()) - sum(bottom.values())
    return {"MA_raw_J": ma_total, "MA_tail_J": tail, "MA_sharp_J": (ma_total + tail) if tail is not None else None,
            "TailOverRaw_pct": (100.0 * tail / ma_total) if tail is not None else None,
            "Deficit_pct": {n: (100.0 * r / ma_total if r is not None else None) for n, r in remainders.items()},
            "FarShare_pct": 100.0 * far / ma_total, "Rings_um": [r["OuterRadius"] for r in rings], "rings_information": information,
            "Ring1Factor": {kind: (prof[kind]["Estimates"][profile.kind_estimator(profile.HEADLINE_ESTIMATOR, kind)] or {}).get("Ring1ResolvedOverModel")
                            for kind in ("top", "bottom")},
            "Rule": profile.HEADLINE_ESTIMATOR + ": top Theory@2 + bottom Fit2-4 (radial_ma_profile.py of the library, imported)"}


def shells_catalog_path(reference_run, shells_dir):
    """<shells_dir>/<basename(reference run)>-shells-catalog.json (the catalog recorded beside the reference config)."""
    return os.path.join(shells_dir, os.path.basename(os.path.normpath(reference_run)) + "-shells-catalog.json")


def edge_out_rows(rows, selector):
    """surface-Q-edge.csv rows passing `selector` -> {interface: {R_um (rounded 9): E_out (J)}}."""
    out = {}
    for r in rows:
        if not selector(r):
            continue
        out.setdefault(int(fnum(r["interface"])), {})[round(fnum(r["R (m)"]) * 1e6, 9)] = fnum(r["E_out (J)"])
    return out


def final_cycle_dir(run):
    """The ES device's final AMR cycle (tabulate_rerun2.cycle_dirs): the iterationNN dirs holding surface-Q-corrected.csv;
    the top-level postpro appended when its error-indicators.csv differs from the last saved iteration's (or none is
    saved). -> (final dir, cycle count)."""
    post = os.path.join(run, "postpro")
    if not os.path.isdir(post):
        raise MissingInput(post)
    cycles = []
    for name in sorted(os.listdir(post), key=lambda s: int(re.sub(r"\D", "", s) or 0)):
        if re.fullmatch(r"iteration\d+", name) and os.path.isfile(os.path.join(post, name, "surface-Q-corrected.csv")):
            cycles.append(os.path.join(post, name))
    if cycles:
        last_ind = os.path.join(cycles[-1], "error-indicators.csv")
        final_ind = os.path.join(post, "error-indicators.csv")
        if os.path.isfile(last_ind) and os.path.isfile(final_ind):
            with open(last_ind) as a, open(final_ind) as b:
                if a.read() == b.read():
                    return cycles[-1], len(cycles)
    if os.path.isfile(os.path.join(post, "surface-Q-corrected.csv")):
        cycles.append(post)
    if not cycles:
        raise MissingInput(os.path.join(post, "surface-Q-corrected.csv"))
    return cycles[-1], len(cycles)


def row_of_source(rows, src):
    out = [r for r in rows if int(fnum(r["i"])) == src]
    if len(out) != 1:
        raise ValueError(f"source {src}: {len(out)} rows (expected 1)")
    return out[0]


def offdiag(row, prefix, src):
    """{j: value} of the `<prefix>[i][j]` columns with j != src."""
    out = {}
    for k, v in row.items():
        m = re.match(re.escape(prefix) + r"\[i\]\[(\d+)\]", k)
        if m and int(m.group(1)) != src:
            out[int(m.group(1))] = fnum(v)
    return out


def es_device_source(final, src, iface, catalog, sp_tails, tail_form="422"):
    """The ES device's per-Type values for source `src` (row i == src of every CSV of the final cycle)."""
    q = read_csv(os.path.join(final, "surface-Q-corrected.csv"))
    row = row_of_source(q, src)
    out = {"source": src, "E_elec_ft_J": fnum(row["E_elec postprocessed fixed-trace (J)"]), "E_elec_raw_J": fnum(row["E_elec raw (J)"])}
    for T, ks in iface.items():
        out[f"{T}_p_ft"] = sum(fnum(row[f"p_surf postprocessed fixed-trace[{k}]"]) for k in ks)
        out[f"{T}_p_raw"] = sum(fnum(row[f"p_surf raw[{k}]"]) for k in ks)
        out[f"{T}_E_ft_J"] = sum(fnum(row[f"E_surf postprocessed fixed-trace[{k}] (J)"]) for k in ks)
        out[f"{T}_E_raw_J"] = sum(fnum(row[f"E_surf raw[{k}] (J)"]) for k in ks)
    edge = read_csv(os.path.join(final, "surface-Q-edge.csv"))
    for T, ks in iface.items():
        out[f"{T}_E_out_J"] = sum(fnum(r["E_out (J)"]) for r in edge if int(fnum(r["i"])) == src and int(fnum(r["interface"])) in ks)
    cc = row_of_source(read_csv(os.path.join(final, "terminal-C-corrected.csv")), src)
    out["C_ft_F"] = fnum(cc[f"C fixed-trace[i][{src}] (F)"])
    out["C_ft_offdiag_F"] = offdiag(cc, "C fixed-trace", src)
    cr = row_of_source(read_csv(os.path.join(final, "terminal-C.csv")), src)
    out["C_raw_F"] = fnum(cr[f"C[i][{src}] (F)"])
    me = read_csv(os.path.join(final, "surface-response-model-energy.csv"))
    F = {}
    for r in me:
        if int(fnum(r["source"])) != src or int(fnum(r["evaluation"])) != 0:
            continue
        F[int(fnum(r["model"]))] = sum(fnum(r[f"fabricated surface energy[{k}] (J)"]) for k in iface["MA"])
    sharp = ma_sharp_convention_b(out["MA_E_ft_J"], F, catalog, sp_tails, tail_form)
    sharp["device_sharp_p"] = sharp["device_sharp_E_J"] / out["E_elec_ft_J"]
    sharp["bracket_p"] = sorted([b / out["E_elec_ft_J"] for b in sharp["bracket_E_J"]])
    out["MA_sharp"] = sharp
    return out


def es_reference_source(ref, src, iface, shells_dir, tail_form="422", rings_information=False):
    """The thick reference's per-Type values for source `src` (the Type sums over `iface`); MA_sharp from the recorded shells."""
    post = os.path.join(ref, "postpro")
    row = row_of_source(read_csv(os.path.join(post, "surface-Q.csv")), src)
    p = {int(k[7:-1]): fnum(v) for k, v in row.items() if k.startswith("p_surf[")}
    drow = row_of_source(read_csv(os.path.join(post, "domain-E.csv")), src)
    e_elec = fnum([v for k, v in drow.items() if k.startswith("E_elec (J)")][0])
    crow = row_of_source(read_csv(os.path.join(post, "terminal-C.csv")), src)
    c_diag = fnum(crow[f"C[i][{src}] (F)"])
    cat_path = shells_catalog_path(ref, shells_dir)
    if not os.path.isfile(cat_path):
        raise MissingInput(cat_path)
    with open(cat_path) as handle:
        cat = json.load(handle)
    eo = edge_out_rows(read_csv(os.path.join(post, "surface-Q-edge.csv")), lambda r: int(fnum(r["i"])) == src)
    sharp = ma_sharp_from_shells(cat, p, e_elec, eo, tuple(iface["MA"]), tail_form, rings_information)
    sharp["reference_depth"], sharp["depth_information"] = reference_depth_information(cat, ref, tail_form)
    return {"source": src, "E_elec_J": e_elec, "SA_p": sum(p[k] for k in iface["SA"]), "MS_p": sum(p[k] for k in iface["MS"]),
            "MA_p": sum(p[k] for k in iface["MA"]), "C_F": c_diag,
            "C_offdiag_F": offdiag(crow, "C", src), "MA_sharp": sharp, "MA_sharp_p": sharp["MA_sharp_J"] / e_elec,
            "shells_catalog": cat_path}


def rows_at_f(rows, F, key="f (GHz)"):
    return [r for r in rows if abs(fnum(r[key]) - F) < FREQ_TOL_GHZ]


def port_v2(rows, F, port=1):
    v = rows_at_f(rows, F)
    if len(v) != 1:
        raise ValueError(f"port-V.csv at {F} GHz: {len(v)} rows (expected 1)")
    return fnum(v[0][f"Re{{V[{port}]}} (V)"]) ** 2 + fnum(v[0][f"Im{{V[{port}]}} (V)"]) ** 2


def mx_device_f(post, F, iface, catalog, sp_tails, p_denominator, with_port=True, excitation=None, tail_form="422"):
    """The Maxwell device's per-Type values at frequency F (the top-level postpro = the final solve). The per-interface
    rows of surface-Q-corrected.csv at F (and excitation `excitation` when the CSV carries an `exc` column and several
    excitations); p = E_surf corrected / denominator (E_norm corrected, or E_norm corrected - E_cap: A1.1.1)."""
    q = read_csv(os.path.join(post, "surface-Q-corrected.csv"))
    rows = rows_at_f(q, F)
    if excitation is not None and rows and "exc" in rows[0]:
        rows = [r for r in rows if int(fnum(r["exc"])) == excitation]
    by_iface = {}
    for r in rows:
        k = int(fnum(r["interface"]))
        if k in by_iface:
            raise ValueError(f"surface-Q-corrected.csv at {F} GHz: interface {k} twice (several excitations? use --excitation)")
        by_iface[k] = r
    needed = {k for ks in iface.values() for k in ks}
    if not needed <= set(by_iface):
        raise ValueError(f"surface-Q-corrected.csv at {F} GHz: interfaces {sorted(needed - set(by_iface))} absent")
    any_row = next(iter(by_iface.values()))
    e_norm_ft = fnum(any_row["E_norm corrected (J)"])
    e_norm_raw = fnum(any_row["E_norm raw (J)"])
    de = rows_at_f(read_csv(os.path.join(post, "domain-E.csv")), F)
    if len(de) != 1:
        raise ValueError(f"domain-E.csv at {F} GHz: {len(de)} rows (expected 1)")
    d = de[0]
    e_cap = fnum(d["E_cap (J)"]) if "E_cap (J)" in d else 0.0
    denom_ft = e_norm_ft - e_cap if p_denominator == "E_norm-minus-E_cap" else e_norm_ft
    denom_raw = e_norm_raw - e_cap if p_denominator == "E_norm-minus-E_cap" else e_norm_raw
    out = {"f_GHz": F, "E_norm_ft_J": e_norm_ft, "E_norm_raw_J": e_norm_raw, "E_elec_J": fnum(d["E_elec (J)"]),
           "E_mag_J": fnum(d["E_mag (J)"]), "E_cap_J": e_cap, "p_denominator": p_denominator, "p_denominator_ft_J": denom_ft}
    out["E_mag_over_E_elec"] = out["E_mag_J"] / out["E_elec_J"] if out["E_elec_J"] else None
    for T, ks in iface.items():
        out[f"{T}_E_ft_J"] = sum(fnum(by_iface[k]["E_surf corrected (J)"]) for k in ks)
        out[f"{T}_E_raw_J"] = sum(fnum(by_iface[k]["E_surf raw (J)"]) for k in ks)
        out[f"{T}_E_ff_J"] = sum(fnum(by_iface[k]["E_surf corrected fixed-flux (J)"]) for k in ks)
        if p_denominator == "E_norm":
            # the binary's own p_surf columns (the printed values: the record form of the stage-E Maxwell grades)
            out[f"{T}_p_ft"] = sum(fnum(by_iface[k]["p_surf corrected"]) for k in ks)
            out[f"{T}_p_raw"] = sum(fnum(by_iface[k]["p_surf raw"]) for k in ks)
            out[f"{T}_p_ff"] = sum(fnum(by_iface[k]["p_surf corrected fixed-flux"]) for k in ks)
        else:
            out[f"{T}_p_ft"] = out[f"{T}_E_ft_J"] / denom_ft
            out[f"{T}_p_raw"] = out[f"{T}_E_raw_J"] / denom_raw
            out[f"{T}_p_ff"] = out[f"{T}_E_ff_J"] / denom_ft
    edge = read_csv(os.path.join(post, "surface-Q-edge.csv"))
    for T, ks in iface.items():
        out[f"{T}_E_out_J"] = sum(fnum(r["E_out (J)"]) for r in rows_at_f(edge, F) if int(fnum(r["interface"])) in ks)
    if with_port:
        out["V2"] = port_v2(read_csv(os.path.join(post, "port-V.csv")), F)
        out["C_eff_ft_F"] = 2.0 * e_norm_ft / out["V2"]
        out["C_eff_raw_F"] = 2.0 * e_norm_raw / out["V2"]
    me = read_csv(os.path.join(post, "surface-response-model-energy.csv"))
    Fm = {}
    for r in rows_at_f(me, F):
        if int(fnum(r["evaluation"])) != 0:
            continue
        if excitation is not None and "source" in r and int(fnum(r["source"])) != excitation:
            continue
        m = int(fnum(r["model"]))
        Fm[m] = Fm.get(m, 0.0) + sum(fnum(r[f"fabricated surface energy[{k}] (J)"]) for k in iface["MA"])
    sharp = ma_sharp_convention_b(out["MA_E_ft_J"], Fm, catalog, sp_tails, tail_form)
    sharp["device_sharp_p"] = sharp["device_sharp_E_J"] / denom_ft
    sharp["bracket_p"] = sorted([b / denom_ft for b in sharp["bracket_E_J"]])
    out["MA_sharp"] = sharp
    return out


def mx_reference_f(ref, F, iface, shells_dir, tail_form="422", rings_information=False):
    """The thick driven reference at F: p_surf[k] of surface-Q.csv, E_elec / E_mag of domain-E.csv, C_eff = 2 E_elec / |V|^2,
    MA_sharp from the recorded shells at F."""
    post = os.path.join(ref, "postpro")
    row = rows_at_f(read_csv(os.path.join(post, "surface-Q.csv")), F)
    if len(row) != 1:
        raise ValueError(f"{ref}: surface-Q.csv at {F} GHz has {len(row)} rows (expected 1)")
    row = row[0]
    p = {int(k[7:-1]): fnum(v) for k, v in row.items() if k.startswith("p_surf[")}
    d = rows_at_f(read_csv(os.path.join(post, "domain-E.csv")), F)[0]
    e_elec, e_mag = fnum(d["E_elec (J)"]), fnum(d["E_mag (J)"])
    V2 = port_v2(read_csv(os.path.join(post, "port-V.csv")), F)
    cat_path = shells_catalog_path(ref, shells_dir)
    if not os.path.isfile(cat_path):
        raise MissingInput(cat_path)
    with open(cat_path) as handle:
        cat = json.load(handle)
    eo = edge_out_rows(read_csv(os.path.join(post, "surface-Q-edge.csv")), lambda r: abs(fnum(r["f (GHz)"]) - F) < FREQ_TOL_GHZ)
    sharp = ma_sharp_from_shells(cat, p, e_elec, eo, tuple(iface["MA"]), tail_form, rings_information)
    sharp["reference_depth"], sharp["depth_information"] = reference_depth_information(cat, ref, tail_form)
    return {"f_GHz": F, "run": ref, "E_elec_J": e_elec, "E_mag_J": e_mag, "E_mag_over_E_elec": e_mag / e_elec, "V2": V2,
            "C_eff_F": 2.0 * e_elec / V2, "SA_p": sum(p[k] for k in iface["SA"]), "MS_p": sum(p[k] for k in iface["MS"]),
            "MA_p": sum(p[k] for k in iface["MA"]), "MA_sharp": sharp, "MA_sharp_p": sharp["MA_sharp_J"] / e_elec,
            "shells_catalog": cat_path}


def reference_from_values(values, member, tail_form="422", rings_information=False):
    """A reference given as values (--reference-values JSON): {"members": {"<source or f>": {"SA_p", "MS_p", "MA_sharp_p",
    "C_F" | "C_eff_F", optional "MA_p", "E_elec_J", "C_offdiag_F", "MA_sharp_rings_um", "MA_sharp_reference_depth"}}} or the same keys
    at the top level for every member. Under tail_form fixed-shell the member's (or top-level) MA_sharp_rings_um = the ring boundaries
    its MA_sharp_p was read on must be FIXED_SHELL_RINGS_UM (absent = not of record): check_reference_rings; and its
    MA_sharp_reference_depth {"r_um", "order"} = the depth of the solve it was read from decides the depth label (absent = not
    recorded = INFORMATION_DEPTH_UNRECORDED; decision 1318): depth_information."""
    members = values.get("members")
    v = members[str(member)] if members else values
    out = dict(v)
    out.setdefault("C_offdiag_F", {})
    out["C_offdiag_F"] = {int(k): x for k, x in out["C_offdiag_F"].items()}
    out["from_values"] = True
    rings = v.get("MA_sharp_rings_um", values.get("MA_sharp_rings_um", []))
    depth = v.get("MA_sharp_reference_depth", values.get("MA_sharp_reference_depth"))
    out["MA_sharp"] = {"Rings_um": rings, "rings_information": check_reference_rings(rings, f"values reference member {member}", tail_form, rings_information),
                       "reference_depth": depth, "depth_information": depth_information(depth, tail_form)}
    return out


def static_limit_fit(freqs_ghz, c_eff):
    """C_eff(f) = C_s (1 + a f^2) fitted by unweighted least squares on f^2 (USER 1230 Q1) -> {C_s, a, residual_pct}.
    The FORM of record (REFREEZE-ROWS-P1.md, R-7 row; review d1230 MINOR-1): unweighted LSQ of C_eff on f^2 over the window's
    frequencies (3 points of record), 2 parameters (C_s, C_s a), per side; x = 100 (C_s,thin / C_s,thick - 1) at the 0.5 % bar."""
    n = len(freqs_ghz)
    if n < 2:
        raise ValueError("the static-limit fit needs >= 2 frequencies")
    xs = [f * f for f in freqs_ghz]
    mx = sum(xs) / n
    my = sum(c_eff) / n
    sxx = sum((x - mx) ** 2 for x in xs)
    if sxx == 0.0:
        raise ValueError("the static-limit fit needs distinct frequencies")
    slope = sum((x - mx) * (y - my) for x, y in zip(xs, c_eff)) / sxx
    c_s = my - slope * mx
    fitted = [c_s + slope * x for x in xs]
    resid = max(abs(pct(fv, y)) for fv, y in zip(fitted, c_eff))
    return {"C_s": c_s, "a_per_GHz2": (slope / c_s) if c_s else None, "fitted": fitted, "max_residual_pct": resid,
            "rule": "unweighted least squares of C_eff on f^2: C_eff(f) = C_s (1 + a f^2); C_s = the static limit"}


def band_check(x, band):
    out = {}
    for T, bt in (("SA", "SA"), ("MS", "MS"), ("MA_sharp", "MA"), ("C", "C")):
        if band is None or bt not in band or x.get(T) is None:
            out[T] = {"x_pt": x.get(T), "band": None, "in": None}
            continue
        lo, hi = band[bt]["low_pt"], band[bt]["high_pt"]
        out[T] = {"x_pt": x[T], "band": [lo, hi], "in": bool(lo - BAND_TOL_PT <= x[T] <= hi + BAND_TOL_PT)}
    return out


def grade_member(dev, ref, band, c_dev_key, c_ref_key):
    """One source / frequency: x per Type, the letters (verdict), the information rows, the off-diagonal C (910)."""
    x = {"SA": pct(dev["SA_p_ft"], ref["SA_p"]), "MS": pct(dev["MS_p_ft"], ref["MS_p"]),
         "MA_sharp": pct(dev["MA_sharp"]["device_sharp_p"], ref["MA_sharp_p"])}
    c_ref = ref.get(c_ref_key)
    x["C"] = pct(dev[c_dev_key], c_ref) if (c_dev_key in dev and c_ref) else None
    br = sorted([pct(dev["MA_sharp"]["bracket_p"][0], ref["MA_sharp_p"]), pct(dev["MA_sharp"]["bracket_p"][1], ref["MA_sharp_p"])])
    letters = {"SA": verdict(x["SA"], D4["SA"]), "MS": verdict(x["MS"], D4["MS"]), "MA_sharp": verdict(x["MA_sharp"], D4["MA_sharp"], br),
               "C": verdict(x["C"], D4["C"]) if x["C"] is not None else "n/a"}
    info = {"SA_raw_pct": pct(dev["SA_p_raw"], ref["SA_p"]), "MS_raw_pct": pct(dev["MS_p_raw"], ref["MS_p"])}
    if ref.get("MA_p"):
        info["MA_raw_pct"] = pct(dev["MA_p_raw"], ref["MA_p"])
        info["MA_raw_ft_pct"] = pct(dev["MA_p_ft"], ref["MA_p"])
        info["MA_raw_ft_letter_information"] = verdict(info["MA_raw_ft_pct"], 5.0)
    if "C_raw_F" in dev and ref.get("C_F"):
        info["C_raw_pct_information"] = pct(dev["C_raw_F"], ref["C_F"])
    if "C_eff_raw_F" in dev and ref.get("C_eff_F"):
        info["C_eff_raw_pct_information"] = pct(dev["C_eff_raw_F"], ref["C_eff_F"])
    if "E_elec_ft_J" in dev and ref.get("E_elec_J"):
        info["E_elec_ft_vs_ref_pct"] = pct(dev["E_elec_ft_J"], ref["E_elec_J"])
        info["E_elec_raw_vs_ref_pct"] = pct(dev["E_elec_raw_J"], ref["E_elec_J"])
    if "E_norm_ft_J" in dev and ref.get("E_elec_J"):
        info["E_norm_ft_vs_ref_E_elec_pct"] = pct(dev["E_norm_ft_J"], ref["E_elec_J"])
        info["E_norm_raw_vs_ref_E_elec_pct"] = pct(dev["E_norm_raw_J"], ref["E_elec_J"])
    if "V2" in dev and ref.get("V2"):
        info["V2_thin_vs_thick_pct"] = pct(dev["V2"], ref["V2"])
    if "SA_p_ff" in dev:
        info["SA_ff_pct"] = pct(dev["SA_p_ff"], ref["SA_p"])
        info["MS_ff_pct"] = pct(dev["MS_p_ff"], ref["MS_p"])
    if "E_mag_over_E_elec" in dev:
        info["E_mag_over_E_elec_thin"] = dev["E_mag_over_E_elec"]
        info["E_mag_over_E_elec_thick"] = ref.get("E_mag_over_E_elec")
    off = {}
    for j, cr in (ref.get("C_offdiag_F") or {}).items():
        cd = (dev.get("C_ft_offdiag_F") or {}).get(j)
        if cd is None or not ref.get("C_F"):
            continue
        off[j] = {"ref_F": cr, "dev_ft_F": cd, "rel_pct": (pct(cd, cr) if cr else None), "abs_over_Cii_pct": 100.0 * abs(cd - cr) / ref["C_F"]}
    return {"x_pct": x, "MA_sharp_bracket_pct": br, "letters": letters, "information": info, "C_offdiag": off, "band": band_check(x, band)}


def aggregate_cell(members, grades, T, rule_text):
    """The fail-closed cell letter over the members (1210 (a) / 1213 (a)): F if any F, else I if any I, else P; the
    worst-|x| member quoted with the range and the counts; in band only if every member is."""
    xs = {m: grades[m]["x_pct"][T] for m in members if grades[m]["x_pct"].get(T) is not None}
    if not xs:
        return {"letter": "n/a", "x_pct": None, "rule": rule_text, "n_members": 0}
    letters = {m: grades[m]["letters"][T] for m in xs}
    ls = list(letters.values())
    letter = "FAIL" if "FAIL" in ls else ("INCONCLUSIVE (ft inside, bracket outside)" if any(l.startswith("INCONCLUSIVE") for l in ls) else "PASS")
    worst = max(xs, key=lambda m: abs(xs[m]))
    bands = [grades[m]["band"][T]["in"] for m in xs]
    return {"letter": letter, "x_pct": xs[worst], "worst_member": worst, "x_range_pct": [min(xs.values()), max(xs.values())],
            "letters_by_member": {str(m): letters[m] for m in xs}, "n_members": len(xs), "n_pass": sum(l == "PASS" for l in ls),
            "n_fail": sum(l == "FAIL" for l in ls), "in_band": (all(bands) if None not in bands else None),
            "n_in_band": sum(bool(b) for b in bands), "rule": rule_text}


def static_limit_cell(freqs, per, band):
    """The Q1 C cell of a ported Maxwell line: C_s thin vs C_s thick from the static-limit fits of both sides."""
    thin = static_limit_fit(freqs, [per[F]["device"]["C_eff_ft_F"] for F in freqs])
    thick = static_limit_fit(freqs, [per[F]["reference"]["C_eff_F"] for F in freqs])
    x = pct(thin["C_s"], thick["C_s"])
    inb = band_check({"C": x}, band)["C"]
    return {"letter": verdict(x, D4["C"]), "x_pct": x, "thin_fit": thin, "thick_fit": thick, "in_band": inb["in"], "band": inb["band"],
            "C_eff_by_f_pct": {str(F): per[F]["grade"]["x_pct"]["C"] for F in freqs},
            "rule": "USER 1230 Q1: the fitted STATIC-LIMIT C over the frequencies is the graded cell on a ported line (C_eff(f) beside)"}


def information_cell(cells, per):
    """The fixed-shell form's INFORMATION labels on the MA_sharp cell (not counted; the member letters stay as read): a member's
    reference MA shells not of record under --reference-rings-information (decision 1310 Q-D), or a member's reference coarser than
    r5 p5 (decision 1314); both labels prefixed when both apply."""
    labels = []
    for key in ("rings_information", "depth_information"):
        found = sorted({m["reference"].get("MA_sharp", {}).get(key) or "" for m in per.values()} - {""})
        labels += found
    if labels:
        cells["MA_sharp"] = dict(cells["MA_sharp"], counted=False, letter=f"{'; '.join(labels)}: {cells['MA_sharp']['letter']}",
                                 reference_rings_um=sorted({tuple(m["reference"]["MA_sharp"].get("Rings_um", [])) for m in per.values()}),
                                 reference_depth=sorted({json.dumps(m["reference"]["MA_sharp"].get("reference_depth")) for m in per.values()}))


def grade_es(a, catalog, cat_src, sp_tails, band, values):
    final, ncyc = final_cycle_dir(a.device)
    iface, _ = interfaces_by_type(os.path.join(a.device, "postpro", "config_resolved.json"))
    q = read_csv(os.path.join(final, "surface-Q-corrected.csv"))
    sources = sorted(int(fnum(r["i"])) for r in q)
    if a.sources:
        sources = [s for s in sources if s in a.sources]
    if values is None:
        qr = read_csv(os.path.join(a.reference, "postpro", "surface-Q.csv"))
        ref_sources = sorted(int(fnum(r["i"])) for r in qr)
        if sorted(int(fnum(r["i"])) for r in q) != ref_sources:
            raise ValueError(f"source sets differ: device {sources} vs reference {ref_sources}")
    per = {}
    for s in sources:
        d = es_device_source(final, s, iface, catalog, sp_tails, a.tail_form)
        r = (reference_from_values(values, s, a.tail_form, a.reference_rings_information) if values is not None
             else es_reference_source(a.reference, s, iface, a.shells_dir, a.tail_form, a.reference_rings_information))
        per[s] = {"device": d, "reference": r, "grade": grade_member(d, r, band, "C_ft_F", "C_F")}
    rule = "every source must PASS (fail-closed; decision 1210 (a)); the worst-|x| source quoted" if len(sources) > 1 else "single source"
    cells = {T: aggregate_cell(sources, {s: per[s]["grade"] for s in sources}, T, rule) for T in TYPES}
    for T in TYPES:
        cells[T]["counted"] = True
    information_cell(cells, per)
    return {"interfaces": iface, "device_final_dir": final, "device_cycles": ncyc, "members": sources, "member_kind": "source",
            "per_member": {str(s): v for s, v in per.items()}, "cells": cells}


def parse_reference_at(items):
    out = {}
    for it in items or []:
        f, _, run = it.partition("=")
        out[float(f)] = run
    return out


def grade_mx(a, catalog, cat_src, sp_tails, band, values):
    post = os.path.join(a.device, "postpro")
    iface, _ = interfaces_by_type(os.path.join(post, "config_resolved.json"))
    freqs = list(a.frequencies)
    refs_at = parse_reference_at(a.reference_at)
    # the port voltage is read whenever a C_eff reference exists (a reference run always; a values file carrying C_eff_F)
    if values is None:
        with_port = True
    else:
        members = values.get("members") or {"": values}
        with_port = any("C_eff_F" in m for m in members.values())
    per = {}
    for F in freqs:
        d = mx_device_f(post, F, iface, catalog, sp_tails, a.p_denominator, with_port=with_port, excitation=a.excitation, tail_form=a.tail_form)
        if values is not None:
            r = reference_from_values(values, F, a.tail_form, a.reference_rings_information)
        else:
            ref_run = refs_at.get(F, a.reference)
            if ref_run is None:
                raise MissingInput(f"--reference or --reference-at {F}=<run> for {F} GHz")
            r = mx_reference_f(ref_run, F, iface, a.shells_dir, a.tail_form, a.reference_rings_information)
        per[F] = {"device": d, "reference": r, "grade": grade_member(d, r, band, "C_eff_ft_F", "C_eff_F")}
    rule = "every frequency must PASS (fail-closed; decision 1213 (a)); the worst-|x| frequency quoted" if len(freqs) > 1 else "single frequency"
    cells = {T: aggregate_cell(freqs, {F: per[F]["grade"] for F in freqs}, T, rule) for T in ("SA", "MS", "MA_sharp")}
    for T in ("SA", "MS", "MA_sharp"):
        cells[T]["counted"] = True
    information_cell(cells, per)
    direct = aggregate_cell(freqs, {F: per[F]["grade"] for F in freqs}, "C", rule + "; C_eff(f) = 2 E_norm_ft / |V|^2 (decision 1214 DIRECT form)")
    if a.c_form == "direct":
        cells["C"] = dict(direct, counted=True, c_form="direct")
    elif a.c_form == "static-limit":
        if direct["n_members"] < 2:
            raise ValueError("--c-form static-limit needs >= 2 frequencies with a port on both sides")
        cells["C"] = dict(static_limit_cell(freqs, per, band), counted=True, c_form="static-limit", direct_form_beside=direct)
    else:
        cells["C"] = dict(direct, counted=False, c_form="information", letter=f"INFORMATION (F-40): {direct['letter']}")
    return {"interfaces": iface, "device_final_dir": post, "members": freqs, "member_kind": "f_GHz", "p_denominator": a.p_denominator,
            "per_member": {str(F): v for F, v in per.items()}, "cells": cells}


def summary_line(label, out):
    parts = []
    for T, c in out["cells"].items():
        if c.get("x_pct") is None:
            parts.append(f"{T} n/a")
            continue
        rng = f" [{c['x_range_pct'][0]:+.3f}, {c['x_range_pct'][1]:+.3f}]" if c.get("x_range_pct") and c.get("n_members", 1) > 1 else ""
        band = "in" if c.get("in_band") else ("OUT" if c.get("in_band") is False else "-")
        parts.append(f"{T} {c['x_pct']:+.3f} %{rng} {c['letter']} band {band}")
    return f"{label}: {out['member_kind']} {len(out['members'])}; " + "; ".join(parts)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--form", required=True, choices=("es", "maxwell"))
    parser.add_argument("--device", required=True, help="the thin CORRECTED run directory (postpro/ inside)")
    parser.add_argument("--reference", help="the thick reference run directory (every source / frequency)")
    parser.add_argument("--reference-at", action="append", metavar="F=RUN", help="Maxwell: the thick reference run at frequency F (repeatable)")
    parser.add_argument("--reference-values", help="a reference given as values (JSON; the Richardson references' form)")
    parser.add_argument("--shells-dir", help="directory of the reference runs' <run>-shells-catalog.json (the MA shells of record)")
    parser.add_argument("--library", help="the process library (the spatial clusters' own MA tails)")
    parser.add_argument("--model-catalog", help="the run's model catalog JSON when postpro/palace.json was not kept")
    parser.add_argument("--band", help="the filed coverage-report.json (block_C = the predicted band per Type)")
    parser.add_argument("--frequencies", nargs="+", type=float, default=[5.2, 5.7, 6.2], help="Maxwell: the graded frequencies (GHz)")
    parser.add_argument("--excitation", type=int, help="Maxwell: the excitation whose rows are graded when the CSVs carry several")
    parser.add_argument("--sources", nargs="+", type=int, help="ES: grade these sources only (default: every source)")
    parser.add_argument("--p-denominator", choices=P_DENOMINATORS, default="E_norm",
                        help="Maxwell: p = E_surf / E_norm (the binary's p_surf) or / (E_norm - E_cap) (the junction-driven form of the cards, A1.1.1)")
    parser.add_argument("--c-form", choices=C_FORMS, help="Maxwell (required): the C cell's form -- static-limit (USER 1230 Q1), direct (1214), information (F-40)")
    parser.add_argument("--tail-form", choices=TAIL_FORMS, default="422",
                        help="the MA_sharp tail form: 422 (the class tails of record; the default until the P4 re-freeze) or fixed-shell (FIXED_SHELL_TAILS on the device, "
                             "the reference's rings required = FIXED_SHELL_RINGS_UM; USER 1308 / decision 1310)")
    parser.add_argument("--reference-rings-information", action="store_true",
                        help="fixed-shell: grade a reference recorded with another ring set and label its MA_sharp cell INFORMATION instead of refusing (1310 Q-D)")
    parser.add_argument("--label", help="the window's name in the output")
    parser.add_argument("--out", required=True, help="the grade JSON")
    a = parser.parse_args(argv)
    if a.form == "maxwell" and a.c_form is None:
        parser.error("--c-form is required for --form maxwell (static-limit | direct | information)")
    if a.reference is None and a.reference_values is None and not a.reference_at:
        parser.error("one of --reference / --reference-at / --reference-values is required")
    if a.reference_values is None and a.shells_dir is None:
        parser.error("--shells-dir is required with a reference run (the MA shells catalog of record)")
    try:
        values = None
        if a.reference_values:
            if not os.path.isfile(a.reference_values):
                raise MissingInput(a.reference_values)
            with open(a.reference_values) as handle:
                values = json.load(handle)
        post = os.path.join(a.device, "postpro")
        catalog, cat_src = model_catalog(post, a.model_catalog)
        tails = spatial_tails(a.library)
        band = None
        if a.band:
            with open(a.band) as handle:
                band = json.load(handle)["block_C"]
        out = grade_es(a, catalog, cat_src, tails, band, values) if a.form == "es" else grade_mx(a, catalog, cat_src, tails, band, values)
    except MissingInput as e:
        print(f"grade_cells.py: missing input {e}", file=sys.stderr)
        return 2
    except ReferenceRingsNotOfRecord as e:
        print(f"grade_cells.py: {e}", file=sys.stderr)
        return 3
    label = a.label or os.path.basename(os.path.normpath(a.device))
    out.update({"tool": "grade_cells.py", "window": label, "form": a.form, "D4": D4, "tail_form": a.tail_form,
                "fixed_shell": ({"rings_um": list(FIXED_SHELL_RINGS_UM), "tails": {k: list(v) for k, v in FIXED_SHELL_TAILS.items()},
                                 "reference_rings_information": a.reference_rings_information} if a.tail_form == "fixed-shell" else None),
                "scorer": "the 994 / 995 tabulators' letters (verdict; MA_sharp (B) tails by --tail-form; the shells estimator) per source / frequency, fail-closed (1210 / 1213)",
                "device_run": a.device, "reference": a.reference or a.reference_at or a.reference_values, "shells_dir": a.shells_dir,
                "library": a.library, "model_catalog_source": cat_src, "models": {str(k): v for k, v in catalog.items()},
                "band_file": a.band, "band": ({T: [band[T]["low_pt"], band[T]["high_pt"]] for T in ("SA", "MS", "MA", "C") if T in band} if band else None),
                "c_form": a.c_form if a.form == "maxwell" else "diagonal C fixed-trace[i][i] (ES)"})
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    with open(a.out, "w") as handle:
        json.dump(out, handle, indent=1, default=str)
    print(summary_line(label, out))
    return 0


if __name__ == "__main__":
    sys.exit(main())
