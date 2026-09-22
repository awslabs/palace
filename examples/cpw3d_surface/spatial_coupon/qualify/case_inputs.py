# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""The Palace run inputs of a coupon from the case's own frozen sources (supervisor
decision 52): a case registered by `coupon-library build` is runnable without a
reference campaign.

Everything comes from the manifest case's source directory and the identity mesh:

* process-library.json - the fabrication (substrate permittivity, interface layer
  thickness / permittivity), the matching radius, the model's edges (slot, conductor,
  process normal) and its interfaces in coupon-index order (Interfaces[].Coupon is the
  interface index -> type map the response matrices are labeled with);
* basis-contract.json + trace-vertices.csv + trace-triangles.csv - the trace basis in
  the process-library frame, placed in the mesh frame by the bound process frame
  (trace_basis.load_trace_basis); every source trace is REGENERATED from it exactly as
  the producer wrote it (generate_spatial_response.write_surface_trace: the hat of
  basis vertex k, the conductor lift of every conductor but the first, the zero trace),
  so a case carries no per-source files; the contract's ZeroTraceIndices and
  ConductorStates bind the source count;
* the identity mesh's $PhysicalNames - the boundary attributes the config may name
  (generate_spatial_response.make_config filters the producer's candidate attributes by
  them, as the producer did on the reference mesh);
* the production recipe's PhysicsRun block (Order, LinearTol) - recorded in the
  manifest with its calibration provenance; every case runs at these values.

The derived config is the producer's (make_config) with Model.Mesh, Problem.Output and
the trace directory substituted; the qualify command compares it with the reference's
own config when a reference exists (reference_campaign) and fails closed on any
difference outside those paths and the recipe's Order / Tol.
"""
import csv
import hashlib
import json
from pathlib import Path
import sys
import tomllib

import numpy as np

HERE = Path(__file__).resolve().parent
TOOLS = HERE.parent
for path in (str(HERE), str(TOOLS)):
    if path not in sys.path:
        sys.path.insert(0, path)
import generate_spatial_response as producer  # noqa: E402
import trace_basis  # noqa: E402
from mixed_mesh import SHELL_LABEL_STRIDE  # noqa: E402

PHYSICS_RUN_KEY = "PhysicsRun"
JOB_POLICY_KEY = "JobPolicy"
JOB_POLICY_MODES = ("speed", "frugal", "fixed")
REDUCER_BLOCK_SIZE_KEY = "ReducerBlockSize"
TRACES_DIRECTORY = "traces"
ZERO_TRACE_FILE = "zero-trace.csv"
# Paths a derived config carries that differ from a reference config by construction.
PATH_FIELDS = ("Model.Mesh", "Problem.Output", "Boundaries.PrescribedPotential[].DataFile")


def sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


class CaseInputError(ValueError):
    """A fail-closed derivation (recorded per case, never a crash of the run)."""


def physics_run_parameters(manifest, case=None):
    """The recipe's Order and Linear.Tol (ProductionRecipe.PhysicsRun); a labeled
    calibration manifest takes them from its Calibration.ProductionManifest.  A
    calibration case declaring Calibration.PhysicsRun.LinearTol runs at that tolerance
    when the manifest's Calibration.GateDeviations.LinearTol labels it (Production equal
    to the recipe value, the case named, ProductionUse FORBIDDEN; fail closed otherwise)."""
    recipe = manifest.get("ProductionRecipe")
    if recipe is None and "Calibration" in manifest:
        production = Path(manifest["Path"]).parent / manifest["Calibration"]["ProductionManifest"]
        recipe = json.loads(production.read_text()).get("ProductionRecipe")
    block = (recipe or {}).get(PHYSICS_RUN_KEY)
    if (not isinstance(block, dict) or not isinstance(block.get("Order"), int) or block["Order"] < 1 or
            not isinstance(block.get("LinearTol"), float) or not 0.0 < block["LinearTol"] < 1.0 or
            not isinstance(block.get("Provenance"), str) or not block["Provenance"]):
        raise CaseInputError("the production recipe must bind PhysicsRun {Order (int >= 1), LinearTol (0 < float < 1), "
                             "Provenance}: no run parameters for the coupon")
    parameters = {"Order": int(block["Order"]), "LinearTol": float(block["LinearTol"]), "Provenance": block["Provenance"]}
    job_policy = block.get(JOB_POLICY_KEY)
    if job_policy is not None:
        # The manifest-recorded default of the per-coupon source split (decision 61b);
        # the command line (--job-policy / --fixed-jobs) overrides it.
        mode = job_policy.get("Mode") if isinstance(job_policy, dict) else None
        fixed = job_policy.get("FixedJobs") if isinstance(job_policy, dict) else None
        if (mode not in JOB_POLICY_MODES or (mode == "fixed") != (isinstance(fixed, int) and not isinstance(fixed, bool) and fixed >= 1)
                or not isinstance(job_policy.get("Rule"), str) or not job_policy["Rule"]):
            raise CaseInputError(f"PhysicsRun.JobPolicy must be {{Mode in {JOB_POLICY_MODES}, FixedJobs (int >= 1, fixed only), "
                                 f"Rule}}, not {job_policy!r}")
        parameters["JobPolicy"] = {"Mode": mode, "FixedJobs": fixed if mode == "fixed" else None, "Rule": job_policy["Rule"]}
    reducer_block_size = block.get(REDUCER_BLOCK_SIZE_KEY)
    if reducer_block_size is not None:
        # The manifest-recorded default of the reducer's PALACE_RESPONSE_BLOCK_SIZE
        # (decision 62(1)); --reducer-block-size on the command line overrides it.
        if (not isinstance(reducer_block_size, dict) or not isinstance(reducer_block_size.get("Value"), int) or
                isinstance(reducer_block_size.get("Value"), bool) or reducer_block_size["Value"] < 1 or
                not isinstance(reducer_block_size.get("Rule"), str) or not reducer_block_size["Rule"]):
            raise CaseInputError(f"PhysicsRun.ReducerBlockSize must be {{Value (int >= 1), Rule}}, not {reducer_block_size!r}")
        parameters["ReducerBlockSize"] = reducer_block_size["Value"]
    deviation_block = ((case or {}).get("Calibration") or {}).get(PHYSICS_RUN_KEY)
    if deviation_block is not None:
        deviation = ((manifest.get("Calibration") or {}).get("GateDeviations") or {}).get("LinearTol")
        tolerance = deviation_block.get("LinearTol") if isinstance(deviation_block, dict) else None
        if (not isinstance(deviation_block, dict) or set(deviation_block) != {"LinearTol"} or not isinstance(tolerance, float) or
                not 0.0 < tolerance < 1.0 or not isinstance(deviation, dict) or deviation.get("Production") != parameters["LinearTol"] or
                deviation.get("Calibration") != tolerance or case.get("Id") not in (deviation.get("Cases") or []) or
                "FORBIDDEN" not in str(deviation.get("ProductionUse", ""))):
            raise CaseInputError(f"{case.get('Id')} declares Calibration.PhysicsRun.LinearTol without a labeled calibration-only "
                                 f"deviation of the recipe value {parameters['LinearTol']} (Calibration.GateDeviations.LinearTol)")
        parameters["LinearTol"] = tolerance
        parameters["Deviation"] = {"LinearTol": {"Production": deviation["Production"], "Calibration": tolerance,
                                                 "Reason": deviation.get("Reason"), "ProductionUse": deviation["ProductionUse"]}}
    return parameters


def interface_types(config):
    """Interface index -> type of a Palace config's postprocessed dielectric interfaces
    (the labels of the surface response matrix rows)."""
    entries = config.get("Boundaries", {}).get("Postprocessing", {}).get("Dielectric", [])
    types = {}
    for entry in entries:
        if "Index" not in entry or "Type" not in entry:
            raise CaseInputError(f"dielectric interface entry without Index / Type: {entry}")
        if int(entry["Index"]) in types:
            raise CaseInputError(f"repeated dielectric interface index {entry['Index']}")
        types[int(entry["Index"])] = str(entry["Type"])
    return types


def read_signature_rows(path):
    with Path(path).open(newline="") as stream:
        return list(csv.DictReader(stream))


def load_case_sources(directory, files):
    """The frozen source files of a case (role -> path) with their digests verified
    against the manifest binding."""
    directory = Path(directory)
    paths = {}
    for role, entry in files.items():
        if "Name" not in entry:
            continue
        path = directory / entry["Name"]
        if not path.is_file():
            raise CaseInputError(f"the bound source {entry['Name']} ({role}) is missing under {directory}")
        actual = sha256(path)
        if actual != entry["SHA256"]:
            raise CaseInputError(f"{entry['Name']} ({role}) SHA256 {actual} differs from the manifest binding {entry['SHA256']}")
        paths[role] = path
    for role in ("Process", "ProcessLibrary", "BasisContract", "TraceVertices", "TraceTriangles", "Boundary", "Signature"):
        if role not in paths:
            raise CaseInputError(f"the case binds no {role}: no run inputs can be derived")
    return paths


def regenerate_traces(basis, labels, out_dir):
    """basis-NNNN.csv for every basis vertex, conductor-N.csv (the lift of conductor N,
    every conductor but the first) and zero-trace.csv under out_dir, in the mesh frame,
    in the producer's format; returns (basis traces, conductor traces, zero trace)."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    for stale in out_dir.glob("*.csv"):
        stale.unlink()
    points, triangles, index_of = basis["Points"], basis["Triangles"], basis["Basis"]
    count = int(index_of.max())
    if sorted(index_of[index_of > 0].tolist()) != list(range(1, count + 1)):
        raise CaseInputError("the trace basis indices are not 1..N")
    traces = []
    for k in range(1, count + 1):
        values = np.zeros(len(points))
        values[index_of == k] = 1.0
        path = out_dir / f"basis-{k:04d}.csv"
        producer.write_surface_trace(path, points, triangles, values)
        traces.append(path)
    conductors = {}
    for conductor in sorted(set(int(label) for label in labels) - {0, 1}):
        values = (labels == conductor).astype(float)
        path = out_dir / f"conductor-{conductor}.csv"
        producer.write_surface_trace(path, points, triangles, values)
        conductors[conductor] = path
    zero = out_dir / ZERO_TRACE_FILE
    producer.write_surface_trace(zero, points, triangles, np.zeros(len(points)))
    return traces, conductors, zero


def local_edges(model, frame):
    """The model's edges with ProcessNormal in the mesh frame (frame @ canonical, -0.0
    normalized to 0.0), the slot / conductor unchanged: what make_config reads."""
    edges = []
    for edge in model["Edges"]:
        normal = frame @ np.asarray(edge["ProcessNormal"], dtype=float)
        edges.append({**edge, "ProcessNormal": [float(value) + 0.0 for value in normal]})
    return edges


def derive(case, directory, *, mesh_path, physics_run, out_dir, mesh=None, output_root="<output>", traces_remote=None,
           radial_shells=None):
    """The run inputs of a case: (config, record).  `mesh` is the value written into
    Model.Mesh (default: the local mesh path), `traces_remote` the directory written
    into every DataFile (default: the local trace directory); Problem.Output is
    `output_root`.  The traces are regenerated under out_dir/traces.  `radial_shells`
    (the census of a radial-shell relabel) derives the base config against the parent
    MA labels and expands it (expand_radial_shells); the record then carries the base
    config under BaseConfig (the reference comparison view) and the shell map under
    RadialShells."""
    files = case["Source"]["Files"]
    paths = load_case_sources(directory, files)
    process = tomllib.loads(paths["Process"].read_text())
    library = json.loads(paths["ProcessLibrary"].read_text())
    contract = json.loads(paths["BasisContract"].read_text())
    if len(library.get("Models", [])) != 1:
        raise CaseInputError(f"the case's process library binds {len(library.get('Models', []))} models, not one")
    model = library["Models"][0]
    if model["Name"] != contract["Model"]:
        raise CaseInputError(f"the basis contract binds model {contract['Model']}, the library {model['Name']}")
    radius = float(library["MatchingRadius"])
    if abs(radius - float(process["Radius"])) > 1e-12 * max(1.0, radius):
        raise CaseInputError(f"process.toml Radius {process['Radius']} differs from the library MatchingRadius {radius}")
    fabrication = library["Fabrication"]
    for key, toml_key in (("MetalThickness", "MetalThickness"), ("OveretchDepth", "Overetch"),
                          ("SidewallAngleDegrees", "SidewallAngle")):
        if abs(float(fabrication[key]) - float(process[toml_key])) > 1e-12:
            raise CaseInputError(f"process.toml {toml_key} {process[toml_key]} differs from the library Fabrication {key} {fabrication[key]}")
    basis = trace_basis.load_trace_basis(paths["BasisContract"], paths["TraceVertices"], paths["TraceTriangles"],
                                         paths["ProcessLibrary"])
    labels = np.array([int(row["conductor"]) for row in read_signature_rows(paths["TraceVertices"])])
    edges = local_edges(model, basis["Frame"])
    conductors = sorted({int(edge["Conductor"]) for edge in edges})
    if conductors != list(range(1, len(conductors) + 1)):
        raise CaseInputError(f"the model's conductors are not 1..N: {conductors}")
    terminal_conductors = conductors[1:]
    if int(contract.get("ConductorStates", len(terminal_conductors))) != len(terminal_conductors):
        raise CaseInputError(f"the contract binds {contract.get('ConductorStates')} conductor states, the model "
                             f"{len(terminal_conductors)} conductors beyond the first")
    if set(labels.tolist()) - {0} - set(conductors):
        raise CaseInputError(f"trace-vertices conductor labels {sorted(set(labels.tolist()))} outside the model's conductors")
    trace_dir = Path(out_dir) / TRACES_DIRECTORY
    traces, conductor_traces, zero_trace = regenerate_traces(basis, labels, trace_dir)
    interfaces = producer.model_interfaces(model)
    declared = {int(entry["Coupon"]): str(entry["Type"]) for entry in model["Interfaces"] if "Coupon" in entry}
    enumerated = {index: entry["Type"] for index, entry in enumerate(interfaces, start=1)}
    if declared and declared != enumerated:
        raise CaseInputError(f"the library's Interfaces[].Coupon indices {declared} differ from the producer's "
                             f"enumeration {enumerated}")
    layers = {name: (float(layer["Thickness"]), float(layer["Permittivity"]))
              for name, layer in fabrication["InterfaceLayers"].items()}
    available = producer.mesh_boundary_attributes(Path(mesh_path))
    shell_parents = set(shell_labels_by_parent(radial_shells)) if radial_shells is not None else set()
    config_available = available
    if radial_shells is None and any(int(a) >= SHELL_LABEL_STRIDE for a in available):
        raise CaseInputError(f"the mesh carries radial-shell labels {sorted(a for a in available if int(a) >= SHELL_LABEL_STRIDE)} "
                             f"but the build record binds no radial-shell census (RadialShells / Relabel)")
    if shell_parents:
        shell_labels = {int(shell["Label"]) for shell in radial_shells["Shells"]}
        if shell_parents & available or not shell_labels <= available:
            raise CaseInputError(f"the mesh's boundary attributes do not match the radial-shell census (parents "
                                 f"{sorted(shell_parents)} must be absent, shell labels present)")
        # The base config sees the parent labels the producer knows; the expansion maps them.
        config_available = (available - shell_labels) | shell_parents
    traces_remote = traces_remote or str(trace_dir)
    config = producer.make_config(
        Path(output_root), "run", mesh if mesh is not None else str(mesh_path),
        [f"{traces_remote}/{path.name}" for path in traces], f"{traces_remote}/{zero_trace.name}",
        terminal_conductors, True, physics_run["Order"], radius, float(fabrication["SubstratePermittivity"]),
        layers, interfaces, edges, available_attributes=config_available,
        terminal_traces={conductor: f"{traces_remote}/{path.name}" for conductor, path in conductor_traces.items()})
    config["Problem"]["Output"] = output_root
    config["Solver"]["Linear"]["Tol"] = physics_run["LinearTol"]
    base_config, shell_map = None, None
    if shell_parents:
        base_config = config
        config, shell_map = expand_radial_shells(base_config, radial_shells)
    attribute_check = check_attributes(config, available)
    sources = []
    local_paths = {path.name: path for path in list(traces) + list(conductor_traces.values())}
    for entry in config["Boundaries"]["PrescribedPotential"]:
        name = Path(entry["DataFile"]).name
        path = local_paths[name]
        sources.append({"Index": int(entry["Index"]), "Name": name, "Path": str(path), "SHA256": sha256(path),
                        "Terminal": bool(entry.get("TerminalAttributes"))})
    zero_indices = [int(i) for i in contract.get("ZeroTraceIndices", [])]
    record = {"Origin": "case", "Directory": str(directory), "Model": model["Name"],
              "SourceSHA256": {role: files[role]["SHA256"] for role in files},
              "PhysicsRun": dict(physics_run), "Interfaces": interface_types(config),
              "InterfaceTypes": sorted(set(interface_types(config).values())),
              "Sources": sources, "ZeroTraceIndices": zero_indices,
              "Traces": {"Directory": str(trace_dir), "Rule": "regenerated from the bound trace basis in the mesh frame "
                         "(trace_basis.load_trace_basis; producer format); coordinates within the contract's "
                         "FrameFitResidual of the producer's own trace files",
                         "FrameFitResidual": contract.get("FrameFitResidual"),
                         "ZeroTrace": str(zero_trace), "ZeroTraceSHA256": sha256(zero_trace),
                         "ContractOutputSourceSHA256Identical": contract_digests_identical(contract, sources)},
              "PlanViewBoundary": str(paths["Boundary"]),
              "RetainedEtch": str(paths["RetainedEtch"]) if "RetainedEtch" in paths else None,
              "Signature": str(paths["Signature"]), "MeshAttributes": sorted(available), "AttributeCheck": attribute_check}
    if base_config is not None:
        record["BaseConfig"] = base_config
        record["BaseInterfaces"] = interface_types(base_config)
        record["RadialShells"] = {"Rule": "the base config (derived against the parent MA labels; the reference comparison "
                                          "view) expanded by expand_radial_shells: one MA Dielectric entry per shell ordinal, "
                                          "the parent labels replaced by the shell labels in every attribute list",
                                  "Census": radial_shells.get("Mesh"), "RingRadii": radial_shells.get("RingRadii"),
                                  "Interfaces": {str(index): value for index, value in sorted(shell_map.items())}}
    return config, record


def contract_digests_identical(contract, sources):
    """Whether the regenerated traces equal the producer's recorded digests
    (OutputSourceSHA256) file by file; None when the contract records none."""
    recorded = contract.get("OutputSourceSHA256")
    if not recorded:
        return None
    by_name = {Path(key).name: value for key, value in recorded.items()}
    compared = [source for source in sources if source["Name"] in by_name]
    if not compared:
        return None
    return all(by_name[source["Name"]] == source["SHA256"] for source in compared)


def shell_labels_by_parent(shells):
    """Parent MA attribute -> [(ordinal, label)] of a radial-shell census
    (relabel_radial_ma_shells.py identity.msh.radial-shells.json), ordinals ascending."""
    by_parent = {}
    for shell in shells["Shells"]:
        by_parent.setdefault(int(shell["Parent"]), []).append((int(shell["Ordinal"]), int(shell["Label"])))
    return {parent: sorted(items) for parent, items in by_parent.items()}


def expand_radial_shells(config, shells):
    """The config of a radial-shell relabel (supervisor decision 56) from the base
    config derived against the parent labels: every attribute list naming a parent MA
    label names its shell labels instead (Ground, TerminalAttributes, edge lists), and
    every Dielectric entry naming parent labels becomes one entry per shell ordinal -
    Attributes = the shell labels of that ordinal over the entry's parents, the same
    Type / layer / edge settings, fresh indices after the base entries' maximum - so a
    participation of the type sums the shells and the per-shell matrices are recorded.
    Returns (config, {index: shell record}) with the base entry the shells replace."""
    by_parent = shell_labels_by_parent(shells)
    by_label = {int(shell["Label"]): shell for shell in shells["Shells"]}
    parents = set(by_parent)

    def expand(attributes):
        out = []
        for attribute in attributes:
            out.extend([label for _, label in by_parent[int(attribute)]] if int(attribute) in parents else [attribute])
        return out

    expanded = json.loads(json.dumps(config))
    boundaries = expanded["Boundaries"]
    if "Ground" in boundaries:
        boundaries["Ground"]["Attributes"] = expand(boundaries["Ground"]["Attributes"])
    for entry in boundaries["PrescribedPotential"]:
        if "TerminalAttributes" in entry:
            entry["TerminalAttributes"] = expand(entry["TerminalAttributes"])
    dielectric = boundaries.get("Postprocessing", {}).get("Dielectric", [])
    next_index = max(int(entry["Index"]) for entry in dielectric) + 1 if dielectric else 1
    entries, shell_map = [], {}
    for entry in dielectric:
        entry["EdgeAttributes"] = expand(entry.get("EdgeAttributes", []))
        if "EdgeExcludeAttributes" in entry:
            entry["EdgeExcludeAttributes"] = expand(entry["EdgeExcludeAttributes"])
        own = [int(a) for a in entry["Attributes"] if int(a) in parents]
        if not own:
            entries.append(entry)
            continue
        if len(own) != len(entry["Attributes"]):
            raise CaseInputError(f"dielectric interface {entry['Index']} mixes relabeled parents {own} with other attributes "
                                 f"{entry['Attributes']}: the shell expansion needs a pure MA entry")
        ordinals = sorted({ordinal for parent in own for ordinal, _ in by_parent[parent]})
        for ordinal in ordinals:
            labels = [label for parent in own for o, label in by_parent[parent] if o == ordinal]
            shell_entry = dict(entry, Index=next_index, Attributes=labels)
            entries.append(shell_entry)
            records = [by_label[label] for label in labels]
            shell_map[next_index] = {"BaseIndex": int(entry["Index"]), "Type": entry["Type"], "Ordinal": ordinal,
                                     "Kind": records[0]["Kind"], "Ring": records[0]["Ring"],
                                     "InnerRadius": records[0]["InnerRadius"], "OuterRadius": records[0]["OuterRadius"],
                                     "Labels": labels, "Parents": [int(r["Parent"]) for r in records],
                                     "Area": sum(float(r["Area"]) for r in records)}
            next_index += 1
    # The base order is kept (the shells take the place of the entry they replace).
    boundaries["Postprocessing"]["Dielectric"] = entries
    return expanded, shell_map


def config_attributes(config):
    """Every boundary attribute a config names, by field."""
    named = {"Ground": list(config["Boundaries"].get("Ground", {}).get("Attributes", []))}
    named["PrescribedPotential"] = sorted({a for e in config["Boundaries"]["PrescribedPotential"] for a in e.get("Attributes", [])})
    named["TerminalAttributes"] = sorted({a for e in config["Boundaries"]["PrescribedPotential"] for a in e.get("TerminalAttributes", [])})
    for entry in config["Boundaries"].get("Postprocessing", {}).get("Dielectric", []):
        named[f"Dielectric[{entry['Index']}].Attributes"] = list(entry.get("Attributes", []))
        named[f"Dielectric[{entry['Index']}].EdgeAttributes"] = list(entry.get("EdgeAttributes", []))
        named[f"Dielectric[{entry['Index']}].EdgeExcludeAttributes"] = list(entry.get("EdgeExcludeAttributes", []))
    return named


def check_attributes(config, available):
    """Every attribute the config names must be a boundary attribute of the mesh
    ($PhysicalNames); fail closed otherwise (before any submission)."""
    named = config_attributes(config)
    missing = {field: sorted(set(values) - set(available)) for field, values in named.items() if set(values) - set(available)}
    if missing:
        raise CaseInputError(f"the config names boundary attributes the mesh does not define: {missing} "
                             f"(mesh $PhysicalNames {sorted(available)})")
    return {"Passed": True, "Named": named, "MeshBoundaryAttributes": sorted(available)}


def strip_paths(config):
    """A config with its path fields (PATH_FIELDS) removed: the comparison view."""
    view = json.loads(json.dumps(config))
    view["Model"].pop("Mesh", None)
    view["Problem"].pop("Output", None)
    for entry in view["Boundaries"]["PrescribedPotential"]:
        entry["DataFile"] = Path(entry["DataFile"]).name
    return view


def config_differences(derived, reference, *, ignore_solver=()):
    """Field paths where the derived and the reference configs differ outside
    PATH_FIELDS (DataFile compared by name) and the named Solver fields."""
    a, b = strip_paths(derived), strip_paths(reference)
    for field in ignore_solver:
        section, key = field.split(".") if "." in field else (None, field)
        for view in (a, b):
            target = view["Solver"][section] if section else view["Solver"]
            target.pop(key, None)
    return _differences(a, b, "")


def _differences(a, b, path):
    out = []
    if isinstance(a, dict) and isinstance(b, dict):
        for key in sorted(set(a) | set(b)):
            if key not in a or key not in b:
                out.append(f"{path}/{key}: only in {'derived' if key in a else 'reference'}")
            else:
                out += _differences(a[key], b[key], f"{path}/{key}")
    elif isinstance(a, list) and isinstance(b, list):
        if len(a) != len(b):
            out.append(f"{path}: {len(a)} vs {len(b)} entries")
        else:
            for index, (x, y) in enumerate(zip(a, b)):
                out += _differences(x, y, f"{path}[{index}]")
    elif a != b:
        out.append(f"{path}: {a!r} vs {b!r}")
    return out
