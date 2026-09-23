# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""The Palace-loadable Version-3 process library qualify writes (decision 65;
qualify/qualify_library.py process_library_entries / preflight_process_library,
qualify/write_process_library.py): the header from the sources (fail closed on a
disagreement), every file copied under models/<slug>/, ThinMatrix null + NotLoadable
without thin matrices, the geometry-only preflight variant, the basis points derived
from the trace vertices, and the merge path rewriting a previous Version-1 library."""
import json
from pathlib import Path
import shutil
import sys
import tempfile
import unittest

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE / "qualify"))
import qualify_library  # noqa: E402
import write_process_library  # noqa: E402

DEVICE_SOURCE = Path("/tmp/coupon-device-transmon-delaunay-20260921/device/sources/spatial-4-edge-e51d7380245e")
FABRICATION = {"LengthUnit": "um", "MetalThickness": 0.1, "OveretchDepth": 0.05,
               "InterfaceLayers": {"SA": {"Thickness": 0.002, "Permittivity": 4.0},
                                   "MS": {"Thickness": 0.002, "Permittivity": 11.47},
                                   "MA": {"Thickness": 0.002, "Permittivity": 10.0}}}
TRACE_VERTICES = ("vertex,x,y,z,basis,conductor\n"
                  "1,1.0000000000000000e+00,0.0000000000000000e+00,0.0000000000000000e+00,2,0\n"
                  "2,0.0000000000000000e+00,1.0000000000000000e+00,0.0000000000000000e+00,0,1\n"
                  "3,0.0000000000000000e+00,0.0000000000000000e+00,1.0000000000000000e+00,1,0\n")


def source_model(name, *, thin=True, support=True):
    model = {"Name": name, "Topology": "SpatialEdgeCluster",
             "FabricatedMatrix": "generated/fabricated/domain-response-matrix.csv",
             "ThinMatrix": "generated/thin/domain-response-matrix.csv",
             "FabricatedSurfaceMatrix": "generated/fabricated/surface-response-matrix.csv",
             "ThinSurfaceMatrix": "generated/thin/surface-response-matrix.csv",
             "BasisPoints": "basis-points.csv",
             "TraceMesh": {"Vertices": "trace-vertices.csv", "Triangles": "trace-triangles.csv"},
             "Edges": [{"Point": [0.0, 0.0, 0.0], "GapDirection": [1.0, 0.0, 0.0], "ProcessNormal": [0.0, 0.0, 1.0],
                        "Interval": [-1.0, 1.0], "Conductor": 1, "InterfaceSlot": 0}],
             "SupportPoints": [[0.0, 0.0, 0.0]] * 8}
    if not support:
        del model["SupportPoints"]
    return model, thin


def write_source(directory, name, *, thin=True, support=True, header=None, basis_points=True):
    """A source directory as the device adapter writes it (thin matrices only when `thin`)."""
    directory.mkdir(parents=True, exist_ok=True)
    model, _ = source_model(name, support=support)
    library = {"Version": 3, "TraceLiftVersion": 3, "Name": name, "MatchingRadius": 2.0, "ExhaustiveSpatialClosure": True,
               "Fabrication": FABRICATION, "Models": [model]}
    library.update(header or {})
    (directory / "process-library.json").write_text(json.dumps(library, indent=2) + "\n")
    (directory / "trace-vertices.csv").write_text(TRACE_VERTICES)
    (directory / "trace-triangles.csv").write_text("triangle,vertex_i,vertex_j,vertex_k\n1,1,2,3\n")
    if basis_points:
        (directory / "basis-points.csv").write_text("x,y,z\n0,0,1\n1,0,0\n")
    if thin:
        (directory / "generated" / "thin").mkdir(parents=True, exist_ok=True)
        (directory / "generated" / "thin" / "domain-response-matrix.csv").write_text("basis_i,basis_j,Q_ij (J)\n1,1,0.5\n")
        (directory / "generated" / "thin" / "surface-response-matrix.csv").write_text("interface,basis_i,basis_j,Q_ij (J)\n1,1,1,0.1\n")
    return directory / "process-library.json"


def write_reducer(directory):
    directory.mkdir(parents=True, exist_ok=True)
    (directory / "domain-response-matrix.csv").write_text("basis_i,basis_j,Q_ij (J)\n1,1,1.5\n")
    (directory / "surface-response-matrix.csv").write_text("interface,basis_i,basis_j,Q_ij (J)\n1,1,1,0.3\n")


class ProcessLibraryWriterTest(unittest.TestCase):
    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp())

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)

    def case(self, case_id, name, *, thin, root, header=None, support=True, basis_points=True):
        """A record / context / manifest triple as run_qualify hands them to the writer."""
        source = write_source(self.tmp / "sources" / case_id, name, thin=thin, header=header, support=support,
                              basis_points=basis_points)
        reducer = root / case_id / "results" / "main" / f"{case_id}-p4" / "reducer"
        write_reducer(reducer)
        manifest = {"RepositoryRoot": ".", "Cases": [{"Id": case_id, "Source": {"Directory": str(source.parent),
                                                                                  "Files": {"ProcessLibrary": {"Name": "process-library.json"}}}}]}
        record = {"Case": case_id, "Root": str(root / case_id), "MainReducer": str(reducer),
                  "Mesh": {"Local": "/nowhere/identity.msh", "SHA256": "0" * 64},
                  "Qualification": {"Verdict": "PendingQualification", "Path": str(root / case_id / "qualification.json"),
                                    "ReferenceAnchor": None}}
        context = {"layout": [{"Role": "main", "Prefix": f"{case_id}-p4", "Order": 4}], "ma_tail": None}
        return record, context, manifest

    def write(self, cases, root, merge_into=None):
        entries = {}
        for case in cases:
            for entry in case[2]["Cases"]:
                entries.setdefault(entry["Id"], entry)
        manifest = {"RepositoryRoot": ".", "Cases": list(entries.values())}
        return qualify_library.process_library_entries([case[0] for case in cases], {case[0]["Case"]: case[1] for case in cases},
                                                       manifest_path=self.tmp / "manifest.json", manifest=manifest, root=root,
                                                       merge_into=merge_into)

    def test_header_copies_and_loadable_model(self):
        root = self.tmp / "root"
        library = self.write([self.case("case-a", "model_a", thin=True, root=root)], root)
        self.assertEqual(library["Version"], 3)
        self.assertEqual(library["MatchingRadius"], 2.0)
        self.assertEqual(library["TraceLiftVersion"], 3)
        self.assertEqual(library["Fabrication"], FABRICATION)
        self.assertNotIn("PreflightOnly", library)
        model = library["Models"][0]
        self.assertIsNone(model["NotLoadable"])
        self.assertTrue(library["Loadable"]["Palace"])
        for field, name in qualify_library.MODEL_FILE_NAMES.items():
            self.assertEqual(model[field], f"models/model-a/{name}")
            self.assertTrue((root / model[field]).is_file(), field)
        self.assertEqual(model["TraceMesh"], {"Vertices": "models/model-a/trace-vertices.csv",
                                              "Triangles": "models/model-a/trace-triangles.csv"})
        self.assertEqual((root / model["FabricatedMatrix"]).read_text(), "basis_i,basis_j,Q_ij (J)\n1,1,1.5\n")
        self.assertEqual((root / model["ThinMatrix"]).read_text(), "basis_i,basis_j,Q_ij (J)\n1,1,0.5\n")
        self.assertEqual(model["SupportPoints"], [[0.0, 0.0, 0.0]] * 8)
        self.assertFalse(model["LibraryQualified"])
        self.assertEqual(model["MA"]["MA_sharp"], None)
        self.assertEqual(qualify_library.preflight_process_library(library)["Models"][0]["ThinMatrix"], model["ThinMatrix"])

    def test_missing_thin_matrices_are_null_and_not_loadable_with_a_preflight_variant(self):
        root = self.tmp / "root"
        library = self.write([self.case("case-a", "model_a", thin=False, root=root)], root)
        model = library["Models"][0]
        self.assertIsNone(model["ThinMatrix"])
        self.assertIsNone(model["ThinSurfaceMatrix"])
        self.assertEqual(model["NotLoadable"]["Missing"], ["ThinMatrix", "ThinSurfaceMatrix"])
        self.assertIn("generated/thin/domain-response-matrix.csv", model["NotLoadable"]["Reason"])
        self.assertEqual(library["Loadable"], {"Palace": False, "Models": 0, "NotLoadable": ["model_a"],
                                               "Preflight": "process-library-preflight.json"})
        preflight = qualify_library.preflight_process_library(library)
        self.assertTrue(preflight["PreflightOnly"])
        self.assertEqual(preflight["Models"][0]["ThinMatrix"], "models/model-a/thin-domain-response-matrix.csv")
        self.assertEqual(preflight["Models"][0]["ThinSurfaceMatrix"], "models/model-a/thin-surface-response-matrix.csv")
        self.assertFalse((root / preflight["Models"][0]["ThinMatrix"]).exists())
        self.assertIsNotNone(preflight["Models"][0]["NotLoadable"])
        # The honest file is untouched by the variant.
        self.assertIsNone(library["Models"][0]["ThinMatrix"])

    def thin_case(self, fabricated, *, root, name, verdict="PendingQualification"):
        """The thin pair of a fabricated case (decision 66): Kind thin, its own reducer, the
        recorded cutoff; the manifest case names FabricatedCase."""
        record, context, manifest = fabricated
        case_id = record["Case"] + "-thin"
        reducer = root / case_id / "results" / "main" / f"{case_id}-p4" / "reducer"
        reducer.mkdir(parents=True, exist_ok=True)
        (reducer / "domain-response-matrix.csv").write_text("basis_i,basis_j,Q_ij (J)\n1,1,0.25\n")
        (reducer / "surface-response-matrix.csv").write_text("interface,basis_i,basis_j,Q_ij (J)\n1,1,1,0.05\n")
        thin_manifest = json.loads(json.dumps(manifest))
        thin_case = json.loads(json.dumps(manifest["Cases"][0]))
        thin_case.update({"Id": case_id, "Kind": "thin", "FabricatedCase": record["Case"]})
        thin_manifest["Cases"].append(thin_case)   # the manifest carries the pair
        thin_record = {"Case": case_id, "Kind": "thin", "FabricatedCase": record["Case"], "ThinCutoff": 0.002,
                       "Root": str(root / case_id), "MainReducer": str(reducer),
                       "Mesh": {"Local": "/nowhere/identity-thin.msh", "SHA256": "1" * 64},
                       "Qualification": {"Verdict": verdict, "Path": str(root / case_id / "qualification.json"), "ReferenceAnchor": None}}
        thin_context = {"layout": [{"Role": "main", "Prefix": f"{case_id}-p4", "Order": 4}], "ma_tail": None}
        return thin_record, thin_context, thin_manifest

    def test_thin_case_supplies_the_thin_matrices_of_its_fabricated_model(self):
        """Decision 66: a thin case is not a model; its reducer matrices become the ThinMatrix /
        ThinSurfaceMatrix of its fabricated model (copied, loadable), with the recorded cutoff;
        an unpaired thin case is recorded, not an entry."""
        root = self.tmp / "root"
        fabricated = self.case("case-a", "model_a", thin=False, root=root)
        thin = self.thin_case(fabricated, root=root, name="model_a")
        library = self.write([fabricated, thin], root)
        self.assertEqual([model["Name"] for model in library["Models"]], ["model_a"])
        model = library["Models"][0]
        self.assertIsNone(model["NotLoadable"])
        self.assertTrue(library["Loadable"]["Palace"])
        self.assertEqual(model["ThinMatrix"], "models/model-a/thin-domain-response-matrix.csv")
        self.assertEqual((root / model["ThinMatrix"]).read_text(), "basis_i,basis_j,Q_ij (J)\n1,1,0.25\n")
        self.assertEqual((root / model["ThinSurfaceMatrix"]).read_text(), "interface,basis_i,basis_j,Q_ij (J)\n1,1,1,0.05\n")
        self.assertEqual((model["ThinCase"], model["ThinCutoff"], model["ThinQualification"]["Order"]), ("case-a-thin", 0.002, 4))
        self.assertEqual(model["ThinCutoffRule"], qualify_library.THIN_CUTOFF_RULE)
        self.assertEqual(library["Thin"]["Paired"], ["case-a"])
        self.assertEqual(library["Thin"]["Unpaired"], {})
        # The thin case alone: no entry, recorded unpaired.
        alone = self.write([thin], root / "alone")
        self.assertEqual(alone["Models"], [])
        self.assertEqual(alone["Thin"]["Unpaired"], {"case-a": "case-a-thin"})
        # Merge: the thin run follows the fabricated run - the kept model takes the thin matrices.
        first = root / "first"
        previous = self.write([self.case("case-a", "model_a", thin=False, root=first)], first)
        (first / "process-library.json").write_text(json.dumps(previous, indent=2) + "\n")
        self.assertIsNotNone(previous["Models"][0]["NotLoadable"])
        second = root / "second"
        merged = self.write([self.thin_case(self.case("case-a", "model_a", thin=False, root=first), root=second, name="model_a")],
                            second, merge_into=first / "process-library.json")
        self.assertEqual(merged["MergedFrom"]["Kept"], ["model_a"])
        model = merged["Models"][0]
        self.assertIsNone(model["NotLoadable"])
        self.assertEqual(model["ThinCase"], "case-a-thin")
        self.assertEqual((second / model["ThinMatrix"]).read_text(), "basis_i,basis_j,Q_ij (J)\n1,1,0.25\n")
        self.assertEqual((second / model["FabricatedMatrix"]).read_text(), "basis_i,basis_j,Q_ij (J)\n1,1,1.5\n")
        self.assertTrue(merged["Loadable"]["Palace"])

    def test_shelled_surface_matrix_is_collapsed_to_the_model_interfaces(self):
        """A radial-shell run labels the MA shells by their own interface indices; the library
        copy sums them into the BaseIndex (the value Palace reads), per (edge, R, basis pair)."""
        root = self.tmp / "root"
        record, context, manifest = self.case("case-a", "model_a", thin=True, root=root)
        header = ("interface,     edge,                      R (m),  basis_i,  basis_j,                   Q_ij (J),"
                  "            Q_ij normal (J),        Q_ij tangential (J),             Q_total_ij (J),"
                  "      Q_total_ij normal (J),  Q_total_ij tangential (J)\n")
        def row(interface, edge, i, j, q):
            return (f" {interface:.2e}, {edge:.2e},        +2.000000000000e-06, {i:.2e}, {j:.2e},        {q:+.12e},"
                    f"        {q:+.12e},        +0.000000000000e+00,        {2 * q:+.12e},        {2 * q:+.12e},        +0.000000000000e+00\n")
        Path(record["MainReducer"], "surface-response-matrix.csv").write_text(
            header + row(2, 1, 1, 1, 0.5) + row(4, 1, 1, 1, 1.0) + row(5, 1, 1, 1, 2.0) + row(4, 2, 1, 1, 8.0) + row(3, 1, 1, 1, 0.25))
        record["Inputs"] = {"RadialShells": {"Interfaces": {"4": {"BaseIndex": 1, "Type": "MA", "Ordinal": 1},
                                                             "5": {"BaseIndex": 1, "Type": "MA", "Ordinal": 2}}}}
        library = self.write([(record, context, manifest)], root)
        model = library["Models"][0]
        self.assertEqual(model["FabricatedSurfaceMatrix"], "models/model-a/fabricated-surface-response-matrix.csv")
        self.assertTrue(model["FabricatedSurfaceMatrixShelled"].endswith("reducer/surface-response-matrix.csv"))
        self.assertEqual(library["CollapsedSurfaceMatrices"]["model_a"]["ShellIndices"], [4, 5])
        self.assertEqual(library["CollapsedSurfaceMatrices"]["model_a"]["CollapsedRows"], 4)
        import csv
        with open(root / model["FabricatedSurfaceMatrix"], newline="") as stream:
            rows = [[cell.strip() for cell in line] for line in csv.reader(stream)]
        self.assertEqual(rows[0][0], "interface")
        by_key = {(int(float(r[0])), int(float(r[1]))): float(r[8]) for r in rows[1:]}
        self.assertEqual(set(by_key), {(2, 1), (1, 1), (1, 2), (3, 1)})
        self.assertAlmostEqual(by_key[(1, 1)], 2 * (1.0 + 2.0))
        self.assertAlmostEqual(by_key[(1, 2)], 16.0)
        self.assertAlmostEqual(by_key[(2, 1)], 1.0)

    def test_missing_support_points_are_not_loadable(self):
        root = self.tmp / "root"
        library = self.write([self.case("case-a", "model_a", thin=True, root=root, support=False)], root)
        self.assertEqual(library["Models"][0]["NotLoadable"]["Reason"], qualify_library.NOT_LOADABLE_SUPPORT)
        self.assertEqual(library["Models"][0]["NotLoadable"]["Missing"], [])

    def test_disagreeing_source_headers_stop_the_writer(self):
        root = self.tmp / "root"
        cases = [self.case("case-a", "model_a", thin=True, root=root),
                 self.case("case-b", "model_b", thin=True, root=root, header={"MatchingRadius": 2.5})]
        with self.assertRaisesRegex(ValueError, "disagree on MatchingRadius"):
            self.write(cases, root)
        other_layers = json.loads(json.dumps(FABRICATION))
        other_layers["InterfaceLayers"]["MS"]["Permittivity"] = 11.45
        cases = [self.case("case-a", "model_a", thin=True, root=root),
                 self.case("case-b", "model_b", thin=True, root=root, header={"Fabrication": other_layers})]
        with self.assertRaisesRegex(ValueError, "disagree on Fabrication"):
            self.write(cases, root)
        cases = [self.case("case-a", "model_a", thin=True, root=root),
                 self.case("case-b", "model_b", thin=True, root=root, header={"Version": 2})]
        with self.assertRaisesRegex(ValueError, "Version 2, not 3"):
            self.write(cases, root)

    def test_basis_points_are_derived_from_the_trace_vertices_when_absent(self):
        root = self.tmp / "root"
        library = self.write([self.case("case-a", "model_a", thin=True, root=root, basis_points=False)], root)
        model = library["Models"][0]
        self.assertEqual(model["BasisPointsRule"], qualify_library.BASIS_POINTS_FROM_TRACE)
        self.assertEqual((root / model["BasisPoints"]).read_text(),
                         "x,y,z\n0.0000000000000000e+00,0.0000000000000000e+00,1.0000000000000000e+00\n"
                         "1.0000000000000000e+00,0.0000000000000000e+00,0.0000000000000000e+00\n")

    @unittest.skipUnless((DEVICE_SOURCE / "basis-points.csv").is_file(), "the device source directory is needed")
    def test_derived_basis_points_reproduce_the_producer_file(self):
        destination = self.tmp / "basis-points.csv"
        qualify_library.write_basis_points_from_trace_vertices(DEVICE_SOURCE / "trace-vertices.csv", destination)
        self.assertEqual(destination.read_bytes(), (DEVICE_SOURCE / "basis-points.csv").read_bytes())

    def test_merge_rewrites_a_previous_version_1_library(self):
        """A previous run's Version-1 file (root-relative fabricated matrices, source-relative
        thin / basis / trace paths, no header): its models are kept, resolved through the
        previous root / the case root / the source directory, copied into this root."""
        previous_root = self.tmp / "previous"
        source = write_source(self.tmp / "sources" / "case-p", "model_p", thin=False)
        write_reducer(previous_root / "case-p" / "results" / "main" / "case-p-p4" / "reducer")
        model = json.loads(source.read_text())["Models"][0]
        model.update({"FabricatedMatrix": "case-p/results/main/case-p-p4/reducer/domain-response-matrix.csv",
                      "FabricatedSurfaceMatrix": "case-p/results/main/case-p-p4/reducer/surface-response-matrix.csv",
                      "Qualification": {"Verdict": "PendingQualification", "Record": str(previous_root / "case-p" / "qualification.json")},
                      "LibraryQualified": False, "SourceProcessLibrary": {"Path": str(source), "SHA256": "x"}})
        previous = previous_root / "process-library.json"
        previous.write_text(json.dumps({"Version": 1, "Root": str(previous_root), "Models": [model]}))
        before = previous.read_bytes()
        root = self.tmp / "root"
        library = self.write([self.case("case-a", "model_a", thin=True, root=root)], root, merge_into=previous)
        self.assertEqual([model["Name"] for model in library["Models"]], ["model_p", "model_a"])
        self.assertEqual(library["MergedFrom"]["Kept"], ["model_p"])
        kept = library["Models"][0]
        self.assertEqual(kept["FabricatedMatrix"], "models/model-p/fabricated-domain-response-matrix.csv")
        self.assertTrue((root / kept["FabricatedMatrix"]).is_file())
        self.assertEqual(kept["BasisPoints"], "models/model-p/basis-points.csv")
        self.assertIsNone(kept["ThinMatrix"])
        self.assertEqual(library["Loadable"]["NotLoadable"], ["model_p"])
        self.assertEqual(previous.read_bytes(), before)
        root2 = self.tmp / "root2"
        root2.mkdir()
        write_process_library.main(["--previous", str(previous), "--root", str(root2)])
        rewritten = json.loads((root2 / "process-library.json").read_text())
        self.assertEqual(rewritten["Version"], 3)
        self.assertEqual([model["Name"] for model in rewritten["Models"]], ["model_p"])
        preflight = json.loads((root2 / "process-library-preflight.json").read_text())
        self.assertTrue(preflight["PreflightOnly"])
        self.assertEqual(preflight["Models"][0]["ThinMatrix"], "models/model-p/thin-domain-response-matrix.csv")

    def test_run_pairs_a_finished_thin_run_with_the_previous_library(self):
        """write_process_library --run: the thin qualify run's library-qualification.json (its
        recorded stage layout and reducer, the manifest through its build record) supplies the
        thin matrices of the previous fabricated library's model without any job or fetch."""
        first = self.tmp / "first"
        fabricated = self.case("case-a", "model_a", thin=False, root=first)
        previous = self.write([fabricated], first)
        (first / "process-library.json").write_text(json.dumps(previous, indent=2) + "\n")
        self.assertEqual(previous["Loadable"]["NotLoadable"], ["model_a"])
        second = self.tmp / "second"
        thin_record, thin_context, thin_manifest = self.thin_case(fabricated, root=second, name="model_a")
        manifest_path = self.tmp / "thin-manifest.json"
        manifest_path.write_text(json.dumps(thin_manifest, indent=2) + "\n")
        build_record = self.tmp / "library-build-thin.json"
        build_record.write_text(json.dumps({"Library": {"Manifest": {"Path": str(manifest_path)}}}) + "\n")
        run_record = second / "library-qualification.json"
        run_record.parent.mkdir(parents=True, exist_ok=True)
        run_record.write_text(json.dumps({"BuildRecord": {"Path": str(build_record)},
                                          "Cases": [thin_record | {"Stages": thin_context["layout"],
                                                                   "MATail": {"Applied": False}},
                                                    {"Case": "case-b-thin", "Stages": None, "Qualification": None}]}) + "\n")
        root = self.tmp / "root"
        write_process_library.main(["--previous", str(first / "process-library.json"), "--root", str(root),
                                    "--run", str(run_record)])
        library = json.loads((root / "process-library.json").read_text())
        self.assertEqual([model["Name"] for model in library["Models"]], ["model_a"])
        model = library["Models"][0]
        self.assertIsNone(model["NotLoadable"])
        self.assertTrue(library["Loadable"]["Palace"])
        self.assertEqual((model["ThinCase"], model["ThinCutoff"], model["ThinQualification"]["Order"]), ("case-a-thin", 0.002, 4))
        self.assertEqual((root / model["ThinMatrix"]).read_text(), "basis_i,basis_j,Q_ij (J)\n1,1,0.25\n")
        self.assertEqual((root / model["FabricatedMatrix"]).read_text(), "basis_i,basis_j,Q_ij (J)\n1,1,1.5\n")
        self.assertEqual(library["MergedFrom"]["Kept"], ["model_a"])
        self.assertEqual(library["Thin"]["Paired"], ["case-a"])
        # --run repeats (a root relaunched for other cases); the same case twice is refused.
        other = second / "library-qualification-other.json"
        other.write_text(json.dumps({"BuildRecord": {"Path": str(build_record)},
                                     "Cases": [{"Case": "case-c-thin", "Stages": None, "Qualification": None}]}) + "\n")
        write_process_library.main(["--previous", str(first / "process-library.json"), "--root", str(self.tmp / "root3"),
                                    "--run", str(run_record), "--run", str(other)])
        self.assertEqual(json.loads((self.tmp / "root3" / "process-library.json").read_text())["Thin"]["Paired"], ["case-a"])
        with self.assertRaisesRegex(ValueError, "both analyzed \\['case-a-thin'\\]"):
            write_process_library.main(["--previous", str(first / "process-library.json"), "--root", str(self.tmp / "root4"),
                                        "--run", str(run_record), "--run", str(run_record)])


if __name__ == "__main__":
    unittest.main()
