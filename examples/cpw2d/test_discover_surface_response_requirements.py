#!/usr/bin/env python3

import json
import os
import stat
import sys
import tempfile
import textwrap
import unittest
from pathlib import Path

CPW2D = Path(__file__).parent
SPATIAL = CPW2D.parent / "cpw3d_surface" / "spatial_coupon"
sys.path.insert(0, str(CPW2D))
sys.path.insert(0, str(SPATIAL))
import cluster_signature_geometry  # noqa: E402
import device_coupons  # noqa: E402
import discover_surface_response_requirements as DISCOVERY  # noqa: E402
import generate_spatial_response as SPATIAL_RESPONSE  # noqa: E402
import signature_library  # noqa: E402  (surface_response_identification, on sys.path through DISCOVERY)

# The PRODUCTION MirrorFormedContract requirement records of the rerun-3 R3 registration jobs'
# pass-01 closures (osc003-O3 d7c875318447, sct002-S2p 613e6f656498; the M7 binary of record on
# main c95769cd28, library s2-r1p9-v3-c, R = 1.9 um): the keys whose closure stalled at pass-02
# on the signature-only placeholder (decision 599).
MIRROR_FORMED_FIXTURE = CPW2D / "testdata" / "mirror-formed-contract-requirements.json"
# The operator's fail-closed Note on a signature-only model (surfaceresponseoperator.cpp, the
# mirror-formed status loop; decision 557 (4)).
SIGNATURE_ONLY_MODEL_NOTE = ("mirror-formed cluster configuration matched by a model without Signature-keyed Edges: "
                             "the real / image split cannot be mapped onto its Edges (fail closed to the raw reading, "
                             "decision 557)")


def mirror_formed_fixture():
    fixture = json.loads(MIRROR_FORMED_FIXTURE.read_text())
    return float(fixture["MatchingRadius"]), fixture["Requirements"]


class DiscoverSurfaceResponseRequirementsTest(unittest.TestCase):
    def requirement(self, topology, geometry, edges=None):
        if edges is not None:
            geometry = {**geometry, "Edges": edges}
        return {
            "Topology": topology,
            "Geometry": geometry,
            "BoundaryCondition": {"Type": "PEC"},
            "Interfaces": [
                {"Slot": 0, "Target": 1, "Type": "SA"},
                {"Slot": 0, "Target": 2, "Type": "MS"},
            ],
        }

    def test_single_and_multiple_conductor_placeholders(self):
        single = self.requirement(
            "SpatialEdgeCluster",
            {"EdgeCount": 2},
            [
                {
                    "Conductor": 1,
                    "Point": [0.0, 0.0, 0.0],
                    "GapDirection": [1.0, 0.0, 0.0],
                    "ProcessNormal": [0.0, 1.0, 0.0],
                    "Interval": [-1.0, 1.0],
                    "InterfaceSlot": 0,
                    "BoundaryCondition": {"Type": "PEC"},
                },
                {
                    "Conductor": 1,
                    "Point": [1.0, 0.0, 0.0],
                    "GapDirection": [-1.0, 0.0, 0.0],
                    "ProcessNormal": [0.0, 1.0, 0.0],
                    "Interval": [-1.0, 1.0],
                    "InterfaceSlot": 0,
                    "BoundaryCondition": {"Type": "PEC"},
                },
            ],
        )
        single["Geometry"]["PlanViewBoundary"] = [
            {
                "Conductor": 1,
                "Segments": [[[0, 0, 0], [1, 0, 0]]],
                "ContinuationSegments": [],
            }
        ]
        single["Geometry"]["MaskRegularization"] = {"Version": 1}
        _, model = DISCOVERY.placeholder_model(
            single,
            2.0,
            {"MetalThickness": 0.1, "OveretchDepth": 0.05},
        )
        self.assertIn("Reference", model)
        self.assertNotIn("ConductorReferences", model)
        self.assertEqual(
            model["PlanViewBoundary"], single["Geometry"]["PlanViewBoundary"]
        )
        self.assertEqual(
            model["MaskRegularization"], single["Geometry"]["MaskRegularization"]
        )
        self.assertEqual(len(model["SupportPoints"]), 8)

        different = self.requirement(
            "DifferentConductorGap", {"EdgeCount": 2, "Separation": 1.0}
        )
        _, model = DISCOVERY.placeholder_model(different, 2.0)
        self.assertEqual(len(model["ConductorReferences"]), 2)
        self.assertNotIn("Reference", model)

    def test_virtual_and_production_spatial_support_are_identical(self):
        direction = 2.0**-0.5
        geometry = {
            "EdgeCount": 2,
            "Edges": [
                {
                    "Conductor": 1,
                    "Point": [0.123456789, 0.0, -0.987654321],
                    "GapDirection": [direction, 0.0, direction],
                    "ProcessNormal": [0.0, 1.0, 0.0],
                    "Interval": [-2.0, 0.0],
                    "InterfaceSlot": 0,
                    "BoundaryCondition": {"Type": "PEC"},
                },
                {
                    "Conductor": 1,
                    "Point": [1.234567891, 0.0, 0.345678912],
                    "GapDirection": [-direction, 0.0, -direction],
                    "ProcessNormal": [0.0, 1.0, 0.0],
                    "Interval": [0.0, 2.0],
                    "InterfaceSlot": 0,
                    "BoundaryCondition": {"Type": "PEC"},
                },
            ],
        }
        fabrication = {"MetalThickness": 0.1, "OveretchDepth": 0.05}
        virtual = DISCOVERY.spatial_support_points(geometry, 2.0, fabrication)
        coupon = {
            "Topology": "SpatialEdgeCluster",
            "Geometry": geometry,
            "BoundaryCondition": {"Type": "PEC"},
        }
        frame, edges, _ = SPATIAL_RESPONSE.normalize_geometry(coupon, 2.0)
        lower, upper = SPATIAL_RESPONSE.coupon_bounds(edges, 2.0, 0.1, 0.05)
        production = SPATIAL_RESPONSE.matching_support_points(
            lower, upper, frame, 2.0
        )
        self.assertEqual(virtual, production)

    # Two version-2 Missing requirements of a fake device: an IsolatedEdge (Hash aaaa...) and a
    # 3-edge ParallelEdgeCluster (Hash bbbb...); the fake Palace below reads every pass's
    # library and reports a requirement Exact once a placeholder named after its Hash exists.
    FAKE_REQUIREMENTS = [
        {"Topology": "IsolatedEdge", "Hash": "a" * 64, "Count": 4, "Instances": 2, "TotalEdgeLength": 10.0,
         "BoundaryCondition": {"Type": "PEC"}, "Interfaces": [{"Slot": 0, "Target": 1, "Type": "SA"}],
         "Signature": {"Type": "IsolatedEdge", "Interfaces": ["SA"], "Law": "{\"Type\":\"PEC\"}"}},
        {"Topology": "ParallelEdgeCluster", "Hash": "b" * 64, "Count": 7, "Instances": 2, "TotalEdgeLength": 8.2,
         "BoundaryCondition": {"Type": "PEC"}, "Interfaces": [{"Slot": 0, "Target": 1, "Type": "SA"}],
         "Geometry": {"EdgeCount": 3},
         "Signature": {"Type": "ParallelEdgeCluster", "Edges": [
             {"Conductor": 1, "GapSide": -1, "Interfaces": ["SA"], "Law": "{\"Type\":\"PEC\"}", "OffsetOverR": 0.0},
             {"Conductor": 1, "GapSide": 1, "Interfaces": ["SA"], "Law": "{\"Type\":\"PEC\"}", "OffsetOverR": 1.052632},
             {"Conductor": 1, "GapSide": -1, "Interfaces": ["SA"], "Law": "{\"Type\":\"PEC\"}", "OffsetOverR": 2.105263}]}},
    ]

    def fake_device(self, tmp):
        """A device config + source library + an executable fake Palace (geometry preflight only)."""
        requirements = tmp / "fake-requirements.json"
        requirements.write_text(json.dumps(self.FAKE_REQUIREMENTS))
        library = tmp / "process-library.json"
        library.write_text(json.dumps({"Version": 3, "Name": "fake-seed", "MatchingRadius": 1.9, "Models": []}))
        config = tmp / "device.json"
        config.write_text(json.dumps({
            "Problem": {"Type": "Electrostatic", "Output": str(tmp / "postpro")},
            "Solver": {"Electrostatic": {"ResponseCorrection": {"Library": str(library), "UnmatchedPolicy": "Warn"}}}}))
        palace = tmp / "fake-palace.py"
        palace.write_text(textwrap.dedent(f"""\
            #!{sys.executable}
            import json, pathlib, sys
            assert sys.argv[1] == "--surface-response-preflight"
            config = json.load(open(sys.argv[2]))
            response = config["Solver"]["Electrostatic"]["ResponseCorrection"]
            library = json.load(open(response["Library"]))
            names = {{model["Name"] for model in library["Models"]}}
            requirements = json.load(open({str(requirements)!r}))
            counts = {{"Exact": 0, "Interpolated": 0, "Missing": 0}}
            lengths = {{"Exact": 0.0, "Interpolated": 0.0, "Missing": 0.0}}
            for requirement in requirements:
                placeholder = "__preflight_placeholder_" + requirement["Hash"][:16]
                if placeholder in names:
                    requirement["Status"] = "Exact"
                    requirement["SelectedModels"] = [{{"Name": placeholder, "Weight": 1.0}}]
                else:
                    requirement["Status"] = "Missing"
                counts[requirement["Status"]] += requirement["Count"]
                lengths[requirement["Status"]] += requirement["TotalEdgeLength"]
            output = pathlib.Path(config["Problem"]["Output"]); output.mkdir(parents=True, exist_ok=True)
            manifest = {{"Version": 2, "Complete": counts["Missing"] == 0, "Requirements": requirements,
                        "Library": {{"Name": library["Name"], "Path": response["Library"],
                                    "MatchingRadius": library["MatchingRadius"]}},
                        "Summary": {{"Counts": counts, "TotalEdgeLengths": lengths}}}}
            (output / "surface-response-requirements.json").write_text(json.dumps(manifest))
            """))
        palace.chmod(palace.stat().st_mode | stat.S_IXUSR)
        return config, palace

    def test_omitted_requirement_gets_no_placeholder_and_stays_missing(self):
        """`--omit-requirement HASH_PREFIX`: the matching Missing requirement gets no placeholder in
        any pass, does not stall the closure, stays Missing in the final manifest (counted in Summary.Missing, Complete false) and is listed under
        OmittedRequirements of the manifest and the closure history; the other requirement
        closes as before. A prefix matching nothing fails closed."""
        with tempfile.TemporaryDirectory(prefix="discover-omit-") as tmp:
            tmp = Path(tmp)
            config, palace = self.fake_device(tmp)
            # Without an omission the fake device closes in two passes with two placeholders.
            manifest = DISCOVERY.discover(config, tmp / "full", palace)
            self.assertEqual(manifest["Summary"]["Counts"], {"Exact": 0, "Interpolated": 0, "Missing": 11})
            self.assertEqual(manifest["OmittedRequirements"], [])
            history = json.loads((tmp / "full" / "closure-history.json").read_text())
            self.assertEqual(history["PlaceholderCount"], 2)
            self.assertEqual(len(history["Passes"]), 2)
            self.assertEqual(sorted(p["Topology"] for p in history["Passes"][0]["AddedPlaceholders"]),
                             ["IsolatedEdge", "ParallelEdgeCluster"])
            # The 3-edge cluster omitted: one placeholder, the omitted requirement Missing throughout.
            manifest = DISCOVERY.discover(config, tmp / "omit", palace, omit_requirements=["bbbbbbbbbbbb"])
            self.assertFalse(manifest["Complete"])
            self.assertEqual(manifest["Summary"]["Counts"], {"Exact": 0, "Interpolated": 0, "Missing": 11})
            by_hash = {requirement["Hash"]: requirement for requirement in manifest["Requirements"]}
            self.assertEqual(by_hash["a" * 64]["Status"], "Missing")   # the restored placeholder
            self.assertEqual(by_hash["a" * 64]["Reason"], "Missing from source library after exhaustive geometry discovery")
            self.assertEqual(by_hash["b" * 64]["Status"], "Missing")   # never had a placeholder
            self.assertNotIn("Reason", by_hash["b" * 64])
            self.assertEqual(manifest["OmittedRequirements"],
                             [{"Hash": "b" * 64, "Prefix": "bbbbbbbbbbbb", "Topology": "ParallelEdgeCluster", "Count": 7,
                               "Instances": 2, "TotalEdgeLength": 8.2}])
            history = json.loads((tmp / "omit" / "closure-history.json").read_text())
            self.assertEqual(history["PlaceholderCount"], 1)
            self.assertEqual(history["OmittedRequirements"], manifest["OmittedRequirements"])
            self.assertEqual(len(history["Passes"]), 2)
            self.assertEqual([p["Topology"] for p in history["Passes"][0]["AddedPlaceholders"]], ["IsolatedEdge"])
            self.assertEqual(history["Passes"][1]["AddedPlaceholders"], [])
            for pass_dir in ("pass-01", "pass-02"):
                library = json.loads((tmp / "omit" / pass_dir / "process-library.json").read_text())
                self.assertNotIn("__preflight_placeholder_" + "b" * 16, {m["Name"] for m in library["Models"]})
            # A prefix that matches no Missing requirement is a mistake, not a no-op.
            with self.assertRaisesRegex(RuntimeError, "matching no Missing requirement.*cccccccc"):
                DISCOVERY.discover(config, tmp / "typo", palace, omit_requirements=["cccccccc"])
            with self.assertRaises(ValueError):
                DISCOVERY.discover(config, tmp / "empty", palace, omit_requirements=[""])

    def test_restore_source_status_marks_only_placeholders_missing(self):
        manifest = {
            "Complete": True,
            "Library": {"Name": "closure", "Path": "/tmp/closure.json"},
            "Summary": {},
            "Requirements": [
                {
                    "Status": "Exact",
                    "Count": 2,
                    "TotalEdgeLength": 3.0,
                    "SelectedModels": [{"Name": "real", "Weight": 1.0}],
                },
                {
                    "Status": "Exact",
                    "Count": 4,
                    "TotalEdgeLength": 5.0,
                    "SelectedModels": [
                        {"Name": "__preflight_placeholder_deadbeef", "Weight": 1.0}
                    ],
                    "NormalizedLibraryDistance": 0.0,
                },
            ],
        }
        source = {"Name": "source", "__SourcePath": Path("/tmp/source.json")}
        result = DISCOVERY.restore_source_status(
            manifest,
            source,
            {
                "__preflight_placeholder_deadbeef": {
                    "Topology": "IsolatedEdge",
                    "Geometry": {"EdgeCount": 1},
                    "BoundaryCondition": {"Type": "PEC"},
                    "Interfaces": [],
                }
            },
        )
        self.assertFalse(result["Complete"])
        self.assertEqual(
            result["Summary"],
            {
                "Counts": {"Exact": 2, "Interpolated": 0, "Missing": 4},
                "TotalEdgeLengths": {
                    "Exact": 3.0,
                    "Interpolated": 0.0,
                    "Missing": 5.0,
                },
            },
        )
        missing = result["Requirements"][1]
        self.assertEqual(missing["Status"], "Missing")
        self.assertNotIn("SelectedModels", missing)
        self.assertNotIn("NormalizedLibraryDistance", missing)


class MirrorFormedContractPlaceholderTest(unittest.TestCase):
    """Decision 599: the discovery placeholder of a SpatialEdgeCluster requirement carrying a
    MirrorFormedContract carries Signature-keyed Edges with the contract's per-Edge Weight and
    the ENTRY MirrorFormed stamp, built by the generator's own functions (the
    device_coupons.stamp_signature_model path), so the M7 operator places it like the real
    coupon's entry and the closure converges; a requirement without a contract keeps the
    signature-only placeholder; an inconsistent contract fails closed by name."""

    def plain_placeholder(self, requirement, radius):
        """The signature-only placeholder of main (decision 599 predecessor): signature_model +
        the BoundaryLawQualification block, nothing else for a SpatialEdgeCluster."""
        model = signature_library.signature_model(
            requirement["Topology"], requirement["Signature"], radius,
            f"__preflight_placeholder_{requirement['Hash'][:16]}", "__preflight_dummy")
        model["BoundaryLawQualification"] = {"Version": 1, "Status": "Unqualified",
                                             "Calibration": "GeometryDiscoveryOnly", "FrequencyUniversal": False}
        return model

    def test_fixture_records_are_production_contract_requirements(self):
        radius, requirements = mirror_formed_fixture()
        self.assertEqual(radius, 1.9)
        self.assertEqual(sorted(requirements), ["o3-f45-d7c875318447", "s2p-f11-613e6f656498"])
        for name, requirement in requirements.items():
            self.assertEqual(requirement["Status"], "Missing")
            self.assertEqual(requirement["Topology"], "SpatialEdgeCluster")
            self.assertIs(requirement["MirrorFormed"], True)
            self.assertTrue(requirement["Hash"].startswith(name.rsplit("-", 1)[1]))
            self.assertEqual(requirement["Hash"], device_coupons.signature_key_hash(requirement["Signature"]))
            self.assertEqual(requirement["MirrorFormedContract"]["Frame"]["Chirality"], 1)
        self.assertEqual(requirements["o3-f45-d7c875318447"]["MirrorFormedContract"]["RealPortions"], [0, 3])
        self.assertEqual(requirements["s2p-f11-613e6f656498"]["MirrorFormedContract"]["RealPortions"], [0, 1])

    def test_contract_placeholder_carries_the_generator_stamp(self):
        """The placeholder's Edges / MirrorFormed record are the real coupon's entry (R3 dry runs
        of record: O3 18 Edges, Weight 1 x 15 (the arc's chords) / 0 / 0 / 1, RealLengthFraction
        0.49976768508632746; S2p 4 Edges, 1 / 1 / 0 / 0, 0.4999999574369773)."""
        radius, requirements = mirror_formed_fixture()
        expected = {"o3-f45-d7c875318447": ([1.0] * 15 + [0.0, 0.0, 1.0], 0.49976768508632746),
                    "s2p-f11-613e6f656498": ([1.0, 1.0, 0.0, 0.0], 0.4999999574369773)}
        for name, requirement in requirements.items():
            with self.subTest(name=name):
                digest, model = DISCOVERY.placeholder_model(requirement, radius, {"MetalThickness": 0.1})
                self.assertEqual(digest, requirement["Hash"][:16])
                self.assertEqual(model["Name"], f"__preflight_placeholder_{digest}")
                self.assertEqual(model["Topology"], "SpatialEdgeCluster")
                self.assertEqual(model["Signature"], requirement["Signature"])
                weights, fraction = expected[name]
                self.assertEqual([edge["Weight"] for edge in model["Edges"]], weights)
                self.assertEqual(model["MirrorFormed"]["RealLengthFraction"], fraction)
                # The generator's functions: the Edges are model_edges with the Weight of their
                # portion (portions_from_signature order = EdgePortions), the record the ENTRY stamp.
                contract = device_coupons.validate_mirror_formed_contract(requirement)
                portions = [p["Portion"] for p in
                            cluster_signature_geometry.portions_from_signature(requirement["Signature"], radius)]
                record, entry_weights = device_coupons.mirror_formed_entry(contract, portions, radius)
                self.assertEqual(model["MirrorFormed"], record)
                self.assertEqual(model["MirrorFormed"]["EdgePortions"], portions)
                self.assertEqual(model["MirrorFormed"]["RealPortions"], requirement["MirrorFormedContract"]["RealPortions"])
                edges = cluster_signature_geometry.model_edges(requirement, radius)
                self.assertEqual(len(edges), len(model["Edges"]))
                for edge, weight, portion, model_edge in zip(edges, entry_weights, portions, model["Edges"]):
                    self.assertEqual({**edge, "Weight": weight}, model_edge)
                    self.assertEqual(weight, 1.0 if portion in contract["RealPortions"] else 0.0)
                # The interface mappings the Edges' InterfaceSlots bind (every slot used), with
                # the dummy surface matrices the discovery declares beside Interfaces.
                self.assertEqual(model["Interfaces"], DISCOVERY.unique_interfaces(requirement))
                self.assertEqual({edge["InterfaceSlot"] for edge in model["Edges"]},
                                 {entry["Slot"] for entry in model["Interfaces"]})
                self.assertEqual(model["FabricatedSurfaceMatrix"], "__preflight_dummy_fabricated_surface.csv")
                self.assertEqual(model["ThinSurfaceMatrix"], "__preflight_dummy_thin_surface.csv")
                # Nothing else differs from the signature-only placeholder.
                plain = self.plain_placeholder(requirement, radius)
                stamped_keys = {"Edges", "MirrorFormed", "Interfaces", "FabricatedSurfaceMatrix", "ThinSurfaceMatrix"}
                self.assertEqual({k: v for k, v in model.items() if k not in stamped_keys}, plain)

    def test_contract_placeholder_equals_stamp_signature_model(self):
        """One spelling: the B5 stamp (device_coupons.stamp_signature_model on the generator's
        one-model library) writes the same Edges and MirrorFormed record."""
        radius, requirements = mirror_formed_fixture()
        with tempfile.TemporaryDirectory(prefix="discover-stamp-") as tmp:
            for name, requirement in requirements.items():
                with self.subTest(name=name):
                    _, model = DISCOVERY.placeholder_model(requirement, radius)
                    library = Path(tmp) / f"{name}-process-library.json"
                    library.write_text(json.dumps({"Models": [{"Name": "generated", "Topology": "SpatialEdgeCluster"}]}))
                    coupon = {"Id": name, **requirement}
                    contract = device_coupons.validate_mirror_formed_contract(coupon)
                    self.assertTrue(device_coupons.stamp_signature_model(library, coupon, radius, mirror_formed=contract))
                    stamped = json.loads(library.read_text())["Models"][0]
                    self.assertEqual(stamped["Edges"], model["Edges"])
                    self.assertEqual(stamped["MirrorFormed"], model["MirrorFormed"])
                    self.assertEqual(stamped["Signature"], model["Signature"])

    def test_requirement_without_contract_keeps_the_signature_only_placeholder(self):
        """Today's placeholder, byte-identical: no contract (the ordinary key), the flag without
        a contract (a REAL key touched by a plane, decision 562), and a CurvedEdge carrying a
        contract (the 2D family's consumer: no SpatialEdgeCluster Edges exist for it)."""
        radius, requirements = mirror_formed_fixture()
        for name, requirement in requirements.items():
            plain = json.dumps(self.plain_placeholder(requirement, radius), sort_keys=True)
            unflagged = {k: v for k, v in requirement.items() if k not in ("MirrorFormed", "MirrorFormedContract")}
            flag_only = {k: v for k, v in requirement.items() if k != "MirrorFormedContract"}
            for label, record in (("unflagged", unflagged), ("flag-only", flag_only)):
                with self.subTest(name=name, case=label):
                    digest, model = DISCOVERY.placeholder_model(record, radius, {"MetalThickness": 0.1})
                    self.assertEqual(digest, requirement["Hash"][:16])
                    self.assertEqual(json.dumps(model, sort_keys=True), plain)
                    for key in ("Edges", "MirrorFormed", "Interfaces", "FabricatedSurfaceMatrix"):
                        self.assertNotIn(key, model)
        curved = {"Topology": "CurvedEdge", "Hash": "c" * 64, "MirrorFormed": True,
                  "MirrorFormedContract": {"Version": 1, "RealPortions": [], "Planes": [0]},
                  "BoundaryCondition": {"Type": "PEC"}, "Interfaces": [{"Slot": 0, "Target": 1, "Type": "SA"}],
                  "Signature": {"Type": "CurvedEdge", "RadiusOverR": 2.0, "Convexity": 1, "Interfaces": ["SA"],
                                "Law": "{\"Type\":\"PEC\"}"}}
        _, model = DISCOVERY.placeholder_model(curved, radius)
        self.assertNotIn("Edges", model)
        self.assertNotIn("MirrorFormed", model)
        self.assertEqual(model["Interfaces"], DISCOVERY.unique_interfaces(curved))

    def test_inconsistent_contract_fails_closed_by_name(self):
        radius, requirements = mirror_formed_fixture()
        base = requirements["s2p-f11-613e6f656498"]

        def mutated(**changes):
            record = json.loads(json.dumps(base))
            for key, value in changes.items():
                target = record["MirrorFormedContract"] if key.startswith("contract.") else record
                field = key.split(".", 1)[1] if key.startswith("contract.") else key
                if value is None:
                    target.pop(field)
                else:
                    target[field] = value
            return record

        cases = [
            ("names every portion", mutated(**{"contract.RealPortions": [0, 1, 2, 3]})),
            ("must be a non-empty sorted list", mutated(**{"contract.RealPortions": []})),
            ("index outside", mutated(**{"contract.RealPortions": [0, 7]})),
            ("RealLengthOverR .* disagrees", mutated(**{"contract.RealLengthOverR": base["MirrorFormedContract"]["RealLengthOverR"] + 0.1})),
            ("Frame.Chirality", mutated(**{"contract.Frame": {k: v for k, v in base["MirrorFormedContract"]["Frame"].items()
                                                              if k != "Chirality"}})),
            ("Version", mutated(**{"contract.Version": 2})),
            ("Hash .* is not sha256", mutated(Hash="0" * 64)),
            ("needs MirrorFormed: true", mutated(MirrorFormed=None)),
            ("is not an object", mutated(MirrorFormedContract="refused")),
        ]
        for text, record in cases:
            with self.subTest(text=text):
                with self.assertRaisesRegex(device_coupons.DeviceAdapterError, text):
                    DISCOVERY.placeholder_model(record, radius)

    def fake_device(self, tmp, requirements):
        """A device of the given Missing requirements, a seed library and an executable fake Palace
        applying the M7 operator's mirror-formed rule: a contract-bearing requirement is placed
        (Exact) only by a model carrying its Signature AND Edges whose Weights follow the ENTRY
        record (1.0 on the RealPortions' edges, 0.0 elsewhere); a signature-only placeholder
        leaves it Missing with the operator's Note (the pass-02 stall of decision 599)."""
        requirements_path = tmp / "fake-requirements.json"
        requirements_path.write_text(json.dumps(requirements))
        library = tmp / "process-library.json"
        library.write_text(json.dumps({"Version": 3, "Name": "fake-seed", "MatchingRadius": 1.9, "Models": []}))
        config = tmp / "device.json"
        config.write_text(json.dumps({
            "Problem": {"Type": "Electrostatic", "Output": str(tmp / "postpro")},
            "Solver": {"Electrostatic": {"ResponseCorrection": {"Library": str(library), "UnmatchedPolicy": "Warn"}}}}))
        palace = tmp / "fake-palace.py"
        palace.write_text(textwrap.dedent(f"""\
            #!{sys.executable}
            import json, pathlib, sys
            assert sys.argv[1] == "--surface-response-preflight"
            config = json.load(open(sys.argv[2]))
            response = config["Solver"]["Electrostatic"]["ResponseCorrection"]
            library = json.load(open(response["Library"]))
            models = {{model["Name"]: model for model in library["Models"]}}
            requirements = json.load(open({str(requirements_path)!r}))
            counts = {{"Exact": 0, "Interpolated": 0, "Missing": 0}}
            lengths = {{"Exact": 0.0, "Interpolated": 0.0, "Missing": 0.0}}
            for requirement in requirements:
                model = models.get("__preflight_placeholder_" + requirement["Hash"][:16])
                requirement["Status"] = "Missing"
                if model is not None and "MirrorFormedContract" in requirement:
                    if not model.get("Signature") or not model.get("Edges"):
                        requirement["Notes"] = [{SIGNATURE_ONLY_MODEL_NOTE!r}]
                    else:
                        entry = model["MirrorFormed"]
                        assert len(entry["EdgePortions"]) == len(model["Edges"])
                        for edge, portion in zip(model["Edges"], entry["EdgePortions"]):
                            assert edge["Weight"] == (1.0 if portion in entry["RealPortions"] else 0.0), "mis-stamped"
                        requirement["Status"] = "Exact"
                elif model is not None:
                    requirement["Status"] = "Exact"
                if requirement["Status"] == "Exact":
                    requirement["SelectedModels"] = [{{"Name": model["Name"], "Weight": 1.0}}]
                counts[requirement["Status"]] += requirement["Count"]
                lengths[requirement["Status"]] += requirement["TotalEdgeLength"]
            output = pathlib.Path(config["Problem"]["Output"]); output.mkdir(parents=True, exist_ok=True)
            manifest = {{"Version": 2, "Complete": counts["Missing"] == 0, "Requirements": requirements,
                        "Library": {{"Name": library["Name"], "Path": response["Library"],
                                    "MatchingRadius": library["MatchingRadius"]}},
                        "Summary": {{"Counts": counts, "TotalEdgeLengths": lengths}}}}
            (output / "surface-response-requirements.json").write_text(json.dumps(manifest))
            """))
        palace.chmod(palace.stat().st_mode | stat.S_IXUSR)
        return config, palace

    def test_contract_bearing_requirements_close(self):
        """The closure of decision 599's stall: the two production contract keys (plus an ordinary
        key) converge in two passes with Signature-keyed placeholders, the final manifest restores
        them Missing with their contract for the planner. FAILS on main ("Geometry closure
        stalled with missing requirements" at pass-02)."""
        radius, fixture = mirror_formed_fixture()
        ordinary = {"Topology": "IsolatedEdge", "Hash": "a" * 64, "Count": 4, "Instances": 2, "TotalEdgeLength": 10.0,
                    "BoundaryCondition": {"Type": "PEC"}, "Interfaces": [{"Slot": 0, "Target": 1, "Type": "SA"}],
                    "Signature": {"Type": "IsolatedEdge", "Interfaces": ["SA"], "Law": "{\"Type\":\"PEC\"}"}}
        requirements = [fixture["o3-f45-d7c875318447"], fixture["s2p-f11-613e6f656498"], ordinary]
        with tempfile.TemporaryDirectory(prefix="discover-contract-") as tmp:
            tmp = Path(tmp)
            config, palace = self.fake_device(tmp, requirements)
            manifest = DISCOVERY.discover(config, tmp / "closure", palace)
            history = json.loads((tmp / "closure" / "closure-history.json").read_text())
            self.assertEqual(len(history["Passes"]), 2)
            self.assertEqual(history["PlaceholderCount"], 3)
            self.assertEqual(history["Passes"][1]["AddedPlaceholders"], [])
            self.assertEqual(history["Passes"][1]["ProductionSummary"]["Counts"]["Missing"], 0)
            self.assertFalse(manifest["Complete"])
            self.assertEqual(manifest["Summary"]["Counts"]["Missing"], 6)
            by_hash = {requirement["Hash"]: requirement for requirement in manifest["Requirements"]}
            for requirement in requirements[:2]:
                restored = by_hash[requirement["Hash"]]
                self.assertEqual(restored["Status"], "Missing")
                self.assertEqual(restored["Reason"], "Missing from source library after exhaustive geometry discovery")
                self.assertEqual(restored["MirrorFormedContract"], requirement["MirrorFormedContract"])
                self.assertIs(restored["MirrorFormed"], True)
                self.assertEqual(restored["Signature"], requirement["Signature"])
                self.assertNotIn("SelectedModels", restored)
            library = json.loads((tmp / "closure" / "pass-02" / "process-library.json").read_text())
            placeholders = {model["Name"]: model for model in library["Models"]}
            for requirement in requirements[:2]:
                model = placeholders[f"__preflight_placeholder_{requirement['Hash'][:16]}"]
                self.assertEqual(model["MirrorFormed"]["RealPortions"], requirement["MirrorFormedContract"]["RealPortions"])
                self.assertEqual([edge["Weight"] for edge in model["Edges"]],
                                 [1.0 if p in requirement["MirrorFormedContract"]["RealPortions"] else 0.0
                                  for p in model["MirrorFormed"]["EdgePortions"]])
            self.assertNotIn("Edges", placeholders["__preflight_placeholder_" + "a" * 16])

    def test_signature_only_placeholder_stalls_the_fake_operator(self):
        """The fake Palace reproduces the stall of record on a signature-only placeholder (the
        rule the closure test depends on)."""
        radius, fixture = mirror_formed_fixture()
        requirement = fixture["s2p-f11-613e6f656498"]
        with tempfile.TemporaryDirectory(prefix="discover-stall-") as tmp:
            tmp = Path(tmp)
            config, palace = self.fake_device(tmp, [requirement])
            library = tmp / "signature-only.json"
            plain = self.plain_placeholder(requirement, radius)
            library.write_text(json.dumps({"Version": 3, "Name": "signature-only", "MatchingRadius": radius, "Models": [plain]}))
            pass_config = json.loads(config.read_text())
            pass_config["Solver"]["Electrostatic"]["ResponseCorrection"]["Library"] = str(library)
            pass_config["Problem"]["Output"] = str(tmp / "stall")
            (tmp / "stall.json").write_text(json.dumps(pass_config))
            DISCOVERY.run_preflight(palace, tmp / "stall.json", tmp / "stall.log")
            manifest = json.loads((tmp / "stall" / "surface-response-requirements.json").read_text())
            self.assertEqual(manifest["Requirements"][0]["Status"], "Missing")
            self.assertEqual(manifest["Requirements"][0]["Notes"], [SIGNATURE_ONLY_MODEL_NOTE])


if __name__ == "__main__":
    unittest.main()
