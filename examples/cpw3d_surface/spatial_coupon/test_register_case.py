# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0
"""register_case turns a source directory plus an explicit footprint declaration into a
manifest case by content: hash reuse, changed-source versioning with a retired entry,
fail-closed footprint declaration and scope classification (decision 48)."""
import json
from pathlib import Path
import shutil
import tempfile
import unittest

import register_case
from register_case import (FOOTPRINT_BOUND, REGISTRATION_RECORD, STATUS_REGISTERED, STATUS_REUSED,
                           STATUS_UNSUPPORTED, RegistrationError, register, sha256)

HERE = Path(__file__).resolve().parent
PRODUCTION_MANIFEST = HERE / "geometry-independence-suite.json"
# A repository source directory whose contract is derived from its inputs and a census
# (its process library binds a trace basis): the stub probe reproduces its census labels.
SOURCE_FIXTURE = HERE / "testdata" / "two-edge-8dd4bc70f183"


def census_probe(labels, summary_status="built", guard=None):
    """A probe standing in for the labels-only pass: writes the label census the real
    pass would record (from the frozen contract) and the build summary."""
    def probe(probe_manifest, case_id, root, log):
        root.mkdir(parents=True)
        Path(log).write_text("stub probe\n")
        (root / "build-census.json").write_text(json.dumps(
            {"LabelsOnly": True, "InterfaceAreas": [{"Attribute": label, "Name": f"surface_{label}"} for label in labels]}))
        probe.manifests.append(json.loads(Path(probe_manifest).read_text()))
        return {"Case": case_id, "Commit": "stub", "Root": str(root), "Status": summary_status,
                "Stage": "gmsh-build" if guard else "labels-only", "ReturnCode": 0 if guard is None else 1,
                "ScopeGuard": guard, "Message": None if guard is None else f"unsupported class {guard}"}
    probe.manifests = []
    return probe


class RegisterCaseTest(unittest.TestCase):
    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp())
        self.addCleanup(shutil.rmtree, self.tmp, True)
        manifest = json.loads(PRODUCTION_MANIFEST.read_text())
        manifest["RepositoryRoot"] = str((HERE / manifest["RepositoryRoot"]).resolve())
        self.manifest_path = self.tmp / "manifest.json"
        self.manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
        self.source = self.tmp / "source"
        shutil.copytree(SOURCE_FIXTURE, self.source)
        (self.source / "semantic-contract.json").unlink()   # the contract is derived, never copied
        frozen = json.loads((SOURCE_FIXTURE / "semantic-contract.json").read_text())
        self.labels = [item["Attribute"] for item in frozen["BoundaryLabels"]]
        self.frozen = frozen
        self.recipe = "examples/cpw3d_surface/spatial_coupon/testdata/generality-mesh-recipe.json"

    def manifest(self):
        return json.loads(self.manifest_path.read_text())

    def register(self, case_id="registered-copy", work=None, **overrides):
        options = dict(footprint="producer-default", inventory_status="RepositoryAssessmentFixture",
                       manifest_path=self.manifest_path, mesh_recipe=self.recipe,
                       work=work or self.tmp / f"work-{len(list(self.tmp.iterdir()))}",
                       probe=census_probe(self.labels), refreeze_calibration=[])
        options.update(overrides)
        return register(case_id, self.source, **options)

    def test_registers_derives_the_contract_and_reuses_by_content(self):
        before = self.manifest()
        record = self.register()
        self.assertEqual(record["Status"], STATUS_REGISTERED)
        self.assertEqual(record["FixtureVersion"], 1)
        manifest = self.manifest()
        self.assertEqual(len(manifest["Cases"]), len(before["Cases"]) + 1)
        case = manifest["Cases"][-1]
        self.assertEqual(case["Id"], "registered-copy")
        self.assertEqual(case["InventoryStatus"], "RepositoryAssessmentFixture")
        self.assertEqual(case["Source"]["Directory"], str(self.source.resolve()))
        self.assertEqual(case["Source"]["EtchFootprint"], "producer-default")
        self.assertEqual(case["Variants"], before["Cases"][0]["Variants"])
        self.assertEqual(case["TransformComparison"], before["Cases"][0]["TransformComparison"])
        self.assertEqual(case["Source"]["Files"]["MeshRecipe"]["RepositoryPath"], self.recipe)
        # Every source file is hashed; the contract is the two-pass derivation (provisional
        # contract bound by the staging probe, final contract from the probe census).
        for role, entry in case["Source"]["Files"].items():
            path = (self.source / entry["Name"] if "Name" in entry
                    else Path(manifest["RepositoryRoot"]) / entry["RepositoryPath"])
            self.assertEqual(sha256(path), entry["SHA256"], role)
        self.assertEqual(set(case["Source"]["Files"]),
                         {"Signature", "Boundary", "Mask", "Process", "ProcessLibrary", "BasisContract",
                          "TraceVertices", "TraceTriangles", "Provenance", "SemanticContract", "MeshRecipe"})
        contract = json.loads((self.source / "semantic-contract.json").read_text())
        self.assertNotIn("Provisional", contract["Derivation"])
        self.assertEqual(contract["Derivation"]["BuildCensusInterfaceLabels"], sorted(self.labels))
        for key in ("BoundaryLabels", "SemanticCorners", "FeatureTopology", "ProtectedSupports"):
            self.assertEqual(contract[key], self.frozen[key], key)
        probe_manifest = json.loads((Path(record["Work"]) / "probe-manifest.json").read_text())
        probe_case = probe_manifest["Cases"][-1]
        self.assertTrue(probe_case["Source"]["Directory"].endswith("probe-source"))
        provisional = json.loads((Path(record["Work"]) / "probe-source" / "semantic-contract.json").read_text())
        self.assertIn("Provisional", provisional["Derivation"])
        self.assertEqual(record["Scope"]["UnsupportedClasses"], [])
        self.assertIn("TraceBasis", record["Scope"]["ExhibitedClasses"])
        self.assertTrue((Path(record["Work"]) / REGISTRATION_RECORD).is_file())
        # Same content again: reused, manifest byte-identical, no probe.
        written = self.manifest_path.read_text()
        again = self.register(probe=census_probe([]))
        self.assertEqual(again["Status"], STATUS_REUSED)
        self.assertEqual(again["FixtureVersion"], 1)
        self.assertEqual(self.manifest_path.read_text(), written)

    def test_changed_source_is_a_new_version_with_the_old_binding_retired(self):
        first = self.register()
        old_case = next(c for c in self.manifest()["Cases"] if c["Id"] == "registered-copy")
        mask = self.source / "plan-view-mask.csv"
        mask.write_text(mask.read_text() + "\n")
        record = self.register()
        self.assertEqual(record["Status"], STATUS_REGISTERED)
        self.assertEqual((record["FixtureVersion"], record["RetiredVersion"]), (2, 1))
        manifest = self.manifest()
        self.assertEqual([c["Id"] for c in manifest["Cases"]].count("registered-copy"), 1)
        case = next(c for c in manifest["Cases"] if c["Id"] == "registered-copy")
        self.assertEqual(case["FixtureVersion"], 2)
        self.assertEqual(case["Source"]["Files"]["Mask"]["SHA256"], sha256(mask))
        self.assertNotEqual(case["Source"]["Files"]["Mask"]["SHA256"], old_case["Source"]["Files"]["Mask"]["SHA256"])
        retired = manifest["RetiredFixtures"]["Entries"][-1]
        self.assertEqual((retired["Id"], retired["RetiredVersion"]), ("registered-copy", 1))
        self.assertEqual(set(retired["RetiredBinding"]), {"Mask", "SemanticContract"})
        self.assertEqual(retired["RetiredBinding"]["Mask"], old_case["Source"]["Files"]["Mask"])
        self.assertIn("['Mask']", retired["Reason"])
        self.assertIn(str(Path(record["Work"]) / REGISTRATION_RECORD), retired["Evidence"])
        self.assertIn("version 1 retired", case["Provenance"])
        self.assertEqual(first["FixtureVersion"], 1)

    def test_footprint_declaration_is_mandatory_and_must_match_the_directory(self):
        with self.assertRaises(RegistrationError) as context:
            self.register(footprint=None)
        self.assertIn("footprint must be declared", str(context.exception))
        with self.assertRaises(RegistrationError):
            self.register(footprint=FOOTPRINT_BOUND)     # no retained-etch.csv in the directory
        (self.source / "retained-etch.csv").write_text("X,Y\n0,0\n")
        with self.assertRaises(RegistrationError):
            self.register(footprint="producer-default")  # a file exists but is not declared
        self.assertNotIn("registered-copy", [c["Id"] for c in self.manifest()["Cases"]])
        record = json.loads((self.tmp / "work-3" / REGISTRATION_RECORD).read_text())
        self.assertEqual(record["Status"], "failed")
        bound = self.register(footprint=FOOTPRINT_BOUND, work=self.tmp / "work-bound")
        self.assertEqual(bound["Status"], STATUS_REGISTERED)
        case = next(c for c in self.manifest()["Cases"] if c["Id"] == "registered-copy")
        self.assertNotIn("EtchFootprint", case["Source"])
        self.assertEqual(case["Source"]["Files"]["RetainedEtch"]["Name"], "retained-etch.csv")
        self.assertIn("DeviceFootprint", bound["Scope"]["ExhibitedClasses"])

    def test_mesh_recipe_and_inventory_status_fail_closed(self):
        with self.assertRaises(RegistrationError):
            self.register(mesh_recipe=None)               # neither mesh-recipe.json nor --mesh-recipe
        with self.assertRaises(RegistrationError):
            self.register(inventory_status="Calibration")
        shutil.copyfile(HERE / "testdata" / "generality-mesh-recipe.json", self.source / "mesh-recipe.json")
        with self.assertRaises(RegistrationError):
            self.register()                               # both bindings present
        record = self.register(mesh_recipe=None)
        self.assertEqual(record["Status"], STATUS_REGISTERED)
        case = next(c for c in self.manifest()["Cases"] if c["Id"] == "registered-copy")
        self.assertEqual(case["Source"]["Files"]["MeshRecipe"]["Name"], "mesh-recipe.json")

    @unittest.skipIf(shutil.which("julia") is None, "julia is required for the labels-only probe")
    def test_labels_only_probe_reproduces_the_frozen_label_set_before_any_mesh(self):
        """Decision 62(2): the real probe (run_gmsh_only_case.py --labels-only) runs the
        production mesher command up to the CAD-entity labelling and writes the label set;
        the derived contract equals the frozen one (labels, corners, topology, supports);
        no mesh, publication or build record is produced."""
        record = self.register(probe=None)
        self.assertEqual(record["Status"], STATUS_REGISTERED)
        probe = Path(record["ProbeRoot"])
        census = json.loads((probe / "build-census.json").read_text())
        self.assertTrue(census["LabelsOnly"])
        self.assertEqual([row["Attribute"] for row in census["InterfaceAreas"]], sorted(self.labels))
        self.assertTrue(all(row["CADSurfaces"] >= 1 and "Area" not in row for row in census["InterfaceAreas"]))
        self.assertIn("TraceBasis", census["Scope"]["ExhibitedClasses"])
        self.assertEqual(census["SemanticContractSHA256"],
                         sha256(Path(record["Work"]) / "probe" / "canonical-semantic.json"))
        for name in ("gmsh-build.msh", "canonical.msh", "canonical-build.json", "canonical-publish.log"):
            self.assertFalse((probe / name).exists(), name)
        summary = json.loads((probe / "build-summary.json").read_text())
        self.assertEqual((summary["Status"], summary["Stage"]), ("built", "labels-only"))
        report = json.loads((probe / "gmsh-build.log.json").read_text())
        self.assertEqual(report["ReturnCode"], 0)
        self.assertNotIn("Stage", report)   # outside the frozen gmsh-build stage contract: no mesh artifact
        self.assertIn("--labels-only", report["Command"])
        self.assertIn("--prism-tubes", report["Command"])
        contract = json.loads((self.source / "semantic-contract.json").read_text())
        for key in ("BoundaryLabels", "SemanticCorners", "FeatureTopology", "ProtectedSupports"):
            self.assertEqual(contract[key], self.frozen[key], key)
        self.assertEqual(contract["Derivation"]["BuildCensusInterfaceLabels"], sorted(self.labels))
        self.assertIn("labels-only", self.manifest()["Cases"][-1]["Provenance"])

    def test_prepare_then_commit_equals_register_and_commits_serially(self):
        """Decision 62(2): prepare_registration (probe + derivation, manifest read only)
        for several cases, then commit_registration in order - the same manifest entries
        as register(); a prepared reuse / stop is final at commit."""
        before = self.manifest_path.read_text()
        second = self.tmp / "source-2"
        shutil.copytree(self.source, second)
        provenance = json.loads((second / "provenance.json").read_text())
        provenance["Copy"] = "second"
        (second / "provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")
        prepared = [register_case.prepare_registration(case_id, directory, footprint="producer-default",
                                                       inventory_status="RepositoryAssessmentFixture",
                                                       manifest_path=self.manifest_path, mesh_recipe=self.recipe,
                                                       work=self.tmp / f"prepare-{case_id}", probe=census_probe(self.labels))
                    for case_id, directory in (("first-copy", self.source), ("second-copy", second))]
        self.assertEqual(self.manifest_path.read_text(), before)   # nothing written by the preparation
        self.assertTrue(all(not item.done for item in prepared))
        self.assertFalse((self.source / "semantic-contract.json").exists())
        records = [register_case.commit_registration(item, manifest_path=self.manifest_path) for item in prepared]
        self.assertEqual([record["Status"] for record in records], [STATUS_REGISTERED] * 2)
        self.assertEqual([(case["Id"], case["FixtureVersion"]) for case in self.manifest()["Cases"][-2:]],
                         [("first-copy", 1), ("second-copy", 1)])
        reused = register_case.prepare_registration("first-copy", self.source, footprint="producer-default",
                                                    inventory_status="RepositoryAssessmentFixture", manifest_path=self.manifest_path,
                                                    mesh_recipe=self.recipe, work=self.tmp / "prepare-again", probe=census_probe([]))
        self.assertTrue(reused.done)
        self.assertEqual(register_case.commit_registration(reused, manifest_path=self.manifest_path)["Status"], STATUS_REUSED)
        stopped = register_case.prepare_registration(
            "third-copy", self.source, footprint="producer-default", inventory_status="RepositoryAssessmentFixture",
            manifest_path=self.manifest_path, mesh_recipe=self.recipe, work=self.tmp / "prepare-guard",
            probe=census_probe([], summary_status="unsupported-class", guard="NarrowHoles"))
        self.assertTrue(stopped.done)
        self.assertEqual(stopped.record["Status"], STATUS_UNSUPPORTED)
        with self.assertRaises(RegistrationError):
            register_case.commit_registration(stopped, manifest_path=self.manifest_path)

    def test_unsupported_class_is_recorded_distinctly_and_never_registered(self):
        # From the inputs: a rounded process is guarded before any probe.
        process = self.source / "process.toml"
        process.write_text(process.read_text().replace("TopRounding = 0.0", "TopRounding = 0.005"))
        with self.assertRaises(RegistrationError) as context:
            self.register(probe=census_probe([]))
        self.assertIn("unsupported class TopRounding", str(context.exception))
        record = json.loads((self.tmp / "work-2" / REGISTRATION_RECORD).read_text())
        self.assertEqual(record["Status"], STATUS_UNSUPPORTED)
        self.assertEqual(record["StoppedBy"], {"Kind": "ScopeGuard", "Id": "TopRounding", "Stage": "inputs"})
        # From the build: a mesher ScopeGuard stop of the probe.
        process.write_text(process.read_text().replace("TopRounding = 0.005", "TopRounding = 0.0"))
        with self.assertRaises(RegistrationError):
            self.register(probe=census_probe([], summary_status="unsupported-class", guard="NarrowHoles"),
                          work=self.tmp / "work-guard")
        record = json.loads((self.tmp / "work-guard" / REGISTRATION_RECORD).read_text())
        self.assertEqual(record["Status"], STATUS_UNSUPPORTED)
        self.assertEqual(record["StoppedBy"]["Id"], "NarrowHoles")
        self.assertNotIn("registered-copy", [c["Id"] for c in self.manifest()["Cases"]])
        self.assertFalse((self.source / "semantic-contract.json").exists())


if __name__ == "__main__":
    unittest.main()
