#!/usr/bin/env python3

# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

import hashlib
import json
import subprocess
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPT = ROOT / "prepare_graded_library_meshes.py"


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


class PrepareGradedLibraryMeshesTest(unittest.TestCase):
    def test_freezes_retained_contract_and_emits_bounded_gmsh22_plan(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = root / "source"
            source.mkdir()
            (source / "mesh-signature.csv").write_text(
                "Index,Slot,Conductor,Px,Py,Pz,Gx,Gy,Gz,Tx,Ty,Tz,Nz,S0,S1,VertexArm\n"
                "1,0,1,0,0,0,1,0,0,0,1,0,1,-1,1,0\n"
            )
            (source / "plan-view-mask.csv").write_text(
                "Facet,Conductor,Plane,X,Y\n1,1,0,-1,-1\n1,1,0,1,-1\n1,1,0,1,1\n"
            )
            (source / "plan-view-boundary.csv").write_text(
                "Loop,Vertex,Conductor,Plane,Hole,Class,X,Y\n"
                "1,1,1,0,0,Physical,-1,-1\n1,2,1,0,0,Physical,1,-1\n1,3,1,0,0,Physical,1,1\n"
            )
            trace = source / "basis.csv"
            trace.write_text(
                "x,y,z,V,triangle\n-2,-2,-2,1,1\n2,-2,-2,0,1\n2,2,-2,0,1\n"
            )
            mesh = root / "retained.msh"
            mesh.write_bytes(b"retained mesh")
            cases = []
            inputs = {str(mesh): {"SHA256": sha(mesh), "Bytes": mesh.stat().st_size}}
            for kind in ("thin", "fabricated"):
                config_path = source / f"spatial_{kind}.json"
                config = {
                    "Problem": {"Output": str(root / kind)},
                    "Model": {"Mesh": str(mesh), "Refinement": {"MaxIts": 0}},
                    "Domains": {"Materials": []},
                    "Boundaries": {
                        "Ground": {"Attributes": [4001]},
                        "PrescribedPotential": [
                            {"Index": 7, "Attributes": [1], "DataFile": str(trace)},
                            {
                                "Index": 9,
                                "Attributes": [1],
                                "TerminalAttributes": [5001, 6001],
                                "DataFile": str(trace),
                            },
                        ],
                        "Postprocessing": {
                            "Dielectric": [
                                {"Attributes": [3000], "EdgeAttributes": [3100]}
                            ]
                        },
                    },
                    "Solver": {
                        "Order": 5,
                        "Linear": {"Tol": 1e-8, "MaxIts": 500},
                    },
                }
                config_path.write_text(json.dumps(config))
                cases.append(
                    {
                        "Key": f"09-{kind}",
                        "Model": "spatialedgecluster_edgecount-10_example",
                        "Kind": kind,
                        "SourceConfig": str(config_path),
                        "Mesh": str(mesh),
                    }
                )
            reference = root / "manifest.json"
            reference.write_text(json.dumps({"Cases": cases, "Inputs": inputs}))
            output = root / "campaign"
            tools = root / "tools"
            tools.mkdir()
            subprocess.run(
                [
                    "python3",
                    str(SCRIPT),
                    str(reference),
                    str(output),
                    "--tools",
                    str(tools),
                    "--julia-project",
                    str(root / "project"),
                    "--matching-trace", "none",
                ],
                check=True,
                capture_output=True,
                text=True,
            )
            plan = json.loads((output / "mesh-plan.json").read_text())
            self.assertEqual(plan["MeshFormat"], "Gmsh 2.2 binary")
            self.assertEqual(len(plan["Cases"]), 2)
            by_kind = {case["Kind"]: case for case in plan["Cases"]}
            self.assertEqual(by_kind["thin"]["FineSurfaceSizeMicrometres"], 0.002)
            self.assertEqual(by_kind["fabricated"]["FineSurfaceSizeMicrometres"], 0.0005)
            for case in plan["Cases"]:
                self.assertEqual(case["MinimumVolumeSizeMicrometres"], 0.002)
                self.assertEqual(case["SourceContract"]["SourceCount"], 2)
                self.assertEqual(case["SourceContract"]["Linear"]["Tol"], 1e-8)
                self.assertEqual(case["ExpectedInterfaceAttributes"], [3000, 3100, 4001, 5001, 6001])
                self.assertEqual(case["SourceContract"]["Sources"][0]["SHA256"], sha(trace))
            runner = (output / "generate-meshes.sh").read_text()
            self.assertIn("run_bounded_mesher.py", runner)
            self.assertIn("TET_GEOMETRY_ORDER=1", runner)
            self.assertIn("--memory-gib 110", runner)


if __name__ == "__main__":
    unittest.main()
