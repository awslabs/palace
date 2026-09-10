#!/usr/bin/env python3
# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""Small p5 multi-conductor covariance check; does not certify mesh convergence."""
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys

from compare_coupon_probe import compare

HERE = Path(__file__).resolve().parent


def run(args):
    args.output.mkdir(parents=True, exist_ok=False)
    env = dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1")
    # Do not inherit response-worker/reducer switches into an ordinary solve.
    env = {k: v for k, v in env.items() if not k.startswith("PALACE_RESPONSE_")}
    def bounded(command, tag, seconds=60):
        subprocess.run([sys.executable, str(HERE/"run_bounded_mesher.py"),
                        "--seconds", str(seconds), "--memory-gib", "6",
                        "--log", str(args.output/f"{tag}.log"), "--", *command],
                       check=True, env=env, timeout=seconds+15)
    transformed = args.output/"rotated.msh"
    bounded([args.julia, "--startup-file=no", f"--project={args.julia_project}",
             str(HERE/"transform_coupon_mesh.jl"), str(args.mesh), str(transformed), "0.63"],
            "transform", 20)
    with Path(str(args.mesh)+".interface-partition.csv").open() as stream:
        attributes = [int(row["attribute"]) for row in csv.DictReader(stream)]
    conductors = sorted({a % 100 for a in attributes if a//1000 in (5, 6)})
    if len(conductors) < 2:
        raise ValueError("This covariance test requires at least two fabricated conductors")
    dielectrics = [dict(Index=i, Attributes=[a for a in attributes if a//1000==family],
                       Type=name, Thickness=.002, Permittivity=permittivity, LossTan=.001)
                   for i, (family, name, permittivity) in enumerate([(6,"MA",10.),(5,"MS",11.47),(3,"SA",4.)],1)]
    for tag, mesh in (("original",args.mesh),("rotated",transformed)):
        config = {
            "Problem": {"Type":"Electrostatic","Verbose":1,"Output":str(args.output/tag),
                        "OutputFormats":{"Paraview":False,"GridFunction":False}},
            "Model": {"Mesh":str(mesh),"L0":1e-6,"Refinement":{"MaxIts":0}},
            "Domains": {"Materials":[{"Attributes":[1],"Permittivity":11.47},
                                      {"Attributes":[2],"Permittivity":1.}]},
            "Boundaries": {"Ground":{"Attributes":[1]},
                           "Terminal":[{"Index":i,"Attributes":[a for a in attributes if a//1000 in (5,6) and a%100==c]}
                                       for i,c in enumerate(conductors,1)],
                           "Postprocessing":{"Dielectric":dielectrics}},
            "Solver": {"Order":5,"Electrostatic":{"Save":0,"ResponseMatrix":False,"AggregateResponseMatrix":False},
                       "Linear":{"Type":"BoomerAMG","KSPType":"CG","Tol":1e-8,"MaxIts":500,
                                 "EstimatorTol":.5,"EstimatorMaxIts":5,"EstimatorMG":True}},
        }
        path=args.output/f"{tag}.json"
        path.write_text(json.dumps(config,indent=2)+"\n")
        bounded(["mpirun","-n",str(args.ranks),str(args.palace),str(path)],tag)
    report=compare(args.output/"original",args.output/"rotated","fabricated",1e-6,1e-6)
    report["Purpose"]="Rigid-coordinate covariance on identical connectivity, not accuracy qualification"
    report["InputSHA256"]={"Mesh":hashlib.sha256(args.mesh.read_bytes()).hexdigest(),
                          "PalaceBinary":hashlib.sha256(args.palace.read_bytes()).hexdigest()}
    report["CapacitanceMatrixByteIdentical"]=(args.output/"original/terminal-C.csv").read_bytes()==(args.output/"rotated/terminal-C.csv").read_bytes()
    (args.output/"comparison.json").write_text(json.dumps(report,indent=2)+"\n")
    print(json.dumps(report,indent=2))
    return report["Passed"]


if __name__ == "__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mesh",type=Path,required=True)
    parser.add_argument("--output",type=Path,required=True)
    parser.add_argument("--palace",type=Path,required=True)
    parser.add_argument("--julia",default="julia")
    parser.add_argument("--julia-project",type=Path,required=True)
    parser.add_argument("--ranks",type=int,default=2,choices=range(1,7))
    raise SystemExit(0 if run(parser.parse_args()) else 1)
