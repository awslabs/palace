# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0
using Gmsh: gmsh
using SHA
using TOML
include(joinpath(@__DIR__,"frozen_volume_study.jl"))
length(ARGS)==1 || error("volume-study.toml")
study=TOML.parsefile(ARGS[1])
for variant in study["Variants"]
    gmsh.initialize()
    try
        gmsh.option.setNumber("General.Verbosity",0)
        gmsh.open(variant["Mesh"])
        actual=frozen_surface_fingerprint()
        actual==study["Boundary"] || error("Serialized boundary changed: $(variant["Name"])")
        bytes2hex(sha256(read(variant["Mesh"])))==variant["MeshSHA256"] || error("Candidate file changed")
        println(variant["Name"],": serialized boundary identical, ",actual["Triangles"]," triangles")
    finally
        gmsh.finalize()
    end
end
