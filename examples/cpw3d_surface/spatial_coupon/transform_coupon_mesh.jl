# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

# Rigidly transform an existing mesh for numerical covariance tests. No remeshing
# or source-definition changes. This is not a mesh-convergence qualification.
using Gmsh: gmsh
using TOML
length(ARGS)==3 || error("input.msh output.msh angle_radians")
isfile(ARGS[2]) && error("Refuse to overwrite output mesh")
angle=parse(Float64,ARGS[3])
isfinite(angle) || error("Non-finite angle")
gmsh.initialize()
try
    gmsh.open(ARGS[1])
    elements=gmsh.model.mesh.getElements()
    groups=[(dim,tag,gmsh.model.getEntitiesForPhysicalGroup(dim,tag)) for (dim,tag) in gmsh.model.getPhysicalGroups()]
    c,s=cos(angle),sin(angle)
    gmsh.model.mesh.affineTransform([c,-s,0.,0.23, s,c,0.,-0.13, 0.,0.,1.,0.17, 0.,0.,0.,1.])
    elements==gmsh.model.mesh.getElements() || error("Connectivity changed")
    groups==[(dim,tag,gmsh.model.getEntitiesForPhysicalGroup(dim,tag)) for (dim,tag) in gmsh.model.getPhysicalGroups()] || error("Physical groups changed")
    gmsh.option.setNumber("Mesh.MshFileVersion",2.2)
    gmsh.option.setNumber("Mesh.Binary",1)
    gmsh.write(ARGS[2])
    open(ARGS[2]*".transform.toml","w") do f
        TOML.print(f,Dict("Input"=>abspath(ARGS[1]),"AngleRadians"=>angle,
                          "Translation"=>[0.23,-0.13,0.17],"ConnectivityPreserved"=>true))
    end
finally
    gmsh.finalize()
end
