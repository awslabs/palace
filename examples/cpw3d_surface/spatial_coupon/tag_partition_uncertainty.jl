# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

# Diagnostic-only copy: split each interface into certified and unresolved
# elements. No overlapping groups, no coordinate/connectivity changes. A solver
# must include BOTH parts in each original conductor boundary condition.
using Gmsh: gmsh
using DelimitedFiles
using SHA
using TOML
const PARTITION_AUDIT_OFFSET=10000
length(ARGS)==2 || error("input.msh diagnostic-output.msh")
input,output=abspath.(ARGS)
isfile(output) && error("Refuse to overwrite diagnostic mesh")
certificate_path=input*".interface-partition.csv.elements.csv"
seal=TOML.parsefile(input*".partition-certificate.toml")
seal["Version"]==1 || error("Unsupported ownership certificate version")
seal["MeshSHA256"]==bytes2hex(sha256(read(input))) || error("Certificate mesh hash mismatch")
seal["ElementCertificateSHA256"]==bytes2hex(sha256(read(certificate_path))) || error("Element certificate hash mismatch")
data,_=readdlm(certificate_path,',',header=true)
certificates=Dict{Tuple{Int,NTuple{3,UInt64}},Bool}()
for row in axes(data,1)
    attribute=Int(data[row,2]);nodes=Tuple(sort(UInt64.(data[row,4:6])))
    key=(attribute,nodes)
    haskey(certificates,key) && error("Duplicate boundary certificate")
    data[row,3] in (0,1) || error("Invalid certificate flag")
    certificates[key]=Bool(data[row,3])
end
gmsh.initialize()
try
    gmsh.open(input)
    before_nodes=gmsh.model.mesh.getNodes()
    before_volume=gmsh.model.mesh.getElements(3)
    groups=Dict{Int,Vector{Tuple{Int,UInt64,Vector{UInt64}}}}()
    old_entities=Set{Int32}();old_groups=Tuple{Int32,Int32}[]
    for (_,attribute) in gmsh.model.getPhysicalGroups(2)
        attribute==1 && continue
        3000<=attribute<7000 || error("Unexpected physical attribute")
        push!(old_groups,(Int32(2),Int32(attribute)))
        for entity in gmsh.model.getEntitiesForPhysicalGroup(2,attribute)
            entity in old_entities && error("Overlapping physical assignments")
            push!(old_entities,entity)
            types,tags,connect=gmsh.model.mesh.getElements(2,entity)
            for (type,etags,nodes) in zip(types,tags,connect)
                _,_,_,nnode,_,primary=gmsh.model.mesh.getElementProperties(type)
                primary==3 || error("Expected triangular interface")
                for (i,element) in enumerate(etags)
                    conn=collect(nodes[(i-1)*nnode+1:i*nnode])
                    key=(Int(attribute),Tuple(sort(conn[1:3])))
                    haskey(certificates,key) || error("Missing/mismatched boundary certificate")
                    certified=pop!(certificates,key)
                    target=Int(attribute)+(certified ? 0 : PARTITION_AUDIT_OFFSET)
                    push!(get!(groups,target,Tuple{Int,UInt64,Vector{UInt64}}[]),(Int(type),element,conn))
                end
            end
        end
    end
    isempty(certificates) || error("Unused boundary certificates")
    for entity in old_entities;gmsh.model.mesh.removeElements(2,entity);end
    gmsh.model.removePhysicalGroups(old_groups)
    for attribute in sort!(collect(keys(groups)))
        entity=gmsh.model.addDiscreteEntity(2)
        for type in unique(r[1] for r in groups[attribute])
            records=[r for r in groups[attribute] if r[1]==type]
            gmsh.model.mesh.addElementsByType(entity,type,[r[2] for r in records],reduce(vcat,[r[3] for r in records]))
        end
        gmsh.model.addPhysicalGroup(2,[entity],attribute,"partition_audit_$attribute")
    end
    before_nodes==gmsh.model.mesh.getNodes() || error("Node coordinates/IDs changed")
    before_volume==gmsh.model.mesh.getElements(3) || error("Volume elements changed")
    gmsh.option.setNumber("Mesh.MshFileVersion",2.2)
    gmsh.option.setNumber("Mesh.Binary",1)
    gmsh.write(output)
    open(output*".attributes.csv","w") do f
        println(f,"attribute,canonical_attribute,certified,elements")
        for attribute in sort!(collect(keys(groups)))
            println(f,"$attribute,$(mod(attribute,PARTITION_AUDIT_OFFSET)),$(Int(attribute<PARTITION_AUDIT_OFFSET)),$(length(groups[attribute]))")
        end
    end
    open(output*".audit.toml","w") do f
        TOML.print(f,Dict("DiagnosticOnly"=>true,"Input"=>input,
                          "InputSHA256"=>bytes2hex(sha256(read(input))),
                          "ConnectivityPreserved"=>true,"Offset"=>PARTITION_AUDIT_OFFSET))
    end
finally
    gmsh.finalize()
end
