import gmsh, sys
gmsh.initialize(); gmsh.option.setNumber("General.Terminal",0)
gmsh.model.add("m")
occ=gmsh.model.occ
L=5.0
box=occ.addBox(-L,-L,-L,2*L,2*L,2*L)
sq=occ.addRectangle(-L,-L,0,2*L,2*L)
disk=occ.addDisk(0,0,0,1.5,1.5)
out,_=occ.fragment([(3,box)],[(2,sq),(2,disk)])
occ.synchronize()
vols=[t for d,t in gmsh.model.getEntities(3)]
surfs=gmsh.model.getEntities(2)
film=[];hole=[];walls=[]
for d,t in surfs:
    com=occ.getCenterOfMass(d,t)
    bb=gmsh.model.getBoundingBox(d,t)
    if abs(bb[2])<1e-4 and abs(bb[5])<1e-4:
        # on z=0 plane
        if bb[3]-bb[0] < 3.1: hole.append(t)
        else: film.append(t)
    else: walls.append(t)
print(len(vols),film,hole,len(walls))
gmsh.model.addPhysicalGroup(3,vols,1)
gmsh.model.addPhysicalGroup(2,walls,2)
gmsh.model.addPhysicalGroup(2,film,8)
gmsh.model.addPhysicalGroup(2,hole,9)
gmsh.option.setNumber("Mesh.MeshSizeMax",0.8)
gmsh.option.setNumber("Mesh.MeshSizeMin",0.2)
gmsh.model.mesh.field.add("Distance",1); gmsh.model.mesh.field.setNumbers(1,"SurfacesList",film+hole)
gmsh.model.mesh.field.add("Threshold",2); gmsh.model.mesh.field.setNumber(2,"InField",1)
for k,v in dict(SizeMin=0.25,SizeMax=0.9,DistMin=0.2,DistMax=3).items(): gmsh.model.mesh.field.setNumber(2,k,v)
gmsh.model.mesh.field.setAsBackgroundMesh(2)
gmsh.model.mesh.generate(3)
gmsh.option.setNumber("Mesh.MshFileVersion",2.2)
gmsh.write("review-repro/pec.msh")
gmsh.finalize()
