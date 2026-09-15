#!/usr/bin/env python3
# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0
"""Topology, planar-support and directional-width checks for native MMG scouts."""
import argparse
import json
from pathlib import Path
import meshio
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from edge_volume_metric import segment_distances
from mesh_array_io import read_mesh


def blocks(mesh,kind):
    key='gmsh:physical' if 'gmsh:physical' in mesh.cell_data else 'medit:ref'
    return (np.concatenate([c.data for c in mesh.cells if c.type==kind]),
            np.concatenate([a for c,a in zip(mesh.cells,mesh.cell_data[key]) if c.type==kind]))


def boundary_component_counts(triangles,labels):
    """Connected metal-surface components, without joining through dielectric faces."""
    metal=triangles[np.isin(labels,[5001,6001])]
    if not len(metal):return 0
    edges=np.sort(metal[:,[(0,1),(1,2),(2,0)]].reshape(-1,2),axis=1)
    owners=np.repeat(np.arange(len(metal)),3)
    order=np.lexsort((edges[:,1],edges[:,0]));edges=edges[order];owners=owners[order]
    first=np.r_[0,np.flatnonzero(np.any(edges[1:]!=edges[:-1],axis=1))+1]
    count=np.diff(np.r_[first,len(edges)])
    if np.any(count>2):raise ValueError('Nonmanifold conductor surface')
    shared=first[count==2]
    graph=coo_matrix((np.ones(len(shared)),(owners[shared],owners[shared+1])),shape=(len(metal),len(metal))).tocsr()
    return int(connected_components(graph,directed=False,return_labels=False))


def analyze(mesh):
    t,material=blocks(mesh,'tetra');b,labels=blocks(mesh,'triangle');p=mesh.points
    xyz=p[t];signed=np.einsum('ij,ij->i',np.cross(xyz[:,1]-xyz[:,0],xyz[:,2]-xyz[:,0]),xyz[:,3]-xyz[:,0])/6
    if np.any(signed<=0) or not np.all(np.isfinite(p)):raise ValueError('Invalid tetrahedra/coordinates')
    volume={int(a):float(signed[material==a].sum()) for a in np.unique(material)}
    face=np.sort(t[:,[(0,1,2),(0,1,3),(0,2,3),(1,2,3)]].reshape(-1,3),axis=1)
    mat=np.repeat(material,4);owner=np.repeat(np.arange(len(t)),4)
    order=np.lexsort((face[:,2],face[:,1],face[:,0]));face=face[order];mat=mat[order];owner=owner[order]
    start=np.r_[0,np.flatnonzero(np.any(face[1:]!=face[:-1],axis=1))+1]
    count=np.diff(np.r_[start,len(face)])
    if np.any(count>2):raise ValueError('Nonmanifold volume face')
    low=np.minimum.reduceat(mat,start);high=np.maximum.reduceat(mat,start)
    unique=face[start];boundary=np.sort(b,axis=1)
    dtype=np.dtype([('a',unique.dtype),('b',unique.dtype),('c',unique.dtype)])
    keys=np.ascontiguousarray(unique).view(dtype).ravel();bk=np.ascontiguousarray(boundary).view(dtype).ravel()
    ids=np.searchsorted(keys,bk)
    if np.any(ids>=len(keys)) or not np.array_equal(unique[ids],boundary) or len(np.unique(ids))!=len(ids):
        raise ValueError('Missing/duplicate boundary faces')
    required=np.flatnonzero((count==1)|(low!=high))
    if not np.array_equal(np.sort(ids),required):raise ValueError('Boundary/interface coverage mismatch')
    for attr in np.unique(labels):
        f=ids[labels==attr]
        good = ((count[f]==2)&(low[f]==1)&(high[f]==2)) if attr==3100 else (
            (count[f]==1)&(low[f]==(1 if attr==5001 else 2)) if attr in (5001,6001) else count[f]==1)
        if not np.all(good):raise ValueError('Incorrect material adjacency for '+str(attr))
    adjacency=start[(count==2)&(low==high)]
    u,v=owner[adjacency],owner[adjacency+1]
    graph=coo_matrix((np.ones(len(u)),(u,v)),shape=(len(t),len(t))).tocsr()
    _,components=connected_components(graph,directed=False)
    component_counts={int(a):len(np.unique(components[material==a])) for a in np.unique(material)}
    xyz=p[b];cross=np.cross(xyz[:,1]-xyz[:,0],xyz[:,2]-xyz[:,0]);twice=np.linalg.norm(cross,axis=1)
    if np.any(twice<=0):raise ValueError('Degenerate boundary triangle')
    normals=cross/twice[:,None];pivot=np.argmax(abs(normals),axis=1)
    normals*=np.sign(normals[np.arange(len(normals)),pivot])[:,None]
    planes=np.column_stack((labels,np.round(normals,8),np.round(np.einsum('ij,ij->i',normals,xyz[:,0]),8)))
    planes[planes==0]=0.0  # Canonicalize signed zero for persistent plane keys.
    unique_planes,index=np.unique(planes,axis=0,return_inverse=True)
    plane_areas=np.bincount(index,weights=twice/2)
    return {'Tetrahedra':len(t),'SurfaceTriangles':len(b),'MaterialVolumes':volume,
            'MaterialComponents':component_counts,'MetalSurfaceComponents':boundary_component_counts(b,labels),
            'PlanarPatchAreas':{' '.join(map(str,k)):float(a) for k,a in zip(unique_planes,plane_areas)}},unique_planes


def directional_widths(mesh,recipe):
    t,_=blocks(mesh,'tetra');xyz=mesh.points[t];center=xyz.mean(axis=1)
    segments=recipe['PhysicalSegments'];corners=np.asarray(recipe['TruePhysicalCorners'])
    distance=np.full(len(t),np.inf);which=np.zeros(len(t),int)
    for i,s in enumerate(segments):
        r,_=segment_distances(center,s);mask=r<distance;which[mask]=i;distance[mask]=r[mask]
    dc=np.full(len(t),np.inf)
    for c in corners:dc=np.minimum(dc,np.linalg.norm(center-c,axis=1))
    fine=recipe['NormalSize'];tangent=recipe['TangentialSize'];reports={}
    for cutoff in (1.,3.):
        chosen=np.flatnonzero((distance<cutoff*fine)&(dc>5*tangent))
        widths=[]
        for i,s in enumerate(segments):
            ids=chosen[which[chosen]==i]
            if not len(ids):continue
            _,direction=segment_distances(center[ids],s)
            transverse=np.eye(3)-np.outer(direction,direction)
            centered=xyz[ids]-center[ids,None,:];projected=centered@transverse
            _,axes=np.linalg.eigh(projected.swapaxes(1,2)@projected)
            values=centered@axes[:,:,1:]
            widths.extend(np.column_stack((np.ptp(xyz[ids]@direction,axis=1),np.ptp(values,axis=1))))
        widths=np.asarray(widths)
        reports[str(cutoff)]={'Cells':len(chosen),'DistanceCutoff':cutoff*fine,
                             'WidthsTangentialTransverse1Transverse2':np.percentile(widths,[10,50,90],axis=0).tolist() if len(widths) else []}
    return reports


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('reference',type=Path);p.add_argument('candidate',type=Path);p.add_argument('recipe',type=Path);p.add_argument('output',type=Path)
    a=p.parse_args();reference=read_mesh(a.reference);candidate=read_mesh(a.candidate)
    before,_=analyze(reference);after,_=analyze(candidate)
    for key in ('MaterialVolumes','PlanarPatchAreas'):
        if before[key].keys()!=after[key].keys():raise ValueError('Changed '+key+' supports')
        error=max(abs(after[key][k]/v-1) for k,v in before[key].items())
        if error>1e-8:raise ValueError('Changed '+key+': '+str(error))
    if before['MaterialComponents']!=after['MaterialComponents']:raise ValueError('Material connectivity changed')
    if before['MetalSurfaceComponents']!=after['MetalSurfaceComponents']:raise ValueError('Conductor surface connectivity changed')
    report={'Reference':before,'Candidate':after,'DirectionalWidths':directional_widths(candidate,json.loads(a.recipe.read_text())),
            'Scope':'Conformity, material-component and planar-area checks; not full CAD-footprint or response qualification','LibraryQualified':False}
    a.output.write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))

if __name__=='__main__':main()
