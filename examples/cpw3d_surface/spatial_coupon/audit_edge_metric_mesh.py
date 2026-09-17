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
from edge_volume_metric import (cluster_coplanar_triangles,edge_layer_cells,match_equivalent_planes,
                                segment_distances)
from mesh_array_io import read_mesh
from semantic_mesh_contract import (boundary_adjacency, boundary_attributes,
                                    load_semantic_contract, metric_surface_attributes,
                                    volume_attributes)


def blocks(mesh,kind):
    key='gmsh:physical' if 'gmsh:physical' in mesh.cell_data else 'medit:ref'
    width=3 if kind=='triangle' else 4 if kind=='tetra' else None
    selected=[(c.data[:,:width],a) for c,a in zip(mesh.cells,mesh.cell_data[key])
              if c.type==kind or c.type.startswith(kind)]
    if not selected:raise ValueError('Mesh lacks '+kind+' elements')
    return (np.concatenate([c for c,_ in selected]),np.concatenate([a for _,a in selected]))


def boundary_component_counts(triangles,labels,metal_attributes):
    """Connected metric-surface components, without joining through other faces."""
    metal=triangles[np.isin(labels,list(metal_attributes))]
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


def physical_names(mesh,dimension):
    result={}
    for name,value in mesh.field_data.items():
        tag,dim=map(int,value[:2])
        if dim==dimension:
            if tag in result:raise ValueError('Duplicate physical-name attribute')
            result[tag]=name
    return result


def analyze(mesh,contract,require_material_names=False):
    t,material=blocks(mesh,'tetra');b,labels=blocks(mesh,'triangle');p=mesh.points
    if set(material)!=volume_attributes(contract):
        raise ValueError('Volume materials differ from the frozen semantic contract')
    names=physical_names(mesh,3)
    expected_names={item['Attribute']:item['Material'] for item in contract['VolumeMaterials']}
    if require_material_names and names!=expected_names:
        raise ValueError('Physical volume names differ from the frozen semantic contract')
    if set(labels)!=boundary_attributes(contract):
        raise ValueError('Boundary labels differ from the frozen semantic contract')
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
    expected_adjacency=boundary_adjacency(contract);actual_adjacency={}
    for attr in np.unique(labels):
        observed=set()
        for f in ids[labels==attr]:
            actual={int(low[f])} if count[f]==1 else {int(low[f]),int(high[f])}
            if actual not in expected_adjacency[int(attr)]:
                raise ValueError('Incorrect material adjacency for '+str(attr))
            observed.update(actual)
        actual_adjacency[int(attr)]=sorted(observed)
    adjacency=start[(count==2)&(low==high)]
    u,v=owner[adjacency],owner[adjacency+1]
    graph=coo_matrix((np.ones(len(u)),(u,v)),shape=(len(t),len(t))).tocsr()
    _,components=connected_components(graph,directed=False)
    component_counts={int(a):len(np.unique(components[material==a])) for a in np.unique(material)}
    xyz=p[b];cross=np.cross(xyz[:,1]-xyz[:,0],xyz[:,2]-xyz[:,0]);twice=np.linalg.norm(cross,axis=1)
    if np.any(twice<=0):raise ValueError('Degenerate boundary triangle')
    normals=cross/twice[:,None];pivot=np.argmax(abs(normals),axis=1)
    normals*=np.sign(normals[np.arange(len(normals)),pivot])[:,None]
    # Planar patches are equivalence classes of (label, plane) under the shared
    # tolerance rule, never rounded coordinates: roundoff-level normal noise on
    # one CAD plane must not fragment it into several keys.
    representatives,patch=cluster_coplanar_triangles(labels,normals,xyz)
    plane_areas=np.bincount(patch,weights=twice/2,minlength=len(representatives))
    return {'Tetrahedra':len(t),'SurfaceTriangles':len(b),'MaterialVolumes':volume,
            'PhysicalVolumeNames':names,'BoundaryAdjacency':actual_adjacency,
            'MaterialComponents':component_counts,
            'MetalSurfaceComponents':boundary_component_counts(
                b,labels,metric_surface_attributes(contract)),
            'PlanarPatchEquivalence':'edge_volume_metric.cluster_coplanar_triangles',
            'PlanarPatchAreas':{planar_patch_key(row):float(a)
                                for row,a in zip(representatives,plane_areas)}},(representatives,patch)


def planar_patch_key(representative):
    """Report key of a planar-patch representative: label and its first triangle's plane."""
    row=np.asarray(representative,dtype=float)
    plane=np.r_[row[1:4],np.dot(row[1:4],row[4:7])]
    plane[plane==0]=0.0  # Canonicalize signed zero for persistent plane keys.
    return ' '.join([str(int(row[0])),*map(str,plane)])


def matched_planar_patch_areas(before,after):
    """Pair the planar patches of two analyses of one geometry by plane equivalence.

    Returns [(reference key, reference area, candidate area)] or None when the
    patch sets do not correspond one to one.
    """
    (left,_),(right,_)=before[1],after[1]
    mapping=match_equivalent_planes(left,right)
    if mapping is None:return None
    return [(planar_patch_key(row),before[0]['PlanarPatchAreas'][planar_patch_key(row)],
             after[0]['PlanarPatchAreas'][planar_patch_key(right[j])])
            for row,j in zip(left,mapping)]


def directional_widths(mesh,recipe,layer_spans=None,layer_reach=None):
    """Tangential / transverse extents of the band cells (centroid within 1 or 3
    NormalSize of a physical segment, away from the corners) as 10/50/90 percentiles.

    With a recorded edge layer (supervisor decision 31) the layer cells
    (edge_volume_metric.EDGE_LAYER_CELL_RULE: a vertex within layer_reach, the
    recipe's RequiredReach, of a layer span) are the layer's own statistics
    ('EdgeLayer', same percentiles) and are excluded from the band statistics: the
    band anisotropy design gate judges the metric-driven band, the layer's design
    statement is the bound EdgeLayer aspect rule.
    """
    t,_=blocks(mesh,'tetra');xyz=mesh.points[t];center=xyz.mean(axis=1)
    segments=recipe['PhysicalSegments'];corners=np.asarray(recipe['TruePhysicalCorners'])
    distance=np.full(len(t),np.inf);which=np.zeros(len(t),int)
    for i,s in enumerate(segments):
        r,_=segment_distances(center,s);mask=r<distance;which[mask]=i;distance[mask]=r[mask]
    dc=np.full(len(t),np.inf)
    for c in corners:dc=np.minimum(dc,np.linalg.norm(center-c,axis=1))
    in_layer=(np.zeros(len(t),dtype=bool) if layer_spans is None else
              edge_layer_cells(mesh.points,t,layer_spans,layer_reach))
    fine=recipe['NormalSize'];tangent=recipe['TangentialSize'];reports={}
    def percentiles(ids):
        widths=[]
        for i,s in enumerate(segments):
            selected=ids[which[ids]==i]
            if not len(selected):continue
            _,direction=segment_distances(center[selected],s)
            transverse=np.eye(3)-np.outer(direction,direction)
            centered=xyz[selected]-center[selected,None,:];projected=centered@transverse
            _,axes=np.linalg.eigh(projected.swapaxes(1,2)@projected)
            values=centered@axes[:,:,1:]
            widths.extend(np.column_stack((np.ptp(xyz[selected]@direction,axis=1),np.ptp(values,axis=1))))
        widths=np.asarray(widths)
        return np.percentile(widths,[10,50,90],axis=0).tolist() if len(widths) else []
    for cutoff in (1.,3.):
        band=(distance<cutoff*fine)&(dc>5*tangent)
        chosen=np.flatnonzero(band&~in_layer);layer=np.flatnonzero(band&in_layer)
        reports[str(cutoff)]={'Cells':len(chosen),'DistanceCutoff':cutoff*fine,
                             'WidthsTangentialTransverse1Transverse2':percentiles(chosen),
                             'ExcludedEdgeLayerCells':int(len(layer))}
        if layer_spans is not None:
            reports[str(cutoff)]['EdgeLayer']={'Cells':int(len(layer)),'Reach':float(layer_reach),
                'WidthsTangentialTransverse1Transverse2':percentiles(layer)}
    return reports


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('reference',type=Path);p.add_argument('candidate',type=Path);p.add_argument('recipe',type=Path);p.add_argument('semantic_contract',type=Path);p.add_argument('output',type=Path)
    a=p.parse_args();reference=read_mesh(a.reference);candidate=read_mesh(a.candidate)
    contract=load_semantic_contract(a.semantic_contract)
    before=analyze(reference,contract);after=analyze(candidate,contract)
    if before[0]['MaterialVolumes'].keys()!=after[0]['MaterialVolumes'].keys():
        raise ValueError('Changed MaterialVolumes supports')
    error=max(abs(after[0]['MaterialVolumes'][k]/v-1) for k,v in before[0]['MaterialVolumes'].items())
    if error>1e-8:raise ValueError('Changed MaterialVolumes: '+str(error))
    matched=matched_planar_patch_areas(before,after)
    if matched is None:raise ValueError('Changed PlanarPatchAreas supports')
    error=max(abs(right/left-1) for _,left,right in matched)
    if error>1e-8:raise ValueError('Changed PlanarPatchAreas: '+str(error))
    before,after=before[0],after[0]
    if before['MaterialComponents']!=after['MaterialComponents']:raise ValueError('Material connectivity changed')
    if before['MetalSurfaceComponents']!=after['MetalSurfaceComponents']:raise ValueError('Conductor surface connectivity changed')
    report={'Reference':before,'Candidate':after,'DirectionalWidths':directional_widths(candidate,json.loads(a.recipe.read_text())),
            'Scope':'Conformity, material-component and planar-area checks; not full CAD-footprint or response qualification','LibraryQualified':False}
    a.output.write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))

if __name__=='__main__':main()
