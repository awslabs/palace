#!/usr/bin/env python3
# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0
"""Bounded planar CAD-support restoration after an experimental MMG adaptation.

Native output is retained. Large displacement, missing support, inconsistent
intersection or inversion rejects the candidate. This does NOT replace footprint,
conductor-topology, achieved-resolution or PDE accuracy checks.
"""
import argparse,json,hashlib
from pathlib import Path
import meshio
import numpy as np


def restore(mesh,recipe,maximum_displacement):
    supports={int(k):v for k,v in recipe['PlanarSupports'].items()}
    refs=mesh.cell_data['medit:ref'];node_supports={};seen=set()
    for block,attributes in zip(mesh.cells,refs):
        if block.type!='triangle':continue
        for tri,ref in zip(block.data,attributes):
            if int(ref) not in supports:raise ValueError('Unexpected planar-support reference')
            seen.add(int(ref))
            for node in tri:node_supports.setdefault(int(node),set()).add(int(ref))
    if seen!=set(supports):raise ValueError('Lost planar support')
    groups={}
    for node,ids in node_supports.items():groups.setdefault(tuple(sorted(ids)),[]).append(node)
    points=mesh.points.copy();largest=0.;moved=0
    for ids,nodes in groups.items():
        a=np.array([supports[i]['Normal'] for i in ids]);b=np.array([supports[i]['Offset'] for i in ids])
        # Normals of coincident supports can differ at floating-point roundoff.
        # Do not interpret that noise as a real extra intersection constraint.
        original=points[nodes].copy();correction=(b-original@a.T)@np.linalg.pinv(a,rcond=1e-10).T
        largest=max(largest,float(np.linalg.norm(correction,axis=1).max()))
        if largest>maximum_displacement:raise ValueError(f'CAD correction exceeds bound: {largest} > {maximum_displacement}')
        points[nodes]+=correction;moved+=int(np.sum(np.any(correction!=0,axis=1)))
        if np.max(abs(points[nodes]@a.T-b))>1e-10:raise ValueError('Inconsistent planar intersection')
    cells=[];attributes=[]
    for block,ref in zip(mesh.cells,refs):
        if block.type not in ('tetra','triangle'):continue
        if block.type=='tetra':
            p=points[block.data]
            det=np.einsum('ij,ij->i',np.cross(p[:,1]-p[:,0],p[:,2]-p[:,0]),p[:,3]-p[:,0])
            if np.any(det<=0):raise ValueError('Projection inverted a tetrahedron')
            attributes.append(ref)
        else:attributes.append(np.array([supports[int(r)]['Attribute'] for r in ref]))
        cells.append((block.type,block.data))
    semantic=recipe.get('SemanticContract')
    if not isinstance(semantic,dict):raise ValueError('Restoration recipe lacks semantic contract')
    field_data={item['Material']:np.array([item['Attribute'],3],dtype=int)
                for item in semantic['VolumeMaterials']}
    field_data.update({('matching_surface' if item['Attribute']==1 else f"surface_{item['Attribute']}"):
                       np.array([item['Attribute'],2],dtype=int)
                       for item in semantic['BoundaryLabels']})
    output=meshio.Mesh(points,cells,cell_data={'gmsh:physical':attributes,'gmsh:geometrical':attributes},
                       field_data=field_data)
    return output,{'MaximumCorrectionUm':largest,'CorrectionBoundUm':maximum_displacement,
                   'MovedBoundaryVertices':moved,'PlanarSupports':len(supports),'LibraryQualified':False}


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('input',type=Path);p.add_argument('recipe',type=Path);p.add_argument('output',type=Path)
    bound=p.add_mutually_exclusive_group(required=True)
    bound.add_argument('--max-displacement',type=float)
    bound.add_argument('--max-displacement-over-normal',type=float)
    a=p.parse_args()
    if a.output.exists():raise ValueError('Do not overwrite candidates')
    recipe=json.loads(a.recipe.read_text())
    maximum=(a.max_displacement if a.max_displacement is not None else
             a.max_displacement_over_normal*float(recipe['NormalSize']))
    if not np.isfinite(maximum) or maximum<=0:raise ValueError('Invalid displacement bound')
    mesh=meshio.read(a.input);output,report=restore(mesh,recipe,maximum)
    if a.max_displacement_over_normal is not None:
        report['CorrectionBoundOverNormalSize']=a.max_displacement_over_normal
    meshio.write(a.output,output,file_format='gmsh22',binary=True)
    report['NativeInputSHA256']=hashlib.sha256(a.input.read_bytes()).hexdigest();report['OutputSHA256']=hashlib.sha256(a.output.read_bytes()).hexdigest()
    a.output.with_suffix('.projection.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))

if __name__=='__main__':main()
