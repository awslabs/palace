#!/usr/bin/env python3
# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0
"""Prepare a sharp, planar tetrahedral seed and native tensor array for MMG.

Diagnostic only. Feature extraction from the seed is not a substitute for CAD
provenance: the caller must independently check geometry and topology afterwards.
"""
import argparse
import hashlib
import json
from pathlib import Path
import meshio
import numpy as np
from edge_volume_metric import surface_features,volume_metric,segment_distances
from mesh_array_io import read_mesh,sha
from semantic_mesh_contract import (boundary_attributes, cut_surface_attributes,
                                    load_semantic_contract, simple_sharp_contract,
                                    volume_attributes)


def _near_segment(point, segments, tolerance):
    return any(segment_distances(np.asarray([point]), segment)[0][0] <= tolerance
               for segment in segments)


def validate_transformed_supports(data, semantic_contract, segments, tolerance=1e-8):
    """Bind transformed source supports to the semantic contract and seed features."""
    if (not isinstance(data, dict) or data.get('Version') != 1 or
            data.get('CoordinateSystem') != 'TransformedGlobal3D'):
        raise ValueError('Unsupported transformed support contract')
    transform_hash = data.get('CanonicalTransformSHA256')
    semantic_hash = data.get('SourceSemanticContractSHA256')
    source_hashes = [data.get(name) for name in
                     ('SourceSignatureSHA256', 'SourceBoundarySHA256', 'SourceMaskSHA256')]
    if (not isinstance(transform_hash, str) or len(transform_hash) != 64 or
            any(character not in '0123456789abcdef' for character in transform_hash) or
            not isinstance(semantic_hash, str) or len(semantic_hash) != 64 or
            any(character not in '0123456789abcdef' for character in semantic_hash) or
            any(not isinstance(value, str) or len(value) != 64 or
                any(character not in '0123456789abcdef' for character in value)
                for value in source_hashes)):
        raise ValueError('Transformed support provenance hashes are invalid')
    if (data.get('RigidTransform') != semantic_contract.get('RigidTransform') or
            transform_hash != semantic_contract.get('CanonicalTransformSHA256') or
            semantic_hash != semantic_contract.get('SourceSemanticContractSHA256')):
        raise ValueError('Transformed support provenance differs from semantic contract')
    edges = data.get('Edges')
    if not isinstance(edges, list) or not edges:
        raise ValueError('Transformed support contract has no edges')
    indices = []
    endpoint_classes = np.asarray(
        semantic_contract['SemanticCorners'] +
        semantic_contract['FeatureTopology']['CADSubdivisionEndpoints'] +
        semantic_contract['FeatureTopology']['CutEndpoints'], dtype=float).reshape(-1, 3)
    for edge in edges:
        if not isinstance(edge, dict):
            raise ValueError('Invalid transformed edge support')
        indices.append(edge.get('Index'))
        point = np.asarray(edge.get('Point'), dtype=float)
        gap = np.asarray(edge.get('GapDirection'), dtype=float)
        tangent = np.asarray(edge.get('TangentDirection'), dtype=float)
        normal = np.asarray(edge.get('ProcessNormalDirection'), dtype=float)
        interval = np.asarray(edge.get('Interval'), dtype=float)
        if (point.shape != (3,) or gap.shape != (3,) or tangent.shape != (3,) or
                normal.shape != (3,) or interval.shape != (2,) or
                not np.all(np.isfinite(np.concatenate((point, gap, tangent, normal,
                                                        interval)))) or
                not isinstance(edge.get('Slot'), int) or edge['Slot'] < 0 or
                not isinstance(edge.get('Conductor'), int) or edge['Conductor'] <= 0 or
                interval[0] >= interval[1]):
            raise ValueError('Invalid transformed edge point/vector support')
        frame = np.stack((gap, tangent, normal))
        if not np.allclose(frame @ frame.T, np.eye(3), rtol=0.0, atol=1e-10):
            raise ValueError('Transformed edge support frame is not orthonormal')
        endpoints = point[None, :] + interval[:, None] * tangent[None, :]
        if (any(not _near_segment(endpoint, segments, tolerance) for endpoint in endpoints) or
                any(np.min(np.linalg.norm(endpoint_classes - endpoint, axis=1)) > tolerance
                    for endpoint in endpoints)):
            raise ValueError('Transformed edge support differs from seed-derived features')
    if any(not isinstance(index, int) for index in indices) or len(set(indices)) != len(indices):
        raise ValueError('Transformed edge indices must be unique integers')
    physical_boundary = [item['Point'] for item in data.get('BoundaryVertices', [])
                         if item.get('Class') == 'Physical']
    if (len(physical_boundary) != len(semantic_contract['SemanticCorners']) or
            any(np.min(np.linalg.norm(np.asarray(physical_boundary) - corner, axis=1)) > tolerance
                for corner in semantic_contract['SemanticCorners'])):
        raise ValueError('Transformed physical boundary differs from semantic corners')
    return data


def prepare(mesh,path,normal,tangent,far,protected_distance=0.,far_growth=1.,protect_surface=0.,
            semantic_contract=None, transformed_supports=None):
    if not np.all(np.isfinite([protected_distance,far_growth,protect_surface])) or protected_distance<0 or far_growth<=0 or protect_surface<0:
        raise ValueError('Invalid grading/protection controls')
    path=Path(path)
    if path.exists():raise ValueError('Use a fresh attempt directory')
    if set(c.type for c in mesh.cells)-{'triangle','tetra'}:raise ValueError('Only linear triangle/tet seeds supported')
    refs=mesh.cell_data.get('gmsh:physical',mesh.cell_data.get('medit:ref'))
    if refs is None:raise ValueError('Missing physical references')
    triangles=np.concatenate([c.data for c in mesh.cells if c.type=='triangle'])
    triangle_refs=np.concatenate([r for c,r in zip(mesh.cells,refs) if c.type=='triangle'])
    tetrahedra=np.concatenate([c.data for c in mesh.cells if c.type=='tetra'])
    tetrahedron_refs=np.concatenate([r for c,r in zip(mesh.cells,refs) if c.type=='tetra'])
    if semantic_contract is None:
        raise ValueError('A frozen semantic contract is required')
    if set(triangle_refs)!=boundary_attributes(semantic_contract):
        raise ValueError('Seed boundary labels differ from the frozen semantic contract')
    if set(tetrahedron_refs)!=volume_attributes(semantic_contract):
        raise ValueError('Seed volume materials differ from the frozen semantic contract')
    features,pins,segments,corners=surface_features(
        mesh.points,triangles,triangle_refs,cut_surface_attributes(semantic_contract))
    if transformed_supports is not None:
        validate_transformed_supports(transformed_supports, semantic_contract, segments)
    semantic_corners=np.asarray(semantic_contract['SemanticCorners'],dtype=float).reshape(-1,3)
    # The contract corners are physical plan-view junctions.  Require them to be
    # represented by the seed instead of silently replacing them with CAD
    # subdivision, extrusion, or coupon-cut vertices discovered from triangles.
    nearest=np.min(np.linalg.norm(mesh.points[:,None]-semantic_corners[None],axis=2),axis=0)
    if np.any(nearest>1e-8):raise ValueError('Semantic corner is absent from the seed mesh')
    xyz=mesh.points[triangles]
    normals=np.cross(xyz[:,1]-xyz[:,0],xyz[:,2]-xyz[:,0]);normals/=np.linalg.norm(normals,axis=1)[:,None]
    pivot=np.argmax(abs(normals),axis=1);normals*=np.sign(normals[np.arange(len(normals)),pivot])[:,None]
    support=np.column_stack((triangle_refs,np.round(normals,8),np.round(np.einsum('ij,ij->i',normals,xyz[:,0]),8)))
    support[support==0]=0.
    planes,first,patch=np.unique(support,axis=0,return_index=True,return_inverse=True)
    exact_planes=np.column_stack((triangle_refs[first],normals[first],
                                  np.einsum('ij,ij->i',normals[first],xyz[first,0])))
    # Separate planar supports during adaptation. Restore original physical labels
    # only after independently checking/projecting each support intersection.
    patch_references=(10000+patch).astype(np.int32)
    path.mkdir(parents=True)
    inputs=meshio.Mesh(mesh.points,[('triangle',triangles.astype(np.int32)),
        ('tetra',tetrahedra.astype(np.int32)),('line',features.astype(np.int32))],
        point_data={'medit:ref':np.zeros(len(mesh.points),dtype=np.int32)},
        cell_data={'medit:ref':[patch_references,tetrahedron_refs.astype(np.int32),
                               np.full(len(features),1,dtype=np.int32)]})
    # This installed MMG build rejects meshio's binary version-3 header. ASCII
    # Medit version 2 is unambiguous; output can still use native MMG .meshb.
    meshio.write(path/'seed.mesh',inputs,file_format='medit')
    np.savetxt(path/'pins.txt',pins+1,fmt='%d')
    fixed=np.isin(triangle_refs,list(cut_surface_attributes(semantic_contract)))
    diameter=np.max(np.stack([np.linalg.norm(xyz[:,i]-xyz[:,j],axis=1)
        for i,j in ((0,1),(1,2),(2,0))]),axis=0)
    if protect_surface>0:
        distance=np.full(len(mesh.points),np.inf)
        for segment in segments:distance=np.minimum(distance,segment_distances(mesh.points,segment)[0])
        # Distance is 1-Lipschitz. This lower bound protects every triangle that
        # could intersect the anisotropic edge band, rather than relying only on
        # its centroid.
        fixed |= distance[triangles].min(axis=1)-diameter <= protect_surface
    corner_distance=np.full(len(mesh.points),np.inf)
    for corner in semantic_corners:
        corner_distance=np.minimum(corner_distance,np.linalg.norm(mesh.points-corner,axis=1))
    # Local semantic isotropy requires surface remeshing too.  Permit inserted
    # vertices only in the same tangentially-scaled balls used by the SPD
    # metric; matching and edge-band triangles remain frozen everywhere else.
    fixed &= corner_distance[triangles].min(axis=1)-diameter > tangent
    np.savetxt(path/'fixed-triangles.txt',np.flatnonzero(fixed)+1,fmt='%d')
    with (path/'metric.f64').open('wb') as f:
        for start in range(0,len(mesh.points),100000):
            metric=volume_metric(mesh.points[start:start+100000],segments,corners,normal,tangent,far,
                                 protected_distance=protected_distance,far_growth=far_growth,
                                 isotropic_corners=semantic_corners,isotropy_radius=tangent)
            if np.any(np.linalg.eigvalsh(metric)<=0):raise ValueError('Metric is not SPD')
            # Native C API order, explicitly NOT the Medit .sol file order.
            metric[:,[0,0,0,1,1,2],[0,1,2,1,2,2]].astype('<f8').tofile(f)
    recipe={'Scope':'Native MMG straight/sharp fabricated diagnostic; no solver qualification',
            'NormalSize':normal,'TangentialSize':tangent,'FarSize':far,
            'RadialGrowth':1.,'CornerGrowth':.25,'ProtectedDistance':protected_distance,
            'FarGrowth':far_growth,'SurfaceProtectionRadius':protect_surface,'FixedSurfaceTriangles':int(fixed.sum()),
            'CornerIsotropyRadius':tangent,
            'MetricOrder':['m11','m12','m13','m22','m23','m33'],
            'Nodes':len(mesh.points),'Tetrahedra':len(tetrahedra),'SurfaceTriangles':len(triangles),
            'PreservedFeatureEdges':len(features),'PinnedGeometryVertices':len(pins),
            'PhysicalSegments':segments.tolist(),'TruePhysicalCorners':semantic_corners.tolist(),
            'SurfaceFeatureCorners':corners.tolist(),
            'PlanarSupports':{str(10000+i):{'Attribute':int(row[0]),'Normal':row[1:4].tolist(),'Offset':float(row[4])} for i,row in enumerate(exact_planes)},
            'SemanticContract':semantic_contract,'LibraryQualified':False}
    (path/'recipe.json').write_text(json.dumps(recipe,indent=2)+'\n')
    print(json.dumps({k:v for k,v in recipe.items() if k not in ('PhysicalSegments','TruePhysicalCorners')},indent=2))
    return recipe


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('mesh',type=Path);p.add_argument('output',type=Path)
    p.add_argument('--normal',type=float,required=True);p.add_argument('--tangent',type=float,required=True);p.add_argument('--far',type=float,required=True)
    p.add_argument('--protected-distance',type=float,default=0.);p.add_argument('--far-growth',type=float,default=1.)
    p.add_argument('--protect-surface',type=float,default=0.)
    contract=p.add_mutually_exclusive_group(required=True)
    contract.add_argument('--semantic-contract',type=Path)
    contract.add_argument('--simple-sharp-contract',action='store_true',
                          help='use the explicit historical 1/3100/5001/6001 compatibility contract')
    p.add_argument('--transformed-supports',type=Path,
                   help='required transformed source-support contract with --semantic-contract')
    a=p.parse_args();m=read_mesh(a.mesh)
    if bool(a.semantic_contract) != bool(a.transformed_supports):
        p.error('--semantic-contract and --transformed-supports are required together')
    semantic=(load_semantic_contract(a.semantic_contract) if a.semantic_contract
              else simple_sharp_contract())
    supports=(json.loads(a.transformed_supports.read_text()) if a.transformed_supports else None)
    r=prepare(m,a.output,a.normal,a.tangent,a.far,
        a.protected_distance,a.far_growth,a.protect_surface,semantic,supports)
    r['SeedArtifact']=str(a.mesh.resolve());r['SeedArtifactSHA256']=sha(a.mesh)
    if a.transformed_supports:
        r['TransformedSupportsArtifact']=str(a.transformed_supports.resolve())
        r['TransformedSupportsSHA256']=sha(a.transformed_supports)
    if a.mesh.suffix=='.toml':
        import tomllib
        data=tomllib.loads(a.mesh.read_text())
        r['SeedMesh']=data['SourceMesh'];r['SeedSHA256']=data['SourceSHA256']
        r['ArraySHA256']={name:sha(a.mesh.parent/name) for name in [data['NodeTags'],data['Coordinates']]+[c['Connectivity'] for c in data['Cells']]}
    else:
        r['SeedMesh']=str(a.mesh.resolve());r['SeedSHA256']=sha(a.mesh)
    (a.output/'recipe.json').write_text(json.dumps(r,indent=2)+'\n')

if __name__=='__main__':main()
