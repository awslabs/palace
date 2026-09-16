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
from edge_volume_metric import (COPLANAR_TOLERANCE,cluster_coplanar_triangles,junction_segments,
                                surface_features,volume_metric,segment_distances)
from mesh_array_io import read_mesh,sha
from mesh_stage_contract import footprint_provenance,footprint_segments
from semantic_mesh_contract import (boundary_attributes, cut_surface_attributes,
                                    load_semantic_contract, material_interface_attributes,
                                    simple_sharp_contract, volume_attributes)


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


def cluster_planar_supports(attributes, normals, xyz):
    """Merge numerically noisy representations of the same CAD plane support.

    One plane per (attribute, equivalence class) under the shared
    edge_volume_metric.cluster_coplanar_triangles rule; the audits key planar
    patches with the same function, so a support the metric stage records is a
    patch the audits recognize.  Returns (planes, patch): planes rows are
    [attribute, nx, ny, nz, offset] of each cluster's first triangle.
    """
    representatives, patch = cluster_coplanar_triangles(attributes, normals, xyz)
    offsets = np.einsum('ij,ij->i', representatives[:, 1:4], representatives[:, 4:7])
    planes = np.column_stack((representatives[:, 0], representatives[:, 1:4], offsets))
    return planes, patch.astype(np.int32)


def budget_aware_far_policy(seed_elements, maximum_elements, far_size, far_growth):
    """Coarsen only the far field when the seed consumes the adaptation budget.

    The policy is geometry/model agnostic.  A seed may use at most 35% of the
    final budget without pressure; above that point the far-size ceiling and
    post-protection growth increase continuously with the square root of load.
    Normal, tangent, semantic-corner, trace, and protected-band targets are not
    inputs and therefore cannot be changed by this policy.
    """
    values = (seed_elements, maximum_elements, far_size, far_growth)
    if (not all(np.isfinite(values)) or seed_elements <= 0 or
            maximum_elements <= 0 or far_size <= 0 or far_growth <= 0):
        raise ValueError('Invalid far-field budget controls')
    seed_fraction = float(seed_elements) / float(maximum_elements)
    pressure = max(1.0, np.sqrt(seed_fraction / 0.35))
    return {'Name': 'seed-fraction-far-field-v1',
            'MaximumElements': int(maximum_elements),
            'SeedElements': int(seed_elements), 'SeedBudgetFraction': seed_fraction,
            'Pressure': float(pressure), 'RequestedFarSize': float(far_size),
            'EffectiveFarSize': float(far_size * pressure),
            'RequestedFarGrowth': float(far_growth),
            'EffectiveFarGrowth': float(far_growth * pressure),
            'UnaffectedTargets': ['normal', 'tangent', 'semantic-corner',
                                  'trace', 'protected-band']}


def protected_corner_ball_triangles(points, triangles, semantic_corners, radius):
    """Surface triangles with at least one vertex inside a semantic corner ball.

    The seed's isotropic corner ball at NormalSize is the intended corner
    discretization; MMG's anisotropic adaptation adds value along edges, not at
    corners, so the balls are protected supports like the edge bands.  Any vertex
    inside the ball freezes the triangle, so the frozen set covers the ball.
    Returns the boolean triangle mask and the frozen-triangle count per corner.
    """
    corners = np.asarray(semantic_corners, dtype=float).reshape(-1, 3)
    if not np.isfinite(radius) or radius <= 0:
        raise ValueError('Invalid corner isotropy radius')
    inside = np.zeros((len(triangles), len(corners)), dtype=bool)
    vertices = np.asarray(points)[np.asarray(triangles)]
    for index, corner in enumerate(corners):
        inside[:, index] = np.linalg.norm(vertices - corner, axis=2).min(axis=1) <= radius
    return inside.any(axis=1), inside.sum(axis=0)


def footprint_segment_record(census, semantic_contract, semantic_contract_sha256):
    """Recipe record of the seed's simplified etch footprint edges (decision 18(c)).

    Etch footprint edges are physical dielectric step edges (trench wall/floor
    and wall/surface junctions) and legitimately carry the NormalSize band that
    the seed-derived PhysicalSegments already give them; recording them from the
    bound census lets the audits distinguish a band along a footprint edge from
    a diagonal over-refinement.  No size is attached to them here.
    """
    if (not isinstance(census, dict) or census.get("Version") != 1 or
            census.get("SemanticContractSHA256") != semantic_contract_sha256 or
            census.get("RigidTransform") != semantic_contract.get("RigidTransform")):
        raise ValueError("Seed census does not belong to the bound semantic contract")
    return {"Provenance": footprint_provenance(census),
            "EtchBoundary": census.get("EtchBoundary"),
            "Tolerance": census.get("FootprintCollinearTolerance"),
            "Polygons": len(census["FootprintPolygons"]),
            "Segments": footprint_segments(census),
            "Rule": "simplified etch footprint polygon edges (device retained etch or "
                    "producer default) are physical dielectric step edges: legitimate "
                    "feature segments for the audits alongside the signature edges; "
                    "they carry the NormalSize band through the seed-derived "
                    "PhysicalSegments and receive no size of their own"}


def junction_segment_record(segments, semantic_contract, census=None):
    """Recipe record of the cut-surface/material-interface junction lines.

    The lines where the Dirichlet cut surface meets a dielectric step or
    material interface (trench floor, trench walls, un-etched substrate-vacuum
    plane) are the features the physics pilot found under-resolved (100% of the
    AMR marks within 0.3 um of the cut/trench junction).  They are derived from
    the semantic contract's roles (CutSurfaceRoles and two-material adjacency)
    applied to the seed's non-coplanar shared edges, and they receive exactly the
    PhysicalSegments band law.  When the bound seed census records the CAD
    junction curves, their total length must agree with the mesh-derived
    segments within the shared dimensionless tolerance.
    """
    segments = np.asarray(segments, dtype=float).reshape(-1, 6)
    lengths = np.linalg.norm(segments[:, 3:] - segments[:, :3], axis=1)
    if not len(segments) or not np.all(np.isfinite(segments)) or np.any(lengths <= 0):
        raise ValueError('Junction segments must be finite and nondegenerate')
    total = float(lengths.sum())
    if census is not None:
        curves = census.get('JunctionCurves')
        if (not isinstance(curves, dict) or not isinstance(curves.get('TotalLength'), (int, float)) or
                isinstance(curves.get('TotalLength'), bool) or
                abs(float(curves['TotalLength']) - total) > COPLANAR_TOLERANCE * total):
            raise ValueError('Seed census junction curves differ from the seed-derived junction segments')
    return {'Segments': segments.tolist(), 'Count': int(len(segments)), 'TotalLength': total,
            'CutSurfaceAttributes': sorted(cut_surface_attributes(semantic_contract)),
            'MaterialInterfaceAttributes': sorted(material_interface_attributes(semantic_contract)),
            'Provenance': 'seed shared edges that are geometric features (not coplanar within '
                          'COPLANAR_TOLERANCE) with a cut-surface triangle and a material-'
                          'interface triangle, chained into straight segments; roles from the '
                          'bound semantic contract (CutSurfaceRoles; labels whose '
                          'AdjacentMaterialSets contain two materials), no coordinates or '
                          'label numbers assumed',
            'Rule': 'junction segments receive exactly the PhysicalSegments band law: '
                    'NormalSize transverse band with ProtectedDistance/FarGrowth grading, '
                    'TangentialSize along the line, SurfaceProtectionRadius freeze of the '
                    'incident seed surface; they are aligned features for the trace-diagonal '
                    'audit; PhysicalSegments are unchanged and precede them in the metric '
                    'intersection order'}


def prepare(mesh,path,normal,tangent,far,protected_distance=0.,far_growth=1.,protect_surface=0.,
            semantic_contract=None, transformed_supports=None, maximum_elements=None,
            footprint_census=None, semantic_contract_sha256=None):
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
    if maximum_elements is None:
        maximum_elements = max(len(tetrahedra), 1)
    far_policy=budget_aware_far_policy(
        len(tetrahedra),maximum_elements,far,far_growth)
    effective_far=far_policy['EffectiveFarSize']
    effective_far_growth=far_policy['EffectiveFarGrowth']
    if semantic_contract is None:
        raise ValueError('A frozen semantic contract is required')
    if set(triangle_refs)!=boundary_attributes(semantic_contract):
        raise ValueError('Seed boundary labels differ from the frozen semantic contract')
    if set(tetrahedron_refs)!=volume_attributes(semantic_contract):
        raise ValueError('Seed volume materials differ from the frozen semantic contract')
    features,pins,pin_kinds,segments,corners=surface_features(
        mesh.points,triangles,triangle_refs,cut_surface_attributes(semantic_contract))
    if transformed_supports is not None:
        validate_transformed_supports(transformed_supports, semantic_contract, segments)
    footprint=(None if footprint_census is None else
               footprint_segment_record(footprint_census,semantic_contract,semantic_contract_sha256))
    # The cut/trench junction lines are metric sources like the physical edges;
    # the metric intersection order (physical, then junction) is part of the recipe.
    junctions=junction_segments(mesh.points,triangles,triangle_refs,
                                cut_surface_attributes(semantic_contract),
                                material_interface_attributes(semantic_contract))
    junction_record=junction_segment_record(junctions,semantic_contract,footprint_census)
    band_segments=np.vstack((np.asarray(segments,dtype=float).reshape(-1,6),junctions))
    semantic_corners=np.asarray(semantic_contract['SemanticCorners'],dtype=float).reshape(-1,3)
    # The contract corners are physical plan-view junctions.  Require them to be
    # represented by the seed instead of silently replacing them with CAD
    # subdivision, extrusion, or coupon-cut vertices discovered from triangles.
    nearest=np.min(np.linalg.norm(mesh.points[:,None]-semantic_corners[None],axis=2),axis=0)
    if np.any(nearest>1e-8):raise ValueError('Semantic corner is absent from the seed mesh')
    xyz=mesh.points[triangles]
    normals=np.cross(xyz[:,1]-xyz[:,0],xyz[:,2]-xyz[:,0]);normals/=np.linalg.norm(normals,axis=1)[:,None]
    pivot=np.argmax(abs(normals),axis=1);normals*=np.sign(normals[np.arange(len(normals)),pivot])[:,None]
    exact_planes,patch=cluster_planar_supports(triangle_refs,normals,xyz)
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
        for segment in band_segments:distance=np.minimum(distance,segment_distances(mesh.points,segment)[0])
        # Distance is 1-Lipschitz. This lower bound protects every triangle that
        # could intersect the anisotropic edge band, rather than relying only on
        # its centroid.
        fixed |= distance[triangles].min(axis=1)-diameter <= protect_surface
    # The seed already carries the isotropic corner ball at NormalSize (its census
    # asserts it), so the ball surface is a protected support: MMG adapts the
    # volume inside the ball to the same isotropic metric but cannot re-mesh the
    # corner surface.
    corner_balls,corner_ball_counts=protected_corner_ball_triangles(
        mesh.points,triangles,semantic_corners,tangent)
    fixed |= corner_balls
    np.savetxt(path/'fixed-triangles.txt',np.flatnonzero(fixed)+1,fmt='%d')
    with (path/'metric.f64').open('wb') as f:
        for start in range(0,len(mesh.points),100000):
            metric=volume_metric(mesh.points[start:start+100000],band_segments,corners,normal,tangent,effective_far,
                                 protected_distance=protected_distance,far_growth=effective_far_growth,
                                 isotropic_corners=semantic_corners,isotropy_radius=tangent)
            if np.any(np.linalg.eigvalsh(metric)<=0):raise ValueError('Metric is not SPD')
            # Native C API order, explicitly NOT the Medit .sol file order.
            metric[:,[0,0,0,1,1,2],[0,1,2,1,2,2]].astype('<f8').tofile(f)
    recipe={'Scope':'Native MMG straight/sharp fabricated diagnostic; no solver qualification',
            'NormalSize':normal,'TangentialSize':tangent,'FarSize':effective_far,
            'RequestedFarSize':far,'RadialGrowth':1.,'CornerGrowth':.25,
            'ProtectedDistance':protected_distance,
            'FarGrowth':effective_far_growth,'RequestedFarGrowth':far_growth,
            'FarFieldBudgetPolicy':far_policy,
            'SurfaceProtectionRadius':protect_surface,'FixedSurfaceTriangles':int(fixed.sum()),
            'CornerIsotropyRadius':tangent,
            'ProtectedCornerBalls':{'Radius':tangent,
                'Rule':'every seed surface triangle with a vertex within CornerIsotropyRadius '
                       'of a semantic corner is frozen; the seed corner ball is the corner '
                       'discretization and MMG adapts only the volume inside it',
                'FrozenTriangles':int(corner_balls.sum()),
                'PerCorner':[{'Point':corner.tolist(),'FrozenTriangles':int(count)}
                             for corner,count in zip(semantic_corners,corner_ball_counts)]},
            'MetricOrder':['m11','m12','m13','m22','m23','m33'],
            'Nodes':len(mesh.points),'Tetrahedra':len(tetrahedra),'SurfaceTriangles':len(triangles),
            'PreservedFeatureEdges':len(features),'PinnedGeometryVertices':len(pins),
            'PinnedVertexKinds':{kind:int(np.sum(pin_kinds==kind))
                                 for kind in ('geometric-corner','reference-turn')},
            'PinnedVertices':[{'Point':mesh.points[node].tolist(),'Kind':str(kind)}
                              for node,kind in zip(pins,pin_kinds)],
            'PinPolicy':'turns/junctions of the complete reference-boundary graph are MMG '
                        'required vertices because MMG reconstructs reference boundaries as '
                        'geometric curves; ridges, metric segments and protected bands are '
                        'non-coplanar features only; reference-turn pins are not semantic '
                        'corners or protected supports',
            'PhysicalSegments':segments.tolist(),'TruePhysicalCorners':semantic_corners.tolist(),
            'JunctionSegments':junction_record,
            'BandSegmentOrder':'PhysicalSegments then JunctionSegments',
            'SurfaceFeatureCorners':corners.tolist(),
            'PlanarSupports':{str(10000+i):{'Attribute':int(row[0]),'Normal':row[1:4].tolist(),'Offset':float(row[4])} for i,row in enumerate(exact_planes)},
            'PlanarSupportEquivalence':{'Tolerance':COPLANAR_TOLERANCE,
                'Rule':'edge_volume_metric.cluster_coplanar_triangles: same attribute, sine of '
                       'the normal angle and point offset relative to max(local size, distance) '
                       'within the tolerance; shared with the planar-patch audits'},
            'SemanticContract':semantic_contract,'LibraryQualified':False}
    if footprint is not None:recipe['FootprintSegments']=footprint
    (path/'recipe.json').write_text(json.dumps(recipe,indent=2)+'\n')
    print(json.dumps({k:v for k,v in recipe.items() if k not in ('PhysicalSegments','TruePhysicalCorners','FootprintSegments','JunctionSegments')},indent=2))
    print(json.dumps({'JunctionSegments':{k:v for k,v in junction_record.items() if k!='Segments'}},indent=2))
    return recipe


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('mesh',type=Path);p.add_argument('output',type=Path)
    p.add_argument('--normal',type=float,required=True);p.add_argument('--tangent',type=float,required=True);p.add_argument('--far',type=float,required=True)
    p.add_argument('--protected-distance',type=float,default=0.);p.add_argument('--far-growth',type=float,default=1.)
    p.add_argument('--protect-surface',type=float,default=0.)
    p.add_argument('--maximum-elements',type=int,required=True,
                   help='frozen final-element budget used only by the generic far-field policy')
    contract=p.add_mutually_exclusive_group(required=True)
    contract.add_argument('--semantic-contract',type=Path)
    contract.add_argument('--simple-sharp-contract',action='store_true',
                          help='use the explicit historical 1/3100/5001/6001 compatibility contract')
    p.add_argument('--transformed-supports',type=Path,
                   help='required transformed source-support contract with --semantic-contract')
    p.add_argument('--seed-census',type=Path,
                   help='required seed corner census (simplified footprint polygons) with --semantic-contract')
    a=p.parse_args();m=read_mesh(a.mesh)
    if not (bool(a.semantic_contract) == bool(a.transformed_supports) == bool(a.seed_census)):
        p.error('--semantic-contract, --transformed-supports and --seed-census are required together')
    semantic=(load_semantic_contract(a.semantic_contract) if a.semantic_contract
              else simple_sharp_contract())
    supports=(json.loads(a.transformed_supports.read_text()) if a.transformed_supports else None)
    census=(json.loads(a.seed_census.read_text()) if a.seed_census else None)
    r=prepare(m,a.output,a.normal,a.tangent,a.far,
        a.protected_distance,a.far_growth,a.protect_surface,semantic,supports,
        a.maximum_elements,census,sha(a.semantic_contract) if a.semantic_contract else None)
    r['SeedArtifact']=str(a.mesh.resolve());r['SeedArtifactSHA256']=sha(a.mesh)
    if a.transformed_supports:
        r['TransformedSupportsArtifact']=str(a.transformed_supports.resolve())
        r['TransformedSupportsSHA256']=sha(a.transformed_supports)
        r['TransformedSupports']=supports
    if a.seed_census:
        r['SeedCensusArtifact']=str(a.seed_census.resolve());r['SeedCensusSHA256']=sha(a.seed_census)
    if a.mesh.suffix=='.toml':
        import tomllib
        data=tomllib.loads(a.mesh.read_text())
        r['SeedMesh']=data['SourceMesh'];r['SeedSHA256']=data['SourceSHA256']
        r['ArraySHA256']={name:sha(a.mesh.parent/name) for name in [data['NodeTags'],data['Coordinates']]+[c['Connectivity'] for c in data['Cells']]}
    else:
        r['SeedMesh']=str(a.mesh.resolve());r['SeedSHA256']=sha(a.mesh)
    (a.output/'recipe.json').write_text(json.dumps(r,indent=2)+'\n')

if __name__=='__main__':main()
