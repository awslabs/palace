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
import sys
from pathlib import Path
import meshio
import numpy as np
from edge_volume_metric import (COPLANAR_TOLERANCE,REQUIRED_TETRAHEDRA_RULE,cluster_coplanar_triangles,
                                corner_grading_record,
                                edge_layer_reach,edge_layer_required_reach,EDGE_LAYER_CELL_RULE,intersect_metrics,
                                junction_segments,required_tetrahedra,surface_features,volume_metric,
                                segment_distances)
from mesh_array_io import read_mesh,sha
from mesh_stage_contract import footprint_provenance,footprint_segments
from semantic_mesh_contract import (boundary_attributes, cut_surface_attributes,
                                    load_semantic_contract, material_interface_attributes,
                                    simple_sharp_contract, volume_attributes)
from trace_basis import (basis_statistics, cut_surface_size_report, load_trace_basis,
                         trace_basis_sizes, transform_trace_basis, unique_edges)


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


EDGE_LAYER_RULE=('the seed carries a geometric transverse layer on every face bounding a metal '
                 'edge (rows of size EdgeSize x GrowthRatio^(k-1) below NormalSize, tangentially '
                 'refined with the ridge grid so every layer cell has aspect <= Aspect: a '
                 'tetrahedron corner with three tangential edges has scaled Jacobian (hn/ht)^2, '
                 'which bounds the anisotropy under the scaled-Jacobian gate); the metric '
                 'prescribes the continuous form of the same layers, tangential size capped at '
                 'Aspect x hn(r) blending into the band law, '
                 'along the recorded spans, hn(r) = EdgeSize + (GrowthRatio - 1) r up to '
                 'Reach = (NormalSize - EdgeSize) / (GrowthRatio - 1) where it equals NormalSize, '
                 'then the ordinary band law (RadialGrowth to ProtectedDistance, FarGrowth beyond) '
                 'continues; the spans are intersected after the NormalSize band segments; the '
                 'seed rows lie within the frozen band (LayerThickness <= SurfaceProtectionRadius) '
                 'so MMG cannot alter the layer surface; the adapter hmin is EdgeSize')


def edge_layer_record(census,edge_size,growth_ratio,aspect,normal,protect_surface,
                      semantic_contract,band_segments,tolerance=1e-8):
    """Bind the seed's edge layer to the metric stage and record its spans.

    The bound census must record the same EdgeSize/GrowthRatio/NormalSize, at
    least one layer and one span, a layer thickness inside the surface protection
    radius (the frozen band must cover the whole layer footprint), and every span
    (placed by the contract's rigid transform) must lie on a band segment.
    Returns the recipe record; the spans are in the metric frame.
    """
    layer=census.get('EdgeLayer') if isinstance(census,dict) else None
    if not isinstance(layer,dict):raise ValueError('Seed census records no edge layer')
    if not np.isfinite(aspect) or aspect<1:raise ValueError('EdgeLayerAspect must be at least 1')
    for name,value in (('EdgeSize',edge_size),('GrowthRatio',growth_ratio),('Aspect',aspect),
                       ('NormalSize',normal)):
        if layer.get(name)!=value:
            raise ValueError(f'Seed census edge layer {name} differs from the metric stage')
    offsets=layer.get('RowOffsets');curves=layer.get('Curves')
    if (not isinstance(offsets,list) or not offsets or layer.get('Layers')!=len(offsets) or
            not isinstance(curves,list) or not curves or
            not all(isinstance(value,(int,float)) and not isinstance(value,bool) for value in offsets) or
            layer.get('LayerThickness')!=offsets[-1]):
        raise ValueError('Seed census edge layer rows are incomplete')
    if not offsets[-1]<=protect_surface:
        raise ValueError('Seed edge layer is thicker than the surface protection radius')
    matrix=np.asarray(semantic_contract.get('RigidTransform',
        [1.,0.,0.,0.,0.,1.,0.,0.,0.,0.,1.,0.,0.,0.,0.,1.]),dtype=float).reshape(4,4)
    spans=[];row_nodes=0
    band=np.asarray(band_segments,dtype=float).reshape(-1,6)
    for curve in curves:
        start=np.asarray(curve['Start'],dtype=float)@matrix[:3,:3].T+matrix[:3,3]
        end=np.asarray(curve['End'],dtype=float)@matrix[:3,:3].T+matrix[:3,3]
        if not (_near_segment(start,band,tolerance) and _near_segment(end,band,tolerance)):
            raise ValueError('Seed edge layer span does not lie on a band segment')
        spans.append([*start.tolist(),*end.tolist()])
        row_nodes+=sum(int(face['RowNodes']) for face in curve['Faces'])
    spans=np.asarray(spans,dtype=float).reshape(-1,6)
    total=float(np.linalg.norm(spans[:,3:]-spans[:,:3],axis=1).sum())
    if abs(total-float(layer['TotalSpanLength']))>COPLANAR_TOLERANCE*total:
        raise ValueError('Seed census edge layer span length is inconsistent')
    record={'EdgeSize':float(edge_size),'GrowthRatio':float(growth_ratio),'Aspect':float(aspect),
            'EdgeSizeOverNormalSize':float(edge_size/normal),
            'Reach':edge_layer_reach(normal,edge_size,growth_ratio),
            'Layers':len(offsets),'RowOffsets':[float(value) for value in offsets],
            'LayerThickness':float(offsets[-1]),'RowZigzag':layer.get('RowZigzag'),
            'SeedTangentialSize':layer.get('TangentialSize'),
            'SeedTangentialSubdivision':layer.get('TangentialSubdivision'),
            'SeedRowSubdivisions':layer.get('RowSubdivisions'),
            'SeedRowTangentialSpacings':layer.get('RowTangentialSpacings'),
            'SeedTaperSubdivisions':layer.get('TaperSubdivisions'),
            'SeedCornerTaperOffset':layer.get('CornerTaperOffset'),
            'SeedRidgeNodesAdded':layer.get('RidgeNodesAdded'),
            'Spans':spans.tolist(),'SpanCount':int(len(spans)),'TotalSpanLength':total,
            'SeedRows':layer.get('Rows'),'SeedRowNodes':row_nodes,
            'SeedCurves':[{'Curve':curve['Curve'],'Faces':curve['Faces'],
                           'SpanLength':curve['SpanLength'],'CurveLength':curve['CurveLength']}
                          for curve in curves],
            'MinimumSize':float(edge_size),'Rule':EDGE_LAYER_RULE}
    # The one recorded layer reach: a tetrahedron with a vertex within it belongs to
    # the layer (required region, restorer, audit statistics and layer quality rule).
    record['RequiredReach']=edge_layer_required_reach(record)
    record['LayerCellRule']=EDGE_LAYER_CELL_RULE
    return record


TRACE_BASIS_RULE=('cut-surface element size <= TraceBasisSizeRatio x the shortest edge of the '
                  'basis triangle containing the point (per-triangle rule: the hat of a basis '
                  'vertex varies linearly over the whole triangle, so its support is resolved '
                  'where it varies only when the whole triangle is discretized at that scale); '
                  'the seed applies it on the frozen cut surface, the metric applies '
                  'min(FarSize, size + FarGrowth x distance to the triangle) to the volume so '
                  'the existing far/grading law is kept away from narrow hats; the only new '
                  'parameter is the dimensionless ratio')


def bound_corner_grading(census,corner_size,growth_ratio,normal,radius):
    """Bind the seed's corner grading (census CornerGrading) to the metric stage:
    both prescribe the same CornerSize (0 / absent = uniform NormalSize) with the
    edge-layer growth ratio inside CornerIsotropyRadius.  Returns the recipe record
    or None."""
    grading=census.get('CornerGrading') if isinstance(census,dict) else None
    seed_size=float(grading.get('CornerSize',0.)) if isinstance(grading,dict) else 0.
    if corner_size is None:
        if seed_size>0:raise ValueError('Seed census records a corner grading the metric stage does not bind')
        return None
    if not np.isfinite(corner_size) or not 0<corner_size<normal:
        raise ValueError('CornerSize must be positive and below NormalSize')
    if (not isinstance(grading,dict) or seed_size!=corner_size or
            grading.get('GrowthRatio')!=growth_ratio or grading.get('NormalSize')!=normal or
            grading.get('Radius')!=radius):
        raise ValueError('Seed census corner grading differs from the metric stage')
    record=corner_grading_record(normal,corner_size,growth_ratio,radius)
    if grading.get('ShellRadii')!=record['ShellRadii'] or grading.get('Reach')!=record['Reach']:
        raise ValueError('Seed census corner grading shells differ from the metric law')
    return record


def trace_basis_sizing_record(basis, ratio, census, semantic_contract, exact_planes, points,
                              cut_triangles, far_size, growth, tolerance=1e-8):
    """Bind the trace basis to the seed and record the cut-surface size rule.

    The census must record the same basis (input digests, ratio, mesh-frame
    triangles) - the seed sized the frozen cut surface with it - and every basis
    vertex, placed by the contract's rigid transform, must lie on a cut-surface
    plane of the seed.  Returns (placed basis, recipe record).
    """
    if not np.isfinite(ratio) or ratio<=0:
        raise ValueError('TraceBasisSizeRatio must be a positive finite dimensionless number')
    record=census.get('TraceBasisSizing') if isinstance(census,dict) else None
    source_triangles=np.asarray(basis['Points'])[np.asarray(basis['Triangles'])]
    scale=float(np.max(np.asarray(basis['Upper'])-np.asarray(basis['Lower'])))
    if (not isinstance(record,dict) or record.get('Ratio')!=ratio or
            record.get('InputSHA256')!=basis['InputSHA256'] or
            not isinstance(record.get('MeshFrameTriangles'),list)):
        raise ValueError('Seed census trace basis sizing differs from the bound trace basis')
    recorded=np.asarray(record['MeshFrameTriangles'],dtype=float)
    if (recorded.shape!=source_triangles.shape or
            not np.allclose(recorded,source_triangles,rtol=0.,atol=tolerance*scale)):
        raise ValueError('Seed census trace basis triangles differ from the bound trace basis')
    placement=semantic_contract.get('RigidTransform',[1.,0.,0.,0.,0.,1.,0.,0.,0.,0.,1.,0.,0.,0.,0.,1.])
    placed=transform_trace_basis(basis,placement)
    cut=cut_surface_attributes(semantic_contract)
    planes=exact_planes[np.isin(exact_planes[:,0],list(cut))]
    if not len(planes):raise ValueError('Seed has no cut-surface planar support')
    offsets=np.abs(placed['Points']@planes[:,1:4].T-planes[:,4][None,:])
    if np.any(offsets.min(axis=1)>tolerance*scale):
        raise ValueError('Trace basis vertices do not lie on the seed cut-surface planes')
    statistics=basis_statistics(placed,ratio,far_size)
    report=cut_surface_size_report(points,cut_triangles,placed,ratio,far_size)
    return placed,{'Ratio':float(ratio),'RatioIsDimensionless':True,'Rule':TRACE_BASIS_RULE,
                   'Model':basis['Model'],'Frame':np.asarray(basis['Frame']).tolist(),
                   'InputSHA256':basis['InputSHA256'],
                   'Lower':np.asarray(basis['Lower']).tolist(),'Upper':np.asarray(basis['Upper']).tolist(),
                   **statistics,'FarSize':float(far_size),'Growth':float(growth),
                   'SeedMeshSizeMinimum':record.get('MeshSizeMinimum'),
                   'SeedGradingSlope':record.get('GradingSlope'),
                   'CutSurfaceSize':report,
                   'MeshFrameTriangles':source_triangles.tolist()}


def prepare(mesh,path,normal,tangent,far,protected_distance=0.,far_growth=1.,protect_surface=0.,
            semantic_contract=None, transformed_supports=None, maximum_elements=None,
            footprint_census=None, semantic_contract_sha256=None, trace_basis=None,
            trace_basis_size_ratio=1.0, edge_size=None, edge_growth_ratio=2.0,
            edge_layer_aspect=4.0, corner_size=None):
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
    census_layer=(footprint_census.get('EdgeLayer') if isinstance(footprint_census,dict) else None)
    if edge_size is None:
        if isinstance(census_layer,dict) and census_layer.get('EdgeSize',0)>0:
            raise ValueError('Seed census records an edge layer the metric stage does not bind')
        edge_layer=None
    else:
        if not np.isfinite(edge_size) or not 0<edge_size<normal:
            raise ValueError('EdgeSize must be positive and below NormalSize')
        edge_layer=edge_layer_record(footprint_census,edge_size,edge_growth_ratio,edge_layer_aspect,
                                     normal,protect_surface,semantic_contract,band_segments)
    layer_spans=None if edge_layer is None else np.asarray(edge_layer['Spans'],dtype=float)
    # Corner grading (decision 33): the seed's graded corner balls are bound and the
    # metric prescribes the same law inside CornerIsotropyRadius (= tangent).
    corner_grading=bound_corner_grading(footprint_census,corner_size,edge_growth_ratio,normal,tangent)
    # The far-field budget policy predicts the adapted load from the seed; the
    # seed cells inside the edge layer footprint (within the surface protection
    # radius of a span) are re-meshed to the recorded layer metric and are not
    # far-field load, so they are excluded.  Without a layer nothing changes.
    layer_cells=0
    if layer_spans is not None:
        centroids=mesh.points[tetrahedra].mean(axis=1)
        centroid_distance=np.full(len(tetrahedra),np.inf)
        for span in layer_spans:
            centroid_distance=np.minimum(centroid_distance,segment_distances(centroids,span)[0])
        layer_cells=int(np.sum(centroid_distance<=protect_surface))
    far_policy=budget_aware_far_policy(
        max(len(tetrahedra)-layer_cells,1),maximum_elements,far,far_growth)
    far_policy['SeedElementsInEdgeLayer']=layer_cells
    far_policy['SeedLoadRule']=('seed tetrahedra minus those whose centroid lies within '
                                'SurfaceProtectionRadius of an edge-layer span')
    effective_far=far_policy['EffectiveFarSize']
    effective_far_growth=far_policy['EffectiveFarGrowth']
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
    cut_mask=np.isin(triangle_refs,list(cut_surface_attributes(semantic_contract)))
    if trace_basis is None:
        if isinstance(footprint_census,dict) and footprint_census.get('TraceBasisSizing') is not None:
            raise ValueError('Seed census records a trace basis the metric stage does not bind')
        placed_basis,trace_record=None,None
    else:
        placed_basis,trace_record=trace_basis_sizing_record(
            trace_basis,trace_basis_size_ratio,footprint_census,semantic_contract,exact_planes,
            mesh.points,triangles[cut_mask],effective_far,effective_far_growth)
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
    fixed=cut_mask.copy()
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
    # The seed cells inside the corner balls and the edge layer are MMG required
    # tetrahedra: the seed defines the near-corner/near-edge discretization and MMG
    # adapts only outside them (supervisor decision 30).
    layer_reach=None if edge_layer is None else edge_layer_required_reach(edge_layer)
    required,required_per_corner,required_per_span=required_tetrahedra(
        mesh.points,tetrahedra,semantic_corners,tangent,layer_spans,layer_reach)
    required_indices=np.flatnonzero(required)+1
    if not len(required_indices):raise ValueError('No seed tetrahedron lies inside a corner ball')
    np.savetxt(path/'required-tetrahedra.txt',required_indices,fmt='%d')
    required_record={'Count':int(len(required_indices)),'Rule':REQUIRED_TETRAHEDRA_RULE,
        'CornerRadius':float(tangent),
        'PerCorner':[{'Point':corner.tolist(),'Tetrahedra':int(count)}
                     for corner,count in zip(semantic_corners,required_per_corner)],
        'LayerRequiredReach':layer_reach,
        'PerSpan':[{'Span':span.tolist(),'Tetrahedra':int(count)}
                   for span,count in zip(([] if layer_spans is None else layer_spans),required_per_span)],
        'File':'required-tetrahedra.txt','IndexBase':1}
    with (path/'metric.f64').open('wb') as f:
        for start in range(0,len(mesh.points),100000):
            metric=volume_metric(mesh.points[start:start+100000],band_segments,corners,normal,tangent,effective_far,
                                 protected_distance=protected_distance,far_growth=effective_far_growth,
                                 isotropic_corners=semantic_corners,isotropy_radius=tangent,
                                 edge_layer_segments=layer_spans,edge_size=edge_size,
                                 growth_ratio=edge_growth_ratio,edge_layer_aspect=edge_layer_aspect,
                                 corner_size=corner_size)
            if placed_basis is not None:
                # Trace rule blended into the far/grading law: an isotropic cap where
                # a narrow basis triangle is near; the far size elsewhere (no change).
                sizes=trace_basis_sizes(mesh.points[start:start+100000],placed_basis,
                                        trace_basis_size_ratio,effective_far,effective_far_growth)
                active=sizes<effective_far
                cap=np.eye(3)[None,:,:]/sizes[active,None,None]**2
                metric[active]=intersect_metrics(metric[active],cap)
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
            'RequiredTetrahedra':required_record,
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
            'BandSegmentOrder':('PhysicalSegments then JunctionSegments' if edge_layer is None else
                                'PhysicalSegments then JunctionSegments then EdgeLayer spans'),
            'SurfaceFeatureCorners':corners.tolist(),
            'PlanarSupports':{str(10000+i):{'Attribute':int(row[0]),'Normal':row[1:4].tolist(),'Offset':float(row[4])} for i,row in enumerate(exact_planes)},
            'PlanarSupportEquivalence':{'Tolerance':COPLANAR_TOLERANCE,
                'Rule':'edge_volume_metric.cluster_coplanar_triangles: same attribute, sine of '
                       'the normal angle and point offset relative to max(local size, distance) '
                       'within the tolerance; shared with the planar-patch audits'},
            'SemanticContract':semantic_contract,'LibraryQualified':False}
    if footprint is not None:recipe['FootprintSegments']=footprint
    if edge_layer is not None:recipe['EdgeLayer']=edge_layer
    if corner_grading is not None:recipe['CornerGrading']=corner_grading
    if trace_record is not None:
        recipe['TraceBasisSizing']=trace_record
        # The bound basis edges in the metric frame: a band lying on one of them is
        # source-driven for the trace-diagonal audit (supervisor decision 21).
        edges=unique_edges(placed_basis['Points'],placed_basis['Triangles'])
        recipe['TraceBasisEdges']={'InputSHA256':trace_basis['InputSHA256'],'Count':int(len(edges)),
            'Segments':edges.reshape(-1,6).tolist(),
            'Rule':'unique edges of the bound trace basis (metric frame); the trace-diagonal audit '
                   'treats a line-like band as source-driven only when it lies on one of them '
                   '(direction aligned and both endpoints within 2 x ShortEdgeThreshold of the '
                   'edge segment) and reports such bands separately'}
    (path/'recipe.json').write_text(json.dumps(recipe,indent=2)+'\n')
    print(json.dumps({k:v for k,v in recipe.items() if k not in ('PhysicalSegments','TruePhysicalCorners','FootprintSegments','JunctionSegments','TraceBasisSizing','EdgeLayer','RequiredTetrahedra')},indent=2))
    print(json.dumps({'RequiredTetrahedra':{k:v for k,v in required_record.items() if k not in ('PerSpan','Rule')}},indent=2))
    if edge_layer is not None:
        print(json.dumps({'EdgeLayer':{k:v for k,v in edge_layer.items() if k not in ('Spans','SeedCurves','Rule')}},indent=2))
    if corner_grading is not None:
        print(json.dumps({'CornerGrading':{k:v for k,v in corner_grading.items() if k!='Rule'}},indent=2))
    print(json.dumps({'JunctionSegments':{k:v for k,v in junction_record.items() if k!='Segments'}},indent=2))
    if trace_record is not None:
        print(json.dumps({'TraceBasisSizing':{k:(v if k!='CutSurfaceSize' else
            {kk:vv for kk,vv in v.items() if kk!='BasisTrianglesBelowFarSize'})
            for k,v in trace_record.items() if k not in ('MeshFrameTriangles','RequestedSizes','Rule')}},indent=2))
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
    p.add_argument('--trace-basis-contract',type=Path,help='bound trace basis contract (basis-contract.json)')
    p.add_argument('--trace-vertices',type=Path,help='bound trace basis vertices (trace-vertices.csv)')
    p.add_argument('--trace-triangles',type=Path,help='bound trace basis triangles (trace-triangles.csv)')
    p.add_argument('--process-library',type=Path,help='bound process library (frame of the trace basis)')
    p.add_argument('--edge-size',type=float,
                   help='EdgeSize of the seed edge layer to bind and prescribe (below --normal); '
                        'requires the seed census edge layer; the adapter hmin must equal it')
    p.add_argument('--edge-growth-ratio',type=float,default=2.0,
                   help='GrowthRatio of the geometric edge layer (must equal the seed value)')
    p.add_argument('--edge-layer-aspect',type=float,default=4.0,
                   help='EdgeLayerAspect: tangential/normal size cap inside the layer (must equal '
                        'the seed value; the scaled Jacobian of a layer cell scales as (hn/ht)^2)')
    p.add_argument('--trace-basis-size-ratio',type=float,default=1.0,
                   help='dimensionless TraceBasisSizeRatio of the cut-surface size rule (default 1.0: '
                        'at least one element per basis edge); only with the four trace basis inputs')
    p.add_argument('--corner-size',type=float,
                   help='CornerSize of the seed corner grading to bind and prescribe (below --normal, '
                        'growing by --edge-growth-ratio to --normal inside the corner balls); requires '
                        'the seed census corner grading; absent = uniform NormalSize balls')
    a=p.parse_args();m=read_mesh(a.mesh)
    if not (bool(a.semantic_contract) == bool(a.transformed_supports) == bool(a.seed_census)):
        p.error('--semantic-contract, --transformed-supports and --seed-census are required together')
    basis_paths=(a.trace_basis_contract,a.trace_vertices,a.trace_triangles,a.process_library)
    if any(path is not None for path in basis_paths) and (
            any(path is None for path in basis_paths) or not a.semantic_contract):
        p.error('--trace-basis-contract, --trace-vertices, --trace-triangles and --process-library '
                'are required together with --semantic-contract')
    if '--trace-basis-size-ratio' in sys.argv and basis_paths[0] is None:
        p.error('--trace-basis-size-ratio requires the bound trace basis inputs')
    if a.edge_size is not None and not a.seed_census:
        p.error('--edge-size requires the bound seed census')
    if any(option in sys.argv for option in ('--edge-growth-ratio','--edge-layer-aspect')) and a.edge_size is None:
        p.error('--edge-growth-ratio and --edge-layer-aspect require --edge-size')
    if a.corner_size is not None and not a.seed_census:
        p.error('--corner-size requires the bound seed census')
    basis=(load_trace_basis(*basis_paths) if basis_paths[0] is not None else None)
    semantic=(load_semantic_contract(a.semantic_contract) if a.semantic_contract
              else simple_sharp_contract())
    supports=(json.loads(a.transformed_supports.read_text()) if a.transformed_supports else None)
    census=(json.loads(a.seed_census.read_text()) if a.seed_census else None)
    r=prepare(m,a.output,a.normal,a.tangent,a.far,
        a.protected_distance,a.far_growth,a.protect_surface,semantic,supports,
        a.maximum_elements,census,sha(a.semantic_contract) if a.semantic_contract else None,
        basis,a.trace_basis_size_ratio,a.edge_size,a.edge_growth_ratio,a.edge_layer_aspect,
        a.corner_size)
    r['SeedArtifact']=str(a.mesh.resolve());r['SeedArtifactSHA256']=sha(a.mesh)
    if basis is not None:
        r['TraceBasisSizing']['Inputs']={name:{'Path':str(path.resolve()),'SHA256':sha(path)}
            for name,path in zip(('BasisContract','TraceVertices','TraceTriangles','ProcessLibrary'),basis_paths)}
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
