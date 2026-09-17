#!/usr/bin/env python3
# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0
"""Bounded planar CAD-support restoration after an experimental MMG adaptation.

Native output is retained. Large displacement, missing support, inconsistent
intersection or inversion rejects the candidate. This does NOT replace footprint,
conductor-topology, achieved-resolution or PDE accuracy checks.
"""
import argparse,copy,json,hashlib
from pathlib import Path
import meshio
import numpy as np
from scipy.optimize import least_squares
from scipy.spatial import cKDTree

from edge_volume_metric import recipe_local_normal_size,segment_distances
from transform_coupon_source_contract import validate_rigid_transform

# A vertex saturated on the displacement ball is clamped to exactly the bound;
# recomputing its displacement from the moved coordinates can then exceed the
# bound by a few ulps. This relative allowance covers that floating-point
# roundoff only; it is not a geometric relaxation of the displacement bound.
DISPLACEMENT_ROUNDOFF_TOLERANCE=1e-12


def _within_displacement_bound(displacement,maximum_displacement):
    """maximum_displacement is a scalar or the per-vertex bound of each displacement."""
    displacement=np.asarray(displacement,dtype=float)
    return bool(np.all(np.isfinite(displacement)) and np.all(
        displacement<=np.asarray(maximum_displacement,dtype=float)*(1.+DISPLACEMENT_ROUNDOFF_TOLERANCE)))


def _bound_statistics(values):
    values=np.asarray(values,dtype=float).reshape(-1)
    if not len(values):return {'Count':0}
    return {'Count':int(len(values)),'MinimumUm':float(values.min()),'MedianUm':float(np.median(values)),
            'MaximumUm':float(values.max())}


def local_bound_size(points,recipe):
    """The local length scale of the restoration bounds at every vertex: the recipe's
    prescribed size (band/edge-layer law, corner balls) capped at NormalSize, so
    far-field and band vertices keep their NormalSize-based bounds while edge-layer
    vertices get EdgeSize-based ones."""
    normal=recipe.get('NormalSize')
    if not isinstance(normal,(int,float)) or not np.isfinite(normal) or normal<=0:
        raise ValueError('Restoration recipe lacks a valid NormalSize')
    return np.minimum(float(normal),recipe_local_normal_size(points,recipe))


def local_size_bounds(points,recipe,ratio):
    """Per-vertex displacement bounds ratio x local_bound_size."""
    if not np.isfinite(ratio) or ratio<=0:raise ValueError('Invalid displacement bound ratio')
    return ratio*local_bound_size(points,recipe)


def frozen_edge_layer_vertices(points,node_supports,recipe):
    """Supported (surface) vertices inside the seed edge layer footprint: the layer
    rows are frozen seed surface which the adapter preserved and the repair must
    not move at all.  The footprint is LayerThickness plus one EdgeSize of margin
    around every recorded span."""
    layer=recipe.get('EdgeLayer')
    if not isinstance(layer,dict):return frozenset()
    reach=float(layer['LayerThickness'])+float(layer['EdgeSize'])
    supported=np.fromiter(node_supports,dtype=int,count=len(node_supports))
    if not len(supported):return frozenset()
    distance=np.full(len(supported),np.inf)
    for span in np.asarray(layer['Spans'],dtype=float).reshape(-1,6):
        distance=np.minimum(distance,segment_distances(points[supported],span)[0])
    return frozenset(int(node) for node in supported[distance<=reach])


def _tetra_quality(points,tetrahedra):
    xyz=points[tetrahedra]
    jacobian=np.stack((xyz[:,1]-xyz[:,0],xyz[:,2]-xyz[:,0],
                       xyz[:,3]-xyz[:,0]),axis=2)
    determinant=np.linalg.det(jacobian)
    scaled=determinant/np.prod(np.linalg.norm(jacobian,axis=1),axis=1)
    singular=np.linalg.svd(jacobian,compute_uv=False)
    return scaled,singular[:,0]/singular[:,-1],determinant


def _movement_basis(node,node_supports,supports,fixed_nodes):
    if node in fixed_nodes:return np.empty((3,0))
    ids=node_supports.get(node)
    if not ids:return np.eye(3)
    normals=np.array([supports[i]['Normal'] for i in ids])
    _,singular,vectors=np.linalg.svd(normals)
    rank=int(np.sum(singular>1e-8))
    return vectors[rank:].T


def _bounded_offset(directions,parameters,maximum_displacement):
    """Map arbitrary basis parameters into the Euclidean displacement ball."""
    coordinates=np.asarray(parameters,dtype=float)
    length=float(np.linalg.norm(coordinates))
    if not np.isfinite(length):raise ValueError('Nonfinite quality-repair displacement')
    if length>1.:coordinates=coordinates/length
    offset=maximum_displacement*(directions@coordinates)
    magnitude=float(np.linalg.norm(offset))
    if magnitude>maximum_displacement:
        offset*=maximum_displacement/magnitude
        magnitude=float(np.linalg.norm(offset))
    if not _within_displacement_bound(magnitude,maximum_displacement):
        raise ValueError('Quality-repair displacement exceeds bound')
    return offset


def _transactional_quality_commit(points,candidate,active,tetrahedra,incident,floor):
    """Commit a repair move only when every original incident-cell floor survives."""
    scaled,_,determinant=_tetra_quality(candidate,tetrahedra[incident])
    degradation=float(np.min(scaled-floor))
    tolerance=max(1e-8,1e-6*float(np.max(np.abs(floor))))
    if np.any(determinant<=0) or degradation < -tolerance:
        return False,degradation,float(determinant.min())
    points[active]=candidate[active]
    return True,degradation,float(determinant.min())


def _pinned_vertices(points,recipe,tolerance=1e-10):
    """Vertices the adapter kept as required vertices, matched on native coordinates."""
    pinned=np.asarray([item['Point'] for item in recipe.get('PinnedVertices',[])],
                      dtype=float).reshape(-1,3)
    if not len(pinned):return frozenset()
    distance,index=cKDTree(points).query(pinned)
    if np.any(distance>tolerance):raise ValueError('Pinned vertex is absent after adaptation')
    return frozenset(map(int,index))


def _collapse_corner_ball_vertices(points,tetrahedra,tetrahedron_refs,node_supports,
                                   pinned_nodes,recipe,quality_target,local_size=None):
    """Collapse MMG-inserted free vertices below the adapter's minimum size in the
    semantic corner balls.

    MMG places free interior vertices at about NormalSize / 2 next to a required
    corner vertex, below the hmin it was given, and the resulting corner cells are
    beyond the bounded smoothing. A free vertex (no planar support, not pinned)
    inside CornerIsotropyRadius of a contract corner whose shortest incident edge
    is below the local prescribed size at the vertex (the recipe law: NormalSize
    in the band and the corner balls, EdgeSize-graded in the edge layer) is
    collapsed onto the corner when adjacent, otherwise onto its nearest non-free
    neighbor. A collapse stands only when every
    remapped cavity cell keeps a positive orientation and a scaled Jacobian of at
    least quality_target; a corner's collapses are committed together and rolled
    back if its corner-incident aspect did not improve. No vertex moves and no
    boundary triangle changes. Returns the compacted points, tetrahedra and
    references, the old-to-new vertex map (-1 for removed vertices) and the
    per-corner collapsed-vertex counts.
    """
    minimum_size=recipe.get('NormalSize')
    if not isinstance(minimum_size,(int,float)) or not np.isfinite(minimum_size) or minimum_size<=0:
        raise ValueError('Restoration recipe lacks a valid NormalSize')
    # Without a local size field the threshold is the recipe NormalSize everywhere.
    local_size=(np.full(len(points),float(minimum_size)) if local_size is None else
                np.asarray(local_size,dtype=float))
    if local_size.shape!=(len(points),) or not np.all(np.isfinite(local_size)) or np.any(local_size<=0):
        raise ValueError('Local size must be a positive finite value per vertex')
    corner_radius=recipe.get('CornerIsotropyRadius')
    if not isinstance(corner_radius,(int,float)) or not np.isfinite(corner_radius) or corner_radius<=0:
        raise ValueError('Restoration recipe lacks a valid corner isotropy radius')
    corners=np.asarray(recipe['TruePhysicalCorners'],dtype=float).reshape(-1,3)
    tetrahedra=np.array(tetrahedra,dtype=int,copy=True)
    alive=np.ones(len(tetrahedra),dtype=bool)
    removed=np.zeros(len(points),dtype=bool)
    free=np.ones(len(points),dtype=bool)
    free[list(node_supports)]=False
    free[list(pinned_nodes)]=False
    tree=cKDTree(points)
    collapsed=[]
    for corner in corners:
        corner_node=[int(i) for i in tree.query_ball_point(corner,1e-10)]
        if len(corner_node)!=1:raise ValueError('Semantic corner is absent during corner collapse')
        corner_node=corner_node[0]
        ball=np.asarray(sorted(int(i) for i in tree.query_ball_point(corner,corner_radius)),dtype=int)
        cells=np.flatnonzero(alive&np.any(np.isin(tetrahedra,ball),axis=1))
        # Working copies: the corner's collapses commit together or not at all.
        working=tetrahedra[cells].copy();working_alive=np.ones(len(cells),dtype=bool)
        working_removed=[]
        incident_before=np.any(working==corner_node,axis=1)
        aspect_before=float(_tetra_quality(points,working[incident_before])[1].max())
        candidates=[int(v) for v in ball if v!=corner_node and free[v]]
        candidates.sort(key=lambda v:float(np.linalg.norm(points[v]-corner)))
        for vertex in candidates:
            incident=np.flatnonzero(working_alive&np.any(working==vertex,axis=1))
            if not len(incident):continue
            neighbors=np.unique(working[incident]);neighbors=neighbors[neighbors!=vertex]
            lengths=np.linalg.norm(points[neighbors]-points[vertex],axis=1)
            if float(lengths.min())>=local_size[vertex]:continue
            # Targets in order: the corner itself, then the nearest non-free
            # neighbors, then the corner's own ring vertices at or beyond the
            # minimum size (free vertices that MMG placed correctly).
            targets=[]
            if corner_node in neighbors:targets.append(corner_node)
            order=np.argsort(lengths)
            targets+=[int(n) for n in neighbors[order] if not free[n] and n!=corner_node]
            ring=set(int(n) for n in np.unique(working[working_alive&np.any(working==corner_node,axis=1)]))
            targets+=[int(n) for n in neighbors[order] if free[n] and int(n) in ring and
                      np.linalg.norm(points[n]-corner)>=local_size[n]]
            # Among the valid cavities, commit the one with the best worst cell, and
            # only if it is no worse than the cells it replaces: a collapse onto
            # the corner itself can be valid yet leave a far worse sliver than a
            # collapse onto a nearby ring vertex.
            local_before=float(_tetra_quality(points,working[incident])[1].max())
            best=None
            for target in targets:
                keep=incident[~np.any(working[incident]==target,axis=1)]
                cavity=working[keep].copy();cavity[cavity==vertex]=target
                if not len(cavity):continue
                # A flattened cavity cell has a zero singular value; its aspect is
                # not used because the orientation test rejects it.
                with np.errstate(divide='ignore',invalid='ignore'):
                    scaled,aspects,determinant=_tetra_quality(points,cavity)
                if (np.any(determinant<=0) or float(scaled.min())<quality_target or
                        float(aspects.max())>local_before):continue
                score=(float(aspects.max()),-float(scaled.min()))
                if best is None or score<best[0]:best=(score,keep,cavity)
            if best is None:continue
            _,keep,cavity=best
            working[keep]=cavity
            working_alive[np.setdiff1d(incident,keep)]=False
            working_removed.append(vertex)
        incident_after=working_alive&np.any(working==corner_node,axis=1)
        aspect_after=float(_tetra_quality(points,working[incident_after])[1].max())
        if working_removed and aspect_after<aspect_before:
            tetrahedra[cells]=working;alive[cells]=working_alive
            removed[working_removed]=True
            collapsed.append({'Point':corner.tolist(),'CollapsedVertices':len(working_removed),
                              'AspectBefore':aspect_before,'AspectAfter':aspect_after})
        else:
            collapsed.append({'Point':corner.tolist(),'CollapsedVertices':0,
                              'AspectBefore':aspect_before,'AspectAfter':aspect_before,
                              'RolledBack':len(working_removed)})
    if np.any(np.isin(tetrahedra[alive],np.flatnonzero(removed))):
        raise ValueError('Collapsed vertex survived in the connectivity')
    vertex_map=np.full(len(points),-1,dtype=int)
    vertex_map[~removed]=np.arange(int(np.sum(~removed)))
    compact_points=points[~removed]
    compact_tetrahedra=vertex_map[tetrahedra[alive]]
    compact_refs=np.asarray(tetrahedron_refs)[alive]
    return compact_points,compact_tetrahedra,compact_refs,vertex_map,collapsed


def _quality_repair(points,tetrahedra,node_supports,supports,recipe,minimum_scaled,
                    maximum_corner_aspect,maximum_displacement,pinned_nodes=frozenset(),
                    frozen_nodes=frozenset()):
    """Constrained post-adaptation repair on frozen planar CAD supports.

    Interior vertices may move freely within a bounded ball (maximum_displacement
    is the per-vertex bound, ratio x the local prescribed size). Surface vertices
    remain on their exact support, ridge vertices remain on support
    intersections, and matching-surface, pinned and frozen edge-layer vertices are
    fixed.
    Semantic-corner neighborhoods are selected only from the frozen contract and
    repaired in alternating one-ring/two-ring passes over the vertices inside the
    corner ball; the best chained candidate is committed as one transaction only
    when it satisfies the displacement bound, the incident floors and the
    corner-aspect gate. A neighborhood without movable vertices is a recorded
    rejection, never an exception.
    """
    # A scalar bound applies to every vertex; the restoration passes the local one.
    maximum_displacement=np.broadcast_to(np.asarray(maximum_displacement,dtype=float),
                                         (len(points),)).copy()
    if not (np.isfinite(minimum_scaled) and 0<minimum_scaled<1 and
            np.isfinite(maximum_corner_aspect) and maximum_corner_aspect>1 and
            maximum_displacement.shape==(len(points),) and
            np.all(np.isfinite(maximum_displacement)) and np.all(maximum_displacement>0)):
        raise ValueError('Invalid post-adaptation quality controls')
    corner_radius=recipe.get('CornerIsotropyRadius')
    if not isinstance(corner_radius,(int,float)) or not np.isfinite(corner_radius) or corner_radius<=0:
        raise ValueError('Restoration recipe lacks a valid corner isotropy radius')
    semantic=recipe['SemanticContract']
    cut_roles=set(semantic['CutSurfaceRoles'])
    cut_attributes={item['Attribute'] for item in semantic['BoundaryLabels']
                    if item['Role'] in cut_roles}
    fixed_nodes={node for node,ids in node_supports.items()
                 if any(supports[i]['Attribute'] in cut_attributes for i in ids)}
    fixed_nodes|=set(pinned_nodes)
    fixed_nodes|=set(frozen_nodes)

    original=points.copy();largest_step=0.;corner_before=[];corner_after=[]
    quality_target=2.*minimum_scaled
    corner_target=.95*maximum_corner_aspect
    # Optimizer convergence leaves achieved aspects within roundoff of the
    # target; use the transactional commit's floor tolerance for that comparison.
    corner_target_tolerance=max(1e-8,1e-6*corner_target)
    original_scaled,_,original_determinant=_tetra_quality(original,tetrahedra)
    global_floor=np.minimum(original_scaled,quality_target)
    if np.any(original_determinant<=0):
        raise ValueError('Quality repair received an inverted tetrahedron')
    rejected_components=0;rejected_corners=0;target_missed_gate_satisfied=0

    def movable(active):
        active=np.asarray(active,dtype=int)
        bases=[_movement_basis(int(node),node_supports,supports,fixed_nodes)
               for node in active]
        selected=[i for i,value in enumerate(bases) if value.shape[1]]
        return active[selected],[bases[i] for i in selected]

    def solve(base,cells,active,objective):
        """Bounded least-squares move of the active vertices from a base state.

        Returns None when no active vertex can move; the caller records that as
        a rejection.
        """
        active,bases=movable(active)
        if not len(active):return None
        offsets=np.cumsum([0]+[value.shape[1] for value in bases])
        incident=np.flatnonzero(np.any(np.isin(tetrahedra,active),axis=1))
        # Floors come from the original global mesh, not a preceding pass.
        floor=global_floor[incident]
        # A parameter saturated on the ball by a preceding pass can re-derive a
        # few ulps outside the optimizer bounds; SciPy rejects such an initial
        # point as infeasible, so clip it onto the bounds (no geometric effect).
        initial=np.clip(np.concatenate([
            directions.T@(base[node]-original[node])/maximum_displacement[node]
            for node,directions in zip(active,bases)]),-.75,.75)
        def updated(value):
            candidate=base.copy()
            for i,(node,directions) in enumerate(zip(active,bases)):
                offset=_bounded_offset(
                    directions,value[offsets[i]:offsets[i+1]],maximum_displacement[node])
                candidate[node]=original[node]+offset
            return candidate
        def residual(value):
            candidate=updated(value)
            scaled,aspects,determinant=_tetra_quality(candidate,tetrahedra[incident])
            result=[100.*np.maximum(floor-scaled,0.),
                    1000.*np.maximum(-determinant,0.),.001*value]
            if objective=='corner':
                _,target_aspects,_=_tetra_quality(candidate,tetrahedra[cells])
                result.insert(0,10.*np.maximum(target_aspects-corner_target,0.))
            else:
                target_scaled,_,_=_tetra_quality(candidate,tetrahedra[cells])
                result.insert(0,100.*np.maximum(quality_target-target_scaled,0.))
            return np.concatenate(result)
        result=least_squares(residual,initial,bounds=(-.75,.75),max_nfev=3000,
                             ftol=1e-11,xtol=1e-11,gtol=1e-11)
        return updated(result.x)

    def commit(candidate,cells,objective):
        """Commit a candidate only within the bound, the floors and the corner gate."""
        nonlocal largest_step,rejected_components,rejected_corners,target_missed_gate_satisfied
        def rejected():
            # The shared point array is unchanged. Keep searching other
            # components; strict global gates remain final.
            nonlocal rejected_components,rejected_corners
            if objective=='corner':rejected_corners+=1
            else:rejected_components+=1
            return False
        if candidate is None:return rejected()
        moved=np.flatnonzero(np.any(candidate!=points,axis=1))
        if not len(moved):return rejected()
        cumulative=np.linalg.norm(candidate[moved]-original[moved],axis=1)
        if not _within_displacement_bound(cumulative,maximum_displacement[moved]):
            return rejected()
        if objective=='corner':
            # The corner-aspect gate, not the optimizer target, decides whether a
            # partial improvement may stand; 3.8 remains the objective margin.
            _,candidate_aspects,_=_tetra_quality(candidate,tetrahedra[cells])
            if float(candidate_aspects.max())>maximum_corner_aspect:return rejected()
            if float(candidate_aspects.max())>corner_target+corner_target_tolerance:
                target_missed_gate_satisfied+=1
        incident=np.flatnonzero(np.any(np.isin(tetrahedra,moved),axis=1))
        step=np.linalg.norm(candidate[moved]-points[moved],axis=1)
        accepted,_,_=_transactional_quality_commit(
            points,candidate,moved,tetrahedra,incident,global_floor[incident])
        if not accepted:return rejected()
        largest_step=max(largest_step,float(step.max()))
        return True

    corners=np.asarray(recipe['TruePhysicalCorners'],dtype=float).reshape(-1,3)
    corner_cells=[];corner_outcomes=[]
    for corner in corners:
        incident=np.flatnonzero(np.any(
            np.linalg.norm(points[tetrahedra]-corner,axis=2)<=1e-10,axis=1))
        corner_cells.append(incident)
        if not len(incident):raise ValueError('Semantic corner is absent during quality repair')
        before=float(_tetra_quality(points,tetrahedra[incident])[1].max())
        corner_before.append(before)
        if before<=corner_target:
            corner_outcomes.append({'Passes':0,'Outcome':'target-satisfied'});continue
        # Both rings are the vertices inside the recipe's corner ball, so the
        # two-ring pass is a superset of the one-ring pass.
        def inside_ball(nodes):
            distance=np.linalg.norm(points[nodes]-corner,axis=1)
            return nodes[(distance>1e-10)&(distance<=corner_radius)]
        one_ring=inside_ball(np.unique(tetrahedra[incident]))
        two_ring=inside_ball(np.unique(tetrahedra[np.any(np.isin(tetrahedra,one_ring),axis=1)]))
        # Alternate one-ring and two-ring passes on a chained candidate; only the
        # best candidate is committed, once, so no intermediate or regressed
        # state can stand alone.
        best=None;achieved=before;passes=0
        for active in (one_ring,two_ring,one_ring,two_ring):
            candidate=solve(points if best is None else best,incident,active,'corner')
            if candidate is None:break
            passes+=1
            improved=float(_tetra_quality(candidate,tetrahedra[incident])[1].max())
            converged=improved<=corner_target+corner_target_tolerance
            stalled=improved>=achieved-corner_target_tolerance
            if improved<achieved:best,achieved=candidate,improved
            if converged or stalled:break
        accepted=commit(best,incident,'corner')
        after=float(_tetra_quality(points,tetrahedra[incident])[1].max())
        corner_outcomes.append({'Passes':passes,'AchievedAspect':achieved,
            'Outcome':('rejected' if not accepted else
                       'target-reached' if after<=corner_target+corner_target_tolerance else
                       'target-missed-gate-satisfied')})

    scaled,_,_= _tetra_quality(points,tetrahedra)
    bad=set(map(int,np.flatnonzero(scaled<quality_target)))
    components=[]
    while bad:
        first=bad.pop();component={first};vertices=set(map(int,tetrahedra[first]))
        changed=True
        while changed:
            changed=False
            for cell in list(bad):
                if any(int(node) in vertices for node in tetrahedra[cell]):
                    bad.remove(cell);component.add(cell)
                    vertices.update(map(int,tetrahedra[cell]));changed=True
        components.append(sorted(component))
    for component in components:
        cells=np.asarray(component,dtype=int)
        commit(solve(points,cells,np.unique(tetrahedra[component]),'scaled'),cells,'scaled')
    final_scaled,_,final_determinant=_tetra_quality(points,tetrahedra)
    corner_after=[float(_tetra_quality(points,tetrahedra[incident])[1].max())
                  for incident in corner_cells]
    if any(value>maximum_corner_aspect for value in corner_after):
        raise ValueError(f'Semantic-corner quality repair failed: {max(corner_after)}')
    final_displacement=np.linalg.norm(points-original,axis=1)
    if not _within_displacement_bound(final_displacement,maximum_displacement):
        raise ValueError('Final quality-repair displacement exceeds bound')
    maximum_final=float(final_displacement.max())
    maximum_final_over_bound=float(np.max(final_displacement/maximum_displacement))
    frozen=np.fromiter(frozen_nodes,dtype=int,count=len(frozen_nodes))
    if frozen.size and np.any(points[frozen]!=original[frozen]):
        raise ValueError('Quality repair moved a frozen edge-layer vertex')
    if np.any(final_determinant<=0) or final_scaled.min()<minimum_scaled:
        raise ValueError(f'Post-adaptation minimum scaled Jacobian is {final_scaled.min()}')
    # Plane/intersection constraints are algebraic, but recheck explicitly before
    # publication so optimizer roundoff cannot silently move protected geometry.
    support_error=0.
    for node,ids in node_supports.items():
        for identifier in ids:
            item=supports[identifier]
            support_error=max(support_error,abs(float(
                np.dot(points[node],item['Normal'])-item['Offset'])))
    if support_error>1e-10:raise ValueError('Quality repair moved a protected support')
    fixed=np.fromiter(fixed_nodes,dtype=int,count=len(fixed_nodes))
    if fixed.size and np.any(points[fixed]!=original[fixed]):
        raise ValueError('Quality repair moved a fixed or pinned vertex')
    return {'QualityRepairVertices':int(np.sum(final_displacement>0)),
            'QualityRepairComponents':len(components),
            'RejectedQualityRepairComponents':rejected_components,
            'RejectedCornerRepairs':rejected_corners,
            'CornerRepairsTargetMissedGateSatisfied':target_missed_gate_satisfied,
            'CornerRepairOutcomes':corner_outcomes,
            'PinnedVerticesFixed':len(pinned_nodes),
            'MaximumQualityStepDisplacementUm':largest_step,
            'MaximumFinalQualityDisplacementUm':maximum_final,
            'MaximumFinalQualityDisplacementOverLocalBound':maximum_final_over_bound,
            'QualityDisplacementBoundUm':(float(maximum_displacement[0]) if np.all(
                maximum_displacement==maximum_displacement[0]) else _bound_statistics(maximum_displacement)),
            'FrozenEdgeLayerVertices':len(frozen_nodes),
            'MinimumScaledJacobianBefore':float(_tetra_quality(original,tetrahedra)[0].min()),
            'MinimumScaledJacobianAfter':float(final_scaled.min()),
            'CornerAspectsBefore':corner_before,'CornerAspectsAfter':corner_after,
            'MaximumSupportConstraintError':support_error}


def restore(mesh,recipe,maximum_displacement,minimum_scaled=None,
            maximum_corner_aspect=None,maximum_quality_displacement=None):
    """maximum_displacement / maximum_quality_displacement are absolute bounds in um
    or (ratio, 'local') pairs: ratio x the local prescribed size at every vertex."""
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
    points=mesh.points.copy();largest=0.;largest_over_bound=0.;moved=0
    def vertex_bounds(control):
        if isinstance(control,tuple):
            ratio,kind=control
            if kind!='local':raise ValueError('Unknown displacement bound kind')
            return local_size_bounds(points,recipe,ratio)
        if not np.isfinite(control) or control<=0:raise ValueError('Invalid displacement bound')
        return np.full(len(points),float(control))
    correction_bound=vertex_bounds(maximum_displacement)
    for ids,nodes in groups.items():
        a=np.array([supports[i]['Normal'] for i in ids]);b=np.array([supports[i]['Offset'] for i in ids])
        # Normals of coincident supports can differ at floating-point roundoff.
        # Do not interpret that noise as a real extra intersection constraint.
        original=points[nodes].copy();correction=(b-original@a.T)@np.linalg.pinv(a,rcond=1e-10).T
        magnitude=np.linalg.norm(correction,axis=1)
        largest=max(largest,float(magnitude.max()))
        over=magnitude/correction_bound[nodes]
        largest_over_bound=max(largest_over_bound,float(over.max()))
        if np.any(over>1.):
            worst=int(np.argmax(over))
            raise ValueError(f'CAD correction exceeds bound: {magnitude[worst]} > {correction_bound[nodes][worst]} '
                             f'(local size bound at {points[nodes][worst].tolist()})')
        points[nodes]+=correction;moved+=int(np.sum(np.any(correction!=0,axis=1)))
        if np.max(abs(points[nodes]@a.T-b))>1e-10:raise ValueError('Inconsistent planar intersection')
    tetrahedra=np.concatenate([block.data for block in mesh.cells if block.type=='tetra'])
    tetrahedron_refs=np.concatenate([np.asarray(ref) for block,ref in zip(mesh.cells,refs)
                                     if block.type=='tetra'])
    if np.any(_tetra_quality(points,tetrahedra)[2]<=0):
        raise ValueError('Projection inverted a tetrahedron')
    quality={};vertex_map=np.arange(len(points))
    controls=(minimum_scaled,maximum_corner_aspect,maximum_quality_displacement)
    if any(value is not None for value in controls):
        if any(value is None for value in controls):
            raise ValueError('All post-adaptation quality controls are required together')
        # Pins are matched on the native adapted coordinates, which the adapter
        # preserved exactly; projection may still correct them onto their supports.
        pinned_nodes=_pinned_vertices(mesh.points,recipe)
        if not (np.isfinite(minimum_scaled) and 0<minimum_scaled<1):
            raise ValueError('Invalid post-adaptation quality controls')
        local_size=local_bound_size(points,recipe)
        frozen_nodes=frozen_edge_layer_vertices(points,node_supports,recipe)
        corner_ball=np.zeros(len(points),dtype=bool)
        for corner in np.asarray(recipe['TruePhysicalCorners'],dtype=float).reshape(-1,3):
            corner_ball|=np.linalg.norm(points-corner,axis=1)<=float(recipe['CornerIsotropyRadius'])
        thresholds=_bound_statistics(local_size[corner_ball])
        points,tetrahedra,tetrahedron_refs,vertex_map,collapsed=_collapse_corner_ball_vertices(
            points,tetrahedra,tetrahedron_refs,node_supports,pinned_nodes,recipe,
            2.*minimum_scaled,local_size)
        node_supports={int(vertex_map[node]):ids for node,ids in node_supports.items()}
        pinned_nodes=frozenset(int(vertex_map[node]) for node in pinned_nodes)
        frozen_nodes=frozenset(int(vertex_map[node]) for node in frozen_nodes)
        if -1 in node_supports or -1 in pinned_nodes or -1 in frozen_nodes:
            raise ValueError('Corner collapse removed a supported, pinned or frozen vertex')
        repair_bound=vertex_bounds(maximum_quality_displacement)
        quality=_quality_repair(points,tetrahedra,node_supports,supports,recipe,
                                minimum_scaled,maximum_corner_aspect,
                                repair_bound,pinned_nodes=pinned_nodes,frozen_nodes=frozen_nodes)
        quality.update({'CornerBallCollapses':collapsed,
                        'CollapsedCornerVertices':int(sum(item['CollapsedVertices']
                                                          for item in collapsed)),
                        'CornerAspectsBeforeCollapse':[item['AspectBefore'] for item in collapsed],
                        'CollapseThresholdUm':thresholds,
                        'LocalSizeUm':_bound_statistics(local_size),
                        'BoundRule':'CAD correction, quality repair displacement and corner '
                                    'collapse thresholds are relative to the local prescribed '
                                    'size at each vertex (recipe band/edge-layer law and corner '
                                    'balls, capped at NormalSize), so only edge-layer vertices '
                                    'get tighter bounds; frozen edge-layer surface vertices '
                                    'are fixed'})
    cells=[];attributes=[]
    for block,ref in zip(mesh.cells,refs):
        if block.type!='triangle':continue
        triangles=vertex_map[block.data]
        if np.any(triangles<0):raise ValueError('Corner collapse removed a boundary vertex')
        attributes.append(np.array([supports[int(r)]['Attribute'] for r in ref]))
        cells.append(('triangle',triangles))
    attributes.append(tetrahedron_refs);cells.append(('tetra',tetrahedra))
    semantic=recipe.get('SemanticContract')
    if not isinstance(semantic,dict):raise ValueError('Restoration recipe lacks semantic contract')
    field_data={item['Material']:np.array([item['Attribute'],3],dtype=int)
                for item in semantic['VolumeMaterials']}
    field_data.update({('matching_surface' if item['Attribute']==1 else f"surface_{item['Attribute']}"):
                       np.array([item['Attribute'],2],dtype=int)
                       for item in semantic['BoundaryLabels']})
    output=meshio.Mesh(points,cells,cell_data={'gmsh:physical':attributes,'gmsh:geometrical':attributes},
                       field_data=field_data)
    supported=np.fromiter(node_supports,dtype=int,count=len(node_supports))
    return output,{'MaximumCorrectionUm':largest,'MaximumCorrectionOverLocalBound':largest_over_bound,
                   'CorrectionBoundUm':_bound_statistics(correction_bound[supported]),
                   'MovedBoundaryVertices':moved,'PlanarSupports':len(supports),**quality,
                   'LibraryQualified':False}


def restore_in_source_frame(mesh,recipe,maximum_displacement,minimum_scaled=None,
                            maximum_corner_aspect=None,maximum_quality_displacement=None):
    """Restore in the source-local frame and return local and published meshes.

    MMG may produce a different valid unstructured topology after a rigid source
    transform.  Pulling coordinates and planar constraints back before numerical
    optimization removes global-axis conditioning from restoration.  The same
    proper transform is reapplied only after all local quality/support checks.
    """
    semantic=recipe.get('SemanticContract',{})
    values=semantic.get('RigidTransform')
    if values is None:
        values=[1.,0.,0.,0.,0.,1.,0.,0.,0.,0.,1.,0.,0.,0.,0.,1.]
    matrix=np.asarray(validate_rigid_transform(values),dtype=float)
    rotation=matrix[:3,:3];translation=matrix[:3,3]
    local_mesh=copy.deepcopy(mesh)
    local_mesh.points=(np.asarray(mesh.points)-translation)@rotation
    local_recipe=copy.deepcopy(recipe)
    for name in ('TruePhysicalCorners','SurfaceFeatureCorners'):
        if name in local_recipe:
            points=np.asarray(local_recipe[name],dtype=float).reshape(-1,3)
            local_recipe[name]=((points-translation)@rotation).tolist()
    for item in local_recipe.get('PinnedVertices',[]):
        item['Point']=((np.asarray(item['Point'],dtype=float)-translation)@rotation).tolist()
    if 'PhysicalSegments' in local_recipe:
        segments=np.asarray(local_recipe['PhysicalSegments'],dtype=float).reshape(-1,2,3)
        local_recipe['PhysicalSegments']=((segments-translation)@rotation).tolist()
    if isinstance(local_recipe.get('JunctionSegments'),dict):
        segments=np.asarray(local_recipe['JunctionSegments']['Segments'],dtype=float).reshape(-1,2,3)
        local_recipe['JunctionSegments']['Segments']=((segments-translation)@rotation).reshape(-1,6).tolist()
    if isinstance(local_recipe.get('EdgeLayer'),dict):
        spans=np.asarray(local_recipe['EdgeLayer']['Spans'],dtype=float).reshape(-1,2,3)
        local_recipe['EdgeLayer']['Spans']=((spans-translation)@rotation).reshape(-1,6).tolist()
    for support in local_recipe['PlanarSupports'].values():
        global_normal=np.asarray(support['Normal'],dtype=float)
        support['Normal']=(rotation.T@global_normal).tolist()
        support['Offset']=float(support['Offset']-np.dot(global_normal,translation))
    local_semantic=local_recipe.get('SemanticContract',{})
    for name in ('SemanticCorners',):
        if name in local_semantic:
            points=np.asarray(local_semantic[name],dtype=float).reshape(-1,3)
            local_semantic[name]=((points-translation)@rotation).tolist()
    local_semantic['RigidTransform']=[1.,0.,0.,0.,0.,1.,0.,0.,0.,0.,1.,0.,0.,0.,0.,1.]
    restored,report=restore(local_mesh,local_recipe,maximum_displacement,minimum_scaled,
                            maximum_corner_aspect,maximum_quality_displacement)
    published=copy.deepcopy(restored)
    published.points=np.asarray(restored.points)@rotation.T+translation
    report.update({'RestorationFrame':'SourceLocal',
                   'RigidTransform':[float(value) for value in matrix.reshape(-1)]})
    return restored,published,report


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('input',type=Path);p.add_argument('recipe',type=Path);p.add_argument('output',type=Path)
    bound=p.add_mutually_exclusive_group(required=True)
    bound.add_argument('--max-displacement',type=float)
    bound.add_argument('--max-displacement-over-normal',type=float)
    p.add_argument('--minimum-scaled-jacobian',type=float)
    p.add_argument('--maximum-corner-aspect',type=float)
    p.add_argument('--maximum-quality-displacement-over-normal',type=float)
    p.add_argument('--source-local-output',type=Path,required=True)
    a=p.parse_args()
    if a.output.exists() or a.source_local_output.exists():raise ValueError('Do not overwrite candidates')
    recipe=json.loads(a.recipe.read_text())
    # The -over-normal ratios bound every vertex relative to its LOCAL prescribed
    # size (the recipe law evaluated at the vertex); an absolute bound stays global.
    maximum=(a.max_displacement if a.max_displacement is not None else
             (a.max_displacement_over_normal,'local'))
    if not isinstance(maximum,tuple) and (not np.isfinite(maximum) or maximum<=0):
        raise ValueError('Invalid displacement bound')
    if isinstance(maximum,tuple) and (not np.isfinite(maximum[0]) or maximum[0]<=0):
        raise ValueError('Invalid displacement bound')
    quality_displacement=(None if a.maximum_quality_displacement_over_normal is None else
                          (a.maximum_quality_displacement_over_normal,'local'))
    mesh=meshio.read(a.input);local_output,output,report=restore_in_source_frame(
        mesh,recipe,maximum,a.minimum_scaled_jacobian,a.maximum_corner_aspect,
        quality_displacement)
    if a.max_displacement_over_normal is not None:
        report['CorrectionBoundOverLocalSize']=a.max_displacement_over_normal
    if a.maximum_quality_displacement_over_normal is not None:
        report['QualityDisplacementBoundOverLocalSize']=a.maximum_quality_displacement_over_normal
    meshio.write(a.source_local_output,local_output,file_format='gmsh22',binary=True)
    meshio.write(a.output,output,file_format='gmsh22',binary=True)
    report['NativeInputSHA256']=hashlib.sha256(a.input.read_bytes()).hexdigest()
    report['SourceLocalOutputSHA256']=hashlib.sha256(a.source_local_output.read_bytes()).hexdigest()
    report['OutputSHA256']=hashlib.sha256(a.output.read_bytes()).hexdigest()
    a.output.with_suffix('.projection.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))

if __name__=='__main__':main()
