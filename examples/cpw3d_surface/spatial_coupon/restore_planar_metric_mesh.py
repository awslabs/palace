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

from transform_coupon_source_contract import validate_rigid_transform


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
    if not np.isfinite(magnitude) or magnitude>maximum_displacement:
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


def _quality_repair(points,tetrahedra,node_supports,supports,recipe,minimum_scaled,
                    maximum_corner_aspect,maximum_displacement):
    """Constrained post-adaptation repair on frozen planar CAD supports.

    Interior vertices may move freely within a bounded ball. Surface vertices
    remain on their exact support, ridge vertices remain on support
    intersections, and matching-surface vertices are fixed. Semantic-corner
    neighborhoods are selected only from the frozen contract.
    """
    if not (np.isfinite(minimum_scaled) and 0<minimum_scaled<1 and
            np.isfinite(maximum_corner_aspect) and maximum_corner_aspect>1 and
            np.isfinite(maximum_displacement) and maximum_displacement>0):
        raise ValueError('Invalid post-adaptation quality controls')
    semantic=recipe['SemanticContract']
    cut_roles=set(semantic['CutSurfaceRoles'])
    cut_attributes={item['Attribute'] for item in semantic['BoundaryLabels']
                    if item['Role'] in cut_roles}
    fixed_nodes={node for node,ids in node_supports.items()
                 if any(supports[i]['Attribute'] in cut_attributes for i in ids)}

    original=points.copy();largest_step=0.;corner_before=[];corner_after=[]
    quality_target=2.*minimum_scaled
    corner_target=.95*maximum_corner_aspect
    original_scaled,_,original_determinant=_tetra_quality(original,tetrahedra)
    global_floor=np.minimum(original_scaled,quality_target)
    if np.any(original_determinant<=0):
        raise ValueError('Quality repair received an inverted tetrahedron')
    rejected_components=0

    def optimize(cells,active,objective):
        nonlocal largest_step,rejected_components
        active=np.asarray(active,dtype=int)
        bases=[_movement_basis(int(node),node_supports,supports,fixed_nodes)
               for node in active]
        selected=[i for i,value in enumerate(bases) if value.shape[1]]
        active=active[selected];bases=[bases[i] for i in selected]
        if not len(active):raise ValueError('Quality repair has no movable vertices')
        offsets=np.cumsum([0]+[value.shape[1] for value in bases])
        incident=np.flatnonzero(np.any(np.isin(tetrahedra,active),axis=1))
        # Floors come from the original global mesh, not a preceding pass.
        floor=global_floor[incident]
        local=points[active].copy()
        initial=np.concatenate([
            directions.T@(local[i]-original[node])/maximum_displacement
            for i,(node,directions) in enumerate(zip(active,bases))])
        def updated(value):
            candidate=points.copy()
            for i,(node,directions) in enumerate(zip(active,bases)):
                offset=_bounded_offset(
                    directions,value[offsets[i]:offsets[i+1]],maximum_displacement)
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
        candidate=updated(result.x)
        step=np.linalg.norm(candidate[active]-local,axis=1)
        cumulative=np.linalg.norm(candidate[active]-original[active],axis=1)
        if np.any(~np.isfinite(cumulative)) or np.any(cumulative>maximum_displacement):
            raise ValueError('Final quality-repair displacement exceeds bound')
        accepted,_,_=_transactional_quality_commit(
            points,candidate,active,tetrahedra,incident,floor)
        if not accepted:
            # The shared point array is unchanged. Keep searching other
            # components; strict global gates remain final.
            rejected_components+=1
            selected_scaled,selected_aspects,_=_tetra_quality(points,tetrahedra[cells])
            return float(selected_scaled.min()),float(selected_aspects.max()),False
        largest_step=max(largest_step,float(step.max()))
        selected_scaled,selected_aspects,_=_tetra_quality(points,tetrahedra[cells])
        return float(selected_scaled.min()),float(selected_aspects.max()),True

    corners=np.asarray(recipe['TruePhysicalCorners'],dtype=float).reshape(-1,3)
    corner_cells=[]
    for corner in corners:
        incident=np.flatnonzero(np.any(
            np.linalg.norm(points[tetrahedra]-corner,axis=2)<=1e-10,axis=1))
        corner_cells.append(incident)
        if not len(incident):raise ValueError('Semantic corner is absent during quality repair')
        before=_tetra_quality(points,tetrahedra[incident])[1].max()
        corner_before.append(float(before))
        if before>corner_target:
            vertices=np.unique(tetrahedra[incident])
            vertices=vertices[np.linalg.norm(points[vertices]-corner,axis=1)>1e-10]
            _,after,_=optimize(incident,vertices,'corner')
        else:after=float(before)
        corner_after.append(float(after))

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
        optimize(np.asarray(component,dtype=int),np.unique(tetrahedra[component]),'scaled')
    final_scaled,_,final_determinant=_tetra_quality(points,tetrahedra)
    corner_after=[float(_tetra_quality(points,tetrahedra[incident])[1].max())
                  for incident in corner_cells]
    if any(value>maximum_corner_aspect for value in corner_after):
        raise ValueError(f'Semantic-corner quality repair failed: {max(corner_after)}')
    final_displacement=np.linalg.norm(points-original,axis=1)
    if np.any(~np.isfinite(final_displacement)) or np.any(
            final_displacement>maximum_displacement):
        raise ValueError('Final quality-repair displacement exceeds bound')
    maximum_final=float(final_displacement.max())
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
    return {'QualityRepairVertices':int(np.sum(final_displacement>0)),
            'QualityRepairComponents':len(components),
            'RejectedQualityRepairComponents':rejected_components,
            'MaximumQualityStepDisplacementUm':largest_step,
            'MaximumFinalQualityDisplacementUm':maximum_final,
            'QualityDisplacementBoundUm':maximum_displacement,
            'MinimumScaledJacobianBefore':float(_tetra_quality(original,tetrahedra)[0].min()),
            'MinimumScaledJacobianAfter':float(final_scaled.min()),
            'CornerAspectsBefore':corner_before,'CornerAspectsAfter':corner_after,
            'MaximumSupportConstraintError':support_error}


def restore(mesh,recipe,maximum_displacement,minimum_scaled=None,
            maximum_corner_aspect=None,maximum_quality_displacement=None):
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
    tetrahedra=np.concatenate([block.data for block in mesh.cells if block.type=='tetra'])
    if np.any(_tetra_quality(points,tetrahedra)[2]<=0):
        raise ValueError('Projection inverted a tetrahedron')
    quality={}
    controls=(minimum_scaled,maximum_corner_aspect,maximum_quality_displacement)
    if any(value is not None for value in controls):
        if any(value is None for value in controls):
            raise ValueError('All post-adaptation quality controls are required together')
        quality=_quality_repair(points,tetrahedra,node_supports,supports,recipe,
                                minimum_scaled,maximum_corner_aspect,
                                maximum_quality_displacement)
    cells=[];attributes=[]
    for block,ref in zip(mesh.cells,refs):
        if block.type not in ('tetra','triangle'):continue
        attributes.append(ref if block.type=='tetra' else
                          np.array([supports[int(r)]['Attribute'] for r in ref]))
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
    if 'PhysicalSegments' in local_recipe:
        segments=np.asarray(local_recipe['PhysicalSegments'],dtype=float).reshape(-1,2,3)
        local_recipe['PhysicalSegments']=((segments-translation)@rotation).tolist()
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
    maximum=(a.max_displacement if a.max_displacement is not None else
             a.max_displacement_over_normal*float(recipe['NormalSize']))
    if not np.isfinite(maximum) or maximum<=0:raise ValueError('Invalid displacement bound')
    quality_displacement=(None if a.maximum_quality_displacement_over_normal is None else
                          a.maximum_quality_displacement_over_normal*float(recipe['NormalSize']))
    mesh=meshio.read(a.input);local_output,output,report=restore_in_source_frame(
        mesh,recipe,maximum,a.minimum_scaled_jacobian,a.maximum_corner_aspect,
        quality_displacement)
    if a.max_displacement_over_normal is not None:
        report['CorrectionBoundOverNormalSize']=a.max_displacement_over_normal
    meshio.write(a.source_local_output,local_output,file_format='gmsh22',binary=True)
    meshio.write(a.output,output,file_format='gmsh22',binary=True)
    report['NativeInputSHA256']=hashlib.sha256(a.input.read_bytes()).hexdigest()
    report['SourceLocalOutputSHA256']=hashlib.sha256(a.source_local_output.read_bytes()).hexdigest()
    report['OutputSHA256']=hashlib.sha256(a.output.read_bytes()).hexdigest()
    a.output.with_suffix('.projection.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))

if __name__=='__main__':main()
