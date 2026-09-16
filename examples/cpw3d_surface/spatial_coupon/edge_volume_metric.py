#!/usr/bin/env python3
# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""Experimental straight-edge volume metric; no solver or production defaults.

The two transverse eigenvalues remain fine. Tangential recovery is limited by
true physical corners, not by mesh/CAD subdivision vertices or box cuts.
"""
import numpy as np

# Adjacent surface triangles are coplanar when the sine of the angle between
# their normals and the offset of their vertices from each other's plane,
# relative to the larger triangle diameter, are both within this dimensionless
# tolerance. It is a roundoff-scale bound: measured seed noise on exact planes is
# below 1e-8, and any physical dihedral is orders of magnitude above it. A label
# change between coplanar triangles is a reference boundary, not a geometric
# feature, so it must not become a ridge, pin, metric source, or protected band.
COPLANAR_TOLERANCE=1e-6


def plane_deviation(normal,point,scale,normals,points,scales):
    """Dimensionless deviation of the planes (normals, points) from (normal, point).

    The deviation is the larger of the sine of the angle between the normals
    (orientation-free) and the offset of every point from the reference plane
    relative to max(scale, scales, distance to the reference point): the offset
    is the angle under which the point leaves the plane as seen from the
    reference point, floored by the local size scale, so plane equivalence is
    independent of the coordinate origin and covariant under rigid motion.
    `points` is (m, 3) or (m, k, 3) (offset maximized over k); `scales` is a
    scalar or (m,).  Planes are equivalent when the result is at most
    COPLANAR_TOLERANCE.
    """
    normal=np.asarray(normal,dtype=float).reshape(3);point=np.asarray(point,dtype=float).reshape(3)
    normals=np.asarray(normals,dtype=float).reshape(-1,3)
    points=np.asarray(points,dtype=float)
    if points.ndim==2:points=points[:,None,:]
    if points.shape[0]!=len(normals) or points.shape[-1]!=3:
        raise ValueError('Plane points must be (m, 3) or (m, k, 3) matching the normals')
    scales=np.broadcast_to(np.asarray(scales,dtype=float),(len(normals),))
    if (not np.all(np.isfinite(normal)) or not np.all(np.isfinite(point)) or
            not np.isfinite(scale) or scale<=0 or not np.all(np.isfinite(normals)) or
            not np.all(np.isfinite(points)) or not np.all(np.isfinite(scales)) or
            np.any(scales<=0)):
        raise ValueError('Plane deviation requires finite planes and positive size scales')
    angle=np.linalg.norm(np.cross(normal,normals),axis=1)
    delta=points-point
    reach=np.maximum(np.maximum(scale,scales)[:,None],np.linalg.norm(delta,axis=2))
    return np.maximum(angle,np.max(np.abs(delta@normal)/reach,axis=1))


def cluster_coplanar_triangles(groups,normals,xyz,tolerance=COPLANAR_TOLERANCE):
    """First-fit clustering of triangles into equivalent planes within each group.

    Every triangle joins the earliest cluster of its group whose representative
    (the cluster's first triangle: unit normal, first vertex, diameter) it
    deviates from by at most `tolerance` under plane_deviation; otherwise it
    starts a cluster.  This is the single plane-equivalence rule for CAD-support
    clustering and every planar-patch audit key; rounding coordinates is not.
    Returns (representatives, patch): representatives is (k, 8) rows
    [group, nx, ny, nz, px, py, pz, diameter] in first-appearance order and
    patch the cluster index of every triangle.
    """
    groups=np.asarray(groups);normals=np.asarray(normals,dtype=float).reshape(-1,3)
    xyz=np.asarray(xyz,dtype=float).reshape(-1,3,3)
    if len(groups)!=len(normals) or len(xyz)!=len(normals):
        raise ValueError('Groups, normals and triangle vertices must have one row per triangle')
    if not np.all(np.isfinite(normals)) or not np.all(np.isfinite(xyz)):
        raise ValueError('Triangle planes must be finite')
    diameter=np.max(np.stack([np.linalg.norm(xyz[:,i]-xyz[:,j],axis=1)
                              for i,j in ((0,1),(1,2),(2,0))]),axis=0) if len(xyz) else np.zeros(0)
    if np.any(np.linalg.norm(np.cross(xyz[:,1]-xyz[:,0],xyz[:,2]-xyz[:,0]),axis=1)<=0):
        raise ValueError('Degenerate triangle in plane clustering')
    patch=np.full(len(groups),-1,dtype=np.int64);representatives=[]
    while True:
        remaining=np.flatnonzero(patch<0)
        if not len(remaining):break
        first=int(remaining[0])
        same=remaining[groups[remaining]==groups[first]]
        deviation=plane_deviation(normals[first],xyz[first,0],diameter[first],
                                  normals[same],xyz[same],diameter[same])
        patch[same[deviation<=tolerance]]=len(representatives)
        representatives.append([groups[first],*normals[first],*xyz[first,0],diameter[first]])
    return np.asarray(representatives,dtype=float).reshape(-1,8),patch


def match_equivalent_planes(reference,candidate,tolerance=COPLANAR_TOLERANCE):
    """One-to-one correspondence of cluster representatives by plane equivalence.

    Returns the candidate row index of every reference row when each reference
    plane is equivalent (same group, mutual plane_deviation <= tolerance) to
    exactly one candidate plane and vice versa; otherwise None.
    """
    reference=np.asarray(reference,dtype=float).reshape(-1,8)
    candidate=np.asarray(candidate,dtype=float).reshape(-1,8)
    if len(reference)!=len(candidate):return None
    if not len(reference):return []
    forward=np.zeros((len(reference),len(candidate)))
    backward=np.zeros((len(candidate),len(reference)))
    for i,row in enumerate(reference):
        forward[i]=plane_deviation(row[1:4],row[4:7],row[7],candidate[:,1:4],candidate[:,4:7],candidate[:,7])
    for j,row in enumerate(candidate):
        backward[j]=plane_deviation(row[1:4],row[4:7],row[7],reference[:,1:4],reference[:,4:7],reference[:,7])
    equivalent=((np.maximum(forward,backward.T)<=tolerance)&
                (reference[:,:1]==candidate[:,0][None,:]))
    if np.any(equivalent.sum(axis=1)!=1) or np.any(equivalent.sum(axis=0)!=1):return None
    return [int(np.flatnonzero(row)[0]) for row in equivalent]


def intersect_metrics(a,b):
    """Deterministic SPD intersection dominating both inputs in Loewner order."""
    a,b=np.asarray(a,dtype=float),np.asarray(b,dtype=float)
    if a.shape!=b.shape or a.ndim<2 or a.shape[-2:]!=(3,3):
        raise ValueError('Metric operands must have matching (..., 3, 3) shapes')
    for operand in (a,b):
        if not np.all(np.isfinite(operand)):
            raise ValueError('Nonfinite input metric')
        if not np.allclose(operand,operand.swapaxes(-1,-2),rtol=1e-12,atol=1e-14):
            raise ValueError('Nonsymmetric input metric')
        if not np.all(np.linalg.eigvalsh(operand)>0):
            raise ValueError('Nonpositive input metric')
    values,vectors=np.linalg.eigh(a)
    root=(vectors*np.sqrt(values)[...,None,:])@vectors.swapaxes(-1,-2)
    inverse=(vectors/np.sqrt(values)[...,None,:])@vectors.swapaxes(-1,-2)
    relative=inverse@b@inverse
    v,q=np.linalg.eigh((relative+relative.swapaxes(-1,-2))*.5)
    result=root@((q*np.maximum(v,1.)[...,None,:])@q.swapaxes(-1,-2))@root
    return (result+result.swapaxes(-1,-2))*.5


def segment_distances(points,segment):
    a,b=np.asarray(segment).reshape(2,3);v=b-a;length=np.linalg.norm(v)
    if not np.isfinite(length) or length<=0:raise ValueError('Invalid segment')
    tangent=v/length;delta=np.asarray(points)-a
    axial=np.clip(delta@tangent,0.,length)
    return np.linalg.norm(delta-axial[:,None]*tangent,axis=1),tangent


def volume_metric(points,segments,corners,normal_size,tangent_size,far_size,
                  radial_growth=1.,corner_growth=.25,protected_distance=0.,far_growth=None,
                  isotropic_corners=None,isotropy_radius=None):
    points=np.asarray(points,dtype=float);segments=np.asarray(segments,dtype=float).reshape(-1,6)
    corners=np.asarray(corners,dtype=float).reshape(-1,3)
    isotropic_corners=(corners if isotropic_corners is None else
                       np.asarray(isotropic_corners,dtype=float).reshape(-1,3))
    if isotropy_radius is None:isotropy_radius=0.
    if points.ndim!=2 or points.shape[1]!=3 or not np.all(np.isfinite(points)):
        raise ValueError('Expected finite 3D points')
    if far_growth is None:far_growth=radial_growth
    if not np.all(np.isfinite([normal_size,tangent_size,far_size,radial_growth,corner_growth,
                              protected_distance,far_growth,isotropy_radius])) or not (
            0<normal_size<=tangent_size<=far_size and radial_growth>0 and corner_growth>0 and
            protected_distance>=0 and far_growth>0 and isotropy_radius>=0):
        raise ValueError('Invalid metric controls')
    if (not len(segments) or not np.all(np.isfinite(segments)) or
            not np.all(np.isfinite(corners)) or not np.all(np.isfinite(isotropic_corners))):
        raise ValueError('Invalid physical feature geometry')
    dc=np.full(len(points),np.inf)
    for c in corners:dc=np.minimum(dc,np.linalg.norm(points-c,axis=1))
    metric=np.broadcast_to(np.eye(3)/far_size**2,(len(points),3,3)).copy()
    for segment in segments:  # Input feature order is part of the recorded recipe.
        r,t=segment_distances(points,segment)
        growth=radial_growth*np.minimum(r,protected_distance)+far_growth*np.maximum(0.,r-protected_distance)
        hn=np.minimum(far_size,normal_size+growth)
        ht=np.maximum(hn,np.minimum(np.minimum(far_size,tangent_size+growth),
                                    normal_size+corner_growth*dc))
        active=hn<far_size
        hn,ht=hn[active],ht[active]
        candidate=np.eye(3)[None,:,:]/hn[:,None,None]**2 + (
            1/ht**2-1/hn**2)[:,None,None]*np.outer(t,t)
        metric[active]=intersect_metrics(metric[active],candidate)
    # A semantic junction is isotropic over one tangential target, not only at
    # one metric node.  This gives MMG a physically scaled ball in which all
    # incident cells see the same SPD target.  CAD subdivisions and coupon-cut
    # endpoints are absent because the caller supplies contract junctions only.
    if isotropy_radius>0 and len(isotropic_corners):
        semantic_distance=np.full(len(points),np.inf)
        for corner in isotropic_corners:
            semantic_distance=np.minimum(semantic_distance,np.linalg.norm(points-corner,axis=1))
        active=semantic_distance<=isotropy_radius
        isotropic=np.broadcast_to(np.eye(3)/normal_size**2,(int(active.sum()),3,3))
        metric[active]=intersect_metrics(metric[active],isotropic)
    return metric


def _continues_straight(v,w):
    """Two edges leaving one node continue each other within COPLANAR_TOLERANCE."""
    v=np.asarray(v,dtype=float);w=np.asarray(w,dtype=float)
    return bool(np.dot(v,w)<0 and
                np.linalg.norm(np.cross(v,w))<=COPLANAR_TOLERANCE*np.linalg.norm(v)*np.linalg.norm(w))


def _graph_turns(points,edges):
    """Nodes of an edge graph that are junctions (degree != 2) or direction changes."""
    adjacency={}
    for a,b in np.asarray(edges,dtype=int).reshape(-1,2):
        adjacency.setdefault(int(a),[]).append(int(b));adjacency.setdefault(int(b),[]).append(int(a))
    return {node for node,neighbors in adjacency.items()
            if len(neighbors)!=2 or not _continues_straight(*(points[neighbors]-points[node]))}


def feature_chains(points,edges,lower,upper,tolerance=1e-8,cut_nodes=None):
    """Merge collinear subdivisions; return straight segments and true corners.

    Edges must already exclude boundary-cut/label-only features from the metric
    source graph. Box-cut endpoints can terminate a chain but do not cap h_t.
    """
    points=np.asarray(points);edges=np.asarray(edges,dtype=int).reshape(-1,2)
    adjacency={}
    for i,(a,b) in enumerate(edges):
        if a==b:raise ValueError('Degenerate feature edge')
        adjacency.setdefault(int(a),[]).append(i);adjacency.setdefault(int(b),[]).append(i)
    breaks=_graph_turns(points,edges)
    visited=set();segments=[]
    for start in sorted(breaks):
        for first_edge in adjacency[start]:
            if first_edge in visited:continue
            node=start;edge=first_edge
            while True:
                visited.add(edge)
                other=next(int(n) for n in edges[edge] if n!=node)
                if other in breaks:
                    a,b=sorted((tuple(points[start]),tuple(points[other])))
                    segments.append((*a,*b));break
                next_edges=[i for i in adjacency[other] if i!=edge]
                if len(next_edges)!=1 or next_edges[0] in visited:raise ValueError('Invalid or curved feature loop')
                node,edge=other,next_edges[0]
    if len(visited)!=len(edges):raise ValueError('Closed smooth/curved features unsupported')
    corners=[points[i] for i in sorted(breaks) if (i not in cut_nodes if cut_nodes is not None else
        not np.any((abs(points[i]-lower)<tolerance)|(abs(points[i]-upper)<tolerance)))]
    # Keep graph traversal order (original vertex identities), not lexicographic
    # coordinate order: rotating the geometry must not reorder SPD intersections.
    return np.asarray(segments),np.asarray(corners).reshape(-1,3)


def _shared_edges(points,triangles):
    """Sorted edge table of a triangle set: (pairs, owner, start, is_feature).

    `pairs[order]`/`owner` list every (edge, incident triangle) pair sorted by
    edge, `start` the first row of each distinct edge, and `is_feature` whether
    the incident triangles of that edge are not coplanar within
    COPLANAR_TOLERANCE (the sine of the normal angle is orientation-free and
    the plane offset is scaled by the local triangle size, so the test is
    dimensionless).
    """
    p=np.asarray(points);tri=np.asarray(triangles)
    xyz=p[tri];n=np.cross(xyz[:,1]-xyz[:,0],xyz[:,2]-xyz[:,0]);length=np.linalg.norm(n,axis=1)
    if np.any(length<=0):raise ValueError('Degenerate surface triangle')
    n/=length[:,None]
    diameter=np.max(np.stack([np.linalg.norm(xyz[:,i]-xyz[:,j],axis=1)
                              for i,j in ((0,1),(1,2),(2,0))]),axis=0)
    pairs=np.sort(tri[:,[(0,1),(1,2),(2,0)]].reshape(-1,2),axis=1)
    owner=np.repeat(np.arange(len(tri)),3)
    order=np.lexsort((pairs[:,1],pairs[:,0]));pairs=pairs[order];owner=owner[order]
    start=np.r_[0,np.flatnonzero(np.any(pairs[1:]!=pairs[:-1],axis=1))+1]
    first=np.repeat(owner[start],np.diff(np.r_[start,len(owner)]))
    normal_deviation=np.linalg.norm(np.cross(n[first],n[owner]),axis=1)
    plane_offset=np.max(abs(np.einsum('ij,ikj->ik',n[first],xyz[owner]-xyz[first][:,:1])),axis=1)
    plane_offset/=np.maximum(diameter[first],diameter[owner])
    deviation=np.maximum(normal_deviation,plane_offset)
    is_feature=np.maximum.reduceat(deviation,start)>COPLANAR_TOLERANCE
    return pairs,owner,start,is_feature


def junction_segments(points,triangles,references,cut_references,interface_references,
                      tolerance=1e-8):
    """Straight lines where the cut (Dirichlet) surface meets a material interface.

    A shared edge is a junction edge when it is a geometric feature (incident
    triangles not coplanar within COPLANAR_TOLERANCE) shared by at least one
    cut-surface triangle and one material-interface triangle: the trench floor
    and walls and the un-etched substrate-vacuum plane meeting the coupon box.
    Cut/cut box edges and cut/conductor edges are not junctions.  Collinear
    junction edges are merged into straight segments exactly like the physical
    feature graph; every junction node lies on the cut, so no corner results.
    Returns the (k, 6) segments in graph-traversal order.  Fails closed when the
    seed has no such edge: a coupon whose interfaces reach the box always has
    them, so their absence is a labeling or geometry error.
    """
    p=np.asarray(points);tri=np.asarray(triangles);refs=np.asarray(references)
    cut=set(cut_references);interface=set(interface_references)
    if not cut or not interface or cut&interface:
        raise ValueError('Cut-surface and material-interface labels must be disjoint nonempty sets')
    if not cut<=set(refs) or not interface<=set(refs):
        raise ValueError('Cut-surface and material-interface labels must be seed labels')
    pairs,owner,start,is_feature=_shared_edges(p,tri)
    on_cut=np.logical_or.reduceat(np.isin(refs[owner],list(cut)),start)
    on_interface=np.logical_or.reduceat(np.isin(refs[owner],list(interface)),start)
    edges=pairs[start][is_feature&on_cut&on_interface]
    if not len(edges):
        raise ValueError('The seed has no cut-surface/material-interface junction edge')
    nodes=set(map(int,edges.ravel()))
    segments,_=feature_chains(p,edges,p.min(axis=0),p.max(axis=0),tolerance,cut_nodes=nodes)
    return segments


def surface_features(points,triangles,references,cut_references,tolerance=1e-8):
    """Feature graph of a conforming piecewise-planar, reference-labeled complex.

    A shared edge is a geometric feature only when its incident triangles are
    not coplanar within COPLANAR_TOLERANCE, whatever their references: coplanar
    label seams and coplanar subdivisions never become features. Geometry
    preservation uses ALL such edges; metric sources exclude box cuts. Returns
    (features, pins, pin kinds, segments, corners); pins are turns/junctions of
    the complete reference-boundary graph, kind 'geometric-corner' when the
    feature graph alone pins them and 'reference-turn' otherwise.
    """
    p=np.asarray(points);tri=np.asarray(triangles);refs=np.asarray(references)
    cut_references=set(cut_references)
    if not cut_references or not cut_references <= set(refs):
        raise ValueError('Cut-surface references must be a nonempty label subset')
    pairs,owner,start,is_feature=_shared_edges(p,tri)
    on_matching=np.isin(refs[owner],list(cut_references))
    features=pairs[start][is_feature]
    # Use labeled matching support, not a global-axis bounding box. This remains
    # valid for a rigidly rotated coupon and identifies artificial cut endpoints.
    touches_matching=np.logical_or.reduceat(on_matching,start)[is_feature]
    physical=features[~touches_matching]
    cut_nodes=set(map(int,tri[np.isin(refs,list(cut_references))].ravel()))
    lower=p.min(axis=0);upper=p.max(axis=0)
    segments,corners=feature_chains(p,physical,lower,upper,tolerance,cut_nodes=cut_nodes)
    # Pin junctions/turns of the complete reference-boundary graph, including box
    # corners. MMG reconstructs every reference boundary (ridge or coplanar label
    # seam) as a discrete geometric curve, so with a tight Hausdorff bound an
    # unpinned seam bend is refined without limit. Seam turns are therefore
    # required vertices only; ridges, metric sources and protected bands stay
    # non-coplanar-only, and straight seam vertices remain free.
    geometric_corners=_graph_turns(p,features)
    reference_change=np.maximum.reduceat(refs[owner],start)!=np.minimum.reduceat(refs[owner],start)
    pins=np.array(sorted(_graph_turns(p,pairs[start][is_feature|reference_change])),dtype=int)
    kinds=np.array(['geometric-corner' if node in geometric_corners else 'reference-turn'
                    for node in pins])
    return features,pins,kinds,segments,corners
