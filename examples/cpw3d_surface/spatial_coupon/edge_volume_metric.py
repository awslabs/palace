#!/usr/bin/env python3
# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""Experimental straight-edge volume metric; no solver or production defaults.

The two transverse eigenvalues remain fine. Tangential recovery is limited by
true physical corners, not by mesh/CAD subdivision vertices or box cuts.
"""
import numpy as np


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
                  radial_growth=1.,corner_growth=.25,protected_distance=0.,far_growth=None):
    points=np.asarray(points,dtype=float);segments=np.asarray(segments,dtype=float).reshape(-1,6)
    corners=np.asarray(corners,dtype=float).reshape(-1,3)
    if points.ndim!=2 or points.shape[1]!=3 or not np.all(np.isfinite(points)):
        raise ValueError('Expected finite 3D points')
    if far_growth is None:far_growth=radial_growth
    if not np.all(np.isfinite([normal_size,tangent_size,far_size,radial_growth,corner_growth,
                              protected_distance,far_growth])) or not (
            0<normal_size<=tangent_size<=far_size and radial_growth>0 and corner_growth>0 and
            protected_distance>=0 and far_growth>0):
        raise ValueError('Invalid metric controls')
    if not len(segments) or not np.all(np.isfinite(segments)) or not np.all(np.isfinite(corners)):
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
    return metric


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
    breaks=set()
    for node,inc in adjacency.items():
        if len(inc)!=2:breaks.add(node);continue
        ends=[next(int(n) for n in edges[i] if n!=node) for i in inc]
        v,w=points[ends]-points[node]
        if np.dot(v,w)/np.linalg.norm(v)/np.linalg.norm(w)>-1+1e-10:breaks.add(node)
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


def surface_features(points,triangles,references,tolerance=1e-8):
    """Feature graph of a conforming piecewise-planar, reference-labeled complex.

    Geometry preservation uses ALL patch junctions. Metric sources exclude box
    cuts; same-reference coplanar subdivisions never become features.
    """
    p=np.asarray(points);tri=np.asarray(triangles);refs=np.asarray(references)
    xyz=p[tri];n=np.cross(xyz[:,1]-xyz[:,0],xyz[:,2]-xyz[:,0]);length=np.linalg.norm(n,axis=1)
    if np.any(length<=0):raise ValueError('Degenerate surface triangle')
    n/=length[:,None]
    pivot=np.argmax(abs(n),axis=1);sign=np.sign(n[np.arange(len(n)),pivot]);n*=sign[:,None]
    planes=np.column_stack((refs,np.round(n,8),np.round(np.einsum('ij,ij->i',n,xyz[:,0]),8)))
    _,patch=np.unique(planes,axis=0,return_inverse=True)
    pairs=np.sort(tri[:,[(0,1),(1,2),(2,0)]].reshape(-1,2),axis=1)
    owner=np.repeat(patch,3)
    order=np.lexsort((pairs[:,1],pairs[:,0]));pairs=pairs[order];owner=owner[order]
    on_matching=np.repeat(refs==1,3)[order]
    start=np.r_[0,np.flatnonzero(np.any(pairs[1:]!=pairs[:-1],axis=1))+1]
    is_feature=np.minimum.reduceat(owner,start)!=np.maximum.reduceat(owner,start)
    features=pairs[start][is_feature]
    # Use labeled matching support, not a global-axis bounding box. This remains
    # valid for a rigidly rotated coupon and identifies artificial cut endpoints.
    touches_matching=np.logical_or.reduceat(on_matching,start)[is_feature]
    physical=features[~touches_matching]
    cut_nodes=set(map(int,tri[refs==1].ravel()))
    lower=p.min(axis=0);upper=p.max(axis=0)
    segments,corners=feature_chains(p,physical,lower,upper,tolerance,cut_nodes=cut_nodes)
    # Pin junctions/turns in the complete preservation graph, including box corners.
    adj={}
    for a,b in features:adj.setdefault(int(a),[]).append(int(b));adj.setdefault(int(b),[]).append(int(a))
    pins=[]
    for node,neighbors in adj.items():
        if len(neighbors)!=2:pins.append(node);continue
        v,w=p[neighbors]-p[node]
        if np.dot(v,w)/np.linalg.norm(v)/np.linalg.norm(w)>-1+1e-10:pins.append(node)
    return features,np.array(sorted(pins)),segments,corners
