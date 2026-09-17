#!/usr/bin/env python3
# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0
"""Read ordinary mesh files or the read-only native-Gmsh array export used by scouts."""
import hashlib
from pathlib import Path
import tomllib
import meshio
import numpy as np


def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda:f.read(8*1024*1024),b''):h.update(block)
    return h.hexdigest()


def read_mesh(path):
    path=Path(path)
    if path.suffix=='.meshb':return read_medit_binary(path)
    if path.suffix!='.toml':return meshio.read(path)
    data=tomllib.loads(path.read_text());root=path.parent
    if sha(data['SourceMesh'])!=data['SourceSHA256']:raise ValueError('Native export source mesh changed')
    def binary(name,dtype):
        file=(root/name).resolve()
        if not file.is_relative_to(root.resolve()):raise ValueError('Array path escapes export')
        return np.fromfile(file,dtype=dtype)
    tags=binary(data['NodeTags'],'<u8');xyz=binary(data['Coordinates'],'<f8')
    if len(tags)!=data['NodeCount'] or len(xyz)!=3*len(tags) or not len(tags):raise ValueError('Node count mismatch')
    if len(np.unique(tags))!=len(tags) or tags.min()<1 or tags.max()>10*len(tags):raise ValueError('Invalid/sparse node numbering')
    if not np.all(np.isfinite(xyz)):raise ValueError('Nonfinite coordinates')
    index=np.full(int(tags.max())+1,-1,dtype=np.int64);index[tags]=np.arange(len(tags))
    blocks=[];attributes=[]
    for cell in data['Cells']:
        if cell['Type'] not in (2,4):raise ValueError('Unsupported native cell type')
        width=3 if cell['Type']==2 else 4
        conn=binary(cell['Connectivity'],'<u8')
        if len(conn)!=width*cell['Count'] or np.any(conn>=len(index)):raise ValueError('Connectivity mismatch')
        conn=index[conn].reshape(-1,width)
        if np.any(conn<0):raise ValueError('Unknown vertex')
        blocks.append(('triangle' if width==3 else 'tetra',conn))
        attributes.append(np.full(len(conn),cell['Attribute'],dtype=np.int32))
    return meshio.Mesh(xyz.reshape(-1,3),blocks,cell_data={'gmsh:physical':attributes})


# libMeshb keyword codes of the native binary Medit (.meshb) fields the tools consume.
MEDIT_VERTICES=4;MEDIT_TRIANGLES=6;MEDIT_TETRAHEDRA=8;MEDIT_REQUIRED_TETRAHEDRA=12;MEDIT_END=54


def read_medit_binary(path):
    """Read a native binary Medit mesh (MMG output) into a meshio Mesh.

    Every keyword carries the file position of the next keyword, so fields the
    tools do not consume (corners, ridges, required vertices/edges/triangles,
    normals) are skipped by position and the reader is complete for any keyword
    set.  MMG's RequiredTetrahedra keyword (kept-verbatim seed cells) is returned
    as the 0/1 cell_data 'medit:required' of the tetrahedron block.
    """
    with Path(path).open('rb') as f:
        magic=np.fromfile(f,count=1,dtype='<i4');version=np.fromfile(f,count=1,dtype='<i4')
        if magic.size!=1 or version.size!=1 or magic.item()!=1 or version.item() not in (1,2,3,4):
            raise ValueError('Not a binary Medit mesh')
        version=version.item()
        position_type='<i8' if version>=3 else '<i4'
        real_type='<f4' if version==1 else '<f8'
        if np.fromfile(f,count=1,dtype='<i4').item()!=3:raise ValueError('Binary Medit mesh lacks its dimension')
        np.fromfile(f,count=1,dtype=position_type)
        dimension=np.fromfile(f,count=1,dtype='<i4').item()
        if dimension!=3:raise ValueError('Only 3D binary Medit meshes are supported')
        points=required=None;blocks=[]
        while True:
            code=np.fromfile(f,count=1,dtype='<i4')
            if code.size==0:raise ValueError('Binary Medit mesh ended before its End keyword')
            code=code.item()
            if code==MEDIT_END:break
            following=np.fromfile(f,count=1,dtype=position_type).item()
            if code==MEDIT_VERTICES:
                count=np.fromfile(f,count=1,dtype='<i4').item()
                points=np.fromfile(f,count=count,dtype=np.dtype([('x',real_type,(3,)),('ref','<i4')]))
            elif code in (MEDIT_TRIANGLES,MEDIT_TETRAHEDRA):
                count=np.fromfile(f,count=1,dtype='<i4').item()
                width=3 if code==MEDIT_TRIANGLES else 4
                block=np.fromfile(f,count=count*(width+1),dtype='<i4').reshape(count,width+1)
                blocks.append(('triangle' if code==MEDIT_TRIANGLES else 'tetra',block))
            elif code==MEDIT_REQUIRED_TETRAHEDRA:
                count=np.fromfile(f,count=1,dtype='<i4').item()
                required=np.fromfile(f,count=count,dtype='<i4')
            if following<=0:raise ValueError('Binary Medit keyword lacks its next position')
            f.seek(following)
    kinds=[kind for kind,_ in blocks]
    if points is None or kinds.count('triangle')!=1 or kinds.count('tetra')!=1:
        raise ValueError('Binary Medit mesh lacks vertices, triangles or tetrahedra')
    xyz=np.ascontiguousarray(points['x'],dtype=float)
    if not np.all(np.isfinite(xyz)):raise ValueError('Nonfinite coordinates')
    cells=[];refs=[];flags=[]
    for kind,block in blocks:
        if np.any(block[:,:-1]<1) or np.any(block[:,:-1]>len(xyz)):raise ValueError('Binary Medit connectivity is out of range')
        kept=np.zeros(len(block),dtype=np.int32)
        if kind=='tetra' and required is not None:
            if np.any(required<1) or np.any(required>len(block)):raise ValueError('Required tetrahedron index is out of range')
            kept[required-1]=1
        cells.append((kind,(block[:,:-1]-1).astype(np.int64)));refs.append(block[:,-1].astype(np.int32));flags.append(kept)
    return meshio.Mesh(xyz,cells,point_data={'medit:ref':points['ref'].astype(np.int32)},
                       cell_data={'medit:ref':refs,'medit:required':flags})
