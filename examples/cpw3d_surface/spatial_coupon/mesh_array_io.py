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
