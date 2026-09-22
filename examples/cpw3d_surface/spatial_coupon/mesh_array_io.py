#!/usr/bin/env python3
# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0
"""Read ordinary mesh files or the read-only native-Gmsh array export used by scouts."""
import hashlib
from pathlib import Path
import shlex
import struct
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
    if path.suffix=='.msh':
        mesh=read_msh22_binary(path)
        return meshio.read(path) if mesh is None else mesh
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


# Gmsh MSH 2.2 element types the direct reader accepts (type -> (meshio name, nodes)); the
# linear types are node-ordered alike in Gmsh and meshio (meshio reorders only the
# quadratic types), so a file with another type falls back to meshio.
MSH22_TYPES={1:('line',2),2:('triangle',3),3:('quad',4),4:('tetra',4),5:('hexahedron',8),6:('wedge',6),7:('pyramid',5),15:('vertex',1)}
MSH22_SECTIONS=(b'$MeshFormat',b'$PhysicalNames',b'$Nodes',b'$Elements')


def read_msh22_binary(path):
    """Read a little-endian MSH 2.2 binary file (8-byte reals; the sections $MeshFormat,
    optional $PhysicalNames, $Nodes, $Elements, in that order, nothing else) directly into
    the meshio.Mesh `meshio.read` returns for it: the same points (float64, in node-tag
    order 1..n), the same consecutive-type cell blocks (int32, 0-based), the same
    'gmsh:physical' / 'gmsh:geometrical' cell data (int32 per block; a tag column absent
    when no element carries it), the same field_data ({name: [tag, dimension]}) and empty
    point_data.  Gmsh writes one header per element (type, 1, tag count), so a run of
    elements of one type and tag count is a fixed-stride record array parsed by numpy at
    once instead of meshio's one numpy.fromfile call per element (decision 62 step 3,
    proposal 4).  Returns None for any other layout (the caller falls back to meshio)."""
    data=Path(path).read_bytes()
    if not data.startswith(b'$MeshFormat\n2.2 1 8\n'):return None
    if data[20:24]!=b'\x01\x00\x00\x00' or data[24:40]!=b'\n$EndMeshFormat\n':return None
    cursor=40
    field_data={}
    if data.startswith(b'$PhysicalNames\n',cursor):
        end=data.find(b'$EndPhysicalNames\n',cursor)
        if end<0:return None
        lines=data[cursor+len(b'$PhysicalNames\n'):end].decode().splitlines()
        if not lines or int(lines[0])!=len(lines)-1:return None
        for line in lines[1:]:
            fields=shlex.split(line)
            if len(fields)!=3:return None
            field_data[fields[2]]=np.array([int(fields[1]),int(fields[0])],dtype=int)
        cursor=end+len(b'$EndPhysicalNames\n')
    if not data.startswith(b'$Nodes\n',cursor):return None
    header_end=data.index(b'\n',cursor+len(b'$Nodes\n'))
    node_count=int(data[cursor+len(b'$Nodes\n'):header_end])
    block=header_end+1;nodes_end=block+28*node_count
    if data[nodes_end:nodes_end+len(b'\n$EndNodes\n')]!=b'\n$EndNodes\n':return None
    records=np.frombuffer(data,dtype=np.dtype([('index','<i4'),('x','<f8',(3,))]),count=node_count,offset=block)
    if not np.array_equal(records['index'],np.arange(1,node_count+1,dtype=np.int32)):return None
    points=np.ascontiguousarray(records['x'])
    cursor=nodes_end+len(b'\n$EndNodes\n')
    if not data.startswith(b'$Elements\n',cursor):return None
    header_end=data.index(b'\n',cursor+len(b'$Elements\n'))
    element_count=int(data[cursor+len(b'$Elements\n'):header_end])
    position=header_end+1
    runs=[]  # (meshio type, node array, physical tags or None, geometrical tags or None)
    read=0
    while read<element_count:
        if position+12>len(data):return None
        element_type,count,tag_count=struct.unpack_from('<3i',data,position)
        if element_type not in MSH22_TYPES or count<1 or tag_count<0:return None
        name,node_count_of=MSH22_TYPES[element_type]
        if count!=1:
            # A multi-element header: one record array of this header alone.
            width=1+tag_count+node_count_of
            values=np.frombuffer(data,dtype='<i4',count=count*width,offset=position+12).reshape(count,width)
            runs.append((name,values[:,1+tag_count:],values[:,1:1+tag_count],tag_count));position+=12+4*count*width;read+=count
            continue
        # A run of single-element headers (type, 1, tag_count): fixed-stride records.
        width=3+1+tag_count+node_count_of
        available=min(element_count-read,(len(data)-position)//(4*width))
        values=np.frombuffer(data,dtype='<i4',count=available*width,offset=position).reshape(available,width)
        same=(values[:,0]==element_type)&(values[:,1]==1)&(values[:,2]==tag_count)
        length=int(np.argmin(same)) if not same.all() else available
        if length<1:return None
        values=values[:length]
        runs.append((name,values[:,4+tag_count:],values[:,4:4+tag_count],tag_count));position+=4*width*length;read+=length
    if read!=element_count or data[position:position+len(b'\n$EndElements\n')]!=b'\n$EndElements\n':return None
    if data[position+len(b'\n$EndElements\n'):].strip():return None
    # meshio merges consecutive runs of one type into one block and slices the per-type
    # tag columns per block; a tag column exists only where some element carries it.
    cells=[];physical=[];geometrical=[]
    for name,nodes,tags,tag_count in runs:
        if cells and cells[-1][0]==name:
            cells[-1][1].append(nodes);physical[-1].append(tags[:,0] if tag_count>0 else None)
            geometrical[-1].append(tags[:,1] if tag_count>1 else None)
        else:
            cells.append((name,[nodes]));physical.append([tags[:,0] if tag_count>0 else None])
            geometrical.append([tags[:,1] if tag_count>1 else None])
    blocks=[];cell_data={}
    for (name,parts),phys,geom in zip(cells,physical,geometrical):
        if any(p is None for p in phys)!=all(p is None for p in phys) or any(g is None for g in geom)!=all(g is None for g in geom):
            return None  # mixed tag counts within one block: let meshio decide
        blocks.append((name,(np.concatenate(parts)-1).astype(np.int32)))
        if phys[0] is not None:cell_data.setdefault('gmsh:physical',[]).append(np.concatenate(phys).astype(np.int32))
        if geom[0] is not None:cell_data.setdefault('gmsh:geometrical',[]).append(np.concatenate(geom).astype(np.int32))
    for key in cell_data:
        if len(cell_data[key])!=len(blocks):return None
    return meshio.Mesh(points,blocks,cell_data=cell_data,field_data=field_data)


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
        # Version 4 stores 64-bit counts; this reader reads 32-bit counts and MMG
        # 5.6 never writes version 4, so it is rejected rather than mis-parsed.
        if version==4:raise ValueError('Binary Medit version 4 (64-bit counts) is not supported')
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
