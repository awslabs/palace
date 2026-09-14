#!/usr/bin/env python3
# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""Prepare the complete nodal matching basis for a tractable mesh-policy comparison."""
import argparse
import json
from pathlib import Path
import tomllib

import numpy as np
from box_trace import validate_box_trace
from rebuild_box_coupon_inputs import read_geometry
from generate_spatial_response import (
    build_matching_surface, conductor_at_points, conductor_trace_lifts,
    coupon_bounds, write_basis, write_surface_trace, mesh_boundary_attributes,
)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('geometry',type=Path)
    parser.add_argument('output',type=Path)
    parser.add_argument('--kind',choices=('thin','fabricated'),required=True)
    parser.add_argument('--mesh',action='append',default=[],help='name=mesh.msh')
    parser.add_argument('--trace-only',action='store_true',help='Prepare the complete trace before meshing')
    parser.add_argument('--order',type=int,default=5)
    parser.add_argument('--ring-size',type=int,default=8)
    args=parser.parse_args()
    if not args.mesh and not args.trace_only:parser.error('Supply mesh cases or --trace-only')
    if args.order<1:parser.error('Order must be positive')
    names=[specification.split('=',1)[0] for specification in args.mesh]
    if len(set(names))!=len(names):parser.error('Duplicate mesh case name')
    root=args.geometry.resolve();output=args.output.resolve()
    if output.exists():parser.error('Refuse to reuse a probe directory')
    process=tomllib.loads((root/'process.toml').read_text())
    if process['Units']!='um':parser.error('Expected um geometry')
    if process['SidewallAngle']!=90 or process['TopRounding']!=0 or process['TrenchRounding']!=0:
        parser.error('This probe adapter currently requires sharp vertical geometry')
    edges,facets=read_geometry(root)
    radius=process['Radius'];thickness=process['MetalThickness'];etch=process['Overetch']
    lo,hi=coupon_bounds(edges,radius,thickness,etch)
    levels=[lo[2],hi[2]]
    for edge in edges:
        z=edge['Point'][2];sign=edge['ProcessNormal'][2]
        levels.extend([z,z+sign*thickness,z-sign*etch])
    points,triangles,groups=build_matching_surface(np.asarray([lo,hi]),levels,args.ring_size,
                                                  edges,radius,thickness,process['SidewallAngle'],facets)
    geometry=validate_box_trace(points,triangles)
    labels=conductor_at_points(points,edges,radius,thickness,process['SidewallAngle'],facets)
    active=np.flatnonzero(labels==0)
    conductors=max(e['Conductor'] for e in edges)
    output.mkdir(parents=True)
    paths,zero=write_basis(output,points,triangles,active.tolist())
    sources=[{'Index':i,'Attributes':[1],'DataFile':str(path)} for i,path in enumerate(paths,start=1)]
    for conductor,values in conductor_trace_lifts(points,labels,conductors).items():
        path=output/f'conductor-{conductor}.csv';write_surface_trace(path,points,triangles,values)
        sources.append({'Index':len(sources)+1,'Attributes':[1],'DataFile':str(path),'Conductor':conductor})
    cases=[]
    for specification in args.mesh:
        name,path=specification.split('=',1);mesh=Path(path).resolve()
        if not name.replace('-','').isalnum() or not mesh.exists():parser.error('Invalid mesh case')
        attributes_present=set(mesh_boundary_attributes(mesh))
        families=(5,6) if args.kind=='fabricated' else (4,)
        ground=sorted(a for a in attributes_present if a//1000 in families and a%100==1)
        if not ground:parser.error('Reference conductor has no physical attributes')
        prescribed=[]
        for source in sources:
            record={k:v for k,v in source.items() if k!='Conductor'}
            if 'Conductor' in source:
                record['TerminalAttributes']=sorted(a for a in attributes_present
                    if a//1000 in families and a%100==source['Conductor'])
                if not record['TerminalAttributes']:parser.error('Conductor state has no physical attributes')
            prescribed.append(record)
        post=[]
        if args.kind=='fabricated':
            settings=[('MA',6,.002,10.,.03),('MS',5,.002,11.47,.0003),('SA',3,.002,4.,.002)]
        else:
            settings=[('MA',4,.002,10.,.03),('MS',4,.002,11.47,.0003),('SA',3,.002,4.,.002)]
        for index,(kind,family,t,eps,loss) in enumerate(settings,start=1):
            attributes=sorted(a for a in attributes_present if a//1000==family)
            if not attributes:parser.error('Missing physical interface family')
            post.append({'Index':index,'Attributes':attributes,'Type':kind,'Thickness':t,
                         'Permittivity':eps,'LossTan':loss,'LocalizeEdgeEnergy':True,
                         'EdgeAttributes':sorted(a for a in attributes_present if a//1000==3),
                         'EdgeDistances':[radius],'EdgeFrameNormal':[0.,0.,1.],'EdgeExcludeAttributes':[1]})
        config={'Problem':{'Type':'Electrostatic','Verbose':1,'Output':str(output/name),
                           'OutputFormats':{'Paraview':False,'GridFunction':False}},
                'Model':{'Mesh':str(mesh),'L0':1e-6,'Refinement':{'MaxIts':0}},
                'Domains':{'Materials':[{'Attributes':[1],'Permittivity':11.47},{'Attributes':[2],'Permittivity':1.0}]},
                'Boundaries':{'Ground':{'Attributes':ground},'PrescribedPotential':prescribed,
                              'Postprocessing':{'Dielectric':post}},
                'Solver':{'Order':args.order,'Electrostatic':{'Save':0,'ResponseMatrix':True,'AggregateResponseMatrix':True},
                          'Linear':{'Type':'BoomerAMG','KSPType':'CG','Tol':1e-8,'MaxIts':500,'EstimatorTol':.5,'EstimatorMaxIts':5,'EstimatorMG':True}}}
        config_path=output/(name+'.json');config_path.write_text(json.dumps(config,indent=2)+'\n')
        cases.append({'Name':name,'Mesh':str(mesh),'Config':str(config_path)})
    (output/'manifest.json').write_text(json.dumps({'Scope':'Complete nodal box trace basis and conductor states, not a production library qualification',
        'Sources':len(sources),'TraceGeometry':geometry,'Cases':cases,'Order':args.order,'Tolerance':1e-8},indent=2)+'\n')
    print('Full trace sources:',len(sources))


if __name__=='__main__':main()
