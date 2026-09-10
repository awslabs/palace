# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

# Exact line/circle distances and conservatively bounded tessellation of validated
# conic models for the native size callback. CAD geometry itself is unchanged;
# unsupported nonlinear curves are rejected by the caller.
using LinearAlgebra

struct GradingArc
    center::NTuple{3,Float64}
    u::NTuple{3,Float64}
    v::NTuple{3,Float64}
    normal::NTuple{3,Float64}
    radius::Float64
    angle::Float64
    first::NTuple{3,Float64}
    last::NTuple{3,Float64}
end

function grading_arc(a,quarter,middle,b; tolerance=1e-8)
    # Using quarter and middle points also handles a full circle (a == b).
    a,quarter,middle,b=collect.((a,quarter,middle,b))
    x,y=quarter-a,middle-a
    normal=cross(x,y)
    den=dot(normal,normal)
    den>0 || error("Degenerate circular grading curve")
    center=a+(dot(x,x)*cross(y,normal)+dot(y,y)*cross(normal,x))/(2den)
    radius=norm(a-center)
    u=(a-center)/radius
    normal/=sqrt(den)
    v=cross(normal,u)
    angle=norm(b-a)<=tolerance ? 2pi : mod(atan(dot(b-center,v),dot(b-center,u)),2pi)
    arc=GradingArc(Tuple(center),Tuple(u),Tuple(v),Tuple(normal),radius,angle,Tuple(a),Tuple(b))
    for p in (a,quarter,middle,b)
        grading_arc_distance2(Tuple(p),arc)<=tolerance^2 ||
            error("CAD curve is not a consistent circular arc")
    end
    return arc
end

function grading_arc_distance2(point,arc::GradingArc)
    delta=ntuple(i->point[i]-arc.center[i],3)
    x,y,z=dot(delta,arc.u),dot(delta,arc.v),dot(delta,arc.normal)
    angle=mod(atan(y,x),2pi)
    if angle<=arc.angle
        return (hypot(x,y)-arc.radius)^2+z*z
    end
    return min(sum((point[i]-arc.first[i])^2 for i in 1:3),
               sum((point[i]-arc.last[i])^2 for i in 1:3))
end

struct GradingConic
    segments::Vector{NTuple{6,Float64}}
    distance_error::Float64
end

function grading_conic(parameters,points,fine)
    first,last=extrema(parameters)
    0<last-first<=2pi+1e-8 || error("Unsupported conic parameter interval")
    design=hcat(ones(length(parameters)),cos.(parameters),sin.(parameters))
    cond(design)<1e10 || error("Ill-conditioned conic parameter fit")
    coefficients=design\reduce(vcat,[reshape(collect(p),1,3) for p in points])
    center,u,v=(vec(coefficients[i,:]) for i in 1:3)
    scale=max(1.,norm(u),norm(v))
    all(norm(center+u*cos(t)+v*sin(t)-collect(p))<1e-9scale for (t,p) in zip(parameters,points)) ||
        error("Trimmed curve does not match an analytic conic")
    # For the validated conic model, ||x''|| <= ||[u v]||. Its piecewise-linear
    # interpolation has a Hausdorff error <= ||[u v]|| * delta_parameter^2 / 8.
    # The actual CAD curve is NOT modified by this sizing approximation.
    derivative_bound=opnorm(hcat(u,v))
    target_error=0.01fine
    count=max(2,ceil(Int,(last-first)*sqrt(derivative_bound/(8target_error))))
    count<=10000 || error("Conic grading tessellation budget exceeded")
    parameters2=range(first,last;length=count+1)
    vertices=[Tuple(center+u*cos(t)+v*sin(t)) for t in parameters2]
    segments=[Tuple((vertices[i]...,vertices[i+1]...)) for i in 1:count]
    error_bound=derivative_bound*((last-first)/count)^2/8
    return GradingConic(segments,error_bound),center,u,v
end

function grading_conic_distance2(point,conic::GradingConic)
    distance=sqrt(minimum(grading_segment_distance2(point,s) for s in conic.segments))
    return max(0.,distance-conic.distance_error)^2
end

function grading_segment_distance2(point,s)
    delta=ntuple(i->s[i+3]-s[i],3)
    length2=dot(delta,delta)
    length2>0 || error("Zero-length grading segment")
    t=clamp(sum((point[i]-s[i])*delta[i] for i in 1:3)/length2,0.,1.)
    return sum((point[i]-s[i]-t*delta[i])^2 for i in 1:3)
end
