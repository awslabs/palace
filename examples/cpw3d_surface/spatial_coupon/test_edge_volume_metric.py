# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0
import unittest
import numpy as np
from edge_volume_metric import (COPLANAR_TOLERANCE, cluster_coplanar_triangles, feature_chains,
                                intersect_metrics, match_equivalent_planes, plane_deviation,
                                surface_features, volume_metric)
from prepare_edge_metric_scout import budget_aware_far_policy, protected_corner_ball_triangles

class EdgeVolumeMetricTest(unittest.TestCase):
    def test_two_transverse_directions_and_corner_recovery(self):
        segments=[[-1,0,0,1,0,0]];corners=[[-1,0,0],[1,0,0]]
        p=np.array([[0,0,0],[-1,0,0],[-.99,0,0],[0,2,0.]])
        m=volume_metric(p,segments,corners,.005,.05,.4)
        np.testing.assert_allclose(m[0],np.diag([1/.05**2,1/.005**2,1/.005**2]))
        np.testing.assert_allclose(m[1],np.eye(3)/.005**2)
        np.testing.assert_allclose(m[2],np.diag([1/.0075**2,1/.005**2,1/.005**2]))
        np.testing.assert_allclose(m[3],np.eye(3)/.4**2)

    def test_intersection_is_spd_dominates_inputs_and_idempotent(self):
        a=np.array([[[4.,1.,0.],[1.,2.,0.],[0.,0.,3.]]]);b=np.array([np.diag([1.,8.,2.])])
        m=intersect_metrics(a,b)
        self.assertGreater(np.linalg.eigvalsh(m).min(),0)
        self.assertGreater(np.linalg.eigvalsh(m-a).min(),-1e-12)
        self.assertGreater(np.linalg.eigvalsh(m-b).min(),-1e-12)
        np.testing.assert_allclose(intersect_metrics(m,m),m,atol=1e-12)

    def test_intersection_rejects_invalid_either_operand(self):
        good=np.eye(3)[None,:,:]
        invalid=(np.diag([1.,1.,-1.])[None,:,:],
                 np.array([[[1.,1.,0.],[0.,1.,0.],[0.,0.,1.]]]),
                 np.array([[[1.,0.,0.],[0.,np.nan,0.],[0.,0.,1.]]]))
        for bad in invalid:
            with self.assertRaises(ValueError):intersect_metrics(bad,good)
            with self.assertRaises(ValueError):intersect_metrics(good,bad)
        with self.assertRaises(ValueError):intersect_metrics(good,np.eye(2)[None,:,:])

    def test_rotation_translation_scale_and_duplicate_sources(self):
        p=np.array([[0.,0.,0.],[.1,.02,.03]])
        segments=np.array([[-1.,0,0,1,0,0],[0,-1,0,0,1,0]])
        corners=np.array([[0.,0.,0.],[-1,0,0],[1,0,0],[0,-1,0],[0,1,0]])
        m=volume_metric(p,segments,corners,.005,.05,.4)
        angle=.37;q=np.array([[np.cos(angle),-np.sin(angle),0],[np.sin(angle),np.cos(angle),0],[0,0,1.]])
        moved=volume_metric(p@q.T+[2,3,4],(segments.reshape(-1,3)@q.T+[2,3,4]).reshape(-1,6),corners@q.T+[2,3,4],.005,.05,.4)
        np.testing.assert_allclose(moved,q@m@q.T,rtol=1e-10,atol=1e-8)
        scaled=volume_metric(3*p,3*segments,3*corners,.015,.15,1.2)
        np.testing.assert_allclose(scaled,m/9,rtol=1e-10,atol=1e-9)
        duplicate=volume_metric(p,np.repeat(segments[:1],2,axis=0),corners,.005,.05,.4)
        single=volume_metric(p,segments[:1],corners,.005,.05,.4)
        np.testing.assert_allclose(duplicate,single,rtol=1e-10,atol=1e-9)

    def test_semantic_isotropy_ball_does_not_blanket_subdivisions_or_cuts(self):
        segments=[[-1.,0,0,1,0,0]]
        geometric_corners=[[-1.,0,0],[1.,0,0]]
        semantic_corners=[[-1.,0,0]]
        points=np.array([[-1.,0,0],[-.96,0,0],[0.,0,0],[.96,0,0]])
        metric=volume_metric(points,segments,geometric_corners,.01,.1,.4,
                             isotropic_corners=semantic_corners,isotropy_radius=.1)
        np.testing.assert_allclose(metric[:2],np.broadcast_to(np.eye(3)/.01**2,(2,3,3)))
        np.testing.assert_allclose(metric[2],np.diag([1/.1**2,1/.01**2,1/.01**2]))
        np.testing.assert_allclose(metric[3],np.diag([1/.02**2,1/.01**2,1/.01**2]))

    def test_subdivision_is_not_a_corner_and_cut_endpoints_not_pinned(self):
        p=np.array([[-1.,0,0],[0,0,0],[1,0,0]])
        segments,corners=feature_chains(p,[[0,1],[1,2]],np.array([-1,-2,-2]),np.array([1,2,2]))
        np.testing.assert_allclose(segments,[[-1,0,0,1,0,0]])
        self.assertEqual(len(corners),0)
        m=volume_metric(p,segments,corners,.005,.05,.4)
        np.testing.assert_allclose(m,np.broadcast_to(np.diag([400.,40000.,40000.]),(3,3,3)))
        t=np.vstack((p,[0,1,0]))
        _,corners=feature_chains(t,[[0,1],[1,2],[1,3]],np.full(3,-2),np.full(3,2))
        self.assertTrue(any(np.all(c==0) for c in corners))

    def test_feature_classification_and_tensor_follow_full_3d_rotation(self):
        cube=np.array([[x,y,z] for z in (0.,1.) for y in (0.,1.) for x in (0.,1.)])
        faces=np.array([[0,2,3],[0,3,1],[4,5,7],[4,7,6],[0,1,5],[0,5,4],
                        [2,6,7],[2,7,3],[0,4,6],[0,6,2],[1,3,7],[1,7,5]])
        points=np.vstack((4*cube-2,cube-.5));tri=np.vstack((faces,faces+8))
        refs=np.r_[np.ones(12,dtype=int),np.full(12,6001,dtype=int)]
        _,_,_,s,c=surface_features(points,tri,refs,{1})
        self.assertEqual(len(s),12);self.assertEqual(len(c),8)
        axis=np.array([1.,2.,3.]);axis/=np.linalg.norm(axis);angle=.43
        cross=np.array([[0,-axis[2],axis[1]],[axis[2],0,-axis[0]],[-axis[1],axis[0],0]])
        q=np.eye(3)*np.cos(angle)+(1-np.cos(angle))*np.outer(axis,axis)+np.sin(angle)*cross
        shift=np.array([2.,-1.,3.])
        _,_,_,rs,rc=surface_features(points@q.T+shift,tri,refs,{1})
        self.assertEqual(len(rs),12);self.assertEqual(len(rc),8)
        query=np.array([[0.,-.5,-.5],[.4,-.49,-.48]])
        m=volume_metric(query,s,c,.005,.05,.4)
        rm=volume_metric(query@q.T+shift,rs,rc,.005,.05,.4)
        np.testing.assert_allclose(rm,q@m@q.T,rtol=1e-10,atol=1e-8)

    @staticmethod
    def _boxed_cube(top_center_lift=0.):
        """Unit cube (6001) inside a matching box (1); top face fanned around a center."""
        cube=np.array([[x,y,z] for z in (0.,1.) for y in (0.,1.) for x in (0.,1.)])
        faces=np.array([[0,2,3],[0,3,1],[4,5,7],[4,7,6],[0,1,5],[0,5,4],
                        [2,6,7],[2,7,3],[0,4,6],[0,6,2],[1,3,7],[1,7,5]])
        inner=cube-.5;center=np.array([[0.,0.,.5+top_center_lift]])
        top=np.array([[4,5,8],[5,7,8],[7,6,8],[6,4,8]])  # 8+8: the center point
        points=np.vstack((4*cube-2,inner,center))
        tri=np.vstack((faces,faces[np.r_[0,1,4,5,6,7,8,9,10,11]]+8,top+8))
        refs=np.r_[np.ones(12,dtype=int),np.full(14,6001,dtype=int)]
        return points,tri,refs

    def test_coplanar_label_seam_is_a_reference_boundary_not_a_feature(self):
        points,tri,refs=self._boxed_cube()
        features,pins,kinds,segments,corners=surface_features(points,tri,refs,{1})
        self.assertEqual(len(features),24);self.assertEqual(len(pins),16)
        self.assertEqual(len(segments),12);self.assertEqual(len(corners),8)
        self.assertTrue(np.all(kinds=='geometric-corner'))
        split=refs.copy();split[-2:]=6101  # two of the four coplanar top fans
        seam_features,seam_pins,seam_kinds,seam_segments,seam_corners=surface_features(
            points,tri,split,{1})
        # The seam is a straight line through the center: no ridge, no metric
        # segment, no corner, and the collinear center vertex is not a pin.
        np.testing.assert_array_equal(seam_features,features)
        np.testing.assert_array_equal(seam_segments,segments)
        np.testing.assert_array_equal(seam_corners,corners)
        np.testing.assert_array_equal(seam_pins,pins)
        np.testing.assert_array_equal(seam_kinds,kinds)
        seam_edges={(12,16),(13,16),(14,16),(15,16)}
        self.assertFalse(seam_edges&set(map(tuple,seam_features)));self.assertNotIn(16,seam_pins)

    def test_coplanar_seam_junction_is_a_reference_turn_pin_only(self):
        points,tri,refs=self._boxed_cube()
        base=surface_features(points,tri,refs,{1})
        split=refs.copy();split[-1]=6101  # one fan: the seam bends at the center
        features,pins,kinds,segments,corners=surface_features(points,tri,split,{1})
        np.testing.assert_array_equal(features,base[0])
        np.testing.assert_array_equal(segments,base[3]);np.testing.assert_array_equal(corners,base[4])
        self.assertEqual(len(pins),17);self.assertIn(16,pins)
        self.assertEqual(dict(zip(pins.tolist(),kinds.tolist()))[16],'reference-turn')
        self.assertEqual(int(np.sum(kinds=='reference-turn')),1)
        self.assertEqual(int(np.sum(kinds=='geometric-corner')),16)

    def test_near_coplanar_within_tolerance_is_not_a_feature_but_a_dihedral_is(self):
        base=surface_features(*self._boxed_cube(),{1})
        within=surface_features(*self._boxed_cube(top_center_lift=1e-8),{1})
        for one,two in zip(base,within):np.testing.assert_array_equal(one,two)
        features,pins,kinds,segments,corners=surface_features(
            *self._boxed_cube(top_center_lift=1e-3),{1})
        self.assertEqual(len(features),28);self.assertIn(16,pins);self.assertEqual(len(pins),17)
        self.assertTrue(np.all(kinds=='geometric-corner'))
        self.assertEqual(len(segments),16)
        self.assertTrue(any(np.allclose(c,[0.,0.,.5+1e-3]) for c in corners))
        self.assertEqual(len(corners),9)

    def test_seam_and_pin_classification_follow_full_3d_rotation(self):
        points,tri,refs=self._boxed_cube()
        split=refs.copy();split[-1]=6101
        axis=np.array([-2.,1.,.5]);axis/=np.linalg.norm(axis);angle=1.1
        cross=np.array([[0,-axis[2],axis[1]],[axis[2],0,-axis[0]],[-axis[1],axis[0],0]])
        q=np.eye(3)*np.cos(angle)+(1-np.cos(angle))*np.outer(axis,axis)+np.sin(angle)*cross
        shift=np.array([-3.,2.,7.])
        for labels in (refs,split):
            plain=surface_features(points,tri,labels,{1})
            moved=surface_features(points@q.T+shift,tri,labels,{1})
            np.testing.assert_array_equal(moved[0],plain[0]);np.testing.assert_array_equal(moved[1],plain[1])
            np.testing.assert_array_equal(moved[2],plain[2])
            self.assertEqual(len(moved[3]),len(plain[3]))
            np.testing.assert_allclose(moved[4],plain[4]@q.T+shift,rtol=0,atol=1e-12)

    def test_protected_radial_band_keeps_old_metric(self):
        s=[[-1.,0,0,1,0,0]];corners=[[-1.,0,0],[1,0,0]]
        p=np.array([[0.,0,0],[0,.01,0],[0,.03,0],[0,.1,0],[0,2,0]])
        old=volume_metric(p,s,corners,.0005,.05,.5)
        new=volume_metric(p,s,corners,.0005,.05,1.,protected_distance=.03,far_growth=2.)
        np.testing.assert_allclose(new[:3],old[:3],rtol=1e-12,atol=1e-8)
        self.assertTrue(np.all(np.linalg.eigvalsh(old[3]-new[3])>=-1e-8))
        np.testing.assert_allclose(new[4],np.eye(3))
        with self.assertRaises(ValueError):volume_metric(p,s,corners,.0005,.05,1.,protected_distance=-1)
        with self.assertRaises(ValueError):volume_metric(p,s,corners,.0005,.05,1.,far_growth=0)

    def test_generic_budget_policy_coarsens_only_far_field_for_small_and_loaded_cases(self):
        small=budget_aware_far_policy(100,1000,.4,1.)
        loaded=budget_aware_far_policy(500,1000,.4,1.)
        self.assertEqual(small['Pressure'],1.)
        self.assertGreater(loaded['Pressure'],1.)
        self.assertEqual(small['UnaffectedTargets'],loaded['UnaffectedTargets'])
        points=np.array([[0.,.01,0.],[0.,.1,0.],[0.,2.,0.]])
        segments=[[-1.,0.,0.,1.,0.,0.]]
        corners=[[-1.,0.,0.],[1.,0.,0.]]
        base=volume_metric(points,segments,corners,.005,.05,small['EffectiveFarSize'],
                           protected_distance=.03,
                           far_growth=small['EffectiveFarGrowth'])
        coarse=volume_metric(points,segments,corners,.005,.05,loaded['EffectiveFarSize'],
                             protected_distance=.03,
                             far_growth=loaded['EffectiveFarGrowth'])
        np.testing.assert_allclose(coarse[0],base[0],rtol=0.,atol=1e-10)
        self.assertTrue(np.all(np.linalg.eigvalsh(base[1:]-coarse[1:])>=-1e-8))

    def test_protected_corner_balls_freeze_straddling_triangles_by_any_vertex(self):
        corners=[[0.,0.,0.],[5.,0.,0.]]
        points=np.array([[0.,0,0],[.05,0,0],[0,.05,0],      # inside the first ball
                         [.09,.09,0],[.3,0,0],[.3,.3,0],    # nearest vertex at 0.127: outside
                         [.11,0,0],[.5,0,0],[.5,.5,0],      # entirely outside
                         [4.95,0,0],[5.2,0,0],[5.,.3,0]])   # straddles the second ball
        triangles=np.array([[0,1,2],[3,4,5],[6,7,8],[9,10,11]])
        frozen,counts=protected_corner_ball_triangles(points,triangles,corners,.1)
        # Triangle 1 crosses the ball's bounding box but has no vertex inside;
        # triangle 3 straddles the second ball through its vertex at x = 4.95.
        np.testing.assert_array_equal(frozen,[True,False,False,True])
        np.testing.assert_array_equal(counts,[1,1])
        straddling=np.array([[0,3,5],[6,7,8]])
        frozen,counts=protected_corner_ball_triangles(points,straddling,corners,.1)
        np.testing.assert_array_equal(frozen,[True,False])
        np.testing.assert_array_equal(counts,[1,0])
        angle=.37;q=np.array([[np.cos(angle),-np.sin(angle),0],[np.sin(angle),np.cos(angle),0],[0,0,1.]])
        moved,moved_counts=protected_corner_ball_triangles(points@q.T+[2,3,4],triangles,
                                                           np.asarray(corners)@q.T+[2,3,4],.1)
        np.testing.assert_array_equal(moved,[True,False,False,True])
        np.testing.assert_array_equal(moved_counts,[1,1])
        for radius in (0.,-1.,float('nan')):
            with self.assertRaises(ValueError):protected_corner_ball_triangles(points,triangles,corners,radius)

    @staticmethod
    def _noisy_matching_plane():
        """Seven triangles of one x = 9.8333... box face whose normals carry
        roundoff-level (1e-8) noise, like the ten-edge seed's matching plane."""
        rng=np.random.default_rng(7);x=9.8333333333
        centers=np.array([[0.,0.],[2.,1.],[-3.,0.5],[1.,-2.],[3.,3.],[-2.,-3.],[0.5,2.5]])
        xyz=[]
        for c in centers:
            corners=np.array([[x,c[0]-.1,c[1]-.1],[x,c[0]+.1,c[1]-.1],[x,c[0],c[1]+.1]])
            corners[:,0]+=rng.uniform(-3.5e-10,3.5e-10,3)  # x spread <= 7e-10 as measured
            xyz.append(corners)
        xyz=np.asarray(xyz);n=np.cross(xyz[:,1]-xyz[:,0],xyz[:,2]-xyz[:,0])
        n/=np.linalg.norm(n,axis=1)[:,None];n*=np.sign(n[:,:1])
        return xyz,n

    def test_noisy_matching_plane_keys_collapse_to_one_equivalence_class(self):
        xyz,n=self._noisy_matching_plane()
        # The retired 8-digit rounding fragments the plane into several keys ...
        rounded=np.round(np.column_stack((n,np.einsum('ij,ij->i',n,xyz[:,0]))),8)
        self.assertGreater(len(np.unique(rounded,axis=0)),1)
        # ... whereas the tolerance rule sees one plane and one label class.
        representatives,patch=cluster_coplanar_triangles(np.ones(len(n),int),n,xyz)
        self.assertEqual(len(representatives),1)
        np.testing.assert_array_equal(patch,0)
        self.assertEqual(representatives[0,0],1)
        # Labels are never merged: the same planes under two labels are two classes.
        representatives,patch=cluster_coplanar_triangles([1,1,1,2,2,2,2],n,xyz)
        self.assertEqual(len(representatives),2)
        np.testing.assert_array_equal(patch,[0,0,0,1,1,1,1])

    def test_distinct_planes_stay_distinct_and_the_rule_is_dimensionless(self):
        xyz,n=self._noisy_matching_plane()
        # A 1e-3 rad dihedral (and a 1e-4 rad one) is a different plane; the
        # roundoff class (1e-8) is not.
        for angle in (1e-3,1e-4):
            rotation=np.array([[np.cos(angle),-np.sin(angle),0],[np.sin(angle),np.cos(angle),0],[0,0,1]])
            tilted=xyz.copy();tilted[3:]=(xyz[3:]-xyz[3,0])@rotation.T+xyz[3,0]
            m=np.cross(tilted[:,1]-tilted[:,0],tilted[:,2]-tilted[:,0]);m/=np.linalg.norm(m,axis=1)[:,None]
            representatives,patch=cluster_coplanar_triangles(np.ones(7,int),m,tilted)
            self.assertEqual(len(representatives),2,angle)
            np.testing.assert_array_equal(patch,[0,0,0,1,1,1,1])
        # A parallel plane offset by 1e-3 of the local size is distinct; the
        # same offset scaled with the geometry stays distinct (dimensionless).
        for scale in (1.,1e-3,1e3):
            shifted=xyz*scale;shifted[3:,:,0]+=1e-3*0.2*scale
            representatives,patch=cluster_coplanar_triangles(np.ones(7,int),n,shifted)
            self.assertEqual(len(representatives),2,scale)
            representatives,_=cluster_coplanar_triangles(np.ones(7,int),n,xyz*scale)
            self.assertEqual(len(representatives),1,scale)
        deviation=plane_deviation([1.,0,0],[0.,0,0],.1,[[1.,0,0]],[[[0.,5,5],[1e-7,5,5],[0,5,6]]],.1)
        # Offset 1e-7 seen from 7 units away: an angle of ~1.4e-8, not 1e-7/0.1.
        np.testing.assert_allclose(deviation,[1e-7/np.linalg.norm([1e-7,5,5])],rtol=1e-6)
        with self.assertRaises(ValueError):plane_deviation([1.,0,0],[0.,0,0],0.,[[1.,0,0]],[[0.,0,0]],.1)
        with self.assertRaises(ValueError):cluster_coplanar_triangles([1],[[1.,0,0]],[[[0,0,0],[0,0,0],[0,1,0]]])

    def test_plane_equivalence_is_rotation_covariant_and_matches_bijectively(self):
        xyz,n=self._noisy_matching_plane()
        angle=1e-3;rotation=np.array([[np.cos(angle),-np.sin(angle),0],[np.sin(angle),np.cos(angle),0],[0,0,1]])
        tilted=xyz.copy();tilted[3:]=(xyz[3:]-xyz[3,0])@rotation.T+xyz[3,0]
        m=np.cross(tilted[:,1]-tilted[:,0],tilted[:,2]-tilted[:,0]);m/=np.linalg.norm(m,axis=1)[:,None]
        labels=np.array([1,1,1,1,2,2,2])
        reference,patch=cluster_coplanar_triangles(labels,m,tilted)
        theta,phi=.63,.41
        rz=np.array([[np.cos(theta),-np.sin(theta),0],[np.sin(theta),np.cos(theta),0],[0,0,1]])
        rx=np.array([[1,0,0],[0,np.cos(phi),-np.sin(phi)],[0,np.sin(phi),np.cos(phi)]])
        R=rx@rz;t=np.array([1.2,-0.7,0.9])
        moved=tilted@R.T+t;mn=m@R.T
        rotated,rotated_patch=cluster_coplanar_triangles(labels,mn,moved)
        np.testing.assert_array_equal(rotated_patch,patch)
        self.assertEqual(len(rotated),len(reference))
        np.testing.assert_allclose(rotated[:,1:4],reference[:,1:4]@R.T,atol=1e-14)
        np.testing.assert_allclose(rotated[:,4:7],reference[:,4:7]@R.T+t,atol=1e-14)
        np.testing.assert_allclose(rotated[:,7],reference[:,7],rtol=1e-13)
        # Matching pairs the same planes of two triangulations (here: the same
        # planes seen through other representative triangles, shuffled) ...
        shuffled=np.array([6,4,5,3,0,2,1])
        candidate,_=cluster_coplanar_triangles(labels[shuffled],m[shuffled],tilted[shuffled])
        mapping=match_equivalent_planes(reference,candidate)
        self.assertIsNotNone(mapping)
        for i,j in enumerate(mapping):self.assertEqual(reference[i,0],candidate[j,0])
        self.assertEqual(sorted(mapping),list(range(len(reference))))
        # ... and fails closed on a missing, an extra, a relabeled or a moved plane.
        self.assertIsNone(match_equivalent_planes(reference,candidate[:-1]))
        self.assertIsNone(match_equivalent_planes(reference,np.vstack((candidate,candidate[-1:]))))
        relabeled=candidate.copy();relabeled[:,0]=np.where(relabeled[:,0]==1,2,1)
        self.assertIsNone(match_equivalent_planes(reference,relabeled))
        displaced=candidate.copy();displaced[0,4]+=1e-3
        self.assertIsNone(match_equivalent_planes(reference,displaced))
        self.assertEqual(match_equivalent_planes(np.zeros((0,8)),np.zeros((0,8))),[])
        self.assertEqual(COPLANAR_TOLERANCE,1e-6)

    def test_bad_controls_fail_closed(self):
        for controls in ((0,.1,1),(.1,.01,1),(.1,2,1)):
            with self.assertRaises(ValueError):volume_metric([[0,0,0]],[[-1,0,0,1,0,0]],[],*controls)
        with self.assertRaises(ValueError):volume_metric([[float('nan'),0,0]],[[-1,0,0,1,0,0]],[],.01,.1,1)

if __name__=='__main__':unittest.main()
