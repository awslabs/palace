# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0
import unittest
import numpy as np
from edge_volume_metric import volume_metric, intersect_metrics, feature_chains, surface_features

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
        _,_,s,c=surface_features(points,tri,refs)
        self.assertEqual(len(s),12);self.assertEqual(len(c),8)
        axis=np.array([1.,2.,3.]);axis/=np.linalg.norm(axis);angle=.43
        cross=np.array([[0,-axis[2],axis[1]],[axis[2],0,-axis[0]],[-axis[1],axis[0],0]])
        q=np.eye(3)*np.cos(angle)+(1-np.cos(angle))*np.outer(axis,axis)+np.sin(angle)*cross
        shift=np.array([2.,-1.,3.])
        _,_,rs,rc=surface_features(points@q.T+shift,tri,refs)
        self.assertEqual(len(rs),12);self.assertEqual(len(rc),8)
        query=np.array([[0.,-.5,-.5],[.4,-.49,-.48]])
        m=volume_metric(query,s,c,.005,.05,.4)
        rm=volume_metric(query@q.T+shift,rs,rc,.005,.05,.4)
        np.testing.assert_allclose(rm,q@m@q.T,rtol=1e-10,atol=1e-8)

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

    def test_bad_controls_fail_closed(self):
        for controls in ((0,.1,1),(.1,.01,1),(.1,2,1)):
            with self.assertRaises(ValueError):volume_metric([[0,0,0]],[[-1,0,0,1,0,0]],[],*controls)
        with self.assertRaises(ValueError):volume_metric([[float('nan'),0,0]],[[-1,0,0,1,0,0]],[],.01,.1,1)

if __name__=='__main__':unittest.main()
