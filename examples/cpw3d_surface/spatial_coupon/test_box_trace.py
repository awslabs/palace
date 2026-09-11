# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0
import unittest
import numpy as np
from box_trace import complete_box_trace, validate_box_trace
from generate_spatial_response import build_matching_surface


class BoxTraceTest(unittest.TestCase):
    def test_missing_corner_is_added_without_dropping_old_nodes(self):
        xy=[(-8.,-8.),(2.,-8.),(6.,-7.),(6.,8.),(-8.,8.)]
        points=np.asarray([(x,y,z) for z in (-2.05,-.05,0.,.1,2.1) for x,y in xy])
        new,triangles,groups,report=complete_box_trace(points)
        self.assertEqual(report['AddedVertices'],5)
        self.assertEqual(groups,[6]*5)
        self.assertTrue(np.array_equal(new[report['RetainedVertexMap']],points))
        self.assertEqual(report['OffBoxTriangles'],0)
        self.assertTrue(report['ClosedOrientedSurface'])
        self.assertLess(report['MaximumRelativeFaceAreaDifference'],1e-12)

    def test_non_square_generator_always_includes_all_corners(self):
        bounds=np.asarray([[-8.,-8.,-2.05],[6.,8.,2.1]])
        points,triangles,_=build_matching_surface(bounds,[-2.05,0.,2.1],8,[],2.,.1,90.,[])
        result=validate_box_trace(points,triangles,bounds)
        self.assertTrue(result['ClosedOrientedSurface'])
        for z in (-2.05,0.,2.1):
            self.assertIn((6.,-8.,z),set(map(tuple,points)))

    def test_corner_cut_is_rejected(self):
        points=np.asarray([[0.,0.,0.],[1.,.1,0.],[0.,0.,1.],[1.,1.,1.]])
        with self.assertRaisesRegex(ValueError,'matching-box face'):
            validate_box_trace(points,np.asarray([[0,1,2]]))


if __name__=='__main__':unittest.main()
