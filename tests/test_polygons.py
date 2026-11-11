"""多边形简单性校验、含孔洞点归属与边界原语测试。"""

import unittest

from datapack.polygons import (
    classify_in_region,
    point_in_ring,
    point_segment_distance,
    region_edges,
    validate_region,
    validate_ring,
)
from datapack.tolerances import GeometryError

CONCAVE = [[0, 0], [10, 0], [10, 10], [6, 10], [6, 6], [0, 6]]
HOLE = [[2, 2], [4, 2], [4, 4], [2, 4]]


class RingValidationTests(unittest.TestCase):
    def test_valid_concave_ring_accepted(self):
        points, eps = validate_ring(CONCAVE)
        self.assertEqual(len(points), 6)
        self.assertGreater(eps, 0)

    def test_explicit_closure_is_degenerate(self):
        with self.assertRaises(GeometryError):
            validate_ring([[0, 0], [1, 0], [1, 1], [0, 1], [0, 0]])

    def test_zero_area_ring_rejected(self):
        with self.assertRaises(GeometryError):
            validate_ring([[0, 0], [0, 0], [0, 0], [0, 0]])

    def test_bowtie_self_intersection_rejected(self):
        with self.assertRaises(GeometryError):
            validate_ring([[0, 0], [2, 2], [0, 2], [3, 0]])

    def test_vertex_touch_rejected(self):
        pinch = [[0, 0], [4, 0], [4, 4], [0, 4], [0, 2], [2, 2], [2, 0]]
        with self.assertRaises(GeometryError):
            validate_ring(pinch)

    def test_scaled_tolerance_distinguishes_large_coordinates(self):
        big = [[0, 0], [1e6, 0], [1e6, 1e6], [0, 1e6]]
        region = validate_region({'id': 'big', 'map_version': 'v', 'priority': 1,
                                  'outer': big, 'holes': []})
        # 1e-9 相对尺度落在容差带内视为边界；1e-3 在带内之外视为内部
        self.assertEqual(point_in_ring((5e5, 1e-9), region['outer'], region['eps']),
                         'boundary')
        self.assertEqual(point_in_ring((5e5, 1e-3), region['outer'], region['eps']),
                         'inside')


class RegionValidationTests(unittest.TestCase):
    def test_hole_outside_outer_rejected(self):
        with self.assertRaises(GeometryError):
            validate_region({'id': 'r', 'map_version': 'v', 'priority': 1,
                             'outer': CONCAVE,
                             'holes': [[[1, 7], [5, 7], [5, 9], [1, 9]]]})

    def test_hole_touching_outer_rejected(self):
        with self.assertRaises(GeometryError):
            validate_region({'id': 'r', 'map_version': 'v', 'priority': 1,
                             'outer': CONCAVE,
                             'holes': [[[0, 1], [2, 1], [2, 2], [0, 2]]]})

    def test_touching_holes_rejected(self):
        outer = [[0, 0], [10, 0], [10, 10], [0, 10]]
        with self.assertRaises(GeometryError):
            validate_region({'id': 'r', 'map_version': 'v', 'priority': 1,
                             'outer': outer,
                             'holes': [[[1, 1], [3, 1], [3, 3], [1, 3]],
                                       [[3, 1], [5, 1], [5, 3], [3, 3]]]})

    def test_valid_region_with_hole(self):
        region = validate_region({'id': 'r', 'map_version': 'v', 'priority': 1,
                                  'outer': CONCAVE, 'holes': [HOLE]})
        self.assertEqual(len(region_edges(region)), 10)


class ClassificationTests(unittest.TestCase):
    def setUp(self):
        self.region = validate_region({'id': 'A', 'map_version': 'v', 'priority': 1,
                                       'outer': CONCAVE, 'holes': [HOLE]})

    def test_inside_outside_boundary_and_excluded(self):
        self.assertEqual(classify_in_region((1, 1), self.region), 'inside')
        self.assertEqual(classify_in_region((3, 3), self.region), 'excluded')
        self.assertEqual(classify_in_region((50, 50), self.region), 'outside')
        self.assertEqual(classify_in_region((0, 3), self.region), 'boundary')
        self.assertEqual(classify_in_region((3, 2), self.region), 'boundary')

    def test_notch_is_outside(self):
        # 凹口 x<6 且 y>6 在外环之外
        self.assertEqual(classify_in_region((3, 8), self.region), 'outside')

    def test_point_segment_distance(self):
        self.assertAlmostEqual(point_segment_distance((0, 1), (0, 0), (10, 0)), 1.0)
        self.assertAlmostEqual(point_segment_distance((5, 3), (0, 0), (10, 0)), 3.0)
        self.assertAlmostEqual(point_segment_distance((20, 0), (0, 0), (10, 0)), 10.0)


if __name__ == '__main__':
    unittest.main()
