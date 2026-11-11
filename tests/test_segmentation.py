"""空间切分与守恒验证器测试。

覆盖：凹区与孔洞多次进出、擦边待核带、误差半径膨胀、区域重叠优先级、
共边不重复计长、退化点、反向轴与毫米单位归一、缺变换/冲突链拒绝、
非法区域拒绝以及长度/参数/时间守恒。
"""

import math
import unittest

from datapack.segmentation import (
    STATE_EXCLUDED,
    STATE_INSIDE,
    STATE_OUTSIDE,
    STATE_UNCERTAIN,
    process_bundle,
)

CONCAVE_A = [[0, 0], [10, 0], [10, 10], [6, 10], [6, 6], [0, 6]]
HOLE = [[2, 2], [4, 2], [4, 4], [2, 4]]
SQUARE_B = [[10, 0], [20, 0], [20, 10], [10, 10]]


def region(id_, version, outer, priority, holes=None):
    return {'id': id_, 'map_version': version, 'outer': outer,
            'holes': holes or [], 'priority': priority}


def path(id_, crs, version, points):
    return {'id': id_, 'crs': crs, 'version': version, 'points': points}


def pt(x, y, time, error=0.0):
    return {'x': x, 'y': y, 'time': time, 'error': error}


class SegmentationFixture:
    """共享构造与断言的混入基类，本身不被 unittest 收集。"""

    def setUp(self):
        self.regions = [
            region('A', 'v1', CONCAVE_A, 1, [HOLE]),
            region('B', 'v1', SQUARE_B, 2),
            region('C', 'v1', [[7, 7], [13, 7], [13, 9], [7, 9]], 3),
        ]
        self.transforms = [
            {'id': 'mm-m', 'version': 'v1', 'source': 'millimeter', 'target': 'local',
             'matrix': [0.001, 0, 0, 0, 0.001, 0]},
            {'id': 'local-map', 'version': 'v1', 'source': 'local', 'target': 'map',
             'matrix': [-1, 0, 20, 0, 1, 0]},
        ]

    def run(self, paths, regions=None, transforms=None):
        report = process_bundle(transforms if transforms is not None else self.transforms,
                                regions if regions is not None else self.regions, paths)
        return report

    def accepted(self, report, path_id):
        result = next(item for item in report['paths'] if item['path_id'] == path_id)
        self.assertTrue(result['accepted'], result.get('reason'))
        return result

    def assert_conservation(self, result, expected_length):
        verification = result['verification']
        self.assertTrue(verification['all_ok'], verification)
        length_check = verification['checks'][0]
        self.assertTrue(length_check['ok'])
        self.assertAlmostEqual(length_check['detail']['actual'], expected_length)
        self.assertAlmostEqual(sum(piece['length'] for piece in result['pieces']),
                               expected_length)


class CuttingTests(SegmentationFixture):
    def test_concave_hole_and_multi_region_sequence(self):
        # 直线横穿外部、A 实体、孔洞排除、A 实体、B 实体、外部
        report = self.run([path('p1', 'map', 'v1',
                                [pt(-2, 3, 0), pt(22, 3, 24)])])
        result = self.accepted(report, 'p1')
        states = [(piece['state'], piece['region_id']) for piece in result['pieces']]
        self.assertEqual(states, [
            (STATE_OUTSIDE, None),
            (STATE_INSIDE, 'A'),
            (STATE_EXCLUDED, 'A'),
            (STATE_INSIDE, 'A'),
            (STATE_INSIDE, 'B'),
            (STATE_OUTSIDE, None),
        ])
        self.assert_conservation(result, 24.0)
        # 孔洞段必须显式标 excluded 而非并入 inside
        self.assertIn((STATE_EXCLUDED, 'A'), states)

    def test_repeated_entry_and_exit(self):
        report = self.run([path('re', 'map', 'v1', [
            pt(-2, 5, 0), pt(5, 5, 7), pt(15, 5, 15), pt(5, 5, 25), pt(-2, 5, 32),
        ])])
        result = self.accepted(report, 're')
        inside = [p for p in result['pieces'] if p['state'] == STATE_INSIDE]
        outside = [p for p in result['pieces'] if p['state'] == STATE_OUTSIDE]
        # 一条折线两次进入区域，产生至少两个 inside 片
        self.assertGreaterEqual(len(inside), 2)
        self.assertGreaterEqual(len(outside), 2)
        self.assert_conservation(result, 34.0)

    def test_millimeter_and_reversed_axis_normalization(self):
        # 毫米源 + x 轴反向：(0,1000mm)->(20,1m)，(20000,1000)->(0,1)
        report = self.run([path('mm', 'millimeter', 'v1',
                                [pt(0, 1000, 0, 0.1), pt(20000, 1000, 20, 0.1)])])
        result = self.accepted(report, 'mm')
        self.assertEqual(result['transform_chain'], ['mm-m', 'local-map'])
        first = result['normalized_vertices'][0]
        self.assertAlmostEqual(first['x'], 20.0)
        self.assertAlmostEqual(first['y'], 1.0)
        self.assert_conservation(result, 20.0)

    def test_error_radius_marks_uncertain_band(self):
        # 误差半径 0.1：距任一区域边界 <=0.1 的部分单列为 uncertain
        report = self.run([path('band', 'map', 'v1',
                                [pt(-1, 5, 0, 0.1), pt(1, 5, 2, 0.1)])])
        result = self.accepted(report, 'band')
        uncertain = [p for p in result['pieces'] if p['state'] == STATE_UNCERTAIN]
        self.assertTrue(uncertain)
        # 待核带总长约 0.2（边界 x=0 两侧各 0.1），不强行归入 A
        total_uncertain = sum(p['length'] for p in uncertain)
        self.assertAlmostEqual(total_uncertain, 0.2, places=6)
        for piece in uncertain:
            self.assertIsNone(piece['region_id'])
        self.assert_conservation(result, 2.0)

    def test_segment_error_uses_larger_endpoint_radius(self):
        report = self.run([path('r', 'map', 'v1',
                                [pt(-1, 5, 0, 0.3), pt(1, 5, 2, 0.1)])])
        result = self.accepted(report, 'r')
        self.assertTrue(all(piece['error_radius'] == 0.3 for piece in result['pieces']))

    def test_collinear_grazing_is_uncertain_even_with_zero_error(self):
        # 沿外环边界 y=0 行走：边界属待核带，零误差也不能算 inside/outside
        report = self.run([path('graze', 'map', 'v1',
                                [pt(2, 0, 0, 0.0), pt(8, 0, 6, 0.0)])])
        result = self.accepted(report, 'graze')
        self.assertEqual({p['state'] for p in result['pieces']}, {STATE_UNCERTAIN})
        self.assert_conservation(result, 6.0)

    def test_diagonal_through_shared_corner_not_falsely_collinear(self):
        # 边起点恰在折线延长线上但方向不平行时，不得误判共线
        regions = [region('E', 'v3', [[4, 4], [16, 4], [16, 16], [4, 16]], 1)]
        report = self.run([path('diag', 'map', 'v3',
                                [pt(5, 5, 0), pt(15, 15, 20)])], regions=regions,
                          transforms=[])
        result = self.accepted(report, 'diag')
        self.assertEqual([p['state'] for p in result['pieces']], [STATE_INSIDE])
        self.assertEqual(result['pieces'][0]['region_id'], 'E')

    def test_overlapping_regions_follow_priority(self):
        # C（priority 3）与 A、B 重叠；重叠确定部分必须归优先级更高的 A/B
        report = self.run([path('ov', 'map', 'v1',
                                [pt(8, 8, 0), pt(12, 8, 8), pt(12, 12, 16)])])
        result = self.accepted(report, 'ov')
        owners = {(round(p['start'][0], 3), p['state']): p['region_id']
                  for p in result['pieces'] if p['state'] == STATE_INSIDE}
        # x=8..10 在 A∩C -> A；x=10..12 在 B∩C -> B；任何片不得归属 C
        self.assertNotIn('C', {p['region_id'] for p in result['pieces']})
        segment_zero = [p for p in result['pieces'] if p['segment_index'] == 0]
        self.assertEqual([p['region_id'] for p in segment_zero], ['A', 'B'])

    def test_shared_edge_counted_once(self):
        report = self.run([path('se', 'map', 'v1',
                                [pt(5, 5, 0), pt(15, 5, 10)])])
        result = self.accepted(report, 'se')
        # 共享边 x=10 只是切点：总长 10，A/B 各 5，无零长重复片
        lengths = {(p['region_id']): 0.0 for p in result['pieces']}
        for piece in result['pieces']:
            lengths[piece['region_id']] = lengths.get(piece['region_id'], 0.0) + piece['length']
        self.assertAlmostEqual(lengths['A'], 5.0)
        self.assertAlmostEqual(lengths['B'], 5.0)
        self.assertAlmostEqual(sum(p['length'] for p in result['pieces']), 10.0)
        no_double = next(c for c in result['verification']['checks']
                         if c['name'] == 'no_double_counting')
        self.assertTrue(no_double['ok'])

    def test_single_point_degenerate_polyline(self):
        report = self.run([
            path('sin', 'map', 'v1', [pt(5, 5, 3)]),
            path('sout', 'map', 'v1', [pt(50, 50, 3)]),
            path('sbd', 'map', 'v1', [pt(0, 5, 3)]),
            path('shole', 'map', 'v1', [pt(3, 2, 3)]),
        ])
        by_id = {item['path_id']: item for item in report['paths']}
        self.assertEqual(by_id['sin']['pieces'][0]['state'], STATE_INSIDE)
        self.assertEqual(by_id['sin']['pieces'][0]['length'], 0.0)
        self.assertEqual(by_id['sout']['pieces'][0]['state'], STATE_OUTSIDE)
        # 单点落在外环或孔洞边界上：待核带，不强行归属
        self.assertEqual(by_id['sbd']['pieces'][0]['state'], STATE_UNCERTAIN)
        self.assertEqual(by_id['shole']['pieces'][0]['state'], STATE_UNCERTAIN)
        self.assert_conservation(by_id['sin'], 0.0)

    def test_time_and_parameter_interpolation_preserved(self):
        report = self.run([path('tm', 'map', 'v1',
                                [pt(-2, 3, 0), pt(22, 3, 24)])])
        result = self.accepted(report, 'tm')
        # 片时间单调、首尾覆盖原时间跨度，切点参数落在原线段 [0,1]
        times = [piece['time0'] for piece in result['pieces']]
        self.assertEqual(times, sorted(times))
        self.assertAlmostEqual(result['pieces'][0]['time0'], 0.0)
        self.assertAlmostEqual(result['pieces'][-1]['time1'], 24.0)
        for piece in result['pieces']:
            self.assertGreaterEqual(piece['t0'], 0.0)
            self.assertLessEqual(piece['t1'], 1.0)
            self.assertLessEqual(piece['t0'], piece['t1'])
        time_check = next(c for c in result['verification']['checks']
                          if c['name'] == 'time_conservation')
        self.assertTrue(time_check['ok'])


class RejectionTests(SegmentationFixture):
    def test_missing_transform_chain_rejects_path_only(self):
        report = self.run([
            path('ok', 'map', 'v1', [pt(1, 1, 0), pt(2, 2, 1)]),
            path('bad', 'missing', 'v1', [pt(0, 0, 0)]),
        ])
        good = next(p for p in report['paths'] if p['path_id'] == 'ok')
        bad = next(p for p in report['paths'] if p['path_id'] == 'bad')
        self.assertTrue(good['accepted'])
        self.assertFalse(bad['accepted'])
        self.assertEqual(bad['pieces'], [])
        self.assertIn('无可达变换链', bad['reason'])

    def test_conflicting_chains_reject_path(self):
        transforms = self.transforms + [
            {'id': 'direct', 'version': 'v2', 'source': 'millimeter', 'target': 'map',
             'matrix': [0.001, 0, 1, 0, 0.001, 0]},
            {'id': 'mm2', 'version': 'v2', 'source': 'millimeter', 'target': 'local',
             'matrix': [0.001, 0, 0, 0, 0.001, 0]},
            {'id': 'lm2', 'version': 'v2', 'source': 'local', 'target': 'map',
             'matrix': [1, 0, 0, 0, 1, 0]},
        ]
        report = self.run([path('conflict', 'millimeter', 'v2',
                                [pt(0, 0, 0), pt(1, 1, 1)])], transforms=transforms)
        result = report['paths'][0]
        self.assertFalse(result['accepted'])
        self.assertIn('不一致的变换链', result['reason'])

    def test_invalid_region_rejects_its_map_version(self):
        bad_regions = self.regions + [
            region('bad', 'v1', [[0, 0], [2, 2], [0, 2], [2, 0]], 9),
        ]
        report = self.run([path('p', 'map', 'v1', [pt(1, 1, 0), pt(2, 2, 1)])],
                          regions=bad_regions)
        self.assertTrue(report['region_failures'])
        self.assertFalse(report['paths'][0]['accepted'])

    def test_invalid_region_does_not_block_other_versions(self):
        bad_regions = self.regions + [
            region('bad', 'v9', [[0, 0], [2, 2], [0, 2], [2, 0]], 9),
            region('safe', 'v9', [[0, 0], [1, 0], [1, 1], [0, 1]], 1),
        ]
        report = self.run([path('p', 'map', 'v1', [pt(1, 1, 0), pt(2, 2, 1)])],
                          regions=bad_regions)
        self.assertTrue(report['paths'][0]['accepted'])

    def test_negative_error_rejected(self):
        report = self.run([path('neg', 'map', 'v1', [pt(0, 0, 0, -0.1)])])
        self.assertFalse(report['paths'][0]['accepted'])


class GeometryConservationTests(SegmentationFixture):
    def test_random_polylines_always_conserve_length(self):
        # 合成随机折线：无论切段多碎，片长之和必须等于归一折线长度
        import random
        random.seed(20261111)
        paths = []
        for index in range(12):
            points = [pt(random.uniform(-5, 25), random.uniform(-5, 15),
                         time, random.choice([0.0, 0.05, 0.2]))
                      for time in range(random.randint(2, 6))]
            paths.append(path('rand%d' % index, 'map', 'v1', points))
        report = self.run(paths)
        for result in report['paths']:
            self.assertTrue(result['accepted'], result.get('reason'))
            self.assertTrue(result['verification']['all_ok'])
            expected = sum(
                math.hypot(b['x'] - a['x'], b['y'] - a['y'])
                for a, b in zip(result['normalized_vertices'],
                                result['normalized_vertices'][1:]))
            self.assertAlmostEqual(sum(p['length'] for p in result['pieces']),
                                   expected, places=6)

    def test_piece_segment_references_cover_every_source_segment(self):
        report = self.run([path('ref', 'map', 'v1', [
            pt(-2, 3, 0), pt(22, 3, 24), pt(22, 13, 30),
        ])])
        result = self.accepted(report, 'ref')
        referenced = {piece['segment_index'] for piece in result['pieces']}
        self.assertEqual(referenced, {0, 1})
        for piece in result['pieces']:
            self.assertEqual(piece['path_id'], 'ref')


if __name__ == '__main__':
    unittest.main()
