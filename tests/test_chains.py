"""坐标变换链解析测试：唯一可达链、多链一致、冲突与缺路径拒绝。"""

import unittest

from datapack.chains import TARGET_CRS, TransformCatalog, apply_matrix
from datapack.tolerances import ChainError


def entry(id_, version, source, target, matrix):
    return {'id': id_, 'version': version, 'source': source,
            'target': target, 'matrix': matrix}


MM_TO_LOCAL = [0.001, 0, 0, 0, 0.001, 0]
FLIP_X = [-1, 0, 20, 0, 1, 0]
IDENTITY = [1, 0, 0, 0, 1, 0]


class ChainTests(unittest.TestCase):
    def setUp(self):
        self.entries = [
            entry('mm-m', 'v1', 'millimeter', 'local', MM_TO_LOCAL),
            entry('local-map', 'v1', 'local', TARGET_CRS, FLIP_X),
            entry('direct', 'v2', 'millimeter', TARGET_CRS, [0.001, 0, 1, 0, 0.001, 0]),
            entry('mm-m-v2', 'v2', 'millimeter', 'local', MM_TO_LOCAL),
            entry('local-map-v2', 'v2', 'local', TARGET_CRS, IDENTITY),
        ]
        self.catalog = TransformCatalog(self.entries)

    def test_single_chain_composes_in_order(self):
        matrix, ids = self.catalog.resolve('millimeter', 'v1')
        self.assertEqual(ids, ['mm-m', 'local-map'])
        # 毫米 -> 米 -> 反向 x 轴并平移 20
        self.assertEqual(apply_matrix(matrix, [1000, 2000]), [19, 2])

    def test_multiple_consistent_chains_accepted(self):
        # v3：直达链与两段复合链数值一致
        consistent = self.entries + [
            entry('mm-m-v3', 'v3', 'millimeter', 'local', MM_TO_LOCAL),
            entry('local-map-v3', 'v3', 'local', TARGET_CRS, [1, 0, 5, 0, 1, 5]),
            entry('direct-v3', 'v3', 'millimeter', TARGET_CRS, [0.001, 0, 5, 0, 0.001, 5]),
        ]
        catalog = TransformCatalog(consistent)
        matrix, ids = catalog.resolve('millimeter', 'v3')
        self.assertEqual(ids, ['mm-m-v3', 'local-map-v3'])
        self.assertEqual(apply_matrix(matrix, [0, 0]), [5, 5])

    def test_conflicting_chains_rejected(self):
        # v2：直达链平移 1，复合链平移 0，必须拒绝而非任选
        with self.assertRaises(ChainError):
            self.catalog.resolve('millimeter', 'v2')

    def test_missing_path_rejected(self):
        with self.assertRaises(ChainError):
            self.catalog.resolve('unknown', 'v1')
        with self.assertRaises(ChainError):
            self.catalog.resolve('millimeter', 'v9')

    def test_map_crs_uses_identity(self):
        matrix, ids = self.catalog.resolve(TARGET_CRS, 'v1')
        self.assertEqual(matrix, [1.0, 0.0, 0.0, 0.0, 1.0, 0.0])
        self.assertEqual(ids, [])

    def test_no_unit_guessing(self):
        # 数值量级不参与解析：同样的点在不同源坐标系走不同矩阵
        matrix_v1, _ = self.catalog.resolve('millimeter', 'v1')
        self.assertNotAlmostEqual(apply_matrix(matrix_v1, [1, 0])[0],
                                  apply_matrix(IDENTITY, [1, 0])[0])

    def test_singular_matrix_rejected(self):
        entries = self.entries + [entry('bad', 'v4', 'millimeter', TARGET_CRS,
                                        [0, 0, 0, 0, 0.001, 0])]
        with self.assertRaises(ChainError):
            TransformCatalog(entries)

    def test_reversed_axis_transform(self):
        catalog = TransformCatalog([entry('rev', 'v1', 'src', TARGET_CRS,
                                          [1, 0, 0, 0, -1, 0])])
        matrix, _ = catalog.resolve('src', 'v1')
        self.assertEqual(apply_matrix(matrix, [3, 7]), [3, -7])

    def test_duplicate_edges_conflict_rejected_identical_accepted(self):
        conflicting = [
            entry('a', 'v1', 'src', TARGET_CRS, [1, 0, 0, 0, 1, 0]),
            entry('b', 'v1', 'src', TARGET_CRS, [1, 0, 9, 0, 1, 0]),
        ]
        with self.assertRaises(ChainError):
            TransformCatalog(conflicting).resolve('src', 'v1')
        identical = [
            entry('a', 'v1', 'src', TARGET_CRS, [1, 0, 0, 0, 1, 0]),
            entry('b', 'v1', 'src', TARGET_CRS, [1, 0, 0, 0, 1, 0]),
        ]
        matrix, ids = TransformCatalog(identical).resolve('src', 'v1')
        self.assertEqual(apply_matrix(matrix, [2, 2]), [2, 2])
        self.assertEqual(len(ids), 1)


if __name__ == '__main__':
    unittest.main()
