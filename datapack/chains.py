"""带版本的二维仿射坐标目录与唯一可达变换链解析。

目录条目形如：
    {id, version, source, target, matrix:[a,b,tx,c,d,ty]}
矩阵含义见 docs/domain.md：目标点 (a*x+b*y+tx, c*x+d*y+ty)。

规则（见 docs/boundary.md）：
- 只允许沿 source -> target 有向边行走，按折线指明的坐标系与地图版本取链；
- 单位（毫米/米）与轴方向完全由矩阵编码，不依据坐标数值猜测单位；
- 枚举两个坐标系之间的全部简单路径，各路径复合结果必须在 MATRIX_EPS 内一致，
  否则拒绝该条折线；无可达路径同样拒绝；
- 源坐标系即目标坐标系时唯一链为单位矩阵。
"""

from .domain import compose
from .tolerances import MATRIX_EPS, SINGULAR_DET_EPS, ChainError

TARGET_CRS = 'map'
_IDENTITY = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0]


def _finite_matrix(matrix):
    if len(matrix) != 6 or any(not isinstance(v, (int, float)) for v in matrix):
        raise ChainError('仿射矩阵需要六个有限数')
    values = [float(v) for v in matrix]
    if any(v != v or v in (float('inf'), float('-inf')) for v in values):
        raise ChainError('仿射矩阵含非有限值')
    a, b, _tx, c, d, _ty = values
    if abs(a * d - b * c) <= SINGULAR_DET_EPS:
        raise ChainError('仿射矩阵近奇异，无法定义坐标变换')
    return values


class TransformCatalog:
    """按地图版本分组的有向变换图。"""

    def __init__(self, entries):
        self._graphs = {}
        for entry in entries:
            for field in ('id', 'version', 'source', 'target', 'matrix'):
                if field not in entry or entry[field] is None:
                    raise ChainError('变换条目缺字段：' + field)
            if entry['source'] == entry['target']:
                raise ChainError('变换条目源与目标坐标系相同：' + str(entry['id']))
            matrix = _finite_matrix(entry['matrix'])
            graph = self._graphs.setdefault(entry['version'], {})
            graph.setdefault(entry['source'], []).append((entry['target'], matrix, entry['id']))

    def versions(self):
        return set(self._graphs)

    def resolve(self, source_crs, version):
        """返回 (复合矩阵, 链上条目id序列)；不可唯一确定时抛 ChainError。"""
        if source_crs == TARGET_CRS:
            return list(_IDENTITY), []
        graph = self._graphs.get(version)
        if graph is None or source_crs not in graph:
            raise ChainError('坐标系 %s 在版本 %s 无可达变换链' % (source_crs, version))
        paths = self._enum_paths(graph, source_crs)
        if not paths:
            raise ChainError('坐标系 %s 在版本 %s 无可达变换链' % (source_crs, version))
        reference_matrix, reference_ids = paths[0]
        for matrix, ids in paths[1:]:
            if any(abs(x - y) > MATRIX_EPS for x, y in zip(matrix, reference_matrix)):
                raise ChainError(
                    '坐标系 %s 在版本 %s 存在不一致的变换链：%s 与 %s'
                    % (source_crs, version, list(reference_ids), list(ids)))
        return reference_matrix, list(reference_ids)

    @staticmethod
    def _enum_paths(graph, start):
        """DFS 枚举 start -> TARGET_CRS 的全部简单路径并复合矩阵。"""
        results = []

        def walk(node, matrix, ids, seen):
            if node == TARGET_CRS:
                results.append((matrix, tuple(ids)))
                return
            for nxt, edge_matrix, edge_id in graph.get(node, ()):
                if nxt in seen:
                    continue
                combined = edge_matrix if matrix is None else compose(edge_matrix, matrix)
                walk(nxt, combined, ids + [edge_id], seen | {nxt})

        walk(start, None, [], {start})
        return results


def apply_matrix(matrix, point):
    """对归一坐标点应用矩阵，仅接受有限二维点。"""
    if len(point) != 2 or any(v != v or v in (float('inf'), float('-inf')) for v in point):
        raise ChainError('坐标点必须为二维有限值')
    a, b, tx, c, d, ty = matrix
    x, y = point
    return [a * x + b * y + tx, c * x + d * y + ty]
