"""含孔洞简单多边形的几何原语与合法性校验。

约定（见 docs/domain.md、docs/boundary.md）：
- 环为隐式闭合点列，不依赖顺逆时针确定意义；外环内部减去孔洞内部为区域确定内部；
- 外环边界本身属于“待核带”，孔洞边界同样属于待核带，孔洞内部属于排除区；
- 退化环（重合点、零面积）与自交/边接触环直接判非法，不产生任何有效归属；
- 同一区域的孔必须互不接触且完整落在外环内部；
- 几何容差按环坐标尺度放缩（tolerances.geometry_eps），并随校验后的区域对象携带，
  跨函数显式传递，不使用隐藏的全局默认值。
"""

from .tolerances import GEOMETRY_EPS, geometry_eps, GeometryError


def _point(point):
    if len(point) != 2 or any(not isinstance(v, (int, float)) for v in point):
        raise GeometryError('点必须为二维数值')
    x, y = float(point[0]), float(point[1])
    if x != x or y != y or x in (float('inf'), float('-inf')) or y in (float('inf'), float('-inf')):
        raise GeometryError('点含非有限值')
    return x, y


def ring_scale(points):
    """环坐标量级：包围盒宽高的最大值，至少为 1，避免小坐标下相对容差趋零。"""
    xs = [p[0] for p in points]
    ys = [p[1] for p in points]
    return max(1.0, max(xs) - min(xs), max(ys) - min(ys))


def ring_edges(ring):
    """闭合环的有向边，跳过零长度边。"""
    edges = []
    for i, start in enumerate(ring):
        end = ring[(i + 1) % len(ring)]
        if start != end:
            edges.append((start, end))
    return edges


def _cross(origin, a, b):
    return (a[0] - origin[0]) * (b[1] - origin[1]) - (a[1] - origin[1]) * (b[0] - origin[0])


def on_segment(point, start, end, eps):
    """点是否落在闭线段上（含端点）。叉积阈值按长度平方尺度 eps*边长 判定。"""
    edge_length = abs(complex(end[0] - start[0], end[1] - start[1]))
    cross_limit = eps * max(1.0, edge_length)
    if abs(_cross(point, start, end)) > cross_limit:
        return False
    return (min(start[0], end[0]) - eps <= point[0] <= max(start[0], end[0]) + eps
            and min(start[1], end[1]) - eps <= point[1] <= max(start[1], end[1]) + eps)


def segments_intersect_proper(start1, end1, start2, end2, eps):
    """两线段是否在非端点处相交（规范相交），用于自交检测。"""
    length1 = abs(complex(end1[0] - start1[0], end1[1] - start1[1]))
    length2 = abs(complex(end2[0] - start2[0], end2[1] - start2[1]))
    limit = eps * max(1.0, length1, length2)
    c1 = _cross(start1, end1, start2)
    c2 = _cross(start1, end1, end2)
    c3 = _cross(start2, end2, start1)
    c4 = _cross(start2, end2, end1)
    return ((c1 > limit and c2 < -limit) or (c1 < -limit and c2 > limit)) and \
           ((c3 > limit and c4 < -limit) or (c3 < -limit and c4 > limit))


def validate_ring(ring):
    """校验单个环：点列合法、无重复点、面积非零、非相邻边不接触。

    返回 (浮点点列, 尺度容差)。
    """
    if not isinstance(ring, list) or len(ring) < 4:
        raise GeometryError('环至少需要四个点（隐式闭合）')
    points = [_point(p) for p in ring]
    eps = geometry_eps(ring_scale(points))
    # 首尾重复属于显式闭合，资料约定为隐式闭合，拒绝以免边语义歧义
    for i in range(len(points)):
        for j in range(i + 1, len(points)):
            if abs(points[i][0] - points[j][0]) <= eps and \
               abs(points[i][1] - points[j][1]) <= eps:
                raise GeometryError('环存在重合点（含首尾显式闭合），属退化环')
    edges = ring_edges(points)
    if len(edges) < 3:
        raise GeometryError('环有效边不足，属退化环')
    area2 = sum(a[0] * b[1] - b[0] * a[1] for a, b in zip(points, points[1:] + points[:1]))
    if abs(area2) <= eps * ring_scale(points):
        raise GeometryError('环面积为零，属退化环')
    count = len(edges)
    for i in range(count):
        s1, e1 = edges[i]
        s2, e2 = edges[(i + 1) % count]
        # 相邻边必须恰好共享 e1==s2；零长边已去除，故等长反向即折返重叠
        if e1 != s2 or s1 == e2:
            raise GeometryError('相邻边折返或接触，环非简单')
    for i in range(count):
        s1, e1 = edges[i]
        for j in range(i + 1, count):
            s2, e2 = edges[j]
            adjacent = (j == i + 1) or (i == 0 and j == count - 1)
            if adjacent:
                continue
            if segments_intersect_proper(s1, e1, s2, e2, eps):
                raise GeometryError('非相邻边规范相交，环自交')
            # 非相邻边的任何端点接触/重叠均使环非简单
            if on_segment(s1, s2, e2, eps) or on_segment(e1, s2, e2, eps) \
                    or on_segment(s2, s1, e1, eps) or on_segment(e2, s1, e1, eps):
                raise GeometryError('非相邻边接触，环自交')
    return points, eps


def point_in_ring(point, ring, eps):
    """点相对单环的归属：'boundary' / 'inside' / 'outside'，射线法，与方向无关。"""
    for start, end in ring_edges(ring):
        if on_segment(point, start, end, eps):
            return 'boundary'
    inside = False
    px, py = point
    for start, end in ring_edges(ring):
        x1, y1 = start
        x2, y2 = end
        if (y1 > py) != (y2 > py):
            cross_x = x1 + (py - y1) * (x2 - x1) / (y2 - y1)
            if cross_x > px + eps:
                inside = not inside
    return 'inside' if inside else 'outside'


def classify_in_region(point, region):
    """点在单个含孔洞区域中的确定归属。

    返回 'boundary'（外环或孔环上）、'excluded'（孔洞内部）、
    'inside'（外环内且不在孔内）、'outside'。
    """
    eps = region['eps']
    outer_state = point_in_ring(point, region['outer'], eps)
    if outer_state == 'boundary':
        return 'boundary'
    if outer_state == 'outside':
        return 'outside'
    for hole in region['holes']:
        hole_state = point_in_ring(point, hole, eps)
        if hole_state == 'boundary':
            return 'boundary'
        if hole_state == 'inside':
            return 'excluded'
    return 'inside'


def point_segment_distance(point, start, end, eps=GEOMETRY_EPS):
    """点到闭线段的欧氏距离。"""
    vx, vy = end[0] - start[0], end[1] - start[1]
    wx, wy = point[0] - start[0], point[1] - start[1]
    length2 = vx * vx + vy * vy
    if length2 <= eps * eps:
        return ((point[0] - start[0]) ** 2 + (point[1] - start[1]) ** 2) ** 0.5
    ratio = max(0.0, min(1.0, (wx * vx + wy * vy) / length2))
    closest_x = start[0] + ratio * vx
    closest_y = start[1] + ratio * vy
    return ((point[0] - closest_x) ** 2 + (point[1] - closest_y) ** 2) ** 0.5


def validate_region(entry):
    """校验区域资料并归一化点列，非法时抛 GeometryError。"""
    for field in ('id', 'map_version', 'outer', 'priority'):
        if field not in entry or entry[field] is None:
            raise GeometryError('区域条目缺字段：' + field)
    if not isinstance(entry['priority'], int) or isinstance(entry['priority'], bool):
        raise GeometryError('区域优先级必须为整数：' + str(entry['id']))
    outer, outer_eps = validate_ring(entry['outer'])
    holes = []
    hole_eps = []
    for hole in entry.get('holes', []):
        hole_points, eps = validate_ring(hole)
        holes.append(hole_points)
        hole_eps.append(eps)
    eps = max([outer_eps] + hole_eps)
    for hole in holes:
        # 孔环须整体位于外环内部：每个顶点在内、且两环边互不接触
        for vertex in hole:
            if point_in_ring(vertex, outer, eps) != 'inside':
                raise GeometryError('孔洞越出外环或触及外环：' + str(entry['id']))
        for hs, he in ring_edges(hole):
            for os_, oe in ring_edges(outer):
                if segments_intersect_proper(hs, he, os_, oe, eps) or \
                        on_segment(hs, os_, oe, eps) or on_segment(he, os_, oe, eps) \
                        or on_segment(os_, hs, he, eps) or on_segment(oe, hs, he, eps):
                    raise GeometryError('孔洞与外环接触：' + str(entry['id']))
    # 孔与孔之间不得相交或接触
    for i in range(len(holes)):
        for j in range(i + 1, len(holes)):
            for s1, e1 in ring_edges(holes[i]):
                for s2, e2 in ring_edges(holes[j]):
                    if segments_intersect_proper(s1, e1, s2, e2, eps) or \
                            on_segment(s1, s2, e2, eps) or on_segment(e1, s2, e2, eps) \
                            or on_segment(s2, s1, e1, eps) or on_segment(e2, s1, e1, eps):
                        raise GeometryError('同区域孔洞相互接触：' + str(entry['id']))
    return {
        'id': entry['id'],
        'map_version': entry['map_version'],
        'outer': outer,
        'holes': holes,
        'priority': entry['priority'],
        'eps': eps,
    }


def region_edges(region):
    """区域的全部边界边（外环与孔环），供待核带距离计算。"""
    edges = ring_edges(region['outer'])
    for hole in region['holes']:
        edges.extend(ring_edges(hole))
    return edges
