"""离线空间归一与边界切分。

处理流水线：
1. 用 chains.TransformCatalog 按折线的坐标系与地图版本解析唯一可达变换链，
   缺链或多链不一致的整条折线被拒绝（不产出任何片段）；
2. 折线全部顶点归一到地图坐标（米），误差半径随点保留；
3. 在所有区域边界（外环与孔环）处精确切段，逐段收集穿越参数与距离临界参数，
   段内误差半径取两端点较大者，距任一边界不大于该半径的部分单列为 uncertain，
   其余按点归属记 inside / outside / excluded；
4. 区域重叠的确定部分按 priority 归属（数值小者优先），共边只计一次长度；
5. 守恒验证器核对长度、参数与时间。

不使用包围盒或首点代替裁切；所有容差集中在 tolerances.py。
"""

from .chains import TransformCatalog, apply_matrix
from .domain import segment_length
from .polygons import classify_in_region, region_edges, validate_region
from .tolerances import (
    BISECT_LIMIT,
    GEOMETRY_EPS,
    PARAM_EPS,
    ChainError,
    GeometryError,
    length_tolerance,
)

STATE_INSIDE = 'inside'
STATE_OUTSIDE = 'outside'
STATE_EXCLUDED = 'excluded'
STATE_UNCERTAIN = 'uncertain'


def _dedup_params(params):
    """参数去重并限定在 [0,1]，PARAM_EPS 内视为同一切点。"""
    result = []
    for value in params:
        if value < -PARAM_EPS or value > 1.0 + PARAM_EPS:
            continue
        value = min(1.0, max(0.0, value))
        if all(abs(value - kept) > PARAM_EPS for kept in result):
            result.append(value)
    return sorted(result)


def _merge_intervals(intervals):
    """合并相交或相接的参数区间。"""
    if not intervals:
        return []
    ordered = sorted(intervals)
    merged = [list(ordered[0])]
    for start, end in ordered[1:]:
        if start <= merged[-1][1] + PARAM_EPS:
            merged[-1][1] = max(merged[-1][1], end)
        else:
            merged.append([start, end])
    return [(max(0.0, s), min(1.0, e)) for s, e in merged]


def segment_intersection_parameter(start, end, edge_start, edge_end, eps=GEOMETRY_EPS):
    """线段 P(t)=start+t*(end-start) 与边界线段的交点参数。

    平行（含共线）返回 None；共线覆盖由调用方按整体待核处理。
    """
    vx, vy = end[0] - start[0], end[1] - start[1]
    ex, ey = edge_end[0] - edge_start[0], edge_end[1] - edge_start[1]
    denom = vx * ey - vy * ex
    # denom=v×E 为面积量纲，平行阈值按两边长度乘积与几何容差放缩
    parallel_limit = eps * (vx * vx + vy * vy) ** 0.5 * (ex * ex + ey * ey) ** 0.5
    if abs(denom) <= parallel_limit:
        return None
    dx = edge_start[0] - start[0]
    dy = edge_start[1] - start[1]
    t = (dx * ey - dy * ex) / denom
    u = (dx * vy - dy * vx) / denom
    if -PARAM_EPS <= t <= 1.0 + PARAM_EPS and -PARAM_EPS <= u <= 1.0 + PARAM_EPS:
        return min(1.0, max(0.0, t))
    return None


def _collinear_overlap(start, end, edge_start, edge_end, eps=GEOMETRY_EPS):
    """共线时返回折线上的覆盖参数区间 (t0,t1)，不共线或无覆盖返回 None。"""
    vx, vy = end[0] - start[0], end[1] - start[1]
    length2 = vx * vx + vy * vy
    if length2 <= eps * eps:
        return None
    ex, ey = edge_end[0] - edge_start[0], edge_end[1] - edge_start[1]
    if ex * ex + ey * ey <= eps * eps:
        return None
    # 共线必须方向平行（v×E≈0）且边起点落在折线所在直线上（w×v≈0）
    parallel_cross = abs(vx * ey - vy * ex)
    parallel_limit = eps * (length2 ** 0.5) * ((ex * ex + ey * ey) ** 0.5)
    if parallel_cross > parallel_limit:
        return None
    cross = (edge_start[0] - start[0]) * vy - (edge_start[1] - start[1]) * vx
    cross_limit = eps * (length2 ** 0.5)
    if abs(cross) > cross_limit:
        return None
    params = [0.0, 1.0]
    for point in (edge_start, edge_end):
        t = ((point[0] - start[0]) * vx + (point[1] - start[1]) * vy) / length2
        if -PARAM_EPS <= t <= 1.0 + PARAM_EPS:
            params.append(min(1.0, max(0.0, t)))
    t0, t1 = min(params), max(params)
    return (t0, t1) if t1 - t0 > PARAM_EPS else None


def distance2_quadratic(start, end, edge_start, edge_end):
    """折线点 P(t)=start+t*v 到边界边（无限直线）最小距离平方的二次型 (A,B,C)。

    令 w=start-edge_start、E=edge_end-edge_start，边方向单位向量 n=E/|E|，
    min_u |P(t)-(edge_start+uE)|^2 = (w+tv)^T (I-nn^T) (w+tv)。
    点到“线段”的距离在边端点投影参数之间与该式一致，之外由端点距离接管，
    端点投影已由 _distance_break_params 纳入临界网格。
    """
    vx, vy = end[0] - start[0], end[1] - start[1]
    ex, ey = edge_end[0] - edge_start[0], edge_end[1] - edge_start[1]
    wx, wy = start[0] - edge_start[0], start[1] - edge_start[1]
    e2 = ex * ex + ey * ey
    if e2 <= GEOMETRY_EPS * GEOMETRY_EPS:
        # 退化边退化为点距离
        a = vx * vx + vy * vy
        b = 2.0 * (wx * vx + wy * vy)
        c = wx * wx + wy * wy
        return a, b, c
    nxv, nyv = ex / (e2 ** 0.5), ey / (e2 ** 0.5)

    def project(px, py):
        dot = px * nxv + py * nyv
        return px - dot * nxv, py - dot * nyv

    mwx, mwy = project(wx, wy)
    mvx, mvy = project(vx, vy)
    a = mvx * mvx + mvy * mvy
    b = 2.0 * (mwx * mvx + mwy * mvy)
    c = mwx * mwx + mwy * mwy
    return a, b, c


def _distance2(coefficients, t):
    a, b, c = coefficients
    return a * t * t + b * t + c


def _distance_break_params(start, end, edge_start, edge_end):
    """收集使距离平方分段单调的临界参数：边端点投影与二次顶点。"""
    vx, vy = end[0] - start[0], end[1] - start[1]
    length2 = vx * vx + vy * vy
    params = []
    if length2 > GEOMETRY_EPS * GEOMETRY_EPS:
        for point in (edge_start, edge_end):
            t = ((point[0] - start[0]) * vx + (point[1] - start[1]) * vy) / length2
            if 0.0 < t < 1.0:
                params.append(t)
    a, b, _c = distance2_quadratic(start, end, edge_start, edge_end)
    if a > GEOMETRY_EPS * GEOMETRY_EPS:
        vertex = -b / (2.0 * a)
        if 0.0 < vertex < 1.0:
            params.append(vertex)
    return params


def _uncertain_intervals(start, end, edges, radius):
    """计算该线段上距任一边界 <= radius 的参数区间并集。

    距离平方在临界参数之间单调，用二分法求 |距离|=radius 的根，
    不依赖包围盒，也不以端点状态代替整段。与边界共线的重叠部分属于边界本身
    （待核带），即使误差半径为 0 也单列为 uncertain。
    """
    collinear = []
    edge_quadratics = []
    for edge_start, edge_end in edges:
        overlap = _collinear_overlap(start, end, edge_start, edge_end)
        if overlap is not None:
            collinear.append(overlap)
            continue
        coefficients = distance2_quadratic(start, end, edge_start, edge_end)
        # 每条边用自己的临界点建网格：网格内该边距离平方单调，根至多一个
        grid = _dedup_params([0.0, 1.0]
                             + _distance_break_params(start, end, edge_start, edge_end))
        edge_quadratics.append((coefficients, grid))
    r2 = radius * radius
    hits = []

    def locate(coefficients, lo, hi, low_inside):
        # 在 (lo,hi) 内求距离平方等于 r2 的根；区间内单调，至多一个根
        for _ in range(BISECT_LIMIT):
            if hi - lo <= PARAM_EPS:
                break
            mid = (lo + hi) / 2.0
            mid_inside = _distance2(coefficients, mid) <= r2
            if mid_inside == low_inside:
                lo = mid
            else:
                hi = mid
        return (lo + hi) / 2.0

    if radius > 0.0:
        for coefficients, grid in edge_quadratics:
            for lo, hi in zip(grid, grid[1:]):
                low_inside = _distance2(coefficients, lo) <= r2 + GEOMETRY_EPS
                high_inside = _distance2(coefficients, hi) <= r2 + GEOMETRY_EPS
                if low_inside and high_inside:
                    hits.append((lo, hi))
                elif low_inside != high_inside:
                    root = locate(coefficients, lo, hi, low_inside)
                    hits.append((lo, root) if low_inside else (root, hi))
    # 共线重叠恒为待核带，不受误差半径是否为 0 影响
    hits.extend(collinear)
    return _merge_intervals(hits)


def _point_at(start, end, t):
    return [start[0] + (end[0] - start[0]) * t, start[1] + (end[1] - start[1]) * t]


def _classify(point, regions):
    """确定点归属：重叠内部按 priority（小者优先，同级按 id 稳定排序）。

    边界点（boundary）不参与确定归属，由调用方的待核带逻辑处理；
    没有任何区域给出 inside/excluded 时即为 outside。
    """
    best = None
    for region in regions:
        state = classify_in_region(point, region)
        if state in (STATE_INSIDE, STATE_EXCLUDED):
            rank = (region['priority'], region['id'])
            if best is None or rank < best[0]:
                best = (rank, region['id'], state)
    if best is None:
        return None, STATE_OUTSIDE
    return best[1], best[2]


def _point_state(point, regions):
    """单点的最终归属：任一边界上即 uncertain，否则走确定归属。"""
    for region in regions:
        if classify_in_region(point, region) == 'boundary':
            return None, STATE_UNCERTAIN
    return _classify(point, regions)


def cut_segment(start, end, edges, radius, regions):
    """在一条原线段上完成边界切段与状态归属，返回片字典列表。"""
    length = segment_length(start, end)
    if length <= GEOMETRY_EPS:
        region_id, state = _point_state(start, regions)
        return [{'t0': 0.0, 't1': 1.0, 'region_id': region_id, 'state': state}]

    crossings = [0.0, 1.0]
    for edge_start, edge_end in edges:
        if _collinear_overlap(start, end, edge_start, edge_end) is not None:
            continue
        parameter = segment_intersection_parameter(start, end, edge_start, edge_end)
        if parameter is not None:
            crossings.append(parameter)
    uncertain = _uncertain_intervals(start, end, edges, radius)

    def in_uncertain(t):
        return any(lo - PARAM_EPS <= t <= hi + PARAM_EPS for lo, hi in uncertain)

    # 边界切点与待核带端点共同构成切段网格，每片状态单一，uncertain 单列
    cuts = _dedup_params([0.0, 1.0] + [v for interval in uncertain for v in interval]
                         + crossings[1:-1])
    pieces = []
    for lo, hi in zip(cuts, cuts[1:]):
        if hi - lo <= PARAM_EPS:
            continue
        mid = (lo + hi) / 2.0
        if in_uncertain(mid):
            pieces.append({'t0': lo, 't1': hi, 'region_id': None, 'state': STATE_UNCERTAIN})
        else:
            region_id, state = _classify(_point_at(start, end, mid), regions)
            pieces.append({'t0': lo, 't1': hi, 'region_id': region_id, 'state': state})
    return pieces


def validate_polyline(entry):
    """校验折线资料结构，返回归一化所需的顶点原始信息。"""
    for field in ('id', 'crs', 'version', 'points'):
        if field not in entry or entry[field] is None:
            raise ChainError('折线条目缺字段：' + field)
    points = entry['points']
    if not isinstance(points, list) or not points:
        raise ChainError('折线至少需要一个点：' + str(entry['id']))
    normalized = []
    previous_time = None
    for index, point in enumerate(points):
        for field in ('x', 'y', 'time', 'error'):
            if field not in point or not isinstance(point[field], (int, float)):
                raise ChainError('折线点缺数值字段 %s：%s' % (field, entry['id']))
        x, y, time_value, error = (float(point[k]) for k in ('x', 'y', 'time', 'error'))
        if any(v != v or v in (float('inf'), float('-inf')) for v in (x, y, time_value, error)):
            raise ChainError('折线点含非有限值：' + entry['id'])
        if error < 0.0:
            raise ChainError('定位误差半径不能为负：' + entry['id'])
        if previous_time is not None and time_value < previous_time - GEOMETRY_EPS:
            raise ChainError('折线时间必须单调不减：' + entry['id'])
        previous_time = time_value
        normalized.append({'x': x, 'y': y, 'time': time_value, 'error': error})
    return normalized


def normalize_polyline(points, matrix):
    """应用变换链，返回地图坐标顶点；误差半径按点保留（资料约定 error 已是米）。"""
    vertices = []
    for point in points:
        mapped = apply_matrix(matrix, [point['x'], point['y']])
        vertices.append({
            'x': mapped[0], 'y': mapped[1],
            'time': point['time'], 'error': point['error'],
        })
    return vertices


def _segment_error(vertices, index):
    """段内定位误差半径保守取两端点较大者。"""
    return max(vertices[index]['error'], vertices[index + 1]['error'])


def _enrich(piece, vertices, segment_index, path_id):
    """补齐片的几何长度、端点时间与原段引用。"""
    i = segment_index
    p0, p1 = vertices[i], vertices[i + 1]
    t0, t1 = piece['t0'], piece['t1']
    start = _point_at([p0['x'], p0['y']], [p1['x'], p1['y']], t0)
    end = _point_at([p0['x'], p0['y']], [p1['x'], p1['y']], t1)
    return {
        'path_id': path_id,
        'segment_index': segment_index,
        't0': t0,
        't1': t1,
        'region_id': piece['region_id'],
        'state': piece['state'],
        'start': start,
        'end': end,
        'time0': p0['time'] + (p1['time'] - p0['time']) * t0,
        'time1': p0['time'] + (p1['time'] - p0['time']) * t1,
        'length': segment_length(start, end),
        'error_radius': _segment_error(vertices, segment_index),
    }


def verify_conservation(path_result, vertices):
    """守恒验证：片长之和=归一折线长度；参数与时间闭合；共边不重复计长。"""
    pieces = path_result['pieces']
    expected_length = sum(
        segment_length([vertices[i]['x'], vertices[i]['y']],
                       [vertices[i + 1]['x'], vertices[i + 1]['y']])
        for i in range(len(vertices) - 1))
    actual_length = sum(piece['length'] for piece in pieces)
    checks = []

    def record(name, ok, detail):
        checks.append({'name': name, 'ok': bool(ok), 'detail': detail})

    record('length_conservation',
           abs(actual_length - expected_length) <= length_tolerance(expected_length),
           {'expected': expected_length, 'actual': actual_length,
            'difference': actual_length - expected_length,
            'tolerance': length_tolerance(expected_length)})

    # 每原段的片参数必须从 0 铺满到 1 且不重不漏；每个原段都必须有片
    by_segment = {}
    for piece in pieces:
        by_segment.setdefault(piece['segment_index'], []).append(piece)
    parameter_ok = len(by_segment) == max(1, len(vertices) - 1)
    parameter_detail = {}
    for index, group in by_segment.items():
        ordered = sorted(group, key=lambda item: item['t0'])
        spans = [(item['t0'], item['t1']) for item in ordered]
        ok = (abs(ordered[0]['t0']) <= PARAM_EPS
              and abs(ordered[-1]['t1'] - 1.0) <= PARAM_EPS
              and all(abs(b[0] - a[1]) <= PARAM_EPS for a, b in zip(spans, spans[1:])))
        parameter_ok = parameter_ok and ok
        parameter_detail[index] = spans
    record('parameter_coverage', parameter_ok, parameter_detail)

    # 时间守恒：片时间差之和等于原折线总时间跨度
    expected_span = vertices[-1]['time'] - vertices[0]['time'] if vertices else 0.0
    actual_span = sum(piece['time1'] - piece['time0'] for piece in pieces)
    time_tolerance = length_tolerance(abs(expected_span))
    record('time_conservation',
           abs(actual_span - expected_span) <= time_tolerance,
           {'expected': expected_span, 'actual': actual_span,
            'difference': actual_span - expected_span,
            'tolerance': time_tolerance})

    # 共边不重复计长：同一原段内任意两片的参数区间不得实质重叠
    overlap_ok = True
    overlap_detail = []
    for index, group in by_segment.items():
        ordered = sorted(group, key=lambda item: item['t0'])
        for earlier, later in zip(ordered, ordered[1:]):
            gap = later['t0'] - earlier['t1']
            if gap < -PARAM_EPS:
                overlap_ok = False
                overlap_detail.append({'segment_index': index, 'overlap': -gap})
    record('no_double_counting', overlap_ok, overlap_detail)

    # 状态合法值
    valid_states = {STATE_INSIDE, STATE_OUTSIDE, STATE_EXCLUDED, STATE_UNCERTAIN}
    record('valid_states', all(piece['state'] in valid_states for piece in pieces),
           sorted({piece['state'] for piece in pieces}))
    return {'path_id': path_result['path_id'], 'all_ok': all(item['ok'] for item in checks),
            'checks': checks}


def process_bundle(transform_entries, region_entries, path_entries):
    """离线加工入口：校验资料、解析变换链、切段并做守恒验证。

    单条折线的链失败只拒绝该条；区域几何非法则对应地图版本的全部折线拒绝。
    """
    catalog = TransformCatalog(transform_entries)

    region_failures = []
    invalid_versions = set()
    regions_by_version = {}
    for entry in region_entries:
        try:
            region = validate_region(entry)
        except GeometryError as exc:
            region_failures.append({'region_id': entry.get('id'), 'reason': str(exc)})
            invalid_versions.add(entry.get('map_version'))
            continue
        regions_by_version.setdefault(region['map_version'], []).append(region)

    edges_by_version = {
        version: [edge for region in regions for edge in region_edges(region)]
        for version, regions in regions_by_version.items()
    }

    results = []
    for entry in path_entries:
        path_id = entry.get('id', '<unknown>')
        try:
            points = validate_polyline(entry)
            version = entry['version']
            if version in invalid_versions:
                raise ChainError('地图版本 %s 存在非法区域，拒绝加工' % version)
            matrix, chain_ids = catalog.resolve(entry['crs'], version)
            vertices = normalize_polyline(points, matrix)
            regions = sorted(regions_by_version.get(version, []),
                             key=lambda region: (region['priority'], region['id']))
            edges = edges_by_version.get(version, [])
            pieces = []
            for index in range(len(vertices) - 1):
                p0, p1 = vertices[index], vertices[index + 1]
                radius = _segment_error(vertices, index)
                for piece in cut_segment([p0['x'], p0['y']], [p1['x'], p1['y']],
                                         edges, radius, regions):
                    pieces.append(_enrich(piece, vertices, index, path_id))
            # 单点折线：退化为零长片，给出确定归属；落在边界上归待核带
            if len(vertices) == 1:
                point = [vertices[0]['x'], vertices[0]['y']]
                region_id, state = _point_state(point, regions)
                pieces.append({
                    'path_id': path_id,
                    'segment_index': 0,
                    't0': 0.0, 't1': 1.0,
                    'region_id': region_id,
                    'state': state,
                    'start': list(point),
                    'end': list(point),
                    'time0': vertices[0]['time'],
                    'time1': vertices[0]['time'],
                    'length': 0.0,
                    'error_radius': vertices[0]['error'],
                })
            path_result = {
                'path_id': path_id,
                'accepted': True,
                'map_version': version,
                'crs': entry['crs'],
                'transform_chain': chain_ids,
                'matrix': matrix,
                'normalized_vertices': vertices,
                'pieces': pieces,
            }
            path_result['verification'] = verify_conservation(path_result, vertices)
            results.append(path_result)
        except (ChainError, GeometryError) as exc:
            results.append({'path_id': path_id, 'accepted': False, 'reason': str(exc),
                            'pieces': [], 'verification': None})

    accepted = [item for item in results if item['accepted']]
    return {
        'region_failures': region_failures,
        'paths': results,
        'summary': {
            'total': len(results),
            'accepted': len(accepted),
            'rejected': len(results) - len(accepted),
            # 仅统计被接受的路径；被拒绝路径不产出片段，其 verification 为 None
            'all_verifications_ok': all(
                item['verification'] is not None and item['verification']['all_ok']
                for item in accepted),
        },
    }
