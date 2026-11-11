"""集中配置的数值容差。

所有“允许的数值容差”只能在此定义，边界判定的含义随容差一起在本模块说明，
docs/boundary.md 与本文件保持一致。资料约定的几何绝对容差为 1e-8（米，地图坐标）。
"""

# 点落在直线/线段上、环自交、孔与外环关系等几何判定的绝对容差（地图单位：米）
GEOMETRY_EPS = 1e-8

# 仿射矩阵逐元素比较的绝对容差；用于判定多条变换链是否给出同一变换
MATRIX_EPS = 1e-9

# 几何判定的相对容差：按资料坐标尺度（环包围盒量级）放缩，避免大坐标下 1e-8 过松
GEOMETRY_REL_EPS = 1e-10

# 切段参数 t 去重、区间合并时的参数容差
PARAM_EPS = 1e-10

# 长度守恒校验的绝对容差（米）；另设相对容差覆盖长线
LENGTH_ABS_EPS = 1e-7
LENGTH_REL_EPS = 1e-9

# 二分法求 |边界距离|=误差半径 交点时的参数精度上限
BISECT_LIMIT = 60

# 仿射变换近奇异（行列式绝对值）阈值：行列式为 0 时无法定义反向单位与面积
SINGULAR_DET_EPS = 1e-15


def geometry_eps(scale):
    """按几何尺度返回点在线上、环接触等判定的实际绝对容差（米）。

    scale 为参与判定的坐标量级（如环包围盒的最大边长）；保留 GEOMETRY_EPS
    作为绝对下限，保证小坐标下仍可分辨 1e-8 级差异。
    """
    return max(GEOMETRY_EPS, abs(scale) * GEOMETRY_REL_EPS)


def length_tolerance(reference_length):
    """按参考长度返回守恒校验实际使用的容差。"""
    return max(LENGTH_ABS_EPS, abs(reference_length) * LENGTH_REL_EPS)


class ChainError(ValueError):
    """坐标变换链无法唯一定位（缺路径、多链不一致或矩阵非法）。"""


class GeometryError(ValueError):
    """区域几何非法（退化、自交、孔越界等），该区域不产生有效结论。"""
