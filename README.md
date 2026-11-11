# 矿区空间片段的坐标归一与边界裁切数据包

矿区二维几何和坐标目录资料工具，提供仿射矩阵、线段参数运算，以及离线空间归一与边界裁切组件。
数据文件位于samples，字段和边界见docs/domain.md，边界与容差约定见docs/boundary.md，持久目录接口见docs/storage.md。

## 本地运行

Python 3.11以上，仅依赖标准库，无需真实账号或外部数据服务。

```bash
python -m datapack.cli inspect samples
python -m datapack.cli segment samples
python -m datapack.cli segment samples --output output/report.json
python -m unittest discover -s tests -v
python -m compileall -q datapack tests
```

第一条命令检查已有资料；第二条执行坐标链解析、空间切分与守恒验证并打印报告；
第三条把报告原子写到本地output目录（已排除出Git）；第四条运行读取、存储、领域原语、
变换链、多边形与裁切测试；第五条构建Python字节码。

## 加工规则摘要

- 变换链：按折线声明的坐标系与地图版本枚举全部可达仿射链，多链结果不一致或缺路径即拒绝该条，
  不依据坐标数值猜测毫米/米或轴方向。
- 区域：含孔洞的简单多边形，先做简单性校验；外环边界与孔环边界属待核带，孔洞内部为排除区，
  退化与自交多边形有明确拒绝结果。
- 切段：在所有区域边界处精确切段，保留切点在原线段上的参数与时间插值；误差半径段内取较大者，
  距任一边界不大于该半径的部分单列为uncertain，其余记inside/outside/excluded。
- 归属：重叠确定部分按priority归属，共边不重复计长；一条折线可多次进出同一区域。
- 守恒：每片输出区域、状态、原段引用与长度，片长加总等于归一后原折线长度，
  验证器同时核对参数覆盖、时间守恒与重复计长。

项目只处理合成平面几何，不给真实矿卡下发路线。

## 容器

```bash
docker build -t local-datapack .
docker run --rm --network none local-datapack
docker run --rm --network none local-datapack python -m unittest discover -s tests -v
```

应用没有必须常驻的服务，入口和库可在同一容器离线运行。state和output是本地工作目录，已排除出Git。样本文件只读，修订目录使用独立SQLite文件，旧版本可继续查询。

## 目录

datapack/storage.py提供摘要、原子文件与修订目录；datapack/records.py提供文件读取及重复主键分类；
datapack/domain.py提供基础领域原语；datapack/chains.py、datapack/polygons.py、datapack/segmentation.py
提供坐标链解析、多边形校验与空间裁切；datapack/tolerances.py集中全部数值容差；tests验证这些行为。
数据加工结果不由资料检查命令自动生成。
