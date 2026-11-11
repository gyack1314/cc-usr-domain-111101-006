# 数据约定

使用二维齐次仿射矩阵[a,b,tx,c,d,ty]，目标点为(ax+by+tx,cx+dy+ty)，单位由变换本身编码。地图坐标为米；外环与孔环不依赖顺逆时针确定意义。折线每点time为秒，error为米；几何绝对容差1e-8，可配置。所有多边形均须先做简单性校验。

样本为合成资料，包含正常、缺失、冲突等记录；这些现象用于保留业务数据原貌，读取成功不代表领域结果有效。所有计算原语位于datapack/domain.py；调用者负责组合、约束和最终验证。

## 离线空间归一与裁切

新增组件在合成平面几何上完成坐标归一、边界切段与守恒验证，不向真实设备下发路线：

- datapack/tolerances.py：集中容差与异常类型，边界含义见docs/boundary.md；
- datapack/chains.py：按版本在有向变换图上枚举全部简单路径，多链不一致或缺路径即拒绝该条折线，不依据坐标数值猜测单位；
- datapack/polygons.py：含孔洞简单多边形的简单性校验、点归属（inside/outside/excluded/boundary）与边界距离；
- datapack/segmentation.py：归一、在所有区域边界处精确切段、误差半径待核带（uncertain）单列、重叠区按priority归属、长度/参数/时间守恒验证；
- 输出每片的区域、状态、原折线与原段引用（path_id、segment_index、t0、t1）、端点时间与长度，片长加总等于归一后原折线长度。

CLI入口：`python -m datapack.cli segment samples [--output 路径]`。
