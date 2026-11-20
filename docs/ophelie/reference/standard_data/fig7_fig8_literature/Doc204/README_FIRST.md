# Doc204：法国 Fig.7 / Fig.8 数字化与下一轮任务包

日期：2026-11-01。

## 先回答状态

上一轮只审核了图号，没有生成数字化数据。本包现在交付的是**有可追溯原图和叠图核对的第一版数字化数据**；不是作者原始 CFD 数据，也不是完整的绝对温度/速度场。

直接把本目录和 `204_CURSOR_NEXT_PLAN.md` 一起交给 Cursor；无需它重新做同一批图像提取。

## 文件与用途

| 文件 | 含义 | 可以做什么 | 不可以做什么 |
|---|---|---|---|
| `reference/fig8_labelled_contour_points.csv` | 4928个沿可读0.5–0.9等温线取得的点 | 比较等温线位置、约束温度重建 | 当作4928个独立实验观测；当成Kelvin |
| `reference/fig8_midline_contour_intersections.csv` | 10个中央剖面与已读等温线交点 | 稀疏剖面对比 | 称为完整温度曲线 |
| `reference/fig8_partial_linear_reconstruction.csv` | 101×101网格，观测凸包内线性插值 | 原型后处理、重建示意 | 直接给Kelvin输入器；当成完整文献场 |
| `reference/fig7_direction_samples.csv` | 20处可辨认箭头的人工近似方向 | 主流方向对照 | 当作速度大小、精确角度或完整速度场 |
| `reference/fig7_topology_regions.csv` | 主次涡中心搜索框及上下行区域 | 宽容差的流型定位 | 将搜索框中心当成精确停滞点 |
| `reference/digitization_provenance.json` | 来源、标定、误差和限制 | 自动加载时的检查 | 忽略unknown字段然后硬写温标 |
| `qa/*.png` | 提取点/方向叠回原图、部分重建图 | 人工核对 | 代替计算结果图 |
| `scripts/digitize_figures.py` | 从随包原图重新生成数字化数据 | 可复现与审计 | 运行时悄悄更改数据定义 |

## 来源和标定

- Jacoutot et al., *Numerical modelling of natural convection in molten glass heated by induction*, Chemical Engineering and Processing 47 (2008) 449–455。
- DOI：10.1016/j.cep.2007.01.013。
- Fig.7/8 均在 PDF第5页、期刊页453。原PDF SHA256见provenance。
- 原图为PDF内嵌的黑白栅格图；本包用原始分辨率提取，反转亮度以恢复黑线白底。未使用生成式图像。
- 横轴r/R直接标定；纵轴按图域顶部z/H=1、底部0标定，该图没有数值纵轴刻度，此处采用完整浴高的几何解释。
- Fig.8 为 `Tstar=T/Tmax`，**不是** `(T-Tcold)/(Thot-Tcold)`。
- R=0.25m为已确认几何；H=0.185m来自博士补充几何，绝非从这张图直接读出的米制尺寸。
- 坐标误差是保守人工读图容差，不是统计置信区间。图上没有给出速度标尺，因此速度大小留空。

## 必须保留的限制

1. Fig.8只提取印有数值且能可靠关联的0.5、0.6、0.7、0.8、0.9曲线。没有给无标签外侧线猜数字。
2. 同一曲线被文字截断时，保留两个观测片段；不把缺口填成观测点。
3. 线性重建是辅助产物，原图没有提供线间唯一解；0.9内部热核的温度高于0.9但具体值未恢复，域边界附近也有未覆盖区。
4. 原文“约1500K”不是精确Tmax。`absolute_Tmax_K=null`；侧壁换热工况对应关系也留空，禁止猜h=300。
5. Fig.7角度按归一化图平面读取并取15度步长，误差约20–30度；物理方向列另按R/H换算，不能把显示角度直接当物理速度角度。
6. 两个涡心框是视觉搜索区域，不是对矢量场解零点得到的坐标；不得使用其几何中心作严格误差验收。
7. 原图标注点高度相关。加密采样并没有增加文献信息量。
8. 用重建温度强迫F1后，可以检验流动响应；不能再拿同一温度场宣称热模块得到独立验证。

## 本包不包含什么

没有伪造 `literature_T_rz.csv` 绝对Kelvin全域文件。下一步重建完整F1输入，需要独立的绝对温标和边界/热核补全假设；按任务书保留为 `RECONSTRUCTED_ASSUMED_INPUT`。这不阻塞对现有F1结果自身是否有闭合环流的检查。

## 重生成

```bash
python3 scripts/digitize_figures.py
python3 scripts/validate_reference.py
```

依赖：Python3、numpy、scipy、matplotlib、Pillow。预期4928个等温线点、20个方向、10个中央剖面交点。方向和ROI为人工读取的版本化标注，不是图像自动反演。
