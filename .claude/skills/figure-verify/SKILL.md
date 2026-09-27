---
name: figure-verify
description: 图像双重验证：像素级机器断言（防缓存欺骗人眼）+ 图↔数据逐值核对。触发词：验证图、figure check、图片审查、像素断言、图和数据对不上。
---

# 图像验证（Figure Verification）

## 方法
1. **像素级机器断言**（生成后必做）：
   ```python
   from PIL import Image; import numpy as np
   img = np.asarray(Image.open(f).convert("RGB")).astype(int)
   near = lambda rgb, tol=50: (np.abs(img-np.array(rgb)).sum(axis=2) < tol).sum()
   # 断言：每条系列色像素数 > 阈值；关键区块（深紫中心/亮黄表层）存在
   ```
2. **图↔数据逐值核对**：柱高/点坐标/标注数字 vs 权威 JSON（Read PNG 目检 + 数据复算双轨）；
3. **口径一致性**：色标共享 norm、聚类口径、单位。

## 已知陷阱（实战教训）
- **图片渲染缓存两级欺骗**：VSCode 预览标签页 AND Read/CDN 渲染都可能给旧图——曾经"目检通过"的图实际是 26px 色斑。**规则：像素断言为准，人眼为辅**；
- matplotlib tripcolor 传「全量坐标+子集三角形」会让轴按全域定标（导线缩成色斑）——子集必须重映射节点+显式 xlim/ylim；
- tripcolor 单元均值会稀释峰值（31.7→9.3 A/mm²）——用顶点 gouraud 着色；
- tripcolor 子集三角各自 clim → 同一半径两种颜色——共享 Normalize；
- 插值器 API 差异被裸 except 吞掉后走退化分支（画出了 |A| 假曲线）——删裸 except，加 `assert jr.max()>阈值`；
- 中文标题里的特殊 Unicode（↔）缺字形显示方框——用 ASCII 替代；
- 图更新后同步刷新所有快照副本（web/public/figures 等）。

## 输出格式
每图一行：`文件名 尺寸 色px=[...] ✅/❌` + 数据核对结论。
