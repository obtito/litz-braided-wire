# 千丝万缕 · Every Strand, Every Position

> 三维拓扑编织 Litz 线电磁优化设计 —— 48 小时数学建模黑客松
> 3D Topologically Braided Litz Wire: Coupling Braid Topology to AC Loss

## 概述

200 kHz、20 A 工况下，实心铜线的交流电阻是直流的 **3.65 倍**（Bessel 精确解 3.647 / 有限差分 3.647 / 2D FEM 3.649，三法互证，δ=147.77 μm）。

直觉认为「细丝成束」能救——六角密排下恰恰相反：**不换位的 Litz 束比等截面实心线还差 2.7%**（Rac/Rdc = 4.653 vs 4.5295，N=127 完整六角束，丝级 FEM）。外层丝链磁通少、阻抗低、抢流 4 倍（η_I=438%），环流把细丝的趋肤优势全部吃掉；且丝径再细也无济于事——失衡由束径决定，扫描证实 d→0.17 mm 后比值平台在实心线水平。

但**排布几何本身就是设计变量**：同样的 331 丝 / 0.14 mm / 同铜面积，只把六角密排换成同心圆环（环内方位对称），环内电流失衡即压到 0.3% 以内，Rac/Rdc = **4.403，反超等面积实心 3.6%**（队友 COMSOL 4.4033；我方 skfem 独立复现 4.408，边界收敛后差 0.12%）。不过距单丝理想值 1.001 仍有 4.4 倍空间——这剩余部分只能靠**真实径向换位**消除，正是 Q3/Q4 的战场。

唯一出路是**真实径向换位**：让每根丝沿长度轮流占据所有径向位置。本项目用辫群语言刻画编织拓扑（braid word → 纯辫/传递性/φ-矩均衡三判据，D→0 即完美换位），用「多站切片 + 编织角修正」的 2D FEM 链条计算任意编织方案的损耗，建立**编织拓扑 ↔ 交流损耗**的定量桥梁——并把全部结果放进可交互网站，任何人可以亲手检验。

## 核心数字（全部可溯源至 findings.md）

| 结论 | 数值 | 出处 |
|---|---|---|
| Q1 实心线锚点（三法互证） | 3.64720 / 3.6473 / 3.6487±0.001 | findings.md H+2, H+4 |
| 等截面实心基准（2.523 mm） | 4.5295（Bessel）/ 4.5310（FEM） | H+10 |
| 未绞合 Litz 束（N=127, d=0.2239 mm） | **4.653（refine=3）** | H+10 |
| 未绞合 vs 实心（六角） | **+2.7%（更差）** | H+10/H+13 |
| 圆形环排布（同铜面积） | **4.403（反超实心 3.6%）** | 队友 COMSOL + 我方复现 4.408（H+15） |
| 单丝自趋肤 | 1.0068（我方设计点）/ 1.0010（N=331 设计） | H+7 |
| 电流失衡 η_I | 438–441% | H+10 |

## 方法链

Bessel 解析锚点 → 丝级 2D FEM（gmsh 距离场网格 + scikit-fem，每丝独立导体 + 并联总流约束）→ 多站切片法（Roebel 式周期平均 + 1/cos α 长度修正）→ 小 N 全 3D 周期模型对照 → 拓扑判据 D 与损耗的相关性（进行中）。

## 仓库结构

```
q1_solid/   Q1 实心线：radial_fd.py（FD）、fem2d_skfem.py（FEM）、make_figs.py
q2_litz/    Q2 Litz 束设计与仿真
sim/        公共仿真库（网格/求解/后处理）
data/       全部仿真输出 CSV/JSON（证据链）
paper/      论文（q1.md 已定稿）与 refs/（15 篇核验文献 + PDF）
geo/        几何生成器（绞合/多级/管状/track-and-column）
web/        Q4 成果网站（VitePress + three.js）
roadshow/   路演讲稿 PITCH.md（数字实时取自 findings.md）
BATTLE_PLAN.md / RESEARCH.md / ROADSHOW.md / findings.md
```

## 快速开始

```bash
python3 -m venv .venv && source .venv/bin/activate
pip install numpy scipy matplotlib pandas meshio scikit-fem gmsh
python q1_solid/fem2d_skfem.py        # 复跑 Q1：Rac/Rdc → 3.6487
```

## 状态

- [x] Q1 实心线（三法互证，网格收敛研究，报告 paper/q1.md 定稿）
- [x] Q2 Litz 束（设计定型 N=127/d=0.2239，未绞合基线 + 丝径权衡扫描，反直觉结论）
- [ ] Q3 编织方案调优（直束/单级/多级 S/Z/管状/真三维，A–E 梯度对比）
- [ ] Q4 完美编织（拓扑判据 D、braid word 设计、可制造性）+ 成果网站

文献引用一律取自 [paper/refs/litz.bib](paper/refs/litz.bib)（15 篇 CrossRef 核验）。
