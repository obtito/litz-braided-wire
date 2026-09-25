# 文献检索与核验记录（2026-09-25，H+7）

> 工具：citation-management skill（OpenAlex 检索 + CrossRef 元数据）。检索面：6 路主题查询 + 6 次精确标题，原始 120 条 → 去重 109 → 相关 64 → 精选 15 篇入 [litz.bib](litz.bib)。**论文引用一律从 litz.bib 取，不再手写。**

## 一、存疑引用核验结果（RESEARCH §十，1–3 项全部解决）

| 原存疑项 | 核验结果 |
|---|---|
| Liu 2025 IJAMT 卷期页码 | ✅ **138(11-12): 5877-5890**, 2025-06；作者 Liu Jiale, Gao Yicen, Shan Zhongde, Sun Zheng, Guo Zitong, Zhu Xiangyu；DOI 10.1007/s00170-025-15743-0 |
| ryz.ece.illinois.edu 镜像 ×2 | ✅ 改引 IEEE DOI：Sullivan 1999 → 10.1109/63.750181；Sullivan & Zhang APEC 2014 → 10.1109/apec.2014.6803681；另 Zhang/White/Kassakian APEC 2014 → 10.1109/apec.2014.6803390 |
| Rosskopf（ResearchGate 页） | ✅ IEEE TPEL：10.1109/TPEL.2013.2293847。续作备选：*Litz wire loss performance and optimization for cryogenic windings*（10.1049/elp2.12279） |
| Ferreira（顺带核验） | ✅ 1994 TPEL：10.1109/63.285503；1992 IEE Proc-B：10.1049/ip-b.1992.0003（**真实存在**，非幻觉引用） |

未处理项（政策性，非 DOI 问题）：Pyrhonen 引 Wiley 正式版；ARS README 统计不用；TU Delft/Dowell 只引不依赖细节。

## 二、新入库文献（15 篇，按用途分组）

**经典损耗模型（Q1/Q2 方法链引言）**
| 论文 | 年 | 引用 | DOI |
|---|---|---|---|
| Ferreira, Improved Analytical Modeling of Conductive Losses | 1994 | — | 10.1109/63.285503 |
| Ferreira, AC resistance of round and rectangular litz windings | 1992 | — | 10.1049/ip-b.1992.0003 |
| Tourkhani/Viarouge 系，Accurate analytical model of winding losses in round Litz wire | 2001 | 231 | 10.1109/20.914375 |
| Nan & Sullivan 系，Winding resistance of litz-wire and multi-strand inductors | 2012 | 217 | 10.1049/iet-pel.2010.0359 |

**丝径/根数设计（Q2 决策链）**
| 论文 | 年 | 引用 | DOI |
|---|---|---|---|
| Sullivan, Cost-constrained selection of strand diameter and number | 2001 | 119 | 10.1109/63.911153 |
| Sullivan & Zhang, Simplified Design Method for Litz Wire | 2014 | — | 10.1109/apec.2014.6803681 |
| Analytical Design Methodology for Litz-Wired HF Power Transformers | 2014 | 99 | 10.1109/tie.2014.2351786（候选，未入 bib） |

**绞合/换位效应（Q3 核心）**
| 论文 | 年 | 引用 | DOI |
|---|---|---|---|
| Analytical model for effects of twisting on litz-wire losses | 2014 | 92 | 10.1109/compel.2014.6877187 |
| Umetani 等，Copper Loss of Litz Wire With Multiple Levels of Twisting | 2021 | — | 10.1109/tia.2021.3063993 |
| Rosskopf 等，Inner Skin- and Proximity Effects in Litz Wires | 2014 | — | 10.1109/tpel.2013.2293847 |
| Fast Numerical Power Loss Calculation for HF Litz Wires | 2020 | 68 | 10.1109/tpel.2020.3008564 |

**换位结构 / Roebel / EV（Q4 + why-now 弹药）**
| 论文 | 年 | 引用 | DOI |
|---|---|---|---|
| Transport & magnetization AC losses of ROEBEL assembled coated conductor cables | 2009 | 97 | 10.1088/0953-2048/23/1/014023 |
| Hybrid Transposed Hairpin Winding for EV Traction | 2022 | 68 | 10.1109/tie.2022.3179571 |
| Liu 等，Graph theory in carrier path research for 3D rotary braiding machines | 2025 | — | 10.1007/s00170-025-15743-0 |

## 三、剩余候选（64 篇相关中的其余高引，未入 bib）

完整清单存 `/tmp/litsearch/merged_relevant.json`（临时）。值得后续注意：等效复磁导率 litz 均质化（10.1109/tia.2009.2013594，146 cit——ethz-pes 均质化路线的理论源头）、平面 litz 结构（10.1109/tpel.2004.843022，109 cit）。

> 检索可复现：6 个查询词见本文件历史版本；`search_openalex.py --limit 20`；核验走 `doi_to_bibtex.py`（CrossRef）。
