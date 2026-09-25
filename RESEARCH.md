# 三维拓扑编织 Litz 线黑客松 —— 研究汇总与决策文档

> 生成于 2026-09-25。由 13 个并行研究 agent（9 角度扫描 + 完备性批评）+ 人工甄别整合而成。
> 72 条发现中 64 条已逐一验证（README/文档/论文实际读取）；未验证条目已标注。
> 原始逐条记录见 [研究原始汇总.md](研究原始汇总.md)。

---

## 〇、十条执行决策（TL;DR）

| # | 决策 | 依据 |
|---|------|------|
| 1 | **主求解器：NGSolve**（`pip install ngsolve`，macOS universal2 wheel 已实证，venv py3.12 可直接装） | PyPI 实测 v6.2.2607 有 cp310–312 macOS wheel；复数时谐 H(curl) 边单元原生成MQS 涡流问题；全 Python 驱动，Q3 优化循环在进程内秒级重解 |
| 2 | 备选：**scikit-fem**（纯 Python，全掌控）→ **GetDP+Gmsh**（手册自带完整 2D a-v 涡流算例 `MagDyn_av_2D.pro`）→ **xfemm**（FEMM 复刻，C++ headless）。FEMM 官方仅 Windows（Wine 可跑，Oct2023 版改善兼容） | 各 README/下载页已验证 |
| 3 | **仿真口径：2D 横截面单位长度相量模型**（每丝独立导体 + 并联约束「各丝每米压降相等、总流 20 A」），绝不做 1 m 的 3D 模型 | GetDP 手册文档化的标准做法 |
| 4 | **绞合/编织 → Rac 的方法学 = 「多站切片法」**：沿一个结构周期取 S 个 z 站 → 每站 2D 涡流解 → 按 dz 加权平均 + 股线加长修正 1/cos α（详见 §三） | Roebel 电缆文献标准近似；填补批评缺口① |
| 5 | **Rdc 口径：Rdc = ρ/(A_cu·cos ᾱ)（每股实际长度）**，所有方案对比时同时报告「每米电缆 Rdc、Rac、Rac/Rdc」三个数 | 填补批评缺口②（IEC 60228 思路） |
| 6 | Q2 设计：**双层设计**——「分析最优设计」d=0.05–0.07 mm（Sullivan 准则，N≈1300–2500，只做解析评估）+「仿真可解设计」d=0.2 mm、N=159（β=1.35，单丝趋肤仅 +0.4%，丝级 2D FEM 可解） | Sullivan & Zhang APEC2014（全文已读）+ 自算锚点 |
| 7 | Q3 至少三方案梯度：**简单绞合 → 多级换向成缆（S/Z）→ 真三维编织（track-and-column 4 步法）**；可加管状编织作第 4 对照（预判其换层有限——这本身是个反直觉论点） | Rosskopf：单级绞合"内丝恒内"；管状编织=定半径双螺旋（题目暗示有误导，见审查） |
| 8 | Q4 判据核心：**辫群形式化 + φ-矩相等**：交叉序列=braid word w；完美 ⇔ ①w^p∈P_n（纯辫，各丝归位）②⟨π(w)⟩≤S_n 传递 ③各丝对场基 {r, rcosθ, rsinθ, r², …} 的停留矩相等（离散=Latin 方/循环赛构造），报告标量差异 D | Wikipedia 辫群（全文读）+ Pyrhonen 经典定义 + Umetani 最小捻距条件 |
| 9 | Q4 网站：**Vite + 原生 three.js**（勿用 R3F，48h 内部件越少越好）；`TubeGeometry`+自定义 `Curve.getPoint(t)` 渲染股线，`mergeGeometries` 合批，`clippingPlanes` 做截面动画；**VitePress(MathJax) 内嵌 viewer** 一个仓库同时承载报告+演示；GitHub Pages 部署（记得 `base:'/<repo>/'`） | 全部 API/文档已验证 |
| 10 | 论文引擎与流程：`OthmanAdi/planning-with-files`（迭代日志持久化，恰好满足 Q3「保留调优记录」）+ 自写 FEM 技能；引用一律核对 DOI（存疑清单见 §十二） | skills 角度已验证 |

---

## 一、赛题审查结论（误导点甄别）

（详见会话记录，此处存档关键数字）

- **δ(200 kHz, Cu) = 147.77 μm**；扫频锚点 Rac/Rdc（2 mm 实心线，精确 Bessel 解）：
  50 k→1.966；100 k→2.662；**200 k→3.647**；500 k→5.609；1 M→7.822。L_int(200 kHz)=14.7 nH/m（DC 50 nH/m）。**这就是 Q1 仿真的收敛目标**。
- **陷阱 1**：Q1 的 2 mm 线 J_dc=6.37 A/mm² > 4——J 约束从 Q2 起生效（A_cu ≥ 5 mm²，等效实心 D=2.523 mm，其 Rac/Rdc=4.529）。Q2「等截面对比」应做两个基准（2 mm 与 2.523 mm）并在论文说明口径。
- **陷阱 2**：「绞合永不换层」对单级绞合成立，但商用 Litz 的多级 S/Z 成缆有部分换位；且**经典管状编织（maypole）股线半径同样近似不变**——真换层需要 track-and-column/多层旋转编织。Q3 方案设计勿踩此坑。
- **陷阱 3**：记号冲突（β=d/δ vs braid word β；σ 电导率 vs 生成元 σ）——论文重命名为 φ 和 B。
- **陷阱 4**：「β<1」是保守法则：d=0.2 mm（β=1.35）单丝趋肤因子仅 1.004；最优丝径在 d≈δ 附近（Sullivan）。Q2 做定量权衡而非教条。
- **陷阱 5**：「忽略漆膜厚度」≠几何可重叠：圆堆积上限 π/(2√3)=0.907；束外径 ≈ d·√(N/0.907)；仿真中每丝必须是独立导体。

---

## 二、工具链（macOS 原生、全开源、已验证）

### 2.1 电磁求解器优先级

| 工具 | 安装 | 角色 | 关键事实 |
|---|---|---|---|
| **NGSolve** | `pip install ngsolve` | 主力：2D/3D 复数时谐 A-v 涡流 | macOS universal2 wheel（PyPI 实证）；自带 Netgen 网格；NGS-Py 全脚本化；i-tutorials 有 Maxwell/边单元教程 (ngsolve.org/docu/latest/i-tutorials) |
| **scikit-fem** | `pip install scikit-fem[all]` | 备选：纯 Python 装配，复数双线性型，Nédélec 单元 | v12.0.2 OS-Independent；无现成涡流算例需自写弱式（Poisson/Helmholtz 模板改） |
| **GetDP+Gmsh** | `brew install gmsh`；getdp.info 下载 | 快速验证：手册自带 `MagDyn_av_2D.pro` 完整算例（Current_2D 总流约束、roj2 损耗后处理、U/I 阻抗） | GPL；macOS 有预编译 |
| **xfemm** | 源码 cmake | FEMM 4.2 复刻（planar+axisymmetric，含 AC 涡流），headless `fmesher/fsolver` | crobarcro/xfemm，GPL-3.0，2026 仍活跃 |
| FEMM 4.2+pyfemm | Windows/Wine | 行业标准 2D 磁（若临时拿到 Windows 机器） | 官方无 mac 版（下载页实证）；planar depth=1 m 直出每米 Rac |
| PyPEEC | `pip install pypeec` | 3D 准静磁 PEEC（FFT 加速），Dartmouth/Sullivan 系出品 | 体素几何；对大量细丝昂贵，适合 3D 小 N 验证 |
| ElmerFEM | 源码/Docker | 不推荐（macOS 无官方包） | WhitneyAVSolver 能力足够但安装成本高 |
| PyAEDT | `pip install pyaedt` | 仅当 ANSYS license 出现时启用 | MIT 包装器，license 在 AEDT 侧 |

### 2.2 几何/数学/可视化代码资产

| 仓库 | 许可 | 用途 |
|---|---|---|
| **otvam/litz_wire_losses_twisting** | BSD-2, MATLAB 无工具箱 | **Q3/Q4 核心数学**：绞合=股阻抗矩阵上的置换矩阵 → 并联均流求解 → Bessel 解析丝损；含 straight/twisted/random 三脚本；COMPEL2017 论文 PDF 附带。逻辑可 1 天移植 Python |
| **ethz-pes/litz_wire_losses_fem_matlab** | BSD-2 | 粗网格场积分 (∫J², ∫H²) → Bessel 丝损，网格无关趋肤深度，MHz 级有效——快速代理模型 |
| **ethz-pes/litz_wire_homogenization_comsol_matlab** | BSD-2 | 复数等效 μ/σ 均质化（<2% 误差至 MHz）；需 COMSOL，**算法思想移植到 NGSolve** |
| **OpenMagnetics/cci_coords** | 无(数据) | N 圆入圆最优堆积坐标——Q2 横截面/逐站截面的现成丝心坐标（基于 Magdeburg packing 库） |
| **DeloongZhang/act-braiding-framework** | Apache-2.0 | **「拓扑→事件→轨迹」参考实现**：M×N 载纱器网格、4 步法运动动画、fiberize 生成 3D 纱线几何、Open3D 可视化、JSON 工程；PySide6/Open3D 全 mac 原生；Abaqus 部分可跳过 |
| DeloongZhang/Virtual-Fiber-Braiding-Simulation-ABAQUS | 无 | 4 步法 quasi-fiber 级仿真（Abaqus 依赖，仅作算法参考） |
| louisepb/TexGen | GPL-2.0 | 纺织几何建模（纱线中心线→实体→体网格）；脚本 API 可用 |
| SageMath `sage.groups.braid` | GPL | 辫群现成机器：`BraidGroup(n)`、Tietze 词、`.permutation()`、纯辫方法、`.plot3d()`；重依赖，用 CoCalc/sagecell 或仅引用其 API 设计 |
| 3-manifolds/Spherogram | GPL-2.0+ | `pip install spherogram`；辫闭包/连接不变量做拓扑健全性检查 |
| rexgreenway/braid-visualiser | MIT | `pip install braidvisualiser`；matplotlib 画 2D 辫图（报告插图） |
| jeanluct/braidlab | GPL-3.0 | MATLAB 辫群包（macOS 二进制）；法式/正则形 |
| marinmersenne2357/openscad-maypole | 无(参考) | 完整 maypole 编织机 OpenSCAD 模型（角齿轮/载纱器/动画）——Q4 可制造性讨论素材 |

### 2.3 网站栈（已验证到 API 级别）

- 脚手架：`doinel1a/vite-three-js`（MIT，2026-09 活跃）`git clone → npm i → npm run dev`。
- 股线渲染：每丝一个 `THREE.Curve` 子类（`getPoint(t)` 返回 (R cos(2πn_t t+φᵢ), R sin(...), L t)）→ `new THREE.TubeGeometry(curve, 256, r_strand, 8, false)`；现成曲线类在 `three/addons/curves/CurveExtras.js`（HelixCurve 等已确认存在）；`BufferGeometryUtils.mergeGeometries` 合批；截面剖切 `renderer.localClippingEnabled + material.clippingPlanes`，随滚动动 plane.constant；`OrbitControls` 自带。
- 文档+演示一体：**VitePress**（MIT）+ `markdown-it-mathjax3@^4` + `markdown:{math:true}`（KaTeX 在 VitePress 有已知构建 bug #229，**别用**）；Vue 组件可直接内嵌 three.js viewer——一个仓库同时是论文站和演示站。
- 部署：VitePress 官方 deploy.yml（configure-pages@v4/upload-pages-artifact@v3/deploy-pages@v4）复制即用；**GitHub Pages 必须 `base:'/<repo>/'`** 否则白屏（经典坑）。备选 Cloudflare Pages（每分支预览）。

### 2.4 Claude Skills（AI 工作流）

| 技能 | 安装 | 用途 |
|---|---|---|
| anthropics/skills（官方） | `/plugin marketplace add anthropics/skills` | web-artifacts-builder、webapp-testing（网站）、skill-creator（自写 FEM 技能）、docx/pptx（提交材料） |
| K-Dense-AI/scientific-agent-skills | `gh skill install …`（166 个, MIT） | 多目标优化（Q3）、符号数学、出版级图表、证据可溯写作 |
| viettranx/3dviz-pro-max | 插件 | Three.js 0.180 十步工作流+**帧捕获自检**（agent 检查自己渲出的帧）——正好治"3D 只能口嗨" |
| OthmanAdi/planning-with-files | 插件 | task_plan/findings/progress.md 磁盘持久+每轮注入，抗 /clear/压缩——**Q3 调优日志天然载体** |
| Imbad0202/academic-research-skills | 插件(CC BY-NC) | 研究→写作→审校管线、引用核查门（有中文 README） |

> 负结果（有价值）：**不存在 FEM/电磁专用 Claude skill**；GitHub 上也没有 K-TLitz/拓扑编织的公开代码——Q4 的「编织拓扑生成器」必须自写，这恰是原创空间。

---

## 三、方法学：从编织几何到 Rac/Rdc 的可计算路径（批评缺口①填补）

**核心矛盾**：2D 横截面模型对「绞合/编织」天然失明——一个 z 站的截面无法区分「直丝束」和「螺旋丝束」。保真度阶梯（成本递增）：

**L0 解析丝损模型**（秒级）
Ferreira/Umetani Kelvin 函数丝因子：
- F(g_s) = g_s(ber·bei′−bei·ber′)/(2(ber′²+bei′²))（趋肤）
- K(g_s) = −2(ber·bei′−bei·ber′)/(g_s(ber²+bei²)) − 1（邻近）
- g_s = d_s·√(πfμ₀σ)；`scipy.special.kelvin` 直接算 ber/bei。
用于 Q2 丝径权衡曲线和一切交叉验证。

**L1 置换矩阵法**（Guillod COMPEL2017，分钟级）——Q3 优化内环
丝坐标 → DC R/L 阻抗矩阵 → 一个结构周期内的换位=置换矩阵序列 → 解并联均流（各丝压降相等）→ 每丝内外场 → Bessel 解析丝损。**绞合/编织方案在这里就是「置换序列」**，与 braid word 的 S_n 像同构。近似：R/L 用 DC 值、丝内场均匀。Python 移植量 ~1 天（或直接跑 MATLAB/Octave）。

**L2 多站切片 2D FEM**（每方案 ~分钟级，推荐主路线）
沿结构周期取 S 站（如 S=16），每站做真 2D 复数涡流解（NGSolve/scikit-fem，丝级网格）：
P̄ = (1/P)∫₀ᴾ P₂D(z)dz ≈ Σ_s w_s·P₂D(z_s)，每站自动给出：总损、各丝电流（→η_I）、J 分布。
修正：①丝长因子 1/cos α（阻值与 Rdc 同乘，Rac/Rdc 部分抵消但不完全）；②截面椭圆化（面积/cos α，小角可忽略）；③忽略轴向场分量与曲率——引用 Sullivan & Zhang COMPEL 2014「twisting effects on litz losses」作为修正依据。**切片状态直接复用为 Q4 网站的逐站截面动画帧**。

**L3 全 3D 周期模型**（小时级，仅小 N 验证）
取 1 个结构周期 + 周期边界，丝级网格（N≤30）或均质化复数 μ/σ（任意 N，ETH 算法移植）。NGSolve 3D 或 PyPEEC。用于证明 L2 误差可接受（如 <5%）。

**推荐配置**：Q1 用 L2 单站（实心线即 1 站）；Q2 用 L2（直束，S=1）+ L0 交叉验证；Q3 优化用 L1+L2 组合（L1 扫参数空间、L2 复核决赛方案）；L3 对最优方案做一次小 N 版验证。

---

## 四、Q1：理论与数值锚点（已算好，直接用）

- 精确解（Wikipedia Skin effect 全文已核）：J(r)=J(a)·J₀(kr)/J₀(ka)，k=√(−jωμσ)=(1−j)/δ（e^{+jωt} 约定）；Z_int/m = kρ·J₀(ka)/(2πa·J₁(ka))；DC 极限 L_int=μ₀/8π=50 nH/m；HF 极限 Rac/Rdc→a/(2δ)；98% 电流在 4δ 内。
- 数值（已用 venv+scipy 复算）：δ=147.77 μm；a/δ=6.767；Rac/Rdc=3.647；Rdc=5.488 mΩ/m；P=8.01 W/m。
- **网格收敛规则**（COMSOL 官方指南，已读）：导体内单元尺寸 ≤δ；边界层首层 ~δ、增长率 1.5–2.5、总厚 ~导体尺寸/10；圆周 ≥8 个二阶单元；阻抗 BC 仅当 δ≪尺寸 ~100×（本题 a/δ=6.8 不适用，必须体网格）。收敛性研究=扫「每 δ 单元数」(1,2,3,4) 画 Rac/Rdc 收敛曲线 → 目标 3.647。
- 验证图配方：FEM |J(r)| vs J₀(kr)/J₀(ka) 曲线叠加；J(r) e 折叠长度 vs δ；L_int 50→14.7 nH/m 趋势。
- Bessel 数值：`scipy.special.jv(0, 复数)` 直接支持，无需 Kelvin 库。

---

## 五、Q2：设计准则（公式已核，全文已读）

- **Sullivan & Zhang APEC2014「Simplified Design Method for Litz Wire」**（免费 PDF ryz.ece.illinois.edu + elektrisola.com 镜像）：
  - 推荐丝数 n_e = k·δ²/(b·N_s)，k：单级绞合/真编织=1.0，两三级成缆=1.33，不绞=2.0；
  - 首级最多丝数 n₁,ₘₐₓ = 4δ²/d_s²；每级 4–5 束；
  - 好设计 d_s ≤ δ/2~δ/4；近邻-only 因子 F_R = 1+(πnN_s)²d_s⁶/(192δ⁴b²)；
  - 圆柱均匀场涡流损 P = πℓd⁴/(64ρ)(dB/dt)² → 损耗因子 ∝ n²d_s⁶。
- **Umetani et al. TIA 2021**（冈山大，免费 PDF 已读）：多级绞合全解析模型（Kelvin F/K 闭式 + 有效电阻率 ρ/m，m=捻距比），**式(51/56/57) 给出「最小捻距条件」——束级邻近损耗可忽略的可计算充分判据**（Q4 判据组件）。其参考文献表钉死了整个学派谱系（Dowell1966→Ferreira1992/94→Sullivan1999/2001→Nan&Sullivan→Rosskopf2014→Sullivan&Zhang2014）。
- **Zhang/White/Kassakian MIT APEC2014**（免费 PDF 已读）：P = F(f)·I²R_dc + G(f)·|H|² 分解 + PEEC 细丝法（比 FEM 快数量级，1 MHz 离散误差 <0.1%）；其捻距敏感性扫描就是 Q3 模板。
- **本题设计数（自算）**：A_cu=5 mm²（J=4 A/mm²）；d=0.05 mm→N=2546（n₁max=35，3 级）；d=0.071→N≈1265；d=0.1→N=637；**d=0.2→N=159（仿真甜点）**；d=0.25→N=102。
- 束径估计：铜填充率（含漆膜+成缆）实测范围 ~0.30–0.60（粗略），取 0.5 时 D_bundle≈3.57 mm。横截面丝心坐标直接用 cci_coords。
- 关键叙事（Sullivan 1999 原文）：**绞合只解决方位角覆盖，不解决径向交换**——内丝恒内、外丝恒外（Rosskopf 实测"损耗差异极端"）→ 这就是 Q3/Q4 的物理动机。

---

## 六、Q3：调优闭环设计

- 指标体系（题目表 2 + 我们的补充）：Rac/Rdc（主）、η_I=max|Iᵢ−Ī|/Ī（切片数据直出）、P=I²Rac、换位充分性 D（见 §七）。
- 变量：N、d_s（受等铜截面约束联动）、编织角 α（=arctan(2πR·n_t/L)）、换位周期 P、层级结构、方案类型。
- **方案梯度**（≥3 本质不同 + 控制变量）：
  1. 直束（不绞，下限基线）
  2. 单级简单绞合（方位角覆盖、无径向）
  3. 多级 S/Z 换向成缆（部分径向）
  4. 管状编织（双螺旋族；预期≈单级绞合——**反直觉论点，值得专门一节**）
  5. 真三维编织（track-and-column；径向完全遍历）
- 敏感性：换位周期 Roebel 经验 13–17 倍丝宽起、越短越好直到弯折/制造极限（Goldacker 2014 arXiv:1406.4244 + Wang 2022）；编织角通过 1/cos α 影响 Rdc、通过换位速率影响邻近损耗——存在最优 α 曲线。
- 记录：planning-with-files 的 findings.md 逐行追加（参数→Rac/Rdc→决策），论文的调优日志直接导出。

---

## 七、Q4：完美换位判据 + 拓扑设计 + 可制造性

### 7.1 判据（三层，全部可计算）

**离散层（组合）**：S 站采样得位置矩阵 M∈[N]^S×N。
完美 ⇔ ①每丝行是位置的均匀重排（Latin 方性质；Berger 循环赛=经典构造：每步整体轮转 N−1 槽）；
②φ-矩相等：对与泄漏场多极展开匹配的基 Φ={r, rcosθ, rsinθ, r², r²cos2θ, r²sin2θ,…}，Σ_s Φ(pos_k(s)) 对所有丝 k 相等；
③标量距离 D = max_{i,j}|ΣΦ_i−ΣΦ_j|/max Φ → 0（论文报告的「换位充分性」）。
（经典定义出处：Pyrhonen 教科书"perfect transposition ⇔ all subconductor currents remain the same"；所有标准 litz 损耗模型都建立在此假设上——Nottingham 博士论文原话"all previous models…under the condition of perfect transposition"。）

**拓扑层（辫群）**：交叉序列=braid word w=σ_{i1}^{±1}σ_{i2}^{±1}…∈B_N（关系 σᵢσᵢ₊₁σᵢ=σᵢ₊₁σᵢσᵢ₊₁，|i−j|≥2 交换）。每站截面状态=w 前缀在 B_N→S_n 满同态下的像（σᵢ↦(i,i+1)）。**周期缆完美 ⇔ w^p ∈ P_n（纯辫=各丝回原位）且 ⟨π(w)⟩ ≤ S_n 传递作用（每丝到达每槽）**。工具：SageMath/braidlab 法式判定；Word problem 可解。

**连续层（自洽）**：用 L2 切片得到的实际损耗密度 φ(r,θ) 替换基函数——判据变为 ∫φ(r_k(z),θ_k(z))dz 对所有 k 相等。Umetani 最小捻距条件作解析旁证。

### 7.2 拓扑设计方向

- 构造性方案：**N=环形槽位、周期 pN 站的「轮转+翻转」字**：每站整体轮转（cyclic permutation，等价于绞合成分）+ 每 m 站插入一次「翻面」纯辫（内外层交换，如 ∏σᵢσ_{N−i} 结构）——轮转管方位均匀，翻面管径向遍历；两者乘积的传递性可证。
- 每丝轨迹 (r(z),θ(z),z)：r(z) 在壳层间按停留-面积加权切换 + 平滑过渡（样条），θ 随轮转走——直接喂 L2 切片和 three.js。
- **原创贡献点（研究确认无人做过）**：纺织文献用网格/置换形式体系（Sontag/Ko 甚至只用"braid groups"指相邻纱线组），辫群数学从不碰角齿轮——**给出「braid word ↔ 角齿轮事件序列 ↔ 载纱器槽位置换」的显式同构并论证其可实现性**即可站住创新性。

### 7.3 可制造性（文献全部真实存在）

- Liu et al. IJAMT 2025（DOI 10.1007/s00170-025-15743-0，摘要在读）：编织机床=图（节点=槽位/轨道交点），载纱器路径=连通图遍历，床面分「典型区/切换区」——事件图的直接学术对应。
- Assi et al. IJAMT 2023（DOI 10.1007/s00170-023-11245-z）：图论载纱器碰撞检测=事件图并行边合法性。
- Emonts et al. Textiles 2021（开放获取，ITA RWTH）：现代 3D 旋转编织机每角齿轮独立控制（无机械开关）→ **任意事件图原则上可实现**；六角排列最密、lace braiding 允许双载纱器。
- Li et al. Polymer Composites 2023：无开关控制 3D 六角环旋转编织的矩阵状态机可视化（Euler 旋转矩阵=事件算子）——NumPy 可直接复刻。
- Bilisik TRJ 2013（~300 引）：4 步法/2 步法经典词汇表。Goldacker 2014：Roebel=1914 年就为趋肤效应发明的完全换位（ punched meander）；**英国 MTC 真的在 switch-track 编织机上试过 Litz 线**（the-mtc.org/insights/switchtrack-braiding-litz-wire）——工程可行性实证。
- CTC（连续换位导线）：每节距将一根矩形漆包线翻越整叠——工程化 Latin 方，工业模板。

---

## 八、K-TLitz / CTMA 调查结论

英文中文网页、GitHub API、Semantic Scholar 全部零命中（`q=K-TLitz` total_count=0）。判定：**内部项目代号，无公开足迹；以赛题文本本身为权威规格**。但「拓扑→事件→轨迹」框架并非杜撰——与 Liu 2025 载纱器路径图论研究 1:1 对应，引用之即可合法化我们的形式体系。论文中不要断言 CTMA 是什么。

---

## 九、批评缺口的处理记录

| 缺口 | 处理 |
|---|---|
| ① 2D 对绞合/编织失明、缺几何→Rac 通路 | §三保真度阶梯（L0–L3）填补 |
| ② Rdc 分母错误（未计丝长 1/cos α） | §〇决策 5 + §三 L2 修正项；IEC 60228 思路（绞线直流电阻按实际长度） |
| ③ 无实测基准/行业数值 | 填充率 0.30–0.60（含漆膜）；Roebel 节距 13–17 丝宽；Sullivan 论文内含实测表；MTC Litz 编织实例；其余在 Q2 时再查 Elektrisola/New England Wire 数据表 |

## 十、存疑引用清单（论文引用前必须核对）

1. ~~**Liu et al. 2025 IJAMT 卷期页码**~~ ✅ 已核验（2026-09-25）：138(11-12):5877-5890，见 paper/refs/LITSEARCH.md。
2. ~~`ryz.ece.illinois.edu` 两个 PDF 是个人站镜像~~ ✅ IEEE DOI 已定：10.1109/63.750181 / 10.1109/apec.2014.6803681。
3. ~~Rosskopf 论文走 ResearchGate 页~~ ✅ IEEE TPEL DOI：10.1109/tpel.2013.2293847（Ferreira 1994/1992 亦已核验，真实存在）。
4. Pyrhonen 引用 Archive.org 扫描本——终稿引正式出版版本（Wiley 2nd ed.）。
5. ARS README 的「Zhao 2026 ~147k 幻觉引用」统计未验证——不用。
6. TU Delft Bessel PDF、Ferreira IEEE 页、Dowell 未逐字读全文——只引不依赖细节。

## 十一、负结果（省时间的排除项）

- GitHub 无 K-TLitz/拓扑编织公开代码；无 FEM/电磁 Claude skill；无成熟可 fork 的 3D 辫查看器（只有玩具）。
- ElmerFEM mac 无官方包（编译成本高，弃）；FEMM 无 mac 原生版。
- KaTeX×VitePress 有已知构建 bug（用 MathJax3）；GitHub Pages 忘设 base 白屏。
- R3F(react-three-fiber) 虽好但部件多——48h 内用 vanilla three.js。
