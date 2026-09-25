# Q4 展示网站（web/）

三维拓扑编织 Litz 线的展示站：原理页（Q1/Q2 证据图）+ 编织查看器（three.js 逐丝轨迹、
z 向剖切、同步 2D 截面）+ 数据页（demo 轨迹指标与 φ-矩判据 D）。
vanilla three.js + Vite，无 UI 框架（决策见 RESEARCH.md §十一）。

## 本地开发

```bash
cd web
npm install        # 仅 three + vite
npm run dev        # 开发服务器（base='/'）
npm run build      # 产出 dist/（base='/litz-braided-wire/'）
npm run preview    # 本地预览构建产物
```

## 启用部署（一次性）

1. 仓库改名为 `litz-braided-wire`（或同步修改 `vite.config.js` 的 `base`，
   base 错误会白屏）。
2. 仓库 **Settings → Pages → Source 选 "GitHub Actions"**。
3. 之后 `main` 分支上 `web/**` 的每次 push 自动触发
   `.github/workflows/deploy.yml`：`npm ci → npm run build →
   configure-pages/upload-pages-artifact/deploy-pages`，发布到
   `https://<user>.github.io/litz-braided-wire/`。

## 轨迹数据契约

`public/demo/*.json` 由 `tools/gen_demo.py` 生成（`./.venv/bin/python web/tools/gen_demo.py`）：

```jsonc
{
  "N": 12, "units": "mm",
  "period_len": 24.0,        // 一个换位周期
  "z_total": 48.0,           // 展示长度（周期的整数倍，首尾截面相同）
  "n_samples": 481,          // z 均匀采样（含端点）
  "strand_r": 0.13,          // 丝半径（渲染管半径）
  "r_min": 0.6, "r_max": 1.45,
  "strands": [ { "r": [...], "theta": [...] } ],   // theta 为【展开】弧度（不取模）
  "metrics": { "basis": ["r", "r^2"], "mean_r": [...], "D": {"r": 0.0034, "r^2": 0.0061} }
}
```

- demo1 单绞向：r 恒定 + 截面刚体旋转（对照，D=1.0，无径向换位）。
- demo2 轮转+翻面：θ 均匀轮转 + tanh 展平的周期壳层交换（示意几何，D≈0.006）。

**Q4 正式数据接入**：`q4_perfect` 流水线产出角齿轮事件序列的真实 (r(z), θ(z)) 轨迹后，
按上述格式写同名 JSON（或新增文件并在 `src/main.js` 的 `DEMOS` 列表登记）即可，
查看器与数据页零改动。指标 `D` 同时在 Python 端（判据计算）与数据页（展示）使用同一口径
（RESEARCH.md §7.1 离散层：基 Φ={r, r²}，周期等权采样）。

## 目录

```
web/
├── index.html          # 三页骨架（原理 / 查看器 / 数据，hash 路由）
├── src/main.js         # 路由、查看器接线、数据页表格
├── src/braid.js        # StrandCurve(THREE.Curve) + BraidViewer + 2D 截面面板 + viridis
├── src/style.css       # 克制样式（白底、系统字体）
├── tools/gen_demo.py   # demo 轨迹生成（.venv 的 numpy）
├── public/demo/        # 轨迹 JSON（构建时原样拷入 dist）
└── public/figures/     # 证据图（从 paper/figures/ 拷贝；重新出图后需更新）
```

注意：`public/figures/` 是 `paper/figures/` 的快照拷贝（Vite 只打包 `web/` 内资源），
论文重新出图后请重新拷贝，保证网站与论文数字/图版一致（findings.md 溯源口径）。
