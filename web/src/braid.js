// braid.js —— 轨迹 JSON → three.js 三维编织查看器 + 2D 截面面板
//
// 数据契约（web/tools/gen_demo.py 输出，q4_perfect 真实轨迹遵循同一格式）：
//   { N, period_len, z_total, n_samples, strand_r, r_min, r_max,
//     strands: [{ r: [M], theta: [M] }] }   —— z 均匀采样含端点，theta 已展开（不取模）
//
// 每丝一个 StrandCurve（THREE.Curve 子类）：getPoint(t) 对 r/θ 做 Catmull-Rom
// 标量插值、z 线性，返回 (r cosθ, r sinθ, z)。TubeGeometry 渲染。
import * as THREE from 'three'
import { OrbitControls } from 'three/addons/controls/OrbitControls.js'

// ---------- viridis 离散色（matplotlib viridis LUT 采样锚点，两端含尽） ----------
const VIRIDIS = [
  [0.267, 0.005, 0.329], [0.267, 0.224, 0.514], [0.194, 0.408, 0.557],
  [0.129, 0.569, 0.549], [0.208, 0.718, 0.475], [0.565, 0.843, 0.263],
  [0.993, 0.906, 0.144]
]

export function viridis(t) {
  t = Math.min(1, Math.max(0, t))
  const x = t * (VIRIDIS.length - 1)
  const i = Math.min(VIRIDIS.length - 2, Math.floor(x))
  const f = x - i
  const c = VIRIDIS[i].map((v, k) => v + f * (VIRIDIS[i + 1][k] - v))
  // 锚点为 sRGB 值：必须显式按 sRGB 写入，否则 getHexString() 会再过一次
  // linear→sRGB 转换把颜色调亮（three r152+ 默认开启色彩管理）
  return new THREE.Color().setRGB(c[0], c[1], c[2], THREE.SRGBColorSpace)
}

// Catmull-Rom 标量插值（端点钳位），f∈[0,1] 在 i 与 i+1 之间
function cr1(arr, i, f) {
  const p0 = arr[Math.max(0, i - 1)], p1 = arr[i]
  const p2 = arr[i + 1], p3 = arr[Math.min(arr.length - 1, i + 2)]
  const f2 = f * f, f3 = f2 * f
  return 0.5 * (2 * p1 + (-p0 + p2) * f +
    (2 * p0 - 5 * p1 + 4 * p2 - p3) * f2 +
    (-p0 + 3 * p1 - 3 * p2 + p3) * f3)
}

// ---------- 每丝一条参数曲线 ----------
export class StrandCurve extends THREE.Curve {
  constructor(rArr, thetaArr, zFirst, zLast) {
    super()
    this.r = rArr
    this.th = thetaArr
    this.z0 = zFirst
    this.z1 = zLast
    this.M = rArr.length
  }
  getPoint(t, target = new THREE.Vector3()) {
    const s = Math.min(this.M - 1.000001, Math.max(0, t * (this.M - 1)))
    const i = Math.floor(s), f = s - i
    const r = cr1(this.r, i, f), th = cr1(this.th, i, f)
    const z = this.z0 + (this.z1 - this.z0) * t
    return target.set(r * Math.cos(th), r * Math.sin(th), z)
  }
}

// ---------- 查看器 ----------
export class BraidViewer {
  /**
   * @param container  3D 视图容器（div）
   * @param canvas2d   2D 截面面板（canvas）
   * @param data       轨迹 JSON
   * @param opts       { onZ(zCut), onFrame(zCut) }
   */
  constructor(container, canvas2d, data, opts = {}) {
    this.container = container
    this.canvas2d = canvas2d
    this.data = data
    this.opts = opts
    this.playing = false
    this.playSpeed = 10 // mm/s
    this.zCut = data.z_total

    const rMax = data.r_max + data.strand_r * 2
    this.rMax = rMax

    // --- 场景 ---
    this.scene = new THREE.Scene()
    this.scene.background = new THREE.Color(0xffffff)
    this.camera = new THREE.PerspectiveCamera(
      40, 1, rMax / 50, (data.z_total + 4 * rMax) * 4)
    this.camera.position.set(2.6 * rMax, -2.0 * rMax, data.z_total * 0.62)

    this.renderer = new THREE.WebGLRenderer({ antialias: true })
    this.renderer.setPixelRatio(Math.min(2, window.devicePixelRatio))
    this.renderer.localClippingEnabled = true
    container.appendChild(this.renderer.domElement)

    this.scene.add(new THREE.HemisphereLight(0xffffff, 0x8a929c, 1.1))
    const dir = new THREE.DirectionalLight(0xffffff, 2.4)
    dir.position.set(1.5 * rMax, -2 * rMax, 2 * rMax)
    this.scene.add(dir)

    // --- 剖切平面：法向 -z，保留 z <= constant 的一段 ---
    this.plane = new THREE.Plane(new THREE.Vector3(0, 0, -1), this.zCut)

    // --- 每丝一条 Tube；viridis 离散配色（按丝序号） ---
    this.geometries = []
    this.materials = []
    const z0 = 0, z1 = data.z_total
    const tubular = Math.min(720, Math.max(240, 2 * data.n_samples))
    for (let k = 0; k < data.N; k++) {
      const s = data.strands[k]
      const curve = new StrandCurve(s.r, s.theta, z0, z1)
      const geo = new THREE.TubeGeometry(curve, tubular, data.strand_r, 10, false)
      const mat = new THREE.MeshStandardMaterial({
        color: viridis(k / Math.max(1, data.N - 1)),
        roughness: 0.35, metalness: 0.5,
        clippingPlanes: [this.plane], side: THREE.DoubleSide
      })
      const mesh = new THREE.Mesh(geo, mat)
      mesh.userData.strand = k
      this.geometries.push(geo)
      this.materials.push(mat)
      this.scene.add(mesh)
    }

    // --- 剖切位置指示环（随 zCut 移动的圆环） ---
    const ringPts = []
    for (let a = 0; a <= 64; a++) {
      const t = (a / 64) * Math.PI * 2
      ringPts.push(new THREE.Vector3(1.18 * rMax * Math.cos(t), 1.18 * rMax * Math.sin(t), 0))
    }
    this.ring = new THREE.Line(
      new THREE.BufferGeometry().setFromPoints(ringPts),
      new THREE.LineBasicMaterial({ color: 0xb0b6bd }))
    this.scene.add(this.ring)

    // --- 轨道控制 ---
    this.controls = new OrbitControls(this.camera, this.renderer.domElement)
    this.controls.enableDamping = true
    this.controls.target.set(0, 0, data.z_total / 2)
    this.controls.update()

    this._ro = new ResizeObserver(() => this._resize())
    this._ro.observe(container)
    this._resize()

    this.clock = new THREE.Clock()
    this._raf = 0
    const tick = () => {
      this._raf = requestAnimationFrame(tick)
      const dt = this.clock.getDelta()
      if (this.playing) {
        let z = this.zCut + this.playSpeed * dt
        if (z > data.z_total) z = 0
        this.setZ(z, { fromPlay: true })
      }
      this.controls.update()
      this.renderer.render(this.scene, this.camera)
    }
    tick()

    this.setZ(this.zCut)
  }

  _resize() {
    const w = this.container.clientWidth, h = this.container.clientHeight
    if (!w || !h) return
    this.renderer.setSize(w, h)
    this.camera.aspect = w / h
    this.camera.updateProjectionMatrix()
  }

  /** 移动剖切平面到 z（mm），同步 2D 面板与外部回调 */
  setZ(z, extra = {}) {
    this.zCut = Math.min(this.data.z_total, Math.max(0, z))
    this.plane.constant = this.zCut
    this.ring.position.z = this.zCut
    drawSectionPanel(this.canvas2d, this.data, this.zCut)
    if (this.opts.onZ) this.opts.onZ(this.zCut, extra)
  }

  play() { this.playing = true }
  pause() { this.playing = false }

  dispose() {
    cancelAnimationFrame(this._raf)
    this._ro.disconnect()
    this.controls.dispose()
    for (const g of this.geometries) g.dispose()
    for (const m of this.materials) m.dispose()
    this.renderer.dispose()
    this.renderer.domElement.remove()
  }
}

// ---------- 2D 截面面板：当前 z 站的丝心位置 + 按当前半径着色 ----------
export function drawSectionPanel(canvas, data, zCut) {
  const dpr = Math.min(2, window.devicePixelRatio || 1)
  const cssW = canvas.clientWidth || 300, cssH = canvas.clientHeight || 300
  if (canvas.width !== cssW * dpr) { canvas.width = cssW * dpr; canvas.height = cssH * dpr }
  const ctx = canvas.getContext('2d')
  ctx.setTransform(dpr, 0, 0, dpr, 0, 0)
  ctx.clearRect(0, 0, cssW, cssH)
  ctx.fillStyle = '#fff'
  ctx.fillRect(0, 0, cssW, cssH)

  const pad = 30
  const side = Math.min(cssW, cssH) - pad * 2 - 14
  const cx = (cssW - 24) / 2, cy = cssH / 2
  const rMaxView = data.r_max + data.strand_r * 2.5
  const scale = side / 2 / rMaxView

  // 参考圆：壳层 r_min / r_max / 束外缘
  ctx.strokeStyle = '#d8dce1'
  ctx.setLineDash([4, 4])
  for (const rr of [data.r_min, data.r_max]) {
    ctx.beginPath()
    ctx.arc(cx, cy, rr * scale, 0, Math.PI * 2)
    ctx.stroke()
  }
  ctx.setLineDash([])
  ctx.strokeStyle = '#eef0f2'
  ctx.beginPath()
  ctx.arc(cx, cy, (data.r_max + data.strand_r) * scale, 0, Math.PI * 2)
  ctx.stroke()

  // 当前 z 站的丝心（最近样本）
  const s = Math.min(data.n_samples - 1.000001,
    Math.max(0, (zCut / data.z_total) * (data.n_samples - 1)))
  const i = Math.floor(s), f = s - i
  for (let k = 0; k < data.N; k++) {
    const st = data.strands[k]
    const r = st.r[i] + f * (st.r[i + 1] - st.r[i])
    const th = st.theta[i] + f * (st.theta[i + 1] - st.theta[i])
    const x = cx + r * Math.cos(th) * scale
    const y = cy - r * Math.sin(th) * scale
    const c = viridis((r - data.r_min) / Math.max(1e-9, data.r_max - data.r_min))
    ctx.fillStyle = '#' + c.getHexString()
    ctx.beginPath()
    ctx.arc(x, y, Math.max(2.5, data.strand_r * scale), 0, Math.PI * 2)
    ctx.fill()
    if (data.N <= 24) {
      ctx.fillStyle = '#666'
      ctx.font = '9px system-ui, sans-serif'
      ctx.fillText(String(k), x + 5, y - 5)
    }
  }

  // 读数 + 图例
  ctx.fillStyle = '#444'
  ctx.font = '12px system-ui, sans-serif'
  ctx.fillText(`z = ${zCut.toFixed(1)} mm / ${data.z_total} mm`, 10, cssH - 10)
  ctx.fillText('截面（按当前半径着色）', 10, 16)
  const barW = 10, barH = side, bx = cssW - barW - 8, by = cy - barH / 2
  const grad = ctx.createLinearGradient(0, by + barH, 0, by)
  for (let a = 0; a <= 10; a++) {
    const c = viridis(a / 10)
    grad.addColorStop(a / 10, '#' + c.getHexString())
  }
  ctx.fillStyle = grad
  ctx.fillRect(bx, by, barW, barH)
  ctx.fillStyle = '#666'
  ctx.font = '10px system-ui, sans-serif'
  ctx.fillText(`${data.r_max.toFixed(2)}`, bx - 4, by + 8)
  ctx.fillText(`${data.r_min.toFixed(2)}`, bx - 4, by + barH)
  ctx.save()
  ctx.translate(bx - 4, by + barH / 2)
  ctx.rotate(-Math.PI / 2)
  ctx.fillText('r / mm', -20, 0)
  ctx.restore()
}
