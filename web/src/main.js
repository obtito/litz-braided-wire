// main.js —— 页面路由 + 查看器接线 + 数据页
import './style.css'
import { BraidViewer, viridis } from './braid.js'

const BASE = import.meta.env.BASE_URL
const DEMOS = ['demo3-threering', 'demo4-untwisted', 'demo1-helix', 'demo2-transpose']

const $ = (sel) => document.querySelector(sel)
const pages = ['principle', 'viewer', 'data']

// ---------- 路由（hash，刷新保持页面） ----------
function route() {
  const hash = (location.hash || '#/principle').replace(/^#\//, '')
  const page = pages.includes(hash) ? hash : 'principle'
  for (const p of pages) {
    $(`#page-${p}`).classList.toggle('active', p === page)
    $(`nav a[href="#/${p}"]`).classList.toggle('active', p === page)
  }
  if (page === 'viewer') initViewer()
  if (page === 'data') initData()
  window.scrollTo(0, 0)
}
window.addEventListener('hashchange', route)

// ---------- 查看器页 ----------
let viewer = null
let viewerReady = false

async function initViewer() {
  if (viewerReady) return
  viewerReady = true

  const sel = $('#demo-select')
  const slider = $('#z-slider')
  const playBtn = $('#play-btn')
  const note = $('#demo-note')

  const load = async (name) => {
    const data = await fetch(`${BASE}demo/${name}.json`).then((r) => {
      if (!r.ok) throw new Error(`HTTP ${r.status}`)
      return r.json()
    })
    note.textContent = data.note || ''
    slider.max = data.z_total
    slider.step = data.z_total / 400
    slider.value = data.z_total
    if (viewer) viewer.dispose()
    viewer = new BraidViewer($('#braid3d'), $('#section2d'), data, {
      onZ: (z) => { slider.value = z }
    })
  }

  sel.addEventListener('change', () => load(sel.value))
  slider.addEventListener('input', () => {
    if (!viewer) return
    if (viewer.playing) { viewer.pause(); playBtn.textContent = '▶ 播放' }
    viewer.setZ(parseFloat(slider.value))
  })
  playBtn.addEventListener('click', () => {
    if (!viewer) return
    if (viewer.playing) { viewer.pause(); playBtn.textContent = '▶ 播放' }
    else { viewer.play(); playBtn.textContent = '⏸ 暂停' }
  })

  await load(sel.value)
}

// ---------- 数据页 ----------
let dataReady = false

async function initData() {
  if (dataReady) return
  dataReady = true
  const wrap = $('#data-demos')
  for (const name of DEMOS) {
    const data = await fetch(`${BASE}demo/${name}.json`).then((r) => r.json())
    wrap.appendChild(buildDemoCard(data))
  }
  wrap.appendChild(comparisonNote())
}

function el(tag, cls, html) {
  const e = document.createElement(tag)
  if (cls) e.className = cls
  if (html !== undefined) e.innerHTML = html
  return e
}

function buildDemoCard(d) {
  const card = el('article', 'card')
  card.appendChild(el('h3', null, `${d.title} <code>${d.name}.json</code>`))
  card.appendChild(el('p', 'muted', d.note || ''))

  const m = d.metrics || {}
  const pct = (x) => (x === undefined ? '—' : (100 * x).toFixed(2) + '%')
  const table = el('table')
  table.innerHTML = `
    <tbody>
      <tr><th>丝数 N</th><td>${d.N}</td><th>丝半径</th><td>${d.strand_r} mm</td></tr>
      <tr><th>换位周期</th><td>${d.period_len} mm</td><th>展示长度</th><td>${d.z_total} mm（${d.z_total / d.period_len} 周期）</td></tr>
      <tr><th>壳层 r 范围</th><td>[${d.r_min}, ${d.r_max}] mm</td><th>z 采样数</th><td>${d.n_samples}</td></tr>
      <tr><th>φ-矩失衡 D（r）</th><td>${pct(m.D && m.D['r'])}</td><th>φ-矩失衡 D（r²）</th><td>${pct(m.D && m.D['r^2'])}</td></tr>
    </tbody>`
  card.appendChild(table)

  if (m.mean_r) card.appendChild(sparkline(m.mean_r, d))
  return card
}

function sparkline(meanR, d) {
  const c = el('canvas', 'spark')
  const dpr = Math.min(2, window.devicePixelRatio || 1)
  const W = 560, H = 90
  c.width = W * dpr; c.height = H * dpr
  c.style.width = '100%'; c.style.maxWidth = W + 'px'
  const ctx = c.getContext('2d')
  ctx.setTransform(dpr, 0, 0, dpr, 0, 0)
  const lo = Math.min(...meanR), hi = Math.max(...meanR)
  const pad = 14
  ctx.font = '10px system-ui, sans-serif'
  meanR.forEach((r, k) => {
    const x = pad + (k / Math.max(1, meanR.length - 1)) * (W - 2 * pad)
    const t = hi > lo ? (r - lo) / (hi - lo) : 0.5
    const y = H - 18 - t * (H - 34)
    const col = viridis((r - d.r_min) / Math.max(1e-9, d.r_max - d.r_min))
    ctx.fillStyle = '#' + col.getHexString()
    ctx.beginPath(); ctx.arc(x, y, 3, 0, Math.PI * 2); ctx.fill()
  })
  ctx.fillStyle = '#888'
  ctx.fillText(`每丝周期平均半径 ⟨r⟩：min ${lo.toFixed(3)} / max ${hi.toFixed(3)} mm —— ${d.name === 'demo1-helix' ? '内外层终身固定' : '各丝几乎重合（换位充分）'}`, pad, H - 4)
  return c
}

function comparisonNote() {
  return el('p', 'muted small',
    'D 的定义见 RESEARCH.md §7.1（离散层判据）：基 Φ={r, r²}，' +
    'D = max<sub>i,j</sub>|ΣΦ<sub>i</sub>−ΣΦ<sub>j</sub>| / max<sub>k</sub>|ΣΦ<sub>k</sub>|，一个周期等权采样。' +
    'D→0 即「完美换位」。此处 demo 为示意几何；Q4 正式轨迹（角齿轮事件序列）产出后直接替换 demo JSON。')
}

// ---------- 启动 ----------
route()
