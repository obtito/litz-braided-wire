import { defineConfig } from 'vite'

// base 必须与 GitHub Pages 仓库名一致（<user>.github.io/litz-braided-wire/）
// 忘设 base 会白屏（见 RESEARCH.md §十一）。
export default defineConfig({
  base: '/litz-braided-wire/',
  build: {
    target: 'es2020',
    assetsInlineLimit: 0
  }
})
