---
name: citation-verify
description: 引用与外部事实核验：bib 逐条 CrossRef/DOI 字段核对、专利走 Google Patents、URL 实开。触发词：引用核验、查引用、verify citations、参考文献检查、专利/DOI 核对。
---

# 引用核验（Citation Verification）

## 方法
1. **bib 逐条**：`curl https://api.crossref.org/works/<DOI>` 核对标题/作者/年/卷期页/DOI 五字段；
2. **专利号必须走 Google Patents**（CrossRef 不覆盖）——分别抓 /en 与 /zh 页防缓存误导；
3. **正文非 bib 引用**（arXiv 号、URL、机构页、案例）逐条打开；
4. **引号直引**：找到原文比对限定词（曾有引文丢 "HTS" 限定词扩大论断范围）；
5. **页码/标签错位**：引文标注的具体页码要真的对应那页内容（防"张冠李戴"）。

## 已知陷阱（实战教训）
- **专利号幻觉存活率极高**：曾有一条编造专利号（实为空调装配专利）穿过三轮文档审计，最终被 Google Patents 双语页抓出——**任何未核验的专利号视为编造**；
- arXiv 号抄错一位（1406.4244→1406.6244）仍像真的——用 arXiv API 核题名+作者配对；
- 摘要页数据 ≠ 全文数据（"仅摘要可得"的结论不得写成全文结论）；
- 卷期页的笔误形态：起始页>结束页、issue 写错（29(11)→29(10)）、卷期对调；
- Markdown 正文里裸 `@20`（本意"在 20 A 下"）会被 citeproc 当引用键渲染成 `(20?)` 进 PDF——全文扫 `@数字` 模式。

## 输出格式
```
[entry/位置] 声称 vs 核验事实 → verdict + 修正
```
