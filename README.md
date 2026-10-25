[README.md](https://github.com/user-attachments/files/33929144/README.md)
# 
The partial reproduction of 《Clonal gene signatures predict prognosis in mesothelioma and lung adenocarcinom》
# PRACME 论文复现（Fig.2c / Fig.4a）

复现 *npj Precision Oncology* (2024) 论文两种肿瘤中 **PRACME 克隆性签名**的预后验证：

> Pei L, et al. Clonal gene signatures predict prognosis in mesothelioma and lung adenocarcinoma.
> DOI: 10.1038/s41698-024-00531-y

签名共 **29 个基因**（LASSO 系数取自论文 Supplemental Table 1）。

---

## 材料清单

| 文件 / 文件夹 | 对应图 | 队列 | 平台 | 内容 |
|---|---|---|---|---|
| `Fig2c_PRACME_repro/` | **Fig.2c** | GSE29354（Bott 间皮瘤） | GPL96 | 完整仓库：代码 + 数据 + 图 |
| `Fig2c_KM_.lnk` | Fig.2c | 同上 | 同上 | 快捷方式，指向 `Fig2c_KM_无注释版.R`（单文件版本） |
| `fig4a.qmd` | **Fig.4a** | GSE31210（肺腺癌） | GPL570 | Quarto 文档，渲染即得 HTML 报告 |

> 注：Fig.4b（GSE68465 队列）不在这批材料中。

---

## 怎么跑

### Fig.2c

```bash
cd Fig2c_PRACME_repro
Rscript install_deps.R      # 首次
Rscript run_all.R           # 出图到 output/Fig2c_KM.png
```

需要先下载 Bott 等（*Nat Genet* 2011）的**补充表 18** 放进 `data-raw/`，
文件名保持 `41588_2011_BFng855_MOESM18_ESM.xlsx`（该文件有版权，未随仓库分发）。
详见 `Fig2c_PRACME_repro/data-raw/GET_DATA.md`。

### Fig.4a

在**项目根目录**（该文档假定为 `PRACME-repro/`）下：

```bash
quarto render fig4a.qmd     # 输出 fig4a.html
```

原始数据 `data/raw/GSE31210_series_matrix.txt.gz` 随文档提供。
246 例样本中 **226 例**为肺腺癌肿瘤组织，另 20 例为癌旁正常组织，文档中会先做筛选。

---

## 结果速览

| | 队列 | 可用基因 | n | log-rank p | Cox HR | 结局 |
|---|---|---|---|---|---|---|
| **Fig.2c** | GSE29354 间皮瘤 | **27 / 29** | 50 | 0.0213 | 1.452 | 总生存 |
| **Fig.4a** | GSE31210 肺腺癌 | **29 / 29** | 226 | 见该文档 | 见该文档 | 无复发生存 |

Fig.2c 只有 27 个基因，因为 **B3GALNT2** 与 **FTSJD1** 在 GPL96 上没有探针——
换成 GPL570 的注释包也没用（那 5 个探针在 Bott 的芯片上物理不存在）。这是平台限制，不是口径问题。
完整说明见 `Fig2c_PRACME_repro/README.md` 第五节。

---

## 数据来源

| 数据 | 来源 |
|---|---|
| GSE29354 表达谱 | GEO（GPL96 = HG-U133A） |
| GSE31210 表达谱 | GEO（GPL570 = HG-U133 Plus 2） |
| Bott 队列随访信息 | Bott et al., *Nat Genet* 2011（doi:10.1038/ng.855）补充表 18，需自行下载 |
| 29 基因系数 | 论文 Supplemental Table 1 |
| 作者原始代码 | https://github.com/YPeiLin/PRACME--NPJ-precision-oncology- （路径均为占位符，预处理产物未公开） |

作者公开的脚本是占位模板（`myinf1 = #path to gene expression`），两套预处理流水线均由本人重建，
因此数值可与论文对照方向，无法逐位比对。

---

## 复现结果摘要（Fig.2c）

```
n = 50    events = 48
分组: 18 / 17 / 15
log-rank p = 0.0213
Cox HR(每 SD) = 1.452 (1.084-1.945), p = 0.0124, C-index = 0.651
```

与原论文 Supplementary Table 2（HR 1.61，p 0.013）高度接近。
