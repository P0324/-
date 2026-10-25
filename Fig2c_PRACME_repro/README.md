# Fig.2c 复现：PRACME 签名在 Bott 间皮瘤队列的预后验证

复现论文 **Fig.2c**（总生存 KM 曲线 + log-rank 检验）：

> Clonal gene signatures predict prognosis in mesothelioma and lung adenocarcinoma
> npj Precision Oncology (2024) · DOI: 10.1038/s41698-024-00531-y
> 作者公开代码：https://github.com/YPeiLin/PRACME--NPJ-precision-oncology-

原图定义：用 PRACME 评分把 Bott 队列（GSE29354）等量三分为低/中/高危组，比较三组总生存。

---

## 一、复现结果

| 指标 | 本复现 | 原论文 |
|---|---|---|
| 分析样本 | **50**（事件 48，96.0%） | 未在图内标注 |
| 分组 | 18 / 17 / 15 | 图形一致 |
| log-rank | χ² = 7.700，df = 2，**p = 0.0213** | 图内标注 |
| 中位生存 | 2.122 / 0.838 / 0.753 年 | 图形方向一致 |
| 单因素 Cox（每 1 SD） | HR = 1.452（1.084–1.945），p = 0.0124 | **HR = 1.61，p = 0.013**（Supplementary Table 2） |
| C-index | 0.651 | — |

Cox 的 HR 与 p 高度接近；曲线排序、分组方向与原图一致。

---

## 二、目录结构

```
Fig2c_PRACME_repro/
├── data/
│   ├── PRACME_coefficients.csv   # 29 基因 + LASSO 系数（论文 Supplemental Table 1）
│   ├── gene_probe_map.csv        # GPL96 探针→基因映射（43 探针 / 28 基因，含 kept 标记）
│   └── sample_map.csv            # GSM 编号 ↔ Array sample #（53 行，来自 GEO 元数据）
├── data-raw/
│   ├── rma_expr.rds              # 22,283 探针 × 53 样本 RMA 矩阵（由公开 CEL 生成）
│   ├── bott_clinical.csv         # 随访数据（运行 01 脚本后生成，不入库）
│   └── GET_DATA.md               # 随访数据的获取方式
├── R/
│   ├── 01_prepare_clinical.R     # 整理随访数据
│   └── 02_fig2c_km.R             # 评分 → 三分位 → KM → 出图
├── output/
│   ├── Fig2c_KM.png
│   └── Fig2c_KM.pdf
├── install_deps.R                # 依赖安装
├── run_all.R                     # 一键复现
├── .gitignore
└── LICENSE
```

---

## 三、快速开始

```bash
git clone <你的仓库地址>
cd Fig2c_PRACME_repro
```

**第一步**：下载 Bott 队列的随访数据（论文没有程序化接口，需手动一步）

把 Bott et al. 2016 (Nat Genet) 的**补充表 18** 放到 `data-raw/` 并重命名为
`41588_2011_BFng855_MOESM18_ESM.xlsx`。
获取方式与字段说明见 [`data-raw/GET_DATA.md`](data-raw/GET_DATA.md)。

**第二步**：安装依赖并运行

```r
Rscript install_deps.R
Rscript run_all.R
```

预期屏幕输出：

```
n = 50    events = 48
log-rank p = 0.0213
```

结果图写入 `output/Fig2c_KM.png` 与 `output/Fig2c_KM.pdf`。

> 如果你已经在本仓库同目录下跑过一次，`data-raw/bott_clinical.csv` 存在时不会重新解析 xlsx。

---

## 四、方法链路

```
GSE29354 原始 CEL (GEO, 公开)
  └─ RMA（背景校正 + 分位数归一化 + 汇总）→ 22,283 探针 × 53 样本，log2 尺度
       └─ 按 gene_probe_map.csv 取 27 个签名基因的探针，多探针取均值
            └─ score = Σ x_i · β_i            # β 为 LASSO 系数
                 └─ 按 score 的三分位切点分三组（作者口径：先在全部有评分样本上定切点，再剔除随访缺失）
                      └─ KM + log-rank + Cox
```

两处刻意保留的作者写法：

1. **不除以基因数**。作者代码为 `t(data) %*% mygen / length(mygen)`，而 `mygen` 是单列数据框，`length()` 返回列数 = 1，实际并未除以基因数。除数对所有样本相同，**不影响排序与分组**。
2. **三分位切点先在全体有评分样本上计算**，随访缺失者在建模时才剔除。若先剔除缺失再算切点，p 会从 0.0213 变为 0.0163；两者都远小于 0.05，结论不变。

---

## 五、为什么是 27 个基因而不是 29

签名共 **29 个基因**（`PRACME_coefficients.csv` 完整保留 29 行），但落到 GPL96（HG-U133A）上只有 **27 个**有可用探针：

| 基因 | LASSO 系数 | 情况 |
|---|---|---|
| **B3GALNT2** | 0.23 | Entrez 不在 GPL96 注释中，该平台无对应探针 |
| **FTSJD1** | −0.01 | Entrez 无探针；若按旧符号别名匹配，会误配到 `212380_at`——该探针 Entrez 为 23070，实为 **CMTR1（旧称 FTSJD2）**，是另一个基因 |

两基因 |系数| 合计 0.24，占全部 |系数| 的 4.02%——**B3GALNT2 的 0.23 在 29 个基因里排前段**，这个缺失是有代价的，但数据侧无法补救。

已实证这不是注释包的问题：换成 `hgu133plus2.db`（GPL570，54,675 探针）后，这两个基因的 5 个探针（1553691_at、1561714_a_at、1562391_at、226233_at、222811_at）**在 Bott 数据里一个都不存在**（0/5）；Bott 的芯片物理上只覆盖 GPL570 探针集的 40.7%。这是**硬件限制**，不是口径问题。

> 换注释字典不会变出探针来：字典再厚，卡片不在手上也读不出信号。

---

## 六、无法与原论文逐位对齐的原因

| 环节 | 原论文 | 本仓库 |
|---|---|---|
| 表达矩阵 | MAS5 + BrainArray 自定义 EntrezCDF（**未公开**） | 自做 RMA（公开路径可重跑） |
| 探针聚合 | CDF 定义的基因级探针组 | GPL96 官方注释 + 多探针取均值 |
| 可用基因 | 29（自定义 CDF 定义） | 27（GPL96 平台限制） |
| 归一化 | 输入已归一化，脚本内无归一化步骤 | RMA 自带分位数归一化 |

作者仓库中的 `Figure2.R` 全部为占位模板：

```r
myinf1 = #path to gene expression
myinf2 = #path to clinical infomation
```

预处理产物（`.pcl` 文件）未公开，因此数值无法逐位比对，只能比对结论方向与量级。

---

## 七、数据来源与许可

| 数据 | 来源 | 可访问性 |
|---|---|---|
| 表达（CEL → RMA 矩阵） | GEO [GSE29354](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE29354)，GPL96 = HG-U133A | 公开，随仓库提供 `.rds`（约 8.8 MB） |
| 29 基因系数 | 论文 Supplemental Table 1 | 公开，已随仓库提供 |
| 随访信息 | Bott et al. 2016 Nat Genet 补充表 18 | 需自行下载，**不随仓库分发** |
| 原作者代码 | 作者 GitHub 仓库 | 公开，路径占位 |

- 本仓库**代码**采用 MIT 许可。
- 论文的补充材料（系数表、补充表 18）版权归原出版方，本仓库仅包含复现所必需的最小引用；补充表 18 请从官方来源获取，不在本仓库内再分发。

## 八、引用

```bibtex
@article{pei2024clonal,
  title   = {Clonal gene signatures predict prognosis in mesothelioma and lung adenocarcinoma},
  journal = {npj Precision Oncology},
  year    = {2024},
  doi     = {10.1038/s41698-024-00531-y}
}
```
