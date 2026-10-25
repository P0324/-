# 随访数据的获取方式

本仓库**不包含** Bott 队列的随访数据，需要你自己下载一次。这样做的原因有两个：该文件来自期刊补充材料（版权限制），且必须履行原始出处的可追溯性。

## 需要的文件

**Supplementary Table 18**，出自：

> Bott M, Brevet M, Taylor BS, et al.
> The nuclear deubiquitinase BAP1 is commonly inactivated by somatic mutations and 3p21.1 losses in malignant pleural mesothelioma.
> *Nature Genetics* 2011;43(7):668–672. doi:10.1038/ng.855

下载后重命名为 `41588_2011_BFng855_MOESM18_ESM.xlsx`，放在本目录（`data-raw/`）下。

## 怎么下

1. 打开文章页（https://doi.org/10.1038/ng.855）
2. 进入 **Supplementary information** 区域
3. 找到 Supplementary Table 18（文件名形如 `41588_2011_BFng855_MOESM18_ESM`）
4. 若学校/机构没有订阅，可尝试该文的 PubMed Central 版本或向通讯作者索取

## 文件应有的字段

`Sheet1`，53 行 × 12 列，关键三列：

| 列名 | 含义 |
|---|---|
| `Array sample #` | 芯片样本编号（与 `data/sample_map.csv` 的 `sample_id` 对应，53 个编号全部可匹配） |
| `Time from surgery to last f/u (months)` | 手术到最后一次随访的月数 |
| `Status at last f/u (0 = Alive)` | 0 = 存活，1 = 死亡；空白表示随访信息缺失 |

其余列（Sex、Age、Histology、Stage 等）在本次复现中未使用。

## 缺了会怎样

`R/01_prepare_clinical.R` 会提示缺失并中止。此时 `R/02_fig2c_km.R` 拿不到随访，`n` 与 `p` 都无从计算。

## 那 3 例缺失

补充表中有 3 例的随访时间或状态为空白，脚本会保留它们参与**分位数切点的计算**（这是作者的口径），只在建模时剔除，因此最终分析集是 **50 例**而非 53 例。

> 如果先剔除缺失再算三分位切点，log-rank p 会从 0.0213 变为 0.0163。两种口径结论一致，本仓库采用作者口径。
