## 安装本仓库所需的 R 包
## tidyverse 覆盖 readr / dplyr / ggplot2 等; survival 与 survminer 需单独安装
need <- c("tidyverse", "survival", "survminer", "readxl")
have <- rownames(installed.packages())
miss <- setdiff(need, have)
if (length(miss)) {
  cat("待安装:", paste(miss, collapse = ", "), "\n")
  install.packages(miss, repos = "https://cloud.r-project.org")
} else {
  cat("依赖已齐备:", paste(need, collapse = ", "), "\n")
}

## Bioconductor 包 hgu133a.db 仅在你想自己重新生成探针映射表时才需要,
## 常规复现走 data/gene_probe_map.csv, 不必安装。
if (!requireNamespace("hgu133a.db", quietly = TRUE)) {
  cat("\n如需重新生成基因映射表, 可安装:\n")
  cat('  BiocManager::install("hgu133a.db")\n')
}
