## 一键复现 Fig.2c
## 用法: 在仓库根目录下执行  Rscript run_all.R
##
## 注意: 部分环境会强制重设 R 的工作目录(本机即如此), 因此这里依据脚本自身位置
##       定位仓库根目录并切换过去, 而不是依赖 getwd()。
local({
  if (dir.exists("data") && dir.exists("R")) return()      # 已是仓库根目录
  a <- commandArgs(trailingOnly = FALSE)
  f <- grep("^--file=", a, value = TRUE)
  if (!length(f)) return()
  p <- tryCatch(normalizePath(sub("^--file=", "", f[1])), error = function(e) NA)
  if (is.na(p)) return()
  root <- if (basename(dirname(p)) == "R") dirname(dirname(p)) else dirname(p)
  if (dir.exists(file.path(root, "data"))) setwd(root)
})
cat("工作目录:", getwd(), "\n\n")

source("R/01_prepare_clinical.R")
source("R/02_fig2c_km.R")
