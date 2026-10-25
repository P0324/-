## 整理 Bott 队列的随访数据
## 输入: data-raw/41588_2011_BFng855_MOESM18_ESM.xlsx  (Bott 等 补充表 18, 需自行下载, 见 GET_DATA.md)
## 输出: data-raw/bott_clinical.csv                    ( gsm, sample_id, time, status )

suppressMessages({ library(readxl); library(dplyr); library(readr) })

## 依据脚本自身位置定位仓库根目录(不依赖 getwd, 部分环境会强制重设工作目录)
local({
  if (dir.exists("data") && dir.exists("R")) return()
  a <- commandArgs(trailingOnly = FALSE)
  f <- grep("^--file=", a, value = TRUE)
  if (!length(f)) return()
  p <- tryCatch(normalizePath(sub("^--file=", "", f[1])), error = function(e) NA)
  if (is.na(p)) return()
  root <- if (basename(dirname(p)) == "R") dirname(dirname(p)) else dirname(p)
  if (dir.exists(file.path(root, "data"))) setwd(root)
})

cli_file <- "data-raw/bott_clinical.csv"
xlsx     <- "data-raw/41588_2011_BFng855_MOESM18_ESM.xlsx"

if (file.exists(cli_file)) {
  cat("随访数据已存在, 跳过解析:", cli_file, "\n")
} else {
  if (!file.exists(xlsx))
    stop("未找到 ", xlsx, "\n请按 data-raw/GET_DATA.md 的说明获取该文件后重试。")

  xl <- readxl::read_excel(xlsx, sheet = "Sheet1")
  cli <- suppressWarnings(xl %>%              # 随访缺失的格子为空, 转数值会产生 NA, 属预期
    transmute(
      sample_id = as.character(`Array sample #`),
      time      = as.numeric(`Time from surgery to last f/u (months)`),
      status    = as.numeric(`Status at last f/u (0 = Alive)`)     # 0 = 存活, 1 = 死亡
    ))

  smap <- readr::read_csv("data/sample_map.csv",
                          col_types = "cc", show_col_types = FALSE)
  out <- left_join(smap, cli, by = "sample_id")

  n_match <- sum(out$sample_id %in% cli$sample_id)   # 编号是否配上
  n_fu    <- sum(!is.na(out$time) & !is.na(out$status))
  cat("编号配对:", n_match, "/", nrow(out), "样本\n")
  cat("随访完整:", n_fu, "/", nrow(out),
      "  (其余样本随访缺失, 按作者口径参与切点计算, 建模时剔除)\n")
  if (n_match < nrow(out))
    warning("有样本的 Array sample # 未能配上, 请核对两表的编号格式")

  write_csv(out, cli_file)
  cat("已写出:", cli_file, "\n")
}
