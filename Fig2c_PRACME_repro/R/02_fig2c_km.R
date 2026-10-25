## Fig.2c 主分析: 27 基因评分 -> 三分位 -> KM / log-rank / Cox -> 出图
## 输入: data-raw/rma_expr.rds, data/{PRACME_coefficients,gene_probe_map}.csv, data-raw/bott_clinical.csv
## 输出: output/Fig2c_KM.png / .pdf

suppressMessages({
  library(readr); library(dplyr); library(survival)
  library(survminer); library(ggplot2)
})

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
dir.create("output", showWarnings = FALSE)

rma  <- readRDS("data-raw/rma_expr.rds")
pm   <- read_csv("data/gene_probe_map.csv",       show_col_types = FALSE)
cf   <- read_csv("data/PRACME_coefficients.csv",  show_col_types = FALSE)
clin <- read_csv("data-raw/bott_clinical.csv",    show_col_types = FALSE)

## ---- 1. 取 27 个可用基因 ----------------------------------------------------
## 29 个签名基因中, B3GALNT2 与 FTSJD1 在 GPL96 上无可用探针;
## FTSJD1 若按旧别名匹配会误配到 212380_at (CMTR1/FTSJD2, 另一基因), 必须显式剔除。
kk <- pm %>%
  filter(kept == TRUE, gene %in% cf$gene) %>%
  filter(gene != "FTSJD1")
stopifnot(nrow(kk) == 27)
cat("可用基因:", nrow(kk), "/", nrow(cf), " 缺失:",
    paste(setdiff(cf$gene, kk$gene), collapse = ", "), "\n")

## ---- 2. 评分 score = Σ x_i · β_i ------------------------------------------
E <- rma[kk$probe, , drop = FALSE]
rownames(E) <- kk$gene
w <- cf$coefficient[match(rownames(E), cf$gene)]
score <- as.numeric(t(E) %*% w)
names(score) <- colnames(E)

## ---- 3. 三分位分组(作者口径: 先在全体有评分样本上定切点) -------------------
d <- clin
d$score   <- score[match(d$gsm, names(score))]
d$time_yr <- d$time / 12
q1 <- quantile(d$score, 1/3, na.rm = TRUE)
q2 <- quantile(d$score, 2/3, na.rm = TRUE)

d <- d %>%
  mutate(grp = factor(case_when(score <= q1 ~ "low",
                                score >= q2 ~ "high",
                                TRUE        ~ "intermediate"),
                      levels = c("low", "intermediate", "high"))) %>%
  filter(!is.na(score), !is.na(time_yr), !is.na(status))

cat("n =", nrow(d), "   events =", sum(d$status), "\n")
cat("分组:", paste(table(d$grp), collapse = " / "), "\n")

## ---- 4. 统计检验 ------------------------------------------------------------
fit <- survfit(Surv(time_yr, status) ~ grp, data = d)
sdf <- survdiff(Surv(time_yr, status) ~ grp, data = d)
p   <- 1 - pchisq(sdf$chisq, df = 2)
cat(sprintf("log-rank chi2 = %.3f (df=2), p = %.5f\n", sdf$chisq, p))

cx <- coxph(Surv(time_yr, status) ~ scale(score), data = d)
s  <- summary(cx)
hr <- s$coefficients[1, 2]                 # 已是 exp(coef), 取对数后再算区间
se <- s$coefficients[1, 3]
cat(sprintf("Cox HR(每 SD) = %.3f (95%% CI %.3f-%.3f), p = %.4f, C-index = %.3f\n",
            hr, exp(log(hr) - 1.96 * se), exp(log(hr) + 1.96 * se),
            s$coefficients[1, 5], s$concordance[1]))
print(tapply(d$time_yr, d$grp, median))

cat("\nlog-rank p =", signif(p, 3), "\n")

## ---- 5. 出图 ----------------------------------------------------------------
plab <- if (p < 1e-4) "p < 0.0001" else sprintf("p = %.2g", p)
PAL  <- c("#4D4D4D", "#7CAE00", "#F8766D")

pfig <- ggsurvplot(
  fit, data = d, size = 1, palette = PAL,
  legend.title = "Strata",
  legend.labs  = c("Low-value", "Intermediate-value", "High-value"),
  xlab = "Survival Time (Years)", ylab = "Overall Survival",
  break.time.by = 2.5, xlim = c(0, 10),
  pval = plab, pval.size = 5, pval.coord = c(0, 0.03),
  censor = TRUE, censor.shape = "+", censor.size = 4.5,
  ggtheme = theme_classic(base_size = 13) +
    theme(legend.position = "top",
          legend.title = element_text(size = 13),
          legend.text  = element_text(size = 13)))
pfig$plot <- pfig$plot + guides(colour = guide_legend(reverse = TRUE))

ggsave("output/Fig2c_KM.png", plot = pfig$plot, width = 7.2, height = 5.6, dpi = 300)
try(grDevices::cairo_pdf("output/Fig2c_KM.pdf", width = 7.2, height = 5.6), silent = TRUE)
print(pfig$plot); dev.off()
cat("图已写出: output/Fig2c_KM.png / .pdf\n")
