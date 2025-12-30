# ==============================================================================
# 脚本名称: 02_fig7K_nomogram.R
# 功能: Fig 7K: 构建基于 MFRS 及临床因素的列线图 (Nomogram)
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. 环境准备
# ------------------------------------------------------------------------------
ensure_package <- function(pkg) {
  if (!require(pkg, character.only = TRUE)) {
    message(paste("正在安装", pkg, "..."))
    install.packages(pkg, repos = "https://cloud.r-project.org")
    if (!require(pkg, character.only = TRUE)) {
      stop(paste("无法安装或加载包:", pkg))
    }
  }
}

packages <- c("survival", "survminer", "tidyverse", "data.table", "rms")
for (pkg in packages) {
  ensure_package(pkg)
}

setwd("d:/Works/R_learning/project")

# ------------------------------------------------------------------------------
# 2. 数据读取与预处理 (复用 7J 的逻辑)
# ------------------------------------------------------------------------------
message("正在读取数据...")
expr_df <- fread("processed_data/expression_matrix.csv", data.table = FALSE)
rownames(expr_df) <- expr_df[, 1]
expr_df <- expr_df[, -1]

clin_df <- read.csv("processed_data/clinical_cleaned.csv")

expr_t <- t(expr_df)
expr_t <- as.data.frame(expr_t)
expr_t$case_submitter_id <- rownames(expr_t)

merged_data <- merge(clin_df, expr_t, by = "case_submitter_id")

analysis_data <- merged_data %>%
  filter(!is.na(OS_time) & OS_time > 0 & !is.na(OS_status))

message(paste("分析样本数:", nrow(analysis_data)))

# ------------------------------------------------------------------------------
# 3. 计算 Risk Score 和 Risk Group
# ------------------------------------------------------------------------------
final_genes <- c("SEMA3C", "RHOB", "PPP1R14A", "PCP4", "FHL2", 
                 "DES", "CUX1", "CALD1", "CA12", "ACTN1")

# 构建 Cox 模型计算 Risk Score
multi_formula <- as.formula(paste("Surv(OS_time, OS_status) ~", paste(final_genes, collapse = " + ")))
multi_cox <- coxph(multi_formula, data = analysis_data)

analysis_data$risk_score <- predict(multi_cox, type = "risk", newdata = analysis_data)

# 分组 (High/Low)
n_survivors <- sum(analysis_data$OS_status == 0)
temp_sorted <- sort(analysis_data$risk_score)
cutoff_score <- temp_sorted[n_survivors]
analysis_data$Risk_group <- ifelse(analysis_data$risk_score > cutoff_score, "High MFRS", "Low MFRS")

# ------------------------------------------------------------------------------
# 4. 临床变量处理 (保持与 7J 一致)
# ------------------------------------------------------------------------------
message("正在处理临床变量...")

# 4.1 Age Group (High/Low)
median_age <- median(analysis_data$age_at_index, na.rm = TRUE)
analysis_data$Age_group <- ifelse(analysis_data$age_at_index > median_age, "High Age", "Low Age")

check_col <- function(df, col_name) {
  if (!col_name %in% colnames(df)) {
    warning(paste("警告: 数据中缺少列", col_name))
    return(rep(NA, nrow(df)))
  } else {
    return(df[[col_name]])
  }
}

# 提取变量
nomo_data <- analysis_data %>%
  transmute(
    OS_time = OS_time,
    OS_status = OS_status,
    Risk_group = factor(Risk_group, levels = c("Low MFRS", "High MFRS")),
    Age_group = factor(Age_group, levels = c("Low Age", "High Age")),
    Race = as.factor(check_col(., "race")),
    Grade = as.factor(check_col(., "tumor_grade")),
    # Stage = as.factor(check_col(., "figo_stage")), # 7J 中已移除 Stage
    T_stage = as.factor(check_col(., "ajcc_pathologic_t")),
    M_stage = as.factor(check_col(., "ajcc_pathologic_m")),
    N_stage = as.factor(check_col(., "ajcc_pathologic_n"))
  )

# 显式设置参考水平 (与 7J 保持一致)
if ("M0" %in% levels(nomo_data$M_stage)) nomo_data$M_stage <- relevel(nomo_data$M_stage, ref = "M0")
if ("N0" %in% levels(nomo_data$N_stage)) nomo_data$N_stage <- relevel(nomo_data$N_stage, ref = "N0")
t_levels <- levels(nomo_data$T_stage)
if ("T1" %in% t_levels) {
  nomo_data$T_stage <- relevel(nomo_data$T_stage, ref = "T1")
} else if ("T0" %in% t_levels) {
  nomo_data$T_stage <- relevel(nomo_data$T_stage, ref = "T0")
}

# 去除缺失值
nomo_data_clean <- na.omit(nomo_data)

# ------------------------------------------------------------------------------
# 5. 再次进行稀疏类别过滤 (与 7J 保持一致)
# ------------------------------------------------------------------------------
valid_cols <- c("Risk_group", "Age_group", "Race", "Grade", "T_stage", "M_stage", "N_stage")
keep_cols <- c()

for (col in valid_cols) {
  if (length(unique(nomo_data_clean[[col]])) < 2) next
  
  counts <- table(nomo_data_clean[[col]])
  sparse_levels <- names(counts)[counts < 3]
  
  if (length(sparse_levels) > 0) {
    nomo_data_clean <- nomo_data_clean[!nomo_data_clean[[col]] %in% sparse_levels, ]
    nomo_data_clean[[col]] <- droplevels(nomo_data_clean[[col]])
  }
  
  if (length(unique(nomo_data_clean[[col]])) >= 2) {
    keep_cols <- c(keep_cols, col)
  }
}

message(paste("最终纳入 Nomogram 的变量:", paste(keep_cols, collapse = ", ")))

# ------------------------------------------------------------------------------
# 6. 构建 Nomogram 模型 (使用 rms 包)
# ------------------------------------------------------------------------------
# rms 包要求设置 datadist
dd <- datadist(nomo_data_clean)
options(datadist = "dd")

# 构建初始变量列表
final_cols <- keep_cols

# --- 新增: 预先检查共线性/奇异性 ---
# 使用 survival::coxph 快速检查是否有系数为 NA (表示奇异)
temp_formula <- as.formula(paste("Surv(OS_time, OS_status) ~", paste(final_cols, collapse = " + ")))
pre_fit <- coxph(temp_formula, data = nomo_data_clean)

if (any(is.na(coef(pre_fit)))) {
  na_coefs <- names(coef(pre_fit))[is.na(coef(pre_fit))]
  message(paste("检测到共线性/奇异性，以下系数为 NA:", paste(na_coefs, collapse=", ")))
  
  # 找出对应的原始变量名
  vars_to_remove <- c()
  for (col in final_cols) {
    # 检查该变量的任何 dummy 变量是否在 na_coefs 中
    if (any(grepl(col, na_coefs))) {
      vars_to_remove <- c(vars_to_remove, col)
    }
  }
  vars_to_remove <- unique(vars_to_remove)
  
  if (length(vars_to_remove) > 0) {
    message(paste("为解决共线性，将移除变量:", paste(vars_to_remove, collapse=", ")))
    final_cols <- setdiff(final_cols, vars_to_remove)
  }
}

if (length(final_cols) == 0) stop("没有有效变量可用于构建 Nomogram")

formula_str <- paste("Surv(OS_time, OS_status) ~", paste(final_cols, collapse = " + "))
message(paste("最终 Nomogram 公式:", formula_str))

# 使用 cph (Cox Proportional Hazards Model) 拟合
# x=TRUE, y=TRUE, surv=TRUE 是 nomogram 必须的
fit <- cph(as.formula(formula_str), data = nomo_data_clean, x = TRUE, y = TRUE, surv = TRUE)

# ------------------------------------------------------------------------------
# 7. 绘制 Nomogram
# ------------------------------------------------------------------------------
# 计算生存概率预测 (例如 1年, 3年, 5年)
# 这里的单位是天，所以 365, 1095, 1825
surv_times <- c(365, 1095, 1825)

# 获取生存函数生成器 (rms 包功能)
# Survival(fit) 返回一个函数，该函数接受 (时间, 线性预测值) 并返回生存概率
surv_func <- Survival(fit)

# 构建 nomogram 对象
nom <- nomogram(fit, 
                fun = list(function(x) surv_func(365, x),
                           function(x) surv_func(1095, x),
                           function(x) surv_func(1825, x)),
                funlabel = c("1-Year Survival", "3-Year Survival", "5-Year Survival"),
                lp = FALSE, # 是否显示线性预测值 (Linear Predictor)
                fun.at = c(0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1) # 概率轴刻度
)

# 保存图片
if (!dir.exists("results")) dir.create("results")

pdf("results/Fig7K_Nomogram.pdf", width = 10, height = 8)
plot(nom, xfrac = .25) # xfrac 调整标签位置
dev.off()

png("results/Fig7K_Nomogram.png", width = 1000, height = 800)
plot(nom, xfrac = .25)
dev.off()

message("Fig 7K 完成，已保存至 results/Fig7K_Nomogram.pdf")

# ------------------------------------------------------------------------------
# 8. 替代方案：rms 原生美化版 Nomogram（稳定无兼容问题）
# ------------------------------------------------------------------------------
message("生成 rms 原生美化版 Nomogram...")

# 重新构建 nomogram 对象（优化刻度和标签）
nom_optimized <- nomogram(
  fit, 
  fun = list(
    function(x) surv_func(365, x),
    function(x) surv_func(1095, x),
    function(x) surv_func(1825, x)
  ),
  funlabel = c("1-Year Survival", "3-Year Survival", "5-Year Survival"),
  lp = FALSE,
  fun.at = seq(0.1, 0.9, by = 0.1), # 简化概率刻度
  maxscale = 100, # 总分刻度统一为 0-100（临床列线图标准）
  abbrev = TRUE, # 简化变量名
  vnames = "labels" # 用变量标签而非原始名
)

# 保存美化版列线图（核心：调整绘图参数，模拟“矩形块”风格）
pdf("results/Fig7K_Nomogram_Optimized.pdf", width = 12, height = 9)
plot(
  nom_optimized,
  xfrac = 0.2, # 标签位置优化
  col.grid = gray(c(0.8, 0.95)), # 网格线颜色（模拟矩形块）
  col.lines = "black", # 线条颜色
  lwd = 1.2, # 线条粗细
  cex.axis = 0.9, # 轴刻度字体
  cex.var = 1.1, # 变量标签字体 (rms 使用 cex.var)
  mgp = c(2, 0.5, 0), # 轴标签位置
  main = "Nomogram for Overall Survival Prediction" # 标题
)
dev.off()

png("results/Fig7K_Nomogram_Optimized.png", width = 1200, height = 900, res = 150)
plot(
  nom_optimized, 
  xfrac = 0.2, 
  col.grid = gray(c(0.8, 0.95)), 
  col.lines = "black",
  lwd = 1.2,
  cex.axis = 0.9,
  cex.var = 1.1
)
dev.off()

message("稳定版列线图已保存至 results/Fig7K_Nomogram_Optimized.pdf/.png")
