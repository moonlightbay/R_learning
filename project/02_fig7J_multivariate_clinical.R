# ==============================================================================
# 脚本名称: 02_fig7J_multivariate_clinical.R
# 功能: Fig 7J: 临床因素与风险评分的多因素 Cox 回归分析 (森林图)
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

packages <- c("survival", "survminer", "tidyverse", "data.table")
for (pkg in packages) {
  ensure_package(pkg)
}

setwd("d:/Works/R_learning/project")

# ------------------------------------------------------------------------------
# 2. 数据读取与预处理
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

# 检查基因是否存在
missing_genes <- final_genes[!final_genes %in% colnames(analysis_data)]
if(length(missing_genes) > 0) stop(paste("缺失基因:", paste(missing_genes, collapse=", ")))

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
# 4. 临床变量处理与分组
# ------------------------------------------------------------------------------
message("正在处理临床变量...")

# 4.1 Age Group (High/Low) - 使用中位数分组
median_age <- median(analysis_data$age_at_index, na.rm = TRUE)
analysis_data$Age_group <- ifelse(analysis_data$age_at_index > median_age, "High Age", "Low Age")
message(paste("年龄分组阈值:", median_age))

# 4.2 其他变量重命名与因子化
# 确保变量存在，如果不存在则创建全NA列以防报错，但会警告
check_col <- function(df, col_name) {
  if (!col_name %in% colnames(df)) {
    warning(paste("警告: 数据中缺少列", col_name))
    return(rep(NA, nrow(df)))
  } else {
    return(df[[col_name]])
  }
}

# 提取并重命名变量以便绘图
cox_data <- analysis_data %>%
  transmute(
    OS_time = OS_time,
    OS_status = OS_status,
    Risk_group = factor(Risk_group, levels = c("Low MFRS", "High MFRS")), # Ref = Low
    Age_group = factor(Age_group, levels = c("Low Age", "High Age")),     # Ref = Low
    Race = as.factor(check_col(., "race")),
    Grade = as.factor(check_col(., "tumor_grade")),
    Stage = as.factor(check_col(., "figo_stage")),
    T_stage = as.factor(check_col(., "ajcc_pathologic_t")),
    M_stage = as.factor(check_col(., "ajcc_pathologic_m")),
    N_stage = as.factor(check_col(., "ajcc_pathologic_n"))
  )

# 显式设置参考水平 (Reference Levels)
# 确保 M0, N0, T1/T0 是参考组，否则 R 可能会默认用字母顺序 (如 M0 在 M1 前，但有时会有空白字符干扰)
if ("M0" %in% levels(cox_data$M_stage)) cox_data$M_stage <- relevel(cox_data$M_stage, ref = "M0")
if ("N0" %in% levels(cox_data$N_stage)) cox_data$N_stage <- relevel(cox_data$N_stage, ref = "N0")
# T_stage 参考组: 优先 T1, 其次 T0, 否则默认
t_levels <- levels(cox_data$T_stage)
if ("T1" %in% t_levels) {
  cox_data$T_stage <- relevel(cox_data$T_stage, ref = "T1")
} else if ("T0" %in% t_levels) {
  cox_data$T_stage <- relevel(cox_data$T_stage, ref = "T0")
}

# 简单的缺失值处理：剔除含有 NA 的样本 (多因素回归要求)
# 或者保留 NA 但 Cox 会自动剔除
cox_data_clean <- na.omit(cox_data)
message(paste("多因素回归纳入样本数 (去除临床缺失值后):", nrow(cox_data_clean)))

# ------------------------------------------------------------------------------
# 5. 多因素 Cox 回归 (Clinical + Risk Score)
# ------------------------------------------------------------------------------
message("正在进行多因素 Cox 回归...")

# 构建公式
# 注意：根据数据实际情况，某些变量可能只有一个水平或全是NA，需要排除
# 修正: 移除 Stage (FIGO Stage) 以避免与 T/N/M 分期产生多重共线性 (Multicollinearity)
# Stage 通常由 T, N, M 决定，同时放入模型会导致系数估计极其不稳定 (如 HR 接近 0 或 无穷大)
valid_cols <- c("Risk_group", "Age_group", "Race", "Grade", "T_stage", "M_stage", "N_stage")
message("注意: 已移除 'Stage' 变量以避免与 T/N/M 分期共线性导致模型崩溃。")

# 检查每个变量的水平数，必须 >= 2
keep_cols <- c()
for (col in valid_cols) {
  # 1. 预先检查水平数
  if (length(unique(cox_data_clean[[col]])) < 2) {
    warning(paste("变量", col, "水平数不足 2，已从模型中排除。"))
    next
  }
  
  # 2. 动态剔除稀疏类别的样本 (防止模型不收敛)
  # 打印各水平的样本量，以便核查 (如 NX 是否存在)
  counts <- table(cox_data_clean[[col]])
  message(paste("变量", col, "的水平分布:"))
  print(counts)
  
  # 如果某个类别的样本量极少 (<3)，会导致系数无限大
  sparse_levels <- names(counts)[counts < 3]
  
  if (length(sparse_levels) > 0) {
    warning(paste("变量", col, "存在稀疏类别 (<3样本):", paste(sparse_levels, collapse=", "), "。将剔除这些样本以保证模型收敛。"))
    # 剔除含有稀疏类别的样本
    cox_data_clean <- cox_data_clean[!cox_data_clean[[col]] %in% sparse_levels, ]
    # 重置因子水平
    cox_data_clean[[col]] <- droplevels(cox_data_clean[[col]])
  }
  
  # 3. 再次检查水平数 (剔除样本后可能只剩1个水平)
  if (length(unique(cox_data_clean[[col]])) < 2) {
    warning(paste("变量", col, "剔除稀疏样本后水平不足 2，已从模型中排除。"))
    next
  }
  
  keep_cols <- c(keep_cols, col)
}

if (length(keep_cols) == 0) stop("没有有效的变量可用于多因素回归！")

formula_str <- paste("Surv(OS_time, OS_status) ~", paste(keep_cols, collapse = " + "))
message(paste("回归公式:", formula_str))

final_model <- coxph(as.formula(formula_str), data = cox_data_clean)

# ------------------------------------------------------------------------------
# 6. 结果提取与处理 (替代 ggforest 以避免报错)
# ------------------------------------------------------------------------------
# 提取模型摘要
multi_summary <- summary(final_model)

# 构建结果数据框
cox_results <- data.frame(
  var_name = rownames(multi_summary$coefficients),
  HR = multi_summary$conf.int[, "exp(coef)"],
  lower = multi_summary$conf.int[, "lower .95"],
  upper = multi_summary$conf.int[, "upper .95"],
  p.value = multi_summary$coefficients[, "Pr(>|z|)"]
)

# --- 关键步骤: 处理无限值/极端值 ---
# 标记极端值 (HR 无限大, 或 CI 范围过大导致绘图报错)
# 增加对 NA 的处理 (防止 if 判断报错)
cox_results$is_extreme <- is.infinite(cox_results$HR) | 
                          is.infinite(cox_results$lower) | 
                          is.infinite(cox_results$upper) |
                          is.na(cox_results$HR) | 
                          is.na(cox_results$lower) | 
                          is.na(cox_results$upper) |
                          cox_results$upper > 50  # 设定一个合理的上限

# 过滤掉极端值行
if (sum(cox_results$is_extreme, na.rm = TRUE) > 0) {
  warning("以下变量因 HR 或 CI 过于极端/无限/NA 而被排除在绘图之外:")
  print(cox_results[which(cox_results$is_extreme), ])
}

plot_data <- cox_results[!cox_results$is_extreme & !is.na(cox_results$is_extreme), ]

# 优化变量名显示
# 将 "Risk_groupHigh MFRS" 转换为 "Risk_group: High MFRS" 等更易读的形式
plot_data$display_name <- plot_data$var_name
for (col in keep_cols) {
  plot_data$display_name <- gsub(col, paste0(col, ": "), plot_data$display_name)
}

# --- 新增：添加参考水平 (Reference Levels) ---
# 用户希望看到低风险组、低年龄组等参考组在图中显示 (HR=1)
ref_rows <- list()
for (col in keep_cols) {
  # 获取参考水平 (第一个 level)
  ref_level <- levels(cox_data_clean[[col]])[1]
  
  # 构建显示名称
  display_name <- paste0(col, ": ", ref_level, " (Ref)")
  
  # 创建数据行
  ref_rows[[length(ref_rows) + 1]] <- data.frame(
    var_name = paste0(col, ref_level), # 虚拟ID
    HR = 1,
    lower = NA, # 参考组没有置信区间
    upper = NA,
    p.value = 1, # 参考组 P值为 1 或 NA
    is_extreme = FALSE,
    display_name = display_name
  )
}

# 合并参考行
if (length(ref_rows) > 0) {
  ref_df <- do.call(rbind, ref_rows)
  # 确保列名一致
  common_cols <- intersect(colnames(plot_data), colnames(ref_df))
  plot_data <- rbind(plot_data[, common_cols], ref_df[, common_cols])
}

# ------------------------------------------------------------------------------
# 7. 绘图 (参考 Fig 7C 风格)
# ------------------------------------------------------------------------------
# 排序: 保持 valid_cols 定义的顺序
# 我们需要根据 valid_cols 的顺序来排列 display_name
# 首先提取 display_name 中的变量部分
plot_data$var_prefix <- sub(":.*", "", plot_data$display_name)

# 定义变量的顺序索引
var_order <- match(plot_data$var_prefix, valid_cols)

# 定义水平的顺序 (参考组在前，其他在后)
# 简单的做法是：如果包含 "(Ref)" 则排在同变量的其他组前面
# 但 ggplot 是从下往上画，所以我们希望 (Ref) 在最下面 (即因子顺序靠前)
# 或者我们希望 (Ref) 在最上面？通常是 Ref 在上，其他在下。
# 让我们构建一个排序键：变量索引 + (是否Ref ? 0 : 1)
plot_data$is_ref <- grepl("\\(Ref\\)", plot_data$display_name)
plot_data <- plot_data %>% arrange(var_order, !is_ref, display_name)

# 设置因子水平 (逆序，以便 ggplot 从上往下画)
plot_data$display_name <- factor(plot_data$display_name, levels = rev(plot_data$display_name))

# --- 准备右侧文本数据 ---
plot_data$hr_ci_text <- ifelse(plot_data$is_ref, "Reference",
                               paste0(sprintf("%.2f", plot_data$HR), " (", 
                                      sprintf("%.2f", plot_data$lower), "-", 
                                      sprintf("%.2f", plot_data$upper), ")"))

plot_data$p_val_text <- ifelse(plot_data$is_ref, "-",
                               ifelse(plot_data$p.value < 0.001, "< 0.001", sprintf("%.3f", plot_data$p.value)))

# 定义文本显示的 x 坐标 (根据 log 坐标轴调整)
# 假设最大 HR 上限在 5-10 左右，我们将文本放在右侧
# Log 坐标下: 10^1=10, 10^1.5≈31, 10^2=100
# 调整: 增大坐标值以避免与森林图重合
text_pos_hr <- 40
text_pos_p <- 250

p_7j <- ggplot(plot_data, aes(x = HR, y = display_name)) +
  
  # 误差棒 (参考组 lower/upper 为 NA，不会绘制)
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.3, size = 0.8, color = "forestgreen") +
  
  # 点 (映射 P 值颜色)
  geom_point(aes(color = p.value), size = 3.5, shape = 18) + 
  
  # 参考线
  geom_vline(xintercept = 1, linetype = "dashed", color = "black", alpha = 0.6) +
  
  # --- 添加右侧文本 ---
  # 1. HR (95% CI)
  geom_text(aes(x = text_pos_hr, label = hr_ci_text), hjust = 0, size = 3.5, color = "black") +
  # 2. P Value
  geom_text(aes(x = text_pos_p, label = p_val_text), hjust = 0, size = 3.5, color = "black") +
  
  # 添加列标题 (使用 annotate)
  annotate("text", x = text_pos_hr, y = length(levels(plot_data$display_name)) + 1, 
           label = "HR (95% CI)", hjust = 0, size = 4, fontface = "bold") +
  annotate("text", x = text_pos_p, y = length(levels(plot_data$display_name)) + 1, 
           label = "P Value", hjust = 0, size = 4, fontface = "bold") +

  # 颜色条
  scale_color_gradientn(colors = c("black", "darkblue", "blue", "#008080", "forestgreen", "lightgreen", "#F0FFF0"),
                        name = "P Value", na.value = "grey50") + 
  
  # 坐标轴调整 (扩大 x 轴范围以容纳文本)
  scale_x_continuous(trans = "log10", breaks = c(0.1, 0.5, 1, 2, 5, 10), limits = c(NA, 600)) + 
  
  labs(
    title = "Multivariate Cox Analysis (Clinical Factors)",
    subtitle = paste0("Filtered Model (n=", nrow(cox_data_clean), ")"),
    x = "Hazard Ratio (Log Scale)",
    y = "" # 移除 y 轴标题，因为变量名已经很清楚
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    plot.subtitle = element_text(hjust = 0.5, size = 12),
    axis.text.y = element_text(size = 10, face = "bold", color = "black"),
    axis.text.x = element_text(size = 10),
    legend.position = "top", # 图例放上面，给右侧文本留空间
    panel.grid.major.y = element_line(color = "gray90", linetype = "dotted"),
    plot.margin = margin(10, 10, 10, 10)
  ) +
  # 允许绘图区域外的文本显示 (防止表头被裁剪)
  coord_cartesian(clip = "off")

# 保存
if (!dir.exists("results")) dir.create("results")

pdf("results/Fig7J_Multivariate_Clinical_Forest.pdf", width = 8, height = 6, onefile = FALSE)
print(p_7j)
dev.off()

png("results/Fig7J_Multivariate_Clinical_Forest.png", width = 800, height = 600)
print(p_7j)
dev.off()

message("Fig 7J 完成，已保存至 results/Fig7J_Multivariate_Clinical_Forest.pdf")

