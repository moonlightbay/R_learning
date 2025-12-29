# ==============================================================================
# 脚本名称: 02_fig7F_7G_heatmap_roc.R
# 功能: 
#   1. Fig 7F: 绘制高低风险组 10 个基因表达差异的热图 (Heatmap)
#   2. Fig 7G: 绘制 1、3、5 年生存预测的 ROC 曲线 (Time-dependent ROC)
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. 环境准备
# ------------------------------------------------------------------------------
# 定义检测和安装包的函数
ensure_package <- function(pkg) {
  if (!require(pkg, character.only = TRUE)) {
    message(paste("正在安装", pkg, "..."))
    install.packages(pkg, repos = "https://cloud.r-project.org")
    if (!require(pkg, character.only = TRUE)) {
      stop(paste("无法安装或加载包:", pkg))
    }
  }
}

packages <- c("survival", "survminer", "tidyverse", "data.table", "pheatmap", "timeROC")
for (pkg in packages) {
  ensure_package(pkg)
}

setwd("d:/Works/R_learning/project")

# ------------------------------------------------------------------------------
# 2. 数据读取与预处理 (保持一致)
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
# 3. 模型构建与风险评分计算 (复用逻辑)
# ------------------------------------------------------------------------------
final_genes <- c("SEMA3C", "RHOB", "PPP1R14A", "PCP4", "FHL2", 
                 "DES", "CUX1", "CALD1", "CA12", "ACTN1")

multi_formula <- as.formula(paste("Surv(OS_time, OS_status) ~", paste(final_genes, collapse = " + ")))
multi_cox <- coxph(multi_formula, data = analysis_data)

# 计算风险评分 (type="risk" 对应 Hazard Ratio)
analysis_data$risk_score <- predict(multi_cox, type = "risk", newdata = analysis_data)

# --- 确定分组阈值 (保持与 Fig 7E 一致) ---
# 使用存活人数作为切分点 (根据用户对 Fig 7E 的要求)
n_survivors <- sum(analysis_data$OS_status == 0)
temp_sorted <- sort(analysis_data$risk_score)
cutoff_score <- temp_sorted[n_survivors]

analysis_data$risk_group <- ifelse(analysis_data$risk_score > cutoff_score, "High MFRS", "Low MFRS")
# -----------------------------------------

# ==============================================================================
# 4. Fig 7F: 差异表达热图 (Heatmap)
# ==============================================================================
message("正在绘制 Fig 7F (Heatmap)...")

# 1. 准备绘图数据
# 按风险评分排序样本
heatmap_data <- analysis_data %>%
  arrange(risk_score)

# 确保是纯 data.frame 并且有正确的行名
heatmap_data <- as.data.frame(heatmap_data)
# 使用 case_submitter_id 作为行名 (确保唯一性)
# 如果有重复ID，make.unique 会处理
rownames(heatmap_data) <- make.unique(as.character(heatmap_data$case_submitter_id))

# 提取 10 个基因的表达矩阵 (行=基因, 列=样本)
# 确保只提取数值列
plot_mat <- t(heatmap_data[, final_genes])

# 2. 准备注释条 (Annotation)
# 必须是 data.frame，不能是 tibble
annotation_col <- data.frame(
  Group = factor(heatmap_data$risk_group, levels = c("Low MFRS", "High MFRS"))
)
# 行名必须与 plot_mat 的列名完全一致
rownames(annotation_col) <- colnames(plot_mat)

# 3. 定义颜色
# 列表名必须与 annotation_col 的列名 ("Group") 完全一致
ann_colors <- list(
  Group = c("Low MFRS" = "blue", "High MFRS" = "#E64B35") # 使用玫红色替换紫色
)

# 热图颜色: 蓝色(低) -> 白色 -> 玫红色(高)
heatmap_cols <- colorRampPalette(c("blue", "white", "#E64B35"))(100)

# 4. 绘图并保存
if (!dir.exists("results")) dir.create("results")

tryCatch({
  pdf("results/Fig7F_Heatmap.pdf", width = 8, height = 6)
  pheatmap(plot_mat,
           scale = "row",           # 按行标准化
           cluster_rows = TRUE,     # 聚类基因
           cluster_cols = FALSE,    # 不聚类样本
           show_colnames = FALSE,   # 隐藏样本名
           annotation_col = annotation_col,
           annotation_colors = ann_colors,
           color = heatmap_cols,
           main = "Expression Heatmap of 10 Risk Genes",
           border_color = NA)
  dev.off()
  
  png("results/Fig7F_Heatmap.png", width = 800, height = 600)
  pheatmap(plot_mat,
           scale = "row",
           cluster_rows = TRUE,
           cluster_cols = FALSE,
           show_colnames = FALSE,
           annotation_col = annotation_col,
           annotation_colors = ann_colors,
           color = heatmap_cols,
           main = "Expression Heatmap of 10 Risk Genes",
           border_color = NA)
  dev.off()
  message("Fig 7F 绘制成功")
}, error = function(e) {
  message("Fig 7F 绘制失败: ", e$message)
  if(!is.null(dev.list())) dev.off()
})


# ==============================================================================
# 5. Fig 7G: ROC 曲线 (1, 3, 5 年)
# ==============================================================================
message("正在绘制 Fig 7G (ROC Curves)...")

# 1. 计算 Time-dependent ROC
# 注意: timeROC 需要时间单位一致。这里 OS_time 是天。
# 1年=365, 3年=1095, 5年=1825
times_of_interest <- c(365, 365*3, 365*5)

roc_res <- timeROC(
  T = analysis_data$OS_time,
  delta = analysis_data$OS_status,
  marker = analysis_data$risk_score,
  cause = 1,                # 1 = 死亡
  times = times_of_interest,
  iid = TRUE                # 计算置信区间 (可选，设为TRUE更严谨)
)

# 2. 提取绘图数据 (手动提取以使用 ggplot2 风格)
roc_plot_data <- data.frame()

for (t in times_of_interest) {
  # 提取每个时间点的 TPR 和 FPR
  # timeROC 返回的 TP 和 FP 是矩阵，列对应时间点
  t_str <- as.character(t)
  # 找到对应时间的索引
  idx <- which(roc_res$times == t)
  
  if (length(idx) > 0) {
    temp_df <- data.frame(
      FPR = roc_res$FP[, idx],
      TPR = roc_res$TP[, idx],
      Time = paste0(t/365, "-Year")
    )
    roc_plot_data <- rbind(roc_plot_data, temp_df)
  }
}

# 3. 提取 AUC 值用于图例
auc_1y <- round(roc_res$AUC[1], 2)
auc_3y <- round(roc_res$AUC[2], 2)
auc_5y <- round(roc_res$AUC[3], 2)

# 更新图例标签
roc_plot_data$Time <- factor(roc_plot_data$Time, 
                             levels = c("1-Year", "3-Year", "5-Year"))

# 4. 绘图
p_7g <- ggplot(roc_plot_data, aes(x = FPR, y = TPR, color = Time)) +
  geom_line(size = 1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +
  
  # 自定义颜色 (参考论文常见配色: 红、绿、蓝)
  scale_color_manual(
    values = c("1-Year" = "#E64B35", "3-Year" = "#4DBBD5", "5-Year" = "#00A087"),
    labels = c(paste0("1-Year (AUC = ", auc_1y, ")"),
               paste0("3-Year (AUC = ", auc_3y, ")"),
               paste0("5-Year (AUC = ", auc_5y, ")"))
  ) +
  
  labs(
    title = "ROC Curves for Survival Prediction",
    x = "False Positive Rate (1 - Specificity)",
    y = "True Positive Rate (Sensitivity)",
    color = "Time Points"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = c(0.7, 0.2), # 图例放在右下角
    legend.background = element_rect(fill = "white", color = "black", size = 0.2)
  )

# 保存 Fig 7G
ggsave("results/Fig7G_ROC.pdf", p_7g, width = 6, height = 6)
ggsave("results/Fig7G_ROC.png", p_7g, width = 6, height = 6)

message("分析完成！")
message("Fig 7F 保存至: results/Fig7F_Heatmap.pdf")
message("Fig 7G 保存至: results/Fig7G_ROC.pdf")
