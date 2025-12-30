# ==============================================================================
# 脚本名称: 02_fig7H_pca.R
# 功能: Fig 7H: 绘制高低风险组的 PCA 散点图 (PC1 vs PC2)
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

packages <- c("survival", "tidyverse", "data.table")
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
# 3. 模型构建与风险评分计算 (用于分组)
# ------------------------------------------------------------------------------
final_genes <- c("SEMA3C", "RHOB", "PPP1R14A", "PCP4", "FHL2", 
                 "DES", "CUX1", "CALD1", "CA12", "ACTN1")

multi_formula <- as.formula(paste("Surv(OS_time, OS_status) ~", paste(final_genes, collapse = " + ")))
multi_cox <- coxph(multi_formula, data = analysis_data)

# 计算风险评分
analysis_data$risk_score <- predict(multi_cox, type = "risk", newdata = analysis_data)

# --- 确定分组阈值 ---
n_survivors <- sum(analysis_data$OS_status == 0)
temp_sorted <- sort(analysis_data$risk_score)
cutoff_score <- temp_sorted[n_survivors]

analysis_data$risk_group <- ifelse(analysis_data$risk_score > cutoff_score, "High MFRS", "Low MFRS")

# ==============================================================================
# 4. Fig 7H: PCA 散点图 (PC1 vs PC2)
# ==============================================================================
message("正在绘制 Fig 7H (PCA)...")

# --- 基于全基因组数据 ---
# 1. 识别所有基因列 (排除临床列、ID列和计算产生的列)
non_gene_cols <- c(colnames(clin_df), "risk_score", "risk_group", "case_submitter_id")
gene_cols <- setdiff(colnames(analysis_data), non_gene_cols)

# 2. 提取全基因组表达矩阵
pca_input <- analysis_data[, gene_cols]
message(paste("正在使用全基因组数据进行 PCA，基因数量:", length(gene_cols)))

# 运行 PCA (scale. = FALSE 以复刻论文的大坐标范围)
pca_res <- prcomp(pca_input, scale. = FALSE)

# 提取 PCA 坐标并合并分组信息
pca_data <- as.data.frame(pca_res$x)
pca_data$risk_group <- analysis_data$risk_group

# 绘图
p_7h <- ggplot(pca_data, aes(x = PC1, y = PC2, color = risk_group)) +
  geom_point(size = 2, alpha = 0.7) +
  scale_color_manual(values = c("High MFRS" = "#E64B35", "Low MFRS" = "#4DBBD5")) +
  stat_ellipse(level = 0.95, show.legend = FALSE, linetype = "dashed") +
  labs(
    title = "PCA of Risk Genes",
    x = paste0("PC1 (", round(summary(pca_res)$importance[2, 1] * 100, 1), "%)"),
    y = paste0("PC2 (", round(summary(pca_res)$importance[2, 2] * 100, 1), "%)"),
    color = "Group"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "top"
  )

if (!dir.exists("results")) dir.create("results")
ggsave("results/Fig7H_PCA.pdf", p_7h, width = 6, height = 6)
ggsave("results/Fig7H_PCA.png", p_7h, width = 6, height = 6)

message("Fig 7H 完成，已保存至 results/Fig7H_PCA.pdf")
