# ==============================================================================
# 脚本名称: 02_fig7C_multivariate_cox.R
# 功能: 对最终选定的 10 个基因进行多因素 Cox 回归分析并绘制森林图 (复现 Fig 7C)
# 说明: "Final genes" 通常指多因素模型，因此这里使用 Multivariate Cox Regression
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. 环境准备与包安装
# ------------------------------------------------------------------------------
packages <- c("survival", "survminer", "tidyverse", "data.table")
for (pkg in packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

setwd("d:/Works/R_learning/project")

# ------------------------------------------------------------------------------
# 2. 数据读取与预处理 (同 Fig 7A)
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

message(paste("最终纳入分析的样本数:", nrow(analysis_data)))

# ------------------------------------------------------------------------------
# 3. 定义最终的 10 个基因
# ------------------------------------------------------------------------------
final_genes <- c("SEMA3C", "RHOB", "PPP1R14A", "PCP4", "FHL2", 
                 "DES", "CUX1", "CALD1", "CA12", "ACTN1")

# 检查基因是否存在
missing <- final_genes[!final_genes %in% colnames(analysis_data)]
if(length(missing) > 0) stop(paste("缺失基因:", paste(missing, collapse=", ")))

# ------------------------------------------------------------------------------
# 4. 多因素 Cox 回归 (Multivariate Cox Regression)
# ------------------------------------------------------------------------------
message("正在进行多因素 Cox 回归...")

# 构建多因素公式: Surv(time, status) ~ Gene1 + Gene2 + ...
multi_formula <- as.formula(paste("Surv(OS_time, OS_status) ~", paste(final_genes, collapse = " + ")))

# 运行模型
multi_cox <- coxph(multi_formula, data = analysis_data)
multi_summary <- summary(multi_cox)

# 提取结果
cox_results <- data.frame(
  gene = rownames(multi_summary$coefficients),
  HR = multi_summary$conf.int[, "exp(coef)"],
  lower = multi_summary$conf.int[, "lower .95"],
  upper = multi_summary$conf.int[, "upper .95"],
  p.value = multi_summary$coefficients[, "Pr(>|z|)"]
)

# 去除基因名称中的反引号 (如果有)
cox_results$gene <- gsub("`", "", cox_results$gene)

# ------------------------------------------------------------------------------
# 5. 绘图 (沿用 Fig 7A 的优化风格)
# ------------------------------------------------------------------------------

# 排序：通常按 HR 大小排序，或者按用户指定的列表顺序
# 这里我们按用户提供的 final_genes 顺序逆序排列 (因为 ggplot 从下往上画)
plot_data <- cox_results
plot_data$gene <- factor(plot_data$gene, levels = rev(final_genes))

p <- ggplot(plot_data, aes(x = HR, y = gene)) +
  
  # 绿色误差棒
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.3, size = 0.8, color = "forestgreen") +
  
  # 菱形点映射 P值颜色
  geom_point(aes(color = p.value), size = 3.5, shape = 18) + 
  
  geom_vline(xintercept = 1, linetype = "dashed", color = "black", alpha = 0.6) +
  
  # 颜色条 (黑 -> 深蓝 -> 蓝 -> 蓝绿 -> 绿 -> 淡绿 -> 极淡绿)
  scale_color_gradientn(colors = c("black", "darkblue", "blue", "#008080", "forestgreen", "lightgreen", "#F0FFF0"),
                        name = "P Value") +
  
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.1))) +
  
  labs(
    title = "Multivariate Cox Regression Analysis",
    subtitle = paste0("Final ", nrow(plot_data), " Genes Model"),
    x = "Hazard Ratio (95% CI)",
    y = "Gene Symbol"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    plot.subtitle = element_text(hjust = 0.5, size = 12),
    axis.text.y = element_text(size = 11, face = "bold", color = "black"),
    axis.text.x = element_text(size = 10),
    legend.position = c(0.95, 0.05),
    legend.justification = c(1, 0),
    legend.background = element_rect(fill = "white", color = "black", size = 0.2),
    panel.grid.major.y = element_line(color = "gray90", linetype = "dotted")
  )

# ------------------------------------------------------------------------------
# 6. 保存结果
# ------------------------------------------------------------------------------
if (!dir.exists("results")) dir.create("results")

ggsave("results/Fig7C_Multivariate_ForestPlot.pdf", plot = p, width = 6, height = 6)
ggsave("results/Fig7C_Multivariate_ForestPlot.png", plot = p, width = 6, height = 6)
write.csv(cox_results, "results/multivariate_cox_results.csv", row.names = FALSE)

message("分析完成！结果已保存至 results/Fig7C_Multivariate_ForestPlot.png")
