# ==============================================================================
# 脚本名称: 02_fig7M_expression_diff.R
# 功能: Fig 7M: 8个风险基因在 High/Low MFRS 组间的表达差异 (箱线图+散点图)
#       布局: 2行4列，列优先排序
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

packages <- c("survival", "survminer", "tidyverse", "data.table", "ggpubr")
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
# 3. 计算 Risk Score 和分组
# ------------------------------------------------------------------------------
# 使用所有10个基因构建模型
all_genes <- c("SEMA3C", "RHOB", "PPP1R14A", "PCP4", "FHL2", 
               "DES", "CUX1", "CALD1", "CA12", "ACTN1")

multi_formula <- as.formula(paste("Surv(OS_time, OS_status) ~", paste(all_genes, collapse = " + ")))
multi_cox <- coxph(multi_formula, data = analysis_data)

analysis_data$risk_score <- predict(multi_cox, type = "risk", newdata = analysis_data)

# 分组 (High/Low)
n_survivors <- sum(analysis_data$OS_status == 0)
temp_sorted <- sort(analysis_data$risk_score)
cutoff_score <- temp_sorted[n_survivors]
analysis_data$Risk_group <- ifelse(analysis_data$risk_score > cutoff_score, "High MFRS", "Low MFRS")
# 修改: High MFRS (Red) 在左，Low MFRS (Blue) 在右
analysis_data$Risk_group <- factor(analysis_data$Risk_group, levels = c("High MFRS", "Low MFRS"))

# ------------------------------------------------------------------------------
# 4. 准备绘图数据
# ------------------------------------------------------------------------------
# 目标8个基因，按列优先顺序排列 (Col 1, Col 2, ...)
# 用户列表: ACTN1, FHL2, CA12, PPP1R14A, CALD1, RHOB, CUX1, SEMA3C
# 布局 2行4列，列优先意味着:
# Col 1: ACTN1, FHL2
# Col 2: CA12, PPP1R14A
# Col 3: CALD1, RHOB
# Col 4: CUX1, SEMA3C
target_genes <- c("ACTN1", "FHL2", "CA12", "PPP1R14A", "CALD1", "RHOB", "CUX1", "SEMA3C")

# 检查基因是否存在
missing <- target_genes[!target_genes %in% colnames(analysis_data)]
if(length(missing) > 0) stop(paste("缺失基因:", paste(missing, collapse=", ")))

# 转换为长格式
plot_data <- analysis_data %>%
  select(case_submitter_id, Risk_group, all_of(target_genes)) %>%
  pivot_longer(cols = all_of(target_genes), names_to = "Gene", values_to = "Expression")

# 设置基因因子的水平，以控制绘图顺序
# facet_wrap 默认按水平顺序填充 (Row 1, Row 2...)
# 如果我们要用 dir="v" (列优先)，我们需要确保因子顺序就是 target_genes 的顺序
plot_data$Gene <- factor(plot_data$Gene, levels = target_genes)

# ------------------------------------------------------------------------------
# 5. 绘图
# ------------------------------------------------------------------------------
message("正在生成箱线散点图...")

# 定义颜色
# Low MFRS: Blue Box, Dark Blue Points
# High MFRS: Red Box, Dark Red Points
box_colors <- c("Low MFRS" = "#4169E1", "High MFRS" = "#DC143C") # RoyalBlue, Crimson
point_colors <- c("Low MFRS" = "#00008B", "High MFRS" = "#8B0000") # DarkBlue, DarkRed

p_7m <- ggplot(plot_data, aes(x = Risk_group, y = Expression)) +
  # 散点图 (Jitter)
  geom_jitter(aes(color = Risk_group), width = 0.2, size = 1.5, alpha = 0.6) +
  # 箱线图 (透明度设为 0.5 以便看到后面的点，或者设为 NA 只显示框线)
  # 这里为了美观，通常箱线图在散点图上面，或者散点在箱线图上面
  # 用户描述: "箱线是红/蓝色，散点是深红/深蓝色"
  # 通常做法: 散点在下，箱线在上(无填充或半透明填充)
  geom_boxplot(aes(color = Risk_group, fill = Risk_group), width = 0.5, alpha = 0.3, outlier.shape = NA) +
  
  # 颜色设置
  scale_color_manual(values = point_colors) + # 散点和箱线边框颜色 (如果箱线边框用color映射)
  scale_fill_manual(values = box_colors) +    # 箱线填充颜色
  
  # 统计检验 (Wilcoxon test)
  stat_compare_means(method = "wilcox.test", label = "p.signif", 
                     symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),
                     label.x.npc = "center", vjust = 1.5, size = 5) +
  
  # 分面 (2行4列，列优先)
  facet_wrap(~Gene, ncol = 4, scales = "free_y", dir = "v") +
  
  # 主题设置
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(face = "bold", size = 12),
    axis.title.x = element_blank(), # 移除 x 轴标题 (Risk Group 显而易见)
    axis.text.x = element_text(face = "bold", color = "black"),
    legend.position = "none" # 移除图例，因为 x 轴已经分组
  ) +
  labs(y = "Gene Expression (log2 CPM + 1)")

# 保存
if (!dir.exists("results")) dir.create("results")

pdf("results/Fig7M_Expression_Diff.pdf", width = 12, height = 8)
print(p_7m)
dev.off()

png("results/Fig7M_Expression_Diff.png", width = 1200, height = 800)
print(p_7m)
dev.off()

message("Fig 7M 完成，已保存至 results/Fig7M_Expression_Diff.pdf")
