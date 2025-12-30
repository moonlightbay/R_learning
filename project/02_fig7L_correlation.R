# ==============================================================================
# 脚本名称: 02_fig7L_correlation.R
# 功能: Fig 7L: 10个预后基因、OS 和 MFRS (Risk Score) 的相关性分析
#       (左下: 散点图; 右上: 热图; 对角线: 密度图/变量名)
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

packages <- c("survival", "survminer", "tidyverse", "data.table", "GGally")
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
# 3. 计算 Risk Score (复用模型)
# ------------------------------------------------------------------------------
final_genes <- c("PCP4", "DES", "FHL2", "RHOB", "SEMA3C", 
                 "ACTN1", "CUX1", "CALD1", "CA12", "PPP1R14A")

# 构建 Cox 模型计算 Risk Score
multi_formula <- as.formula(paste("Surv(OS_time, OS_status) ~", paste(final_genes, collapse = " + ")))
multi_cox <- coxph(multi_formula, data = analysis_data)

analysis_data$risk_score <- predict(multi_cox, type = "risk", newdata = analysis_data)

# ------------------------------------------------------------------------------
# 4. 准备绘图数据
# ------------------------------------------------------------------------------
# 变量顺序: OS, 10个基因, Risk
# 注意: OS 通常指 OS_time (连续变量) 用于相关性分析
plot_data <- analysis_data %>%
  select(OS_time, all_of(final_genes), risk_score) %>%
  rename(OS = OS_time, Risk = risk_score)

message("绘图变量:")
print(colnames(plot_data))

# ------------------------------------------------------------------------------
# 5. 定义自定义绘图函数 (GGally)
# ------------------------------------------------------------------------------

# 5.1 右上三角: 热图 (Blue-White-Red) + 数值
upper_fn <- function(data, mapping, ...) {
  # 提取 x 和 y 数据
  x <- eval_data_col(data, mapping$x)
  y <- eval_data_col(data, mapping$y)
  
  # 计算相关系数
  corr <- cor(x, y, use = "complete.obs", method = "pearson")
  
  # 定义颜色映射 (Blue - White - Red)
  # 使用 colorRampPalette 生成渐变色
  col_pal <- colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(200)
  
  # 将相关系数 [-1, 1] 映射到颜色索引 [1, 200]
  # -1 -> 1, 0 -> 100, 1 -> 200
  idx <- round((corr + 1) / 2 * 199) + 1
  fill_color <- col_pal[idx]
  
  # 文本颜色: 深色背景用白色，浅色背景用黑色
  text_color <- ifelse(abs(corr) > 0.6, "white", "black")
  
  # 绘制矩形和文本
  ggplot() +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, fill = fill_color) +
    annotate("text", x = 0.5, y = 0.5, label = sprintf("%.2f", corr), color = text_color, size = 4, fontface = "bold") +
    theme_void() +
    theme(panel.border = element_rect(color = "white", fill = NA, size = 1)) # 白色边框分隔
}

# 5.2 左下三角: 散点图 (黑色点) + 拟合线 (黑色线)
lower_fn <- function(data, mapping, ...) {
  ggplot(data = data, mapping = mapping) +
    geom_point(color = "black", size = 0.5, alpha = 0.5) +
    geom_smooth(method = "lm", color = "black", se = FALSE, size = 0.6) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank()
    )
}

# 5.3 对角线: 密度图 (可选)
diag_fn <- function(data, mapping, ...) {
  ggplot(data = data, mapping = mapping) +
    geom_density(fill = "grey80", color = "black", alpha = 0.7) +
    theme_void() +
    theme(panel.border = element_rect(color = "black", fill = NA))
}

# ------------------------------------------------------------------------------
# 6. 绘制图形
# ------------------------------------------------------------------------------
message("正在生成相关性矩阵图 (这可能需要几秒钟)...")

p_7l <- ggpairs(
  plot_data,
  upper = list(continuous = upper_fn),
  lower = list(continuous = lower_fn),
  diag = list(continuous = "densityDiag"), # 使用默认密度图，或者自定义
  axisLabels = "none", # 隐藏坐标轴刻度，保持整洁
  title = "Correlation between Prognostic Genes, OS and MFRS"
) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    strip.background = element_rect(fill = "white", color = "black"),
    strip.text = element_text(face = "bold", size = 10, color = "black")
  )

# 保存
if (!dir.exists("results")) dir.create("results")

# PDF
pdf("results/Fig7L_Correlation_Matrix.pdf", width = 12, height = 12)
print(p_7l)
dev.off()

# PNG
png("results/Fig7L_Correlation_Matrix.png", width = 1200, height = 1200, res = 100)
print(p_7l)
dev.off()

message("Fig 7L 完成，已保存至 results/Fig7L_Correlation_Matrix.pdf")
