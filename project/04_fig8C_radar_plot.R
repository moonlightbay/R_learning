# ==============================================================================
# 脚本名称: 04_fig8C_radar_plot.R
# 功能: Fig 8C: 10个风险基因在 High/Low MFRS 组间的表达模式 (雷达图/蜘蛛网图)
#       展示 High MFRS (深红) vs Low MFRS (深蓝) 的表达差异
#       使用 fmsb 包绘制传统雷达图
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

packages <- c("survival", "tidyverse", "data.table", "fmsb")
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

# ------------------------------------------------------------------------------
# 3. 计算 Risk Score (MFRS)
# ------------------------------------------------------------------------------
# 10个预后基因
risk_genes <- c("SEMA3C", "RHOB", "PPP1R14A", "PCP4", "FHL2", 
                "DES", "CUX1", "CALD1", "CA12", "ACTN1")

# 检查基因是否存在
missing_genes <- risk_genes[!risk_genes %in% colnames(analysis_data)]
if (length(missing_genes) > 0) {
  stop(paste("缺失基因:", paste(missing_genes, collapse = ", ")))
}

message("正在构建多变量 Cox 模型计算 Risk Score...")
multi_formula <- as.formula(paste("Surv(OS_time, OS_status) ~", paste(risk_genes, collapse = " + ")))
multi_cox <- coxph(multi_formula, data = analysis_data)

# 预测 Risk Score
analysis_data$risk_score <- predict(multi_cox, type = "risk", newdata = analysis_data)

# 分组 (High vs Low MFRS)
cutoff_score <- median(analysis_data$risk_score)
analysis_data$Risk_group <- ifelse(analysis_data$risk_score > cutoff_score, "High MFRS", "Low MFRS")
# 注意: fmsb 绘图时，数据行的顺序决定了图例顺序。
# 我们希望 High MFRS (Red) 和 Low MFRS (Blue)

# ------------------------------------------------------------------------------
# 4. 准备 fmsb 雷达图数据
# ------------------------------------------------------------------------------
# 指定的顺时针顺序
# 注意: 如果 fmsb 默认逆时针绘制，为了达到视觉上的顺时针效果 (PCP4 -> PPP1R14A...)
# 我们需要将除第一个基因外的其他基因倒序排列
ordered_genes <- c("PCP4", "DES", "FHL2", "RHOB", "SEMA3C", 
                   "ACTN1", "CUX1", "CALD1", "CA12", "PPP1R14A")

# 计算每组的平均表达量
# 结果是一个 2行 x 10列 的数据框 (行是组，列是基因)
radar_data_mean <- analysis_data %>%
  group_by(Risk_group) %>%
  summarise(across(all_of(ordered_genes), mean)) %>%
  arrange(desc(Risk_group)) # 确保 High MFRS 在前 (如果 High > Low 字母顺序的话，High 在前)
  # High MFRS, Low MFRS

# 转换为纯数据框，并将组名作为行名
radar_df <- as.data.frame(radar_data_mean[, -1])
rownames(radar_df) <- radar_data_mean$Risk_group

# fmsb 要求数据结构:
# 第1行: Max 值 (每个变量的上限)
# 第2行: Min 值 (每个变量的下限)
# 第3行+: 实际数据

# 确定上下限 (为了美观，可以设为全局最大值+一点余量，和0)
max_val <- max(radar_df) * 1.2
min_val <- 0

# 构建 fmsb 数据框
data_for_plot <- rbind(rep(max_val, 10), rep(min_val, 10), radar_df)

# ------------------------------------------------------------------------------
# 5. 绘制雷达图 (fmsb)
# ------------------------------------------------------------------------------
message("正在绘制 Fig 8C (Radar Plot using fmsb)...")

# 颜色定义 (带透明度)
# High MFRS: Dark Red (#8B0000)
# Low MFRS: Dark Blue (#00008B)
colors_border <- c("#8B0000", "#00008B")
# 添加透明度用于填充 (例如 40% 透明度 -> 66)
colors_fill <- c(rgb(139/255, 0, 0, 0.4), rgb(0, 0, 139/255, 0.4))

# 保存 PDF
if (!dir.exists("results")) dir.create("results")
pdf("results/Fig8C_Radar_Plot.pdf", width = 8, height = 8)

# 绘图参数调整
radarchart(data_for_plot,
           axistype = 1,       # 轴的类型
           pcol = colors_border, # 线条颜色
           pfcol = colors_fill,  # 填充颜色
           plwd = 2,           # 线条宽度
           plty = 1,           # 线条类型 (实线)
           cglcol = "grey",    # 网格线颜色
           cglty = 1,          # 网格线类型 (实线)
           axislabcol = "grey",# 轴标签颜色
           caxislabels = seq(0, round(max_val, 1), length.out = 5), # 轴刻度标签
           cglwd = 0.8,        # 网格线宽度
           vlcex = 1.0,        # 变量标签字体大小
           title = "Expression Patterns of Prognostic Genes"
)

# 添加图例
legend(x = "bottom", legend = rownames(radar_df), horiz = TRUE,
       bty = "n", pch = 20, col = colors_border,
       text.col = "black", cex = 1.2, pt.cex = 2)

dev.off()

# 保存 PNG
png("results/Fig8C_Radar_Plot.png", width = 800, height = 800)
radarchart(data_for_plot,
           axistype = 1,
           pcol = colors_border,
           pfcol = colors_fill,
           plwd = 2,
           plty = 1,
           cglcol = "grey",
           cglty = 1,
           axislabcol = "grey",
           caxislabels = seq(0, round(max_val, 1), length.out = 5),
           cglwd = 0.8,
           vlcex = 1.0,
           title = "Expression Patterns of Prognostic Genes"
)
legend(x = "bottom", legend = rownames(radar_df), horiz = TRUE,
       bty = "n", pch = 20, col = colors_border,
       text.col = "black", cex = 1.2, pt.cex = 2)
dev.off()

message("Fig 8C 完成，已保存至 results/Fig8C_Radar_Plot.pdf")
