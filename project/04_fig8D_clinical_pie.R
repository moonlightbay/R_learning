# ==============================================================================
# 脚本名称: 04_fig8D_clinical_pie.R
# 功能: Fig 8D: 临床特征在 High/Low MFRS 组间的分布 (环形饼图)
#       布局: 2行8列 (行: Status, Age, Race, Stage, Grade, T, M, N; 列: High, Low)
#       颜色: 每列不同色系，深浅表示占比
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

packages <- c("survival", "tidyverse", "data.table", "patchwork", "grid", "cowplot")
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
risk_genes <- c("SEMA3C", "RHOB", "PPP1R14A", "PCP4", "FHL2", 
                "DES", "CUX1", "CALD1", "CA12", "ACTN1")

multi_formula <- as.formula(paste("Surv(OS_time, OS_status) ~", paste(risk_genes, collapse = " + ")))
multi_cox <- coxph(multi_formula, data = analysis_data)
analysis_data$risk_score <- predict(multi_cox, type = "risk", newdata = analysis_data)

cutoff_score <- median(analysis_data$risk_score)
analysis_data$Risk_group <- ifelse(analysis_data$risk_score > cutoff_score, "High", "Low")
analysis_data$Risk_group <- factor(analysis_data$Risk_group, levels = c("High", "Low"))

# ------------------------------------------------------------------------------
# 4. 准备临床变量
# ------------------------------------------------------------------------------
# 目标变量: Status, Age, Race, Stage, Grade, T, M, N
# 映射到数据列名
# Status -> OS_status (0=Alive, 1=Dead)
# Age -> age_at_index (需要分组: >60 vs <=60 或 Median)
# Race -> race
# Stage -> figo_stage
# Grade -> tumor_grade
# T -> ajcc_pathologic_t
# M -> ajcc_pathologic_m
# N -> ajcc_pathologic_n

# 4.1 变量清洗与重命名
plot_df <- analysis_data %>%
  transmute(
    Risk_group = Risk_group,
    Status = ifelse(OS_status == 1, "Dead", "Alive"),
    Age = ifelse(age_at_index > 60, ">60", "<=60"), # 常用分组
    Race = race,
    Stage = figo_stage,
    Grade = tumor_grade,
    T = ajcc_pathologic_t,
    M = ajcc_pathologic_m,
    N = ajcc_pathologic_n
  )

# 简单的缺失值/未知值处理
plot_df[plot_df == ""] <- NA
plot_df[plot_df == "not reported"] <- NA

# 变量列表 (顺序对应列颜色)
vars <- c("Status", "Age", "Race", "Stage", "Grade", "T", "M", "N")

# 颜色定义 (每列一个主色调)
# 淡绿、墨绿、橘红、深红、橙黄、青蓝、红色、蓝灰色
# 对应 vars 的顺序
base_colors <- c(
  "#90EE90", # Status: LightGreen
  "#006400", # Age: DarkGreen
  "#FF4500", # Race: OrangeRed
  "#8B0000", # Stage: DarkRed
  "#FFA500", # Grade: Orange
  "#00BFFF", # T: DeepSkyBlue (青蓝)
  "#FF0000", # M: Red
  "#708090"  # N: SlateGray (蓝灰)
)
names(base_colors) <- vars

# ------------------------------------------------------------------------------
# 5. 绘图函数 (生成单个环形图与图例)
# ------------------------------------------------------------------------------
# 辅助函数: 获取变量的所有类别并生成统一颜色映射
get_palette <- function(data, var_name, base_color) {
  cats <- sort(unique(na.omit(data[[var_name]])))
  n_cats <- length(cats)
  # 生成颜色: 从浅到深
  cols <- colorRampPalette(c("white", base_color))(n_cats + 2)[3:(n_cats + 2)]
  names(cols) <- cats
  return(cols)
}

create_donut_with_legend <- function(data, var_name, group_name, palette) {
  # 提取数据
  sub_data <- data %>%
    filter(Risk_group == group_name) %>%
    select(all_of(var_name)) %>%
    na.omit()
  
  # 统计频数
  counts <- table(sub_data[[1]])
  df <- as.data.frame(counts)
  colnames(df) <- c("Category", "Count")
  
  # 确保包含所有类别 (即使计数为0，为了颜色一致性)
  # 但饼图只画存在的
  
  # 计算比例
  df$Fraction <- df$Count / sum(df$Count)
  df$ymax <- cumsum(df$Fraction)
  df$ymin <- c(0, head(df$ymax, n = -1))
  
  # 绘图
  p <- ggplot(df, aes(ymax = ymax, ymin = ymin, xmax = 4, xmin = 3, fill = Category)) +
    geom_rect(color = "white") +
    coord_polar(theta = "y") +
    xlim(c(2, 4)) + # 空心圆 (环形)
    scale_fill_manual(values = palette) +
    theme_void() +
    theme(legend.position = "none") 
    
  return(p)
}

get_legend_plot <- function(data, var_name, palette) {
  # 创建一个虚拟图来提取图例
  cats <- names(palette)
  dummy_df <- data.frame(Category = factor(cats, levels = cats), Value = 1)
  
  p <- ggplot(dummy_df, aes(x=Category, y=Value, fill=Category)) +
    geom_bar(stat="identity") +
    scale_fill_manual(values = palette) +
    theme_void() +
    theme(
      legend.position = "bottom",
      legend.direction = "vertical",
      legend.title = element_blank(),
      legend.text = element_text(size = 10),
      legend.key.size = unit(0.5, "cm")
    )
    
  leg <- get_legend(p)
  return(as_ggplot(leg)) # 将 grob 转换为 ggplot 对象以便 patchwork 使用
}

# ------------------------------------------------------------------------------
# 6. 循环生成所有子图并拼图
# ------------------------------------------------------------------------------
message("正在生成子图与布局...")

# 存储每一列的组合图
col_plots <- list()

# 循环 8 个变量 (列)
for (i in 1:length(vars)) {
  v <- vars[i]
  base_col <- base_colors[v]
  
  # 1. 统一颜色映射
  palette <- get_palette(plot_df, v, base_col)
  
  # 2. 计算 P 值
  tbl <- table(plot_df[[v]], plot_df$Risk_group)
  if (nrow(tbl) > 1 && ncol(tbl) > 1) {
    test <- chisq.test(tbl)
    p_val <- test$p.value
  } else {
    p_val <- NA
  }
  p_txt <- ifelse(is.na(p_val), "NA", 
                  ifelse(p_val < 0.001, "p < 0.001", sprintf("p = %.3f", p_val)))
  
  # 3. 生成 High 组图
  p_high <- create_donut_with_legend(plot_df, v, "High", palette) +
    labs(title = v) + 
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))
  
  # 4. 生成 Low 组图 (带 P 值)
  p_low <- create_donut_with_legend(plot_df, v, "Low", palette) +
    labs(caption = p_txt) +
    theme(plot.caption = element_text(hjust = 0.5, size = 12, face = "italic", margin = margin(t = 5)))
  
  # 5. 生成图例图
  p_leg <- get_legend_plot(plot_df, v, palette)
  
  # 6. 组合该列 (High / Low / Legend)
  # 使用 patchwork 垂直堆叠
  col_combined <- p_high / p_low / p_leg + 
    plot_layout(heights = c(1, 1, 0.8)) # 调整高度比例，图例留够空间
  
  col_plots[[i]] <- col_combined
}

# ------------------------------------------------------------------------------
# 7. 添加左侧行标签
# ------------------------------------------------------------------------------
# 创建左侧标签列
lbl_high <- ggplot() + annotate("text", x = 0, y = 0, label = "High MFRS", angle = 90, size = 6, fontface = "bold") + theme_void()
lbl_low <- ggplot() + annotate("text", x = 0, y = 0, label = "Low MFRS", angle = 90, size = 6, fontface = "bold") + theme_void()
lbl_empty <- ggplot() + theme_void()

# 组合左侧列
left_col <- lbl_high / lbl_low / lbl_empty + plot_layout(heights = c(1, 1, 0.8))

# ------------------------------------------------------------------------------
# 8. 最终拼图与保存
# ------------------------------------------------------------------------------
# 将左侧列与 8 个数据列水平组合
final_plot <- left_col | col_plots[[1]] | col_plots[[2]] | col_plots[[3]] | 
              col_plots[[4]] | col_plots[[5]] | col_plots[[6]] | 
              col_plots[[7]] | col_plots[[8]]

# 调整整体宽度比例: 左侧标签占 0.5 份，其他各占 1 份
final_plot <- final_plot + plot_layout(widths = c(0.3, rep(1, 8)))

# 保存 (高分辨率)
if (!dir.exists("results")) dir.create("results")

# PDF
pdf("results/Fig8D_Clinical_Pie.pdf", width = 20, height = 8)
print(final_plot)
dev.off()

# PNG (高分辨率)
png("results/Fig8D_Clinical_Pie.png", width = 4000, height = 1600, res = 300)
print(final_plot)
dev.off()

message("Fig 8D 完成，已保存至 results/Fig8D_Clinical_Pie.pdf (高分辨率)")

