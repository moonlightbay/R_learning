# ==============================================================================
# 脚本名称: 02_fig7I_km.R
# 功能: Fig 7I: 绘制高低 SDC1 表达组的 Kaplan-Meier 生存曲线
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

# 引入 survminer
packages <- c("xfun","survival", "survminer", "tidyverse", "data.table")
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

# ==============================================================================
# 3. Fig 7I: Kaplan-Meier 生存曲线 (SDC1 High vs Low)
# ==============================================================================
message("正在绘制 Fig 7I (KM Curve for SDC1)...")

# 检查 SDC1 是否存在
if (!"SDC1" %in% colnames(analysis_data)) {
  stop("错误: 数据中未找到 SDC1 基因表达数据！")
}

# 1. 根据 SDC1 表达量中位数分组
cutoff_sdc1 <- median(analysis_data$SDC1)
analysis_data$SDC1_group <- ifelse(analysis_data$SDC1 > cutoff_sdc1, "High SDC1", "Low SDC1")

# 2. 拟合 KM 模型 (时间单位转换为年)
fit <- survfit(Surv(OS_time/365, OS_status) ~ SDC1_group, data = analysis_data)

# 3. 使用 survminer 绘图
p_7i <- ggsurvplot(
  fit,
  data = analysis_data,
  pval = TRUE,             # 显示 P 值
  pval.method = TRUE,      # 显示 P 值计算方法
  conf.int = FALSE,        # 不显示置信区间 (根据原图风格)
  risk.table = TRUE,       # 显示风险表
  risk.table.col = "strata", # 风险表颜色跟随曲线
  linetype = "solid",      # 线型
  surv.median.line = "hv", # 标注中位生存时间
  ggtheme = theme_classic(), # 经典主题
  palette = c("#E64B35", "#4DBBD5"), # 配色: High=红色, Low=蓝色
  legend.title = "Group",
  legend.labs = c("High SDC1", "Low SDC1"),
  xlab = "Time (Years)",
  ylab = "Survival Probability",
  title = "Kaplan-Meier Survival Analysis (SDC1)",
  break.time.by = 1
)

# 调整标题位置
p_7i$plot <- p_7i$plot + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

if (!dir.exists("results")) dir.create("results")

# 保存 Fig 7I
pdf("results/Fig7I_KM_Curve_SDC1.pdf", width = 6, height = 7, onefile = FALSE)
print(p_7i)
dev.off()

png("results/Fig7I_KM_Curve_SDC1.png", width = 600, height = 700)
print(p_7i)
dev.off()

message("Fig 7I 完成，已保存至 results/Fig7I_KM_Curve_SDC1.pdf")
