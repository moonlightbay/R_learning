# ==============================================================================
# 脚本名称: 03_fig8B_km_mfrs.R
# 功能: Fig 8B: 高低 MFRS (Risk Score) 组的 Kaplan-Meier 生存曲线
#       展示 High MFRS 组生存率显著低于 Low MFRS 组
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

# ------------------------------------------------------------------------------
# 4. 分组 (High vs Low MFRS)
# ------------------------------------------------------------------------------
# 使用中位数作为截断值 (标准做法)
cutoff_score <- median(analysis_data$risk_score)
analysis_data$Risk_group <- ifelse(analysis_data$risk_score > cutoff_score, "High MFRS", "Low MFRS")

# 设置因子水平，High MFRS 为基准或特定顺序
# 通常 High Risk 对应较差生存
analysis_data$Risk_group <- factor(analysis_data$Risk_group, levels = c("Low MFRS", "High MFRS"))

message(paste("High MFRS 样本数:", sum(analysis_data$Risk_group == "High MFRS")))
message(paste("Low MFRS 样本数:", sum(analysis_data$Risk_group == "Low MFRS")))

# ------------------------------------------------------------------------------
# 5. 绘制 Kaplan-Meier 曲线
# ------------------------------------------------------------------------------
message("正在绘制 Fig 8B (KM Curve for MFRS)...")

# 拟合 KM 模型 (时间单位转换为年)
fit <- survfit(Surv(OS_time/365, OS_status) ~ Risk_group, data = analysis_data)

# 绘图
p_8b <- ggsurvplot(
  fit,
  data = analysis_data,
  pval = TRUE,             # 显示 P 值
  pval.method = TRUE,      # 显示 P 值计算方法
  conf.int = FALSE,        # 不显示置信区间 (保持整洁)
  risk.table = TRUE,       # 显示风险表
  risk.table.col = "strata", # 风险表颜色跟随曲线
  linetype = "solid",      # 线型
  surv.median.line = "hv", # 标注中位生存时间
  ggtheme = theme_classic(), # 经典主题
  palette = c("#377EB8", "#E41A1C"), # 配色: Low=蓝色, High=红色 (注意 levels 顺序: Low, High)
  legend.title = "MFRS Group",
  legend.labs = c("Low MFRS", "High MFRS"),
  xlab = "Time (Years)",
  ylab = "Overall Survival",
  title = "Kaplan-Meier Survival Analysis (MFRS)",
  break.time.by = 1
)

# 调整标题位置
p_8b$plot <- p_8b$plot + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16))

if (!dir.exists("results")) dir.create("results")

# 保存 Fig 8B
pdf("results/Fig8B_KM_MFRS.pdf", width = 6, height = 7, onefile = FALSE)
print(p_8b)
dev.off()

png("results/Fig8B_KM_MFRS.png", width = 600, height = 700)
print(p_8b)
dev.off()

message("Fig 8B 完成，已保存至 results/Fig8B_KM_MFRS.pdf")
