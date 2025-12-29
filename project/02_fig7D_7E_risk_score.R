# ==============================================================================
# 脚本名称: 02_fig7D_7E_risk_score.R
# 功能: 
#   1. Fig 7D: 绘制 10 个风险基因的系数柱状图 (Coefficient Bar Plot)
#   2. Fig 7E: 绘制风险评分曲线 (Risk Score Curve) 和生存状态散点图 (Survival Status)
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. 环境准备
# ------------------------------------------------------------------------------
packages <- c("survival", "survminer", "tidyverse", "data.table", "ggpubr")
for (pkg in packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

setwd("d:/Works/R_learning/project")

# ------------------------------------------------------------------------------
# 2. 数据读取与预处理 (与 Fig 7C 保持一致)
# ------------------------------------------------------------------------------
message("正在读取数据...")
# 读取表达矩阵
expr_df <- fread("processed_data/expression_matrix.csv", data.table = FALSE)
rownames(expr_df) <- expr_df[, 1]
expr_df <- expr_df[, -1]

# 读取临床数据
clin_df <- read.csv("processed_data/clinical_cleaned.csv")

# 转置表达矩阵并合并临床数据
expr_t <- t(expr_df)
expr_t <- as.data.frame(expr_t)
expr_t$case_submitter_id <- rownames(expr_t)

merged_data <- merge(clin_df, expr_t, by = "case_submitter_id")

# 过滤无效生存数据
analysis_data <- merged_data %>%
  filter(!is.na(OS_time) & OS_time > 0 & !is.na(OS_status))

message(paste("分析样本数:", nrow(analysis_data)))

# ------------------------------------------------------------------------------
# 3. 多因素 Cox 模型构建
# ------------------------------------------------------------------------------
final_genes <- c("SEMA3C", "RHOB", "PPP1R14A", "PCP4", "FHL2", 
                 "DES", "CUX1", "CALD1", "CA12", "ACTN1")

# 构建公式
multi_formula <- as.formula(paste("Surv(OS_time, OS_status) ~", paste(final_genes, collapse = " + ")))
multi_cox <- coxph(multi_formula, data = analysis_data)

# ==============================================================================
# 4. Fig 7D: 基因系数柱状图 (Coefficient Bar Plot)
# ==============================================================================
message("正在绘制 Fig 7D...")

# 提取系数和P值
multi_summary <- summary(multi_cox)
coef_data <- data.frame(
  Gene = rownames(multi_summary$coefficients),
  Coef = multi_summary$coefficients[, "coef"],
  P_val = multi_summary$coefficients[, "Pr(>|z|)"]
)

# 去除反引号
coef_data$Gene <- gsub("`", "", coef_data$Gene)

# --- 用户特殊要求: P值颜色条范围 0-2.5 ---
# 推测论文使用的是 -log10(P_value)，因为原始P值不可能超过1
# 标记: 用户要求复刻论文的 0-2.5 范围
coef_data$logP <- -log10(coef_data$P_val)
# -----------------------------------------

# 排序: 按系数大小排序
coef_data$Gene <- factor(coef_data$Gene, levels = coef_data$Gene[order(coef_data$Coef)])

# 绘图
p_7d <- ggplot(coef_data, aes(x = Gene, y = Coef)) +
  geom_bar(stat = "identity", aes(fill = logP), width = 0.7) +
  
  # 配色方案参考 Fig 7C (P值颜色映射: 黑->蓝->绿)
  # 设置 limits = c(0, 2.5) 以复刻论文
  scale_fill_gradientn(colors = c("black", "darkblue", "blue", "#008080", "forestgreen", "lightgreen", "#F0FFF0"),
                       name = "-log10(P Value)",
                       limits = c(0, 2.5),
                       oob = scales::squish) + # 超出范围的值会被挤压到边界颜色
  
  labs(
    title = "Coefficients of the 10 Risk Genes",
    x = "Gene Symbol",
    y = "Coefficient (log Hazard Ratio)"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black"), # 竖直放置时X轴标签旋转
    axis.text.y = element_text(color = "black"),
    legend.position = "right"
  )

# 保存 Fig 7D
if (!dir.exists("results")) dir.create("results")
ggsave("results/Fig7D_Coefficients.pdf", p_7d, width = 6, height = 5)
ggsave("results/Fig7D_Coefficients.png", p_7d, width = 6, height = 5)


# ==============================================================================
# 5. Fig 7E: 风险评分曲线与生存状态图 (Risk Score & Survival Status)
# ==============================================================================
message("正在绘制 Fig 7E...")

# 计算风险评分
# 用户指出纵坐标为 0, 5, 10，这通常是 Hazard Ratio (exp(coef))，即 type="risk"
risk_score <- predict(multi_cox, type = "risk", newdata = analysis_data)
analysis_data$risk_score <- risk_score

# --- 用户特殊要求: 竖虚线对应存活患者数 ---
# 论文图中，中间竖虚线对应的横坐标值应该是存活的患者数
# 这意味着高低风险的分组界限可能不是中位数，而是根据存活/死亡人数比例确定的
n_survivors <- sum(analysis_data$OS_status == 0) # 0=Alive
n_total <- nrow(analysis_data)

# 为了让竖线和颜色分界对齐，我们需要找到对应排名的风险评分作为阈值
# 先临时排序找到阈值
temp_sorted <- sort(analysis_data$risk_score)
cutoff_score <- temp_sorted[n_survivors]

# 划分高低风险组 (使用新的阈值)
analysis_data$risk_group <- ifelse(analysis_data$risk_score > cutoff_score, "High MFRS", "Low MFRS")
# -----------------------------------------

# 数据排序 (按风险评分从小到大排序)
plot_data <- analysis_data %>%
  arrange(risk_score) %>%
  mutate(patient_rank = 1:n())

# 转换生存时间为年 (假设原数据为天，且用户觉得单位很大)
plot_data$OS_time_years <- plot_data$OS_time / 365

# --- 子图 1: 风险评分曲线 (Risk Score Curve) ---
p_risk <- ggplot(plot_data, aes(x = patient_rank, y = risk_score)) +
  geom_area(aes(fill = risk_group), alpha = 0.8) + 
  geom_line(color = "black", size = 0.2) +
  geom_point(aes(color = risk_group), size = 0.5) +
  
  scale_fill_manual(values = c("High MFRS" = "#E64B35", "Low MFRS" = "#4DBBD5")) +
  scale_color_manual(values = c("High MFRS" = "#E64B35", "Low MFRS" = "#4DBBD5")) +
  
  # 竖线位置改为存活人数
  geom_vline(xintercept = n_survivors, linetype = "dashed", color = "black") +
  
  labs(y = "Risk Score", x = NULL) +
  theme_classic() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),
    legend.position = "top",
    legend.title = element_blank()
  )

# --- 子图 2: 生存状态散点图 (Survival Status Scatter) ---
plot_data$status_label <- ifelse(plot_data$OS_status == 1, "Dead", "Alive")

p_status <- ggplot(plot_data, aes(x = patient_rank, y = OS_time_years)) +
  geom_point(aes(color = status_label), size = 1.5, alpha = 0.8) +
  
  scale_color_manual(values = c("Dead" = "red", "Alive" = "forestgreen")) +
  
  # 竖线位置改为存活人数
  geom_vline(xintercept = n_survivors, linetype = "dashed", color = "black") +
  
  labs(y = "Survival Time (Years)", x = "Patients (sorted by risk score)") +
  theme_classic() +
  theme(
    legend.position = "top",
    legend.title = element_blank()
  )

# --- 拼图 ---
p_7e <- ggarrange(p_risk, p_status, 
                  ncol = 1, nrow = 2, 
                  heights = c(1, 1), 
                  align = "v")

# 保存 Fig 7E
ggsave("results/Fig7E_RiskScore_Survival.pdf", p_7e, width = 7, height = 8)
ggsave("results/Fig7E_RiskScore_Survival.png", p_7e, width = 7, height = 8)

message("分析完成！")
message("Fig 7D 保存至: results/Fig7D_Coefficients.pdf")
message("Fig 7E 保存至: results/Fig7E_RiskScore_Survival.pdf")
