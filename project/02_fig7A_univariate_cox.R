# ==============================================================================
# 脚本名称: 02_fig7A_univariate_cox.R
# 功能: 批量单因素 Cox 回归分析并绘制森林图 (复现 Fig 7A)
# 时间: 2025-12-29
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. 环境准备与包安装
# ------------------------------------------------------------------------------
# 检查并安装必要的包
packages <- c("survival", "survminer", "tidyverse", "data.table")

for (pkg in packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

# 设置工作目录 (请根据实际情况修改)
setwd("d:/Works/R_learning/project")

# ------------------------------------------------------------------------------
# 2. 数据读取与预处理
# ------------------------------------------------------------------------------
message("正在读取数据...")

# 读取表达矩阵 (行=基因, 列=样本)
# check.names=FALSE 防止 R 自动修改列名中的连字符
expr_df <- fread("processed_data/expression_matrix.csv", data.table = FALSE)
rownames(expr_df) <- expr_df[, 1]
expr_df <- expr_df[, -1]

# 读取临床数据
clin_df <- read.csv("processed_data/clinical_cleaned.csv")

# --- 数据合并 ---
# 1. 转置表达矩阵: 变为 (行=样本, 列=基因)
expr_t <- t(expr_df)
expr_t <- as.data.frame(expr_t)
expr_t$case_submitter_id <- rownames(expr_t)

# 2. 合并临床信息
# 仅保留有生存数据的样本
merged_data <- merge(clin_df, expr_t, by = "case_submitter_id")

# 3. 过滤无效数据
# ------------------------------------------------------------------------------
# 关于生存数据的处理原则:
# 1. OS_time 为 NA: 必须剔除 (无法定位时间点)
# 2. OS_time = 0: 通常剔除 (可能是手术当天死亡或数据缺失，会干扰模型)
# 3. OS_status (0=Alive, 1=Dead): 两者都保留。0 代表删失(Censored)，即"至少活到了这一天"
# ------------------------------------------------------------------------------
analysis_data <- merged_data %>%
  filter(!is.na(OS_time) & OS_time > 0 & !is.na(OS_status))

message(paste("最终纳入分析的样本数:", nrow(analysis_data)))

# ------------------------------------------------------------------------------
# 3. 基因筛选策略 (核心修改)
# ------------------------------------------------------------------------------
message("开始进行 Cox 回归筛选...")

# 策略 A: 使用论文图片中识别的 13 个基因 (精准复现)
custom_genes <- c("SMTN", "SEMA3C", "RHOB", "PPP1R14A", "PDLIM7", "PCP4",
                  "FHL2", "DSTN", "DES", "CUX1", "CALD1", "CA12", "ACTN1")

# 获取所有基因列表 (排除临床列)
all_genes <- colnames(expr_t)[colnames(expr_t) != "case_submitter_id"]

# 验证基因是否存在
missing_genes <- custom_genes[!custom_genes %in% all_genes]
if (length(missing_genes) > 0) {
  warning(paste("以下基因在表达矩阵中未找到 (可能是拼写错误或被过滤):", paste(missing_genes, collapse = ", ")))
}

# 仅保留存在的基因
candidate_genes <- custom_genes[custom_genes %in% all_genes]
message(paste("成功匹配到", length(candidate_genes), "/", length(custom_genes), "个目标基因"))

if (length(candidate_genes) == 0) {
  stop("没有匹配到任何目标基因！请检查基因名拼写或数据源。")
}

# ------------------------------------------------------------------------------
# 4. 批量单因素 Cox 回归
# ------------------------------------------------------------------------------
# 定义 Cox 回归函数
run_cox <- function(gene) {
  # 构建公式
  formula <- as.formula(paste0("Surv(OS_time, OS_status) ~ `", gene, "`"))
  
  tryCatch({
    cox_fit <- coxph(formula, data = analysis_data)
    cox_summary <- summary(cox_fit)
    
    # 提取关键指标
    HR <- cox_summary$conf.int[1]
    HR.confint.lower <- cox_summary$conf.int[3]
    HR.confint.upper <- cox_summary$conf.int[4]
    p.value <- cox_summary$coefficients[5]
    
    return(c(gene = gene, HR = HR, lower = HR.confint.lower, upper = HR.confint.upper, p.value = p.value))
  }, error = function(e) {
    return(NULL)
  })
}

# 运行批量分析
cox_results_list <- lapply(candidate_genes, run_cox)

# 整理结果
cox_results <- do.call(rbind, cox_results_list)
cox_results <- as.data.frame(cox_results)

# 转换数据类型
cox_results$HR <- as.numeric(cox_results$HR)
cox_results$lower <- as.numeric(cox_results$lower)
cox_results$upper <- as.numeric(cox_results$upper)
cox_results$p.value <- as.numeric(cox_results$p.value)

# ------------------------------------------------------------------------------
# 5. 结果筛选与可视化 (Forest Plot)
# ------------------------------------------------------------------------------

# 保留所有计算结果 (用户要求: 不论显著性，展示所有指定基因)
sig_genes <- cox_results %>%
  arrange(match(gene, custom_genes)) # 按原始基因顺序

message(paste("纳入绘图的基因数量:", nrow(sig_genes)))

# 保存显著基因列表 (这对后续步骤很重要)
write.csv(sig_genes, "results/univariate_cox_sig_genes.csv", row.names = FALSE)

# 取前 15 个基因用于绘图 (模拟 Fig 7A)
# 如果是自定义列表，则全部展示
if (!is.null(custom_genes)) {
  plot_genes <- sig_genes
} else {
  plot_genes <- head(sig_genes, 15)
}

# 准备绘图数据
plot_data <- plot_genes %>%
  mutate(
    # 格式化 P 值和 HR 文本
    p_text = ifelse(p.value < 0.001, "< 0.001", sprintf("%.3f", p.value)),
    hr_text = paste0(sprintf("%.2f", HR), " (", sprintf("%.2f", lower), "-", sprintf("%.2f", upper), ")"),
    # 区分风险因子
    type = ifelse(HR > 1, "Risk (HR > 1)", "Protective (HR < 1)")
  )

# 核心修正：强制指定因子水平 (Levels) 以控制绘图顺序
# ggplot2 的 Y 轴默认从下往上 (Level 1 -> Level N)
# 我们希望列表第一个基因 (SMTN) 在最上面，最后一个 (ACTN1) 在最下面
# 所以 levels 应该设为 candidate_genes 的倒序
final_levels <- rev(candidate_genes[candidate_genes %in% plot_data$gene])
plot_data$gene <- factor(plot_data$gene, levels = final_levels)

# 绘制森林图 (优化样式)
# 修改 1: 移除全局 color 映射，只在 point 中映射 P值
p <- ggplot(plot_data, aes(x = HR, y = gene)) +
  
  # 修改 2: 线条颜色固定为绿色 (forestgreen)
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.3, size = 0.8, color = "forestgreen") +
  
  # 修改 3: 只有菱形点映射 P值颜色
  geom_point(aes(color = p.value), size = 3.5, shape = 18) + 
  
  geom_vline(xintercept = 1, linetype = "dashed", color = "black", alpha = 0.6) +
  
  # 修改 4: 移除 geom_text (不再显示 HR/P值 文本)
  
  # 修改 5: 自定义颜色条 (黑 -> 深蓝 -> 蓝 -> 蓝绿 -> 绿 -> 淡绿 -> 极淡绿)
  # P值越小(越显著) -> 黑色; P值越大 -> 淡绿
  scale_color_gradientn(colors = c("black", "darkblue", "blue", "#008080", "forestgreen", "lightgreen", "#F0FFF0"),
                        name = "P Value") +
  
  # 调整 X 轴范围 (不再需要为右侧文本预留大量空间)
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.1))) +
  
  labs(
    title = "Univariate Cox Regression Analysis",
    subtitle = paste0("Analysis of ", nrow(plot_data), " Target Genes"),
    x = "Hazard Ratio (95% CI)",
    y = "mRNA"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    plot.subtitle = element_text(hjust = 0.5, size = 12),
    axis.text.y = element_text(size = 11, face = "bold", color = "black"),
    axis.text.x = element_text(size = 10),
    # 图例位置: 右下角
    legend.position = c(0.95, 0.05),
    legend.justification = c(1, 0),
    legend.background = element_rect(fill = "white", color = "black", size = 0.2),
    panel.grid.major.y = element_line(color = "gray90", linetype = "dotted")
  )

# ------------------------------------------------------------------------------
# 6. 保存结果
# ------------------------------------------------------------------------------
# 创建输出目录
if (!dir.exists("results")) dir.create("results")

# 保存图片
ggsave("results/Fig7A_ForestPlot.pdf", plot = p, width = 6, height = 6)
ggsave("results/Fig7A_ForestPlot.png", plot = p, width = 6, height = 6)

# 保存显著基因列表 (供后续 Lasso 分析使用)
write.csv(sig_genes, "results/univariate_cox_sig_genes.csv", row.names = FALSE)

message("分析完成！结果已保存至 results/ 文件夹。")
