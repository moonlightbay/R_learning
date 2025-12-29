# ==============================================================================
# 脚本名称: 02_fig7B_lasso.R
# 功能: 使用 Lasso 回归筛选基因并绘制 Lambda 图 (复现 Fig 7B)
# 输入: 13 个候选基因 (来自 Fig 7A)
# 输出: Lasso 系数路径图, Lasso CV 误差图, 最终筛选的基因
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. 环境准备
# ------------------------------------------------------------------------------
packages <- c("glmnet", "survival", "survminer", "tidyverse", "data.table", "ggpubr")
for (pkg in packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
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

# 过滤无效生存数据
analysis_data <- merged_data %>%
  filter(!is.na(OS_time) & OS_time > 0 & !is.na(OS_status))

# ------------------------------------------------------------------------------
# 3. 准备 Lasso 输入矩阵
# ------------------------------------------------------------------------------
# 候选基因列表 (来自 Fig 7A 的 13 个基因)
candidate_genes <- c("SMTN", "SEMA3C", "RHOB", "PPP1R14A", "PDLIM7", "PCP4",
                     "FHL2", "DSTN", "DES", "CUX1", "CALD1", "CA12", "ACTN1")

# 检查基因是否存在
exist_genes <- candidate_genes[candidate_genes %in% colnames(analysis_data)]
if (length(exist_genes) < length(candidate_genes)) {
  warning("部分基因未在表达矩阵中找到")
}

# 构建 X 矩阵 (基因表达) 和 Y 对象 (生存信息)
x <- as.matrix(analysis_data[, exist_genes])
y <- Surv(analysis_data$OS_time, analysis_data$OS_status)

# ------------------------------------------------------------------------------
# 4. 运行 Lasso 回归
# ------------------------------------------------------------------------------
message("正在运行 Lasso Cox 回归...")

# 1. 拟合模型 (用于绘制系数路径图)
fit <- glmnet(x, y, family = "cox")

# 2. 交叉验证 (用于寻找最佳 Lambda)
set.seed(123) # 设置种子以保证结果可复现
cv_fit <- cv.glmnet(x, y, family = "cox", nfolds = 10)

# 输出最佳 Lambda
message(paste("Lambda.min:", round(cv_fit$lambda.min, 4)))
message(paste("Lambda.1se:", round(cv_fit$lambda.1se, 4)))

# ------------------------------------------------------------------------------
# 5. 绘图 (Fig 7B 的两张子图) - 使用 ggplot2 重绘
# ------------------------------------------------------------------------------
if (!dir.exists("results")) dir.create("results")

# 准备绘图数据 1: 系数路径
fit_df <- as.data.frame(as.matrix(t(fit$beta)))
fit_df$LogLambda <- log(fit$lambda)
fit_long <- fit_df %>% 
  pivot_longer(cols = -LogLambda, names_to = "Gene", values_to = "Coefficient")

# 准备绘图数据 2: 交叉验证误差
cv_df <- data.frame(
  LogLambda = log(cv_fit$lambda),
  CVM = cv_fit$cvm,
  CVUP = cv_fit$cvup,
  CVLO = cv_fit$cvlo
)

# 绘制图 A: Lasso 系数路径
p1 <- ggplot(fit_long, aes(x = LogLambda, y = Coefficient, color = Gene)) +
  geom_line(size = 0.8) +
  theme_bw() +
  labs(x = "Log Lambda", y = "Coefficients", title = "Lasso Coefficient Profiles") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "right",
    legend.key.size = unit(0.4, "cm")
  )

# 绘制图 B: 交叉验证误差
p2 <- ggplot(cv_df, aes(x = LogLambda, y = CVM)) +
  geom_errorbar(aes(ymin = CVLO, ymax = CVUP), width = 0.1, color = "gray") +
  geom_point(color = "#E64B35", size = 2) + # 使用类似论文的红色
  geom_vline(xintercept = log(cv_fit$lambda.min), linetype = "dashed", color = "black") +
  geom_vline(xintercept = log(cv_fit$lambda.1se), linetype = "dashed", color = "black") +
  annotate("text", x = log(cv_fit$lambda.min), y = max(cv_df$CVUP), 
           label = paste0("min: ", round(cv_fit$lambda.min, 4)), vjust = -1, size = 3) +
  theme_bw() +
  labs(x = "Log Lambda", y = "Partial Likelihood Deviance", title = "Lasso Cross-Validation") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

# 组合图片
p_combined <- ggarrange(p1, p2, ncol = 2, labels = c("A", "B"))

# 保存
ggsave("results/Fig7B_Lasso.pdf", p_combined, width = 12, height = 6)
ggsave("results/Fig7B_Lasso.png", p_combined, width = 12, height = 6)

message("图片已保存至 results/Fig7B_Lasso.pdf")

# ------------------------------------------------------------------------------
# 6. 提取最终筛选的基因
# ------------------------------------------------------------------------------
# 使用 lambda.min 提取系数
coef_min <- coef(cv_fit, s = "lambda.min")
active_index <- which(coef_min != 0)
active_genes <- rownames(coef_min)[active_index]
active_coefs <- coef_min[active_index]

result_df <- data.frame(
  Gene = active_genes,
  Coefficient = active_coefs
)

# 按系数绝对值排序
result_df <- result_df[order(abs(result_df$Coefficient), decreasing = TRUE), ]

message("Lasso 筛选出的基因:")
print(result_df)

# 保存筛选结果
write.csv(result_df, "results/lasso_selected_genes.csv", row.names = FALSE)
