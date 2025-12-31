# ==============================================================================
# 脚本名称: 04_fig9B_volcano.R
# 功能: Fig 9B: High/Low MFRS 组间差异基因火山图 (使用 DESeq2)
# 策略: 
#   1. 读取 TPM 矩阵用于 Risk Score 计算和分组
#   2. 读取 Raw Counts 矩阵用于 DESeq2 差异分析
#   3. 运行 DESeq2
#   4. 绘制火山图 (x=Log2FC, y=-log10(P))
#   5. 标记显著基因 (|Log2FC|>1, p<0.05)
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

# BiocManager 包安装辅助函数
ensure_bioc_package <- function(pkg) {
  if (!require(pkg, character.only = TRUE)) {
    message(paste("正在安装 Bioconductor 包", pkg, "..."))
    if (!require("BiocManager", quietly = TRUE))
        install.packages("BiocManager", repos = "https://cloud.r-project.org")
    BiocManager::install(pkg)
    if (!require(pkg, character.only = TRUE)) {
      stop(paste("无法安装或加载包:", pkg))
    }
  }
}

packages <- c("survival", "tidyverse", "data.table", "survminer", "ggrepel")
for (pkg in packages) {
  ensure_package(pkg)
}

# 安装 DESeq2
ensure_bioc_package("DESeq2")

setwd("d:/Works/R_learning/project")

# ------------------------------------------------------------------------------
# 2. 数据读取与预处理 (Hybrid Strategy)
# ------------------------------------------------------------------------------
message("正在读取数据...")

# 1. 读取 TPM 数据 (用于生存分析/Risk Score)
if (!file.exists("processed_data/expression_matrix.csv")) {
  stop("找不到 processed_data/expression_matrix.csv")
}
message("读取 TPM 矩阵 (用于 Risk Score)...")
tpm_df <- fread("processed_data/expression_matrix.csv", data.table = FALSE)
rownames(tpm_df) <- tpm_df[, 1]
tpm_df <- tpm_df[, -1]

# 2. 读取 Raw Counts 数据 (用于 DESeq2)
if (!file.exists("processed_data/expression_matrix_counts.csv")) {
  stop("找不到 processed_data/expression_matrix_counts.csv")
}
message("读取 Raw Counts 矩阵 (用于 DESeq2)...")
counts_df <- fread("processed_data/expression_matrix_counts.csv", data.table = FALSE)
rownames(counts_df) <- counts_df[, 1]
counts_df <- counts_df[, -1]

# 3. 读取临床数据
clin_df <- read.csv("processed_data/clinical_cleaned.csv")

# 确保样本交集 (TPM, Counts, Clinical 三者交集)
common_samples <- intersect(colnames(tpm_df), colnames(counts_df))
common_samples <- intersect(common_samples, clin_df$case_submitter_id)
message(paste("共有样本数:", length(common_samples)))

# 对齐数据
tpm_data <- tpm_df[, common_samples]
counts_data <- counts_df[, common_samples]
clin_data <- clin_df[clin_df$case_submitter_id %in% common_samples, ]
rownames(clin_data) <- clin_data$case_submitter_id
clin_data <- clin_data[common_samples, ] # 确保顺序一致

# ------------------------------------------------------------------------------
# 3. 计算 Risk Score 和分组 (基于 TPM)
# ------------------------------------------------------------------------------
message("正在计算 Risk Score...")

# 准备 Cox 模型输入数据 (tpm_data 已经是 Log2(TPM+1))
expr_t <- as.data.frame(t(tpm_data))
expr_t$case_submitter_id <- rownames(expr_t)
analysis_data <- merge(clin_data, expr_t, by = "case_submitter_id")

# 过滤无效生存时间
analysis_data <- analysis_data %>%
  filter(!is.na(OS_time) & OS_time > 0 & !is.na(OS_status))

# 关键基因
risk_genes <- c("SEMA3C", "RHOB", "PPP1R14A", "PCP4", "FHL2", 
                "DES", "CUX1", "CALD1", "CA12", "ACTN1")

# 3. Cox 模型计算 Risk Score
multi_formula <- as.formula(paste("Surv(OS_time, OS_status) ~", paste(risk_genes, collapse = " + ")))
multi_cox <- coxph(multi_formula, data = analysis_data)
analysis_data$risk_score <- predict(multi_cox, type = "risk", newdata = analysis_data)

# 4. 分组 (使用最佳截断值)
res.cut <- surv_cutpoint(analysis_data, time = "OS_time", event = "OS_status", variables = "risk_score")
cutoff_score <- res.cut$cutpoint$cutpoint
message(paste("最佳截断值 (Optimal Cutoff):", round(cutoff_score, 4)))

analysis_data$Risk_group <- ifelse(analysis_data$risk_score > cutoff_score, "High", "Low")

# 提取分组信息用于 DESeq2
group_info <- analysis_data[, c("case_submitter_id", "Risk_group")]
rownames(group_info) <- group_info$case_submitter_id

# ------------------------------------------------------------------------------
# 4. 运行 DESeq2 (使用 Raw Counts)
# ------------------------------------------------------------------------------
message("正在运行 DESeq2 (这可能需要几分钟)...")

# 确保 DESeq2 输入样本一致 (Cox模型可能过滤掉了一些无生存数据的样本)
valid_samples <- rownames(group_info)
counts_for_deseq <- counts_data[, valid_samples]
col_data_for_deseq <- group_info[valid_samples, , drop = FALSE]

# 确保因子水平顺序 (Low 为对照组)
col_data_for_deseq$Risk_group <- factor(col_data_for_deseq$Risk_group, levels = c("Low", "High"))

# 创建 DESeqDataSet 对象
dds <- DESeqDataSetFromMatrix(countData = counts_for_deseq,
                              colData = col_data_for_deseq,
                              design = ~ Risk_group)

# 过滤低表达基因
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]

# 运行 DESeq
dds <- DESeq(dds)

# 获取结果 (High vs Low)
res <- results(dds, contrast = c("Risk_group", "High", "Low"))
results_df <- as.data.frame(res)
results_df$Gene <- rownames(results_df)

# ------------------------------------------------------------------------------
# 5. 整理绘图数据 (使用 padj)
# ------------------------------------------------------------------------------
# 定义阈值
logfc_cutoff <- 1
pval_cutoff <- 0.05

# 处理 NA 值
results_df <- results_df %>% filter(!is.na(padj) & !is.na(log2FoldChange))

# 分类 (使用 padj)
results_df$Group <- "Not Significant"
results_df$Group[results_df$padj < pval_cutoff & results_df$log2FoldChange > logfc_cutoff] <- "Up"
results_df$Group[results_df$padj < pval_cutoff & results_df$log2FoldChange < -logfc_cutoff] <- "Down"

# 统计数量
message("差异基因统计 (基于 padj < 0.05):")
print(table(results_df$Group))

# 挑选要标记的基因 (Label)
top_genes <- results_df %>%
  filter(Group != "Not Significant") %>%
  arrange(desc(abs(log2FoldChange))) %>%
  head(10)

results_df$Label <- NA
results_df$Label[results_df$Gene %in% top_genes$Gene] <- results_df$Gene[results_df$Gene %in% top_genes$Gene]

# ------------------------------------------------------------------------------
# 6. 绘制火山图 (Y轴使用 -log10(padj))
# ------------------------------------------------------------------------------
message("正在绘制火山图...")

# 颜色设置
cols <- c("Down" = "forestgreen", "Not Significant" = "grey", "Up" = "firebrick3")

p <- ggplot(results_df, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = Group), alpha = 0.6, size = 1.5) +
  scale_color_manual(values = cols) +
  
  # 阈值线
  geom_vline(xintercept = c(-logfc_cutoff, logfc_cutoff), linetype = "dashed", color = "black", size = 0.5) +
  geom_hline(yintercept = -log10(pval_cutoff), linetype = "dashed", color = "black", size = 0.5) +
  
  # 标签
  geom_text_repel(aes(label = Label), 
                  box.padding = 0.5, 
                  point.padding = 0.3,
                  max.overlaps = 20,
                  size = 3) +
  
  labs(title = "Volcano Plot of DEGs (High vs Low MFRS) - DESeq2",
       x = "Log2 Fold Change",
       y = "-Log10(Adjusted P-value)",
       color = "Expression") +
  
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "top")


# 保存
if (!dir.exists("results")) dir.create("results")

ggsave("results/Fig9B_Volcano.pdf", plot = p, width = 8, height = 7)
ggsave("results/Fig9B_Volcano.png", plot = p, width = 2400, height = 2100, units = "px", dpi = 300)

message("Fig 9B 完成，已保存至 results/Fig9B_Volcano.pdf")
