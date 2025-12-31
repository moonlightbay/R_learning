# ==============================================================================
# 脚本名称: 04_fig9A_heatmap.R
# 功能: Fig 9A: High/Low MFRS 组间差异基因热图 (使用 DESeq2)
# 策略: 
#   1. 读取 TPM 矩阵用于 Risk Score 计算和分组
#   2. 读取 Raw Counts 矩阵用于 DESeq2 差异分析
#   3. 使用 DESeq2 进行差异分析
#   4. 筛选显著差异基因 (padj < 0.05, |Log2FC| > 1)
#   5. 选取最显著的 30 个基因绘制热图
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

packages <- c("survival", "tidyverse", "data.table", "pheatmap", "survminer")
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
group_info <- analysis_data[, c("case_submitter_id", "Risk_group", "risk_score")]
rownames(group_info) <- group_info$case_submitter_id

# ------------------------------------------------------------------------------
# 4. 运行 DESeq2
# ------------------------------------------------------------------------------
message("正在运行 DESeq2 (这可能需要几分钟)...")

# 确保 DESeq2 输入样本一致
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
# 5. 筛选 Top 30 基因
# ------------------------------------------------------------------------------
# 过滤 padj < 0.05 且 |Log2FC| > 1
sig_genes <- results_df %>% 
  filter(!is.na(padj) & padj < 0.05 & abs(log2FoldChange) > 1)

message(paste("满足 padj < 0.05 且 |Log2FC| > 1 的基因数量:", nrow(sig_genes)))

# 选取 Top 30 (按 padj 排序)
if (nrow(sig_genes) > 30) {
  top_30_genes <- sig_genes %>% 
    arrange(padj) %>% 
    head(30)
} else {
  top_30_genes <- sig_genes
  if (nrow(top_30_genes) == 0) {
    stop("没有基因满足筛选条件，请检查数据或放宽阈值。")
  }
}

# 按照 Log2FC 排序 (上调在前，下调在后)
top_30_genes <- top_30_genes %>% arrange(desc(log2FoldChange))

message(paste("最终用于绘图的基因数量:", nrow(top_30_genes)))
print(top_30_genes[, c("Gene", "log2FoldChange", "padj")])

# ------------------------------------------------------------------------------
# 6. 绘制热图
# ------------------------------------------------------------------------------
message("正在绘制热图...")

# 获取标准化后的表达矩阵 (VST 转换，适合热图)
# vst 转换后的数据近似正态分布，且消除了均值-方差依赖性
vst_data <- assay(vst(dds, blind = FALSE))

# 提取 Top 30 基因的表达量
plot_genes <- top_30_genes$Gene
plot_mat <- vst_data[plot_genes, ]

# 排序样本: 按 Risk Group 和 Risk Score 排序
# 注意: col_data_for_deseq 已经包含了 Risk_group 和 risk_score
ord <- order(col_data_for_deseq$Risk_group, col_data_for_deseq$risk_score)
plot_mat <- plot_mat[, ord]
annotation_col <- data.frame(Risk_group = col_data_for_deseq$Risk_group[ord])
rownames(annotation_col) <- colnames(plot_mat)

# 颜色设置
ann_colors <- list(
  Risk_group = c(Low = "forestgreen", High = "firebrick3")
)

# 保存
if (!dir.exists("results")) dir.create("results")

# PDF
pdf("results/Fig9A_Heatmap.pdf", width = 8, height = 10)
pheatmap(plot_mat,
         annotation_col = annotation_col,
         annotation_colors = ann_colors,
         scale = "row",       # 按行(基因)归一化，显示相对表达量
         show_colnames = FALSE, # 样本太多，不显示列名
         show_rownames = TRUE,
         cluster_cols = FALSE, # 不聚类列，因为我们已经按组排序了
         cluster_rows = FALSE, # 不聚类行，保持我们按 LogFC 排序的顺序
         color = colorRampPalette(c("forestgreen", "white", "firebrick3"))(100),
         main = "Expression Profile of Top DEGs (High vs Low MFRS) - DESeq2",
         fontsize_row = 10)
dev.off()

# PNG
png("results/Fig9A_Heatmap.png", width = 2400, height = 3000, res = 300)
pheatmap(plot_mat,
         annotation_col = annotation_col,
         annotation_colors = ann_colors,
         scale = "row",
         show_colnames = FALSE,
         show_rownames = TRUE,
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         color = colorRampPalette(c("forestgreen", "white", "firebrick3"))(100),
         main = "Expression Profile of Top DEGs (High vs Low MFRS) - DESeq2",
         fontsize_row = 10)
dev.off()

message("Fig 9A 完成，已保存至 results/Fig9A_Heatmap.pdf")
