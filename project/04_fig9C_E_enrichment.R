# ==============================================================================
# 脚本名称: 04_fig9C_E_enrichment.R
# 功能: Fig 9C-E: High/Low MFRS 组间差异基因富集分析 (GO-BP, GO-CC, KEGG)
# 策略: (混合策略 Hybrid Strategy)
#   1. 读取 TPM 矩阵用于 Risk Score 计算和分组
#   2. 读取 Raw Counts 矩阵用于 DESeq2 差异分析
#   3. 筛选显著差异基因 (padj < 0.05, |Log2FC| > 1)
#   4. 使用 clusterProfiler 进行 GO (BP, CC) 和 KEGG 富集分析
#   5. 绘制条形图
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

packages <- c("survival", "tidyverse", "data.table", "survminer", "ggplot2")
for (pkg in packages) {
  ensure_package(pkg)
}

# 安装生物信息学相关包
bioc_packages <- c("DESeq2", "clusterProfiler", "org.Hs.eg.db", "enrichplot", "DOSE")
for (pkg in bioc_packages) {
  ensure_bioc_package(pkg)
}

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

# 确保样本交集
common_samples <- intersect(colnames(tpm_df), colnames(counts_df))
common_samples <- intersect(common_samples, clin_df$case_submitter_id)
message(paste("共有样本数:", length(common_samples)))

# 对齐数据
tpm_data <- tpm_df[, common_samples]
counts_data <- counts_df[, common_samples]
clin_data <- clin_df[clin_df$case_submitter_id %in% common_samples, ]
rownames(clin_data) <- clin_data$case_submitter_id
clin_data <- clin_data[common_samples, ]

# ------------------------------------------------------------------------------
# 3. 计算 Risk Score 和分组 (基于 TPM)
# ------------------------------------------------------------------------------
message("正在计算 Risk Score...")

# 准备 Cox 模型输入数据
expr_t <- as.data.frame(t(tpm_data))
expr_t$case_submitter_id <- rownames(expr_t)
analysis_data <- merge(clin_data, expr_t, by = "case_submitter_id")

# 过滤无效生存时间
analysis_data <- analysis_data %>%
  filter(!is.na(OS_time) & OS_time > 0 & !is.na(OS_status))

# 关键基因
risk_genes <- c("SEMA3C", "RHOB", "PPP1R14A", "PCP4", "FHL2", 
                "DES", "CUX1", "CALD1", "CA12", "ACTN1")

# Cox 模型
multi_formula <- as.formula(paste("Surv(OS_time, OS_status) ~", paste(risk_genes, collapse = " + ")))
multi_cox <- coxph(multi_formula, data = analysis_data)
analysis_data$risk_score <- predict(multi_cox, type = "risk", newdata = analysis_data)

# 分组
res.cut <- surv_cutpoint(analysis_data, time = "OS_time", event = "OS_status", variables = "risk_score")
cutoff_score <- res.cut$cutpoint$cutpoint
analysis_data$Risk_group <- ifelse(analysis_data$risk_score > cutoff_score, "High", "Low")

# 提取分组信息
group_info <- analysis_data[, c("case_submitter_id", "Risk_group")]
rownames(group_info) <- group_info$case_submitter_id

# ------------------------------------------------------------------------------
# 4. 运行 DESeq2 获取差异基因
# ------------------------------------------------------------------------------
message("正在运行 DESeq2...")

valid_samples <- rownames(group_info)
counts_for_deseq <- counts_data[, valid_samples]
col_data_for_deseq <- group_info[valid_samples, , drop = FALSE]
col_data_for_deseq$Risk_group <- factor(col_data_for_deseq$Risk_group, levels = c("Low", "High"))

dds <- DESeqDataSetFromMatrix(countData = counts_for_deseq,
                              colData = col_data_for_deseq,
                              design = ~ Risk_group)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]
dds <- DESeq(dds)
res <- results(dds, contrast = c("Risk_group", "High", "Low"))
results_df <- as.data.frame(res)
results_df$Gene <- rownames(results_df)

# 筛选显著差异基因
sig_genes_df <- results_df %>% 
  filter(!is.na(padj) & padj < 0.05 & abs(log2FoldChange) > 1)

message(paste("显著差异基因数量:", nrow(sig_genes_df)))

if (nrow(sig_genes_df) == 0) {
  stop("没有找到显著差异基因，无法进行富集分析。")
}

# ------------------------------------------------------------------------------
# 5. ID 转换 (Symbol -> Entrez ID)
# ------------------------------------------------------------------------------
message("正在进行 ID 转换...")

gene_list <- sig_genes_df$Gene

# 使用 bitr 进行转换
gene_ids <- bitr(gene_list, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")

message(paste("成功转换 ID 数量:", nrow(gene_ids)))

# ------------------------------------------------------------------------------
# 6. 富集分析 (GO-BP, GO-CC, KEGG)
# ------------------------------------------------------------------------------
if (!dir.exists("results")) dir.create("results")

# 定义颜色梯度: 黑(显著/低P) -> 深蓝 -> 蓝绿 -> 淡绿(不显著/高P)
# 对应用户要求: 淡绿- 蓝绿 -深蓝- 黑 (从高P到低P)
custom_cols <- c("black", "darkblue", "#20B2AA", "#90EE90") # Black, DarkBlue, LightSeaGreen, LightGreen

# --- 6.1 GO Biological Process (BP) ---
message("正在进行 GO-BP 富集分析...")
ego_bp <- enrichGO(gene          = gene_ids$ENTREZID,
                   OrgDb         = org.Hs.eg.db,
                   ont           = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 0.05,
                   readable      = TRUE)

if (!is.null(ego_bp) && nrow(ego_bp) > 0) {
  p <- barplot(ego_bp, showCategory = 20, title = "GO Biological Process Enrichment") +
       scale_fill_gradientn(colours = custom_cols) +
       theme(plot.title = element_text(hjust = 0.5))
  
  pdf("results/Fig9C_GO_BP.pdf", width = 10, height = 8)
  print(p)
  dev.off()
  
  png("results/Fig9C_GO_BP.png", width = 3000, height = 2400, res = 300)
  print(p)
  dev.off()
  
  write.csv(as.data.frame(ego_bp), "results/Fig9C_GO_BP_table.csv")
} else {
  message("未发现显著的 GO-BP 富集结果。")
}

# --- 6.2 GO Cellular Component (CC) ---
message("正在进行 GO-CC 富集分析...")
ego_cc <- enrichGO(gene          = gene_ids$ENTREZID,
                   OrgDb         = org.Hs.eg.db,
                   ont           = "CC",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 0.05,
                   readable      = TRUE)

if (!is.null(ego_cc) && nrow(ego_cc) > 0) {
  p <- barplot(ego_cc, showCategory = 20, title = "GO Cellular Component Enrichment") +
       scale_fill_gradientn(colours = custom_cols) +
       theme(plot.title = element_text(hjust = 0.5))
  
  pdf("results/Fig9D_GO_CC.pdf", width = 10, height = 8)
  print(p)
  dev.off()
  
  png("results/Fig9D_GO_CC.png", width = 3000, height = 2400, res = 300)
  print(p)
  dev.off()
  
  write.csv(as.data.frame(ego_cc), "results/Fig9D_GO_CC_table.csv")
} else {
  message("未发现显著的 GO-CC 富集结果。")
}

# --- 6.3 KEGG Pathway ---
message("正在进行 KEGG 富集分析...")
kk <- enrichKEGG(gene         = gene_ids$ENTREZID,
                 organism     = 'hsa',
                 pvalueCutoff = 0.05)

if (!is.null(kk) && nrow(kk) > 0) {
  p <- barplot(kk, showCategory = 20, title = "KEGG Pathway Enrichment") +
       scale_fill_gradientn(colours = custom_cols) +
       theme(plot.title = element_text(hjust = 0.5))
  
  pdf("results/Fig9E_KEGG.pdf", width = 10, height = 8)
  print(p)
  dev.off()
  
  png("results/Fig9E_KEGG.png", width = 3000, height = 2400, res = 300)
  print(p)
  dev.off()
  
  write.csv(as.data.frame(kk), "results/Fig9E_KEGG_table.csv")
} else {
  message("未发现显著的 KEGG 富集结果。")
}

message("Fig 9C-E 富集分析完成。")
