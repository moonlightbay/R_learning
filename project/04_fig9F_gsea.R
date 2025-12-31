# ==============================================================================
# 脚本名称: 04_fig9F_gsea.R
# 功能: Fig 9F: High/Low MFRS 组间 GSEA 富集分析 (GO-BP)
# 策略: (混合策略 Hybrid Strategy)
#   1. 读取 TPM 矩阵用于 Risk Score 计算和分组
#   2. 读取 Raw Counts 矩阵用于 DESeq2 差异分析
#   3. 获取所有基因的 Log2FC 排序列表
#   4. 使用 clusterProfiler 进行 GSEA 分析
#   5. 绘制 GSEA 富集图 (Top Up & Top Down)
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
# 4. 运行 DESeq2 获取所有基因的 Log2FC
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

# ------------------------------------------------------------------------------
# 5. 准备 GSEA 排序列表
# ------------------------------------------------------------------------------
message("正在准备 GSEA 排序列表...")

# 移除 Log2FC 为 NA 的行
results_df <- results_df %>% filter(!is.na(log2FoldChange))

# ID 转换 (Symbol -> Entrez ID)
gene_ids <- bitr(results_df$Gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")

# 合并 Entrez ID 和 Log2FC
# 注意：bitr 可能会导致一对多或多对一，这里简单处理，取第一个匹配
merged_df <- merge(results_df, gene_ids, by.x = "Gene", by.y = "SYMBOL")

# 去重：如果同一个 Entrez ID 对应多个 Symbol，取 Log2FC 绝对值最大的那个
merged_df <- merged_df %>% 
  group_by(ENTREZID) %>% 
  arrange(desc(abs(log2FoldChange))) %>% 
  slice(1) %>% 
  ungroup()

# 创建排序列表
gene_list <- merged_df$log2FoldChange
names(gene_list) <- merged_df$ENTREZID

# 必须按降序排列
gene_list <- sort(gene_list, decreasing = TRUE)

message(paste("用于 GSEA 的基因数量:", length(gene_list)))

# ------------------------------------------------------------------------------
# 6. 运行 GSEA (GO-BP)
# ------------------------------------------------------------------------------
message("正在运行 GSEA (GO-BP)...")

# 设置种子以保证可重复性
set.seed(123)

gse_bp <- gseGO(geneList     = gene_list,
                OrgDb        = org.Hs.eg.db,
                ont          = "BP",
                minGSSize    = 10,
                maxGSSize    = 500,
                pvalueCutoff = 0.05,
                verbose      = FALSE,
                eps          = 0) # 设置 eps=0 以获得更精确的 P 值

if (is.null(gse_bp) || nrow(gse_bp) == 0) {
  stop("未发现显著的 GSEA 富集结果。")
}

message(paste("显著富集的通路数量:", nrow(gse_bp)))

# ------------------------------------------------------------------------------
# 7. 绘制 GSEA 结果图
# ------------------------------------------------------------------------------
if (!dir.exists("results")) dir.create("results")

# 挑选 Top 5 Up (NES > 0) 和 Top 5 Down (NES < 0)
gse_res <- as.data.frame(gse_bp)

top_up <- gse_res %>% 
  filter(NES > 0) %>% 
  arrange(desc(NES)) %>% 
  head(5) %>% 
  pull(ID)

top_down <- gse_res %>% 
  filter(NES < 0) %>% 
  arrange(NES) %>% 
  head(5) %>% 
  pull(ID)

selected_pathways <- c(top_up, top_down)

if (length(selected_pathways) > 0) {
  message("正在绘制 GSEA 详细富集图...")
  
  # 绘制多通路 GSEA 图
  p <- gseaplot2(gse_bp, 
                 geneSetID = selected_pathways, 
                 title = "GSEA Enrichment (Top 5 Up & Down)",
                 pvalue_table = TRUE,
                 ES_geom = "line") # 使用线条展示 ES
  
  pdf("results/Fig9F_GSEA.pdf", width = 12, height = 10)
  print(p)
  dev.off()
  
  png("results/Fig9F_GSEA.png", width = 3600, height = 3000, res = 300)
  print(p)
  dev.off()
  
  # 保存结果表格
  write.csv(gse_res, "results/Fig9F_GSEA_table.csv")
  
  message("Fig 9F 完成，已保存至 results/Fig9F_GSEA.pdf")
} else {
  message("没有找到足够的显著通路进行绘图。")
}
