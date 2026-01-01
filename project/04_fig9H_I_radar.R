# ==============================================================================
# 脚本名称: 04_fig9H_I_radar.R
# 功能: 
#   Fig 9H: High/Low MFRS 组间 Top 20 差异富集通路雷达图
#   Fig 9I: High/Low MFRS 组间 Top 20 免疫细胞浸润差异雷达图
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

packages <- c("survival", "tidyverse", "data.table", "fmsb", "msigdbr")
for (pkg in packages) {
  ensure_package(pkg)
}

bioc_packages <- c("GSVA", "BiocParallel")
for (pkg in bioc_packages) {
  ensure_bioc_package(pkg)
}

setwd("d:/Works/R_learning/project")

# ------------------------------------------------------------------------------
# 2. 数据读取与预处理
# ------------------------------------------------------------------------------
message("正在读取数据...")
tpm_df <- fread("processed_data/expression_matrix.csv", data.table = FALSE)
rownames(tpm_df) <- tpm_df[, 1]
tpm_df <- tpm_df[, -1]

clin_df <- read.csv("processed_data/clinical_cleaned.csv")

common_samples <- intersect(colnames(tpm_df), clin_df$case_submitter_id)
tpm_data <- tpm_df[, common_samples]
clin_data <- clin_df[clin_df$case_submitter_id %in% common_samples, ]
rownames(clin_data) <- clin_data$case_submitter_id
clin_data <- clin_data[common_samples, ]

# ------------------------------------------------------------------------------
# 3. 计算 Risk Score 和分组
# ------------------------------------------------------------------------------
message("正在计算 Risk Score...")
expr_t <- as.data.frame(t(tpm_data))
expr_t$case_submitter_id <- rownames(expr_t)
analysis_data <- merge(clin_data, expr_t, by = "case_submitter_id")
analysis_data <- analysis_data %>% filter(!is.na(OS_time) & OS_time > 0 & !is.na(OS_status))

risk_genes <- c("SEMA3C", "RHOB", "PPP1R14A", "PCP4", "FHL2", 
                "DES", "CUX1", "CALD1", "CA12", "ACTN1")

multi_formula <- as.formula(paste("Surv(OS_time, OS_status) ~", paste(risk_genes, collapse = " + ")))
multi_cox <- coxph(multi_formula, data = analysis_data)
analysis_data$risk_score <- predict(multi_cox, type = "risk", newdata = analysis_data)

# 分组 (High vs Low)
# 使用中位数分组，或者使用 surv_cutpoint (这里为了雷达图对比明显，使用中位数或最佳截断值均可)
# 为了与之前一致，使用 surv_cutpoint
res.cut <- survminer::surv_cutpoint(analysis_data, time = "OS_time", event = "OS_status", variables = "risk_score")
cutoff_score <- res.cut$cutpoint$cutpoint
analysis_data$Risk_group <- ifelse(analysis_data$risk_score > cutoff_score, "High", "Low")

# ------------------------------------------------------------------------------
# 4. 通用雷达图绘制函数
# ------------------------------------------------------------------------------
draw_radar_plot <- function(data_df, title_text, filename) {
  # data_df 结构: 行名为 High/Low, 列名为变量
  
  # 准备 fmsb 数据结构: Max, Min, Data
  # 自动计算范围 (移除边距，使图形顶格)
  max_val <- max(data_df)
  min_val <- min(data_df)
  
  plot_data <- rbind(rep(max_val, ncol(data_df)), 
                     rep(min_val, ncol(data_df)), 
                     data_df)
  
  # 颜色: High (Red), Low (Blue)
  # 确保数据行顺序为: 第1行 High, 第2行 Low
  colors_border <- c("#8B0000", "#00008B")
  colors_fill <- c(rgb(139/255, 0, 0, 0.2), rgb(0, 0, 139/255, 0.2))
  
  if (!dir.exists("results")) dir.create("results")
  
  # PDF
  pdf(paste0("results/", filename, ".pdf"), width = 10, height = 10)
  radarchart(plot_data,
             axistype = 1,
             pcol = colors_border,
             pfcol = colors_fill,
             plwd = 2,
             cglcol = "grey",
             cglty = 1,
             axislabcol = "grey",
             caxislabels = seq(round(min_val, 2), round(max_val, 2), length.out = 5),
             cglwd = 0.8,
             vlcex = 0.8,
             title = title_text
  )
  legend(x = "bottom", legend = rownames(data_df), horiz = TRUE,
         bty = "n", pch = 20, col = colors_border,
         text.col = "black", cex = 1.2, pt.cex = 2)
  dev.off()
  
  # PNG
  png(paste0("results/", filename, ".png"), width = 2000, height = 2000, res = 300)
  radarchart(plot_data,
             axistype = 1,
             pcol = colors_border,
             pfcol = colors_fill,
             plwd = 2,
             cglcol = "grey",
             cglty = 1,
             axislabcol = "grey",
             caxislabels = seq(round(min_val, 2), round(max_val, 2), length.out = 5),
             cglwd = 0.8,
             vlcex = 0.8,
             title = title_text
  )
  legend(x = "bottom", legend = rownames(data_df), horiz = TRUE,
         bty = "n", pch = 20, col = colors_border,
         text.col = "black", cex = 1.2, pt.cex = 2)
  dev.off()
  
  message(paste("已保存:", filename))
}

# ------------------------------------------------------------------------------
# 5. Fig 9H: 差异富集通路雷达图 (GSVA)
# ------------------------------------------------------------------------------
message("=== 开始处理 Fig 9H (GSVA Pathways) ===")

# 5.1 获取基因集 (同 9G)
message("获取基因集...")
local_gmt_path <- "data/c2.all.v2023.2.Hs.symbols.gmt"
gs_list <- tryCatch({
  m_df <- msigdbr(species = "Homo sapiens", category = "C2")
  split(x = m_df$gene_symbol, f = m_df$gs_name)
}, error = function(e) {
  # 搜索 ../data 或 data
  gmt_files <- list.files("../data", pattern = "^c2.*\\.gmt$", full.names = TRUE)
  if (length(gmt_files) == 0) gmt_files <- list.files("data", pattern = "^c2.*\\.gmt$", full.names = TRUE)
  
  if (length(gmt_files) > 0) {
    message(paste("读取本地 GMT:", gmt_files[1]))
    c2_gmt <- clusterProfiler::read.gmt(gmt_files[1])
    split(c2_gmt$gene, c2_gmt$term)
  } else {
    stop("无法获取基因集，请确保已下载 GMT 文件。")
  }
})

# 5.2 运行 GSVA (带优化)
message("运行 GSVA (C2)...")
expr_mat <- as.matrix(tpm_data)
gsva_cache_file <- "processed_data/gsva_c2_res.rds"

# 预定义并行参数 (供后续使用)
if (exists("gsvaParam")) {
  n_cores <- 4
  if (.Platform$OS.type == "windows") {
    bp_param <- SnowParam(workers = n_cores, progressbar = TRUE)
  } else {
    bp_param <- MulticoreParam(workers = n_cores, progressbar = TRUE)
  }
}

if (file.exists(gsva_cache_file)) {
  message(paste("读取缓存的 GSVA 结果:", gsva_cache_file))
  gsva_res <- readRDS(gsva_cache_file)
} else {
  if (exists("gsvaParam")) {
    gspar <- gsvaParam(exprData = expr_mat, geneSets = gs_list, 
                       kcdf = "Gaussian", minSize = 10, maxSize = 500)
    gsva_res <- gsva(gspar, verbose = TRUE, BPPARAM = bp_param)
  } else {
    gsva_res <- gsva(expr_mat, gs_list, kcdf = "Gaussian", verbose = TRUE, 
                     parallel.sz = 4, min.sz=10, max.sz=500)
  }
  saveRDS(gsva_res, gsva_cache_file)
}

# 5.3 筛选 Top 20 差异/相关通路
message("筛选 Top 20 通路...")
# 关键修正: 确保 analysis_data 行名为样本ID，否则后续索引会产生 NA
rownames(analysis_data) <- analysis_data$case_submitter_id

risk_vec <- analysis_data$risk_score
names(risk_vec) <- analysis_data$case_submitter_id
common_samps <- intersect(colnames(gsva_res), names(risk_vec))
gsva_mat <- gsva_res[, common_samps]
risk_vec <- risk_vec[common_samps]

# 计算相关性
cor_res <- apply(gsva_mat, 1, function(x) cor(x, risk_vec, method = "pearson"))
cor_res <- cor_res[!is.na(cor_res)] # 去除可能产生的 NA
top_pathways <- names(sort(abs(cor_res), decreasing = TRUE)[1:20]) # 选绝对值最大的20个

# 5.4 准备绘图数据
plot_mat <- gsva_mat[top_pathways, ]
group_info <- analysis_data[common_samps, "Risk_group"]

# 计算均值
radar_data_h <- data.frame(t(plot_mat)) %>%
  mutate(Group = group_info) %>%
  group_by(Group) %>%
  summarise(across(everything(), mean)) %>%
  arrange(Group) # High, Low (字母顺序 H 在 L 前)

radar_df_h <- as.data.frame(radar_data_h[, -1])
rownames(radar_df_h) <- radar_data_h$Group

# 简化通路名称 (用户要求保留完整名称)
# colnames(radar_df_h) <- gsub("^KEGG_|^REACTOME_|^BIOCARTA_", "", colnames(radar_df_h))
# colnames(radar_df_h) <- substr(colnames(radar_df_h), 1, 15)

# 5.5 绘图
draw_radar_plot(radar_df_h, "Differential Enriched Pathways (High vs Low Risk)", "Fig9H_Enrichment_Radar")


# ------------------------------------------------------------------------------
# 6. Fig 9I: 免疫细胞浸润雷达图 (ssGSEA)
# ------------------------------------------------------------------------------
message("=== 开始处理 Fig 9I (Immune Infiltration) ===")

# 6.1 定义免疫细胞基因集 (24种常见免疫细胞, Bindea et al.)
immune_genes_list <- list(
  "B cells" = c("CD19", "MS4A1", "CD79A", "CD79B", "BLK"),
  "T cells" = c("CD3D", "CD3E", "CD3G", "CD247"),
  "CD8 T cells" = c("CD8A", "CD8B"),
  "Cytotoxic cells" = c("GZMA", "GZMB", "PRF1", "GNLY"),
  "DC" = c("CCL13", "CD209", "HSD11B1"),
  "Eosinophils" = c("IL5RA", "PRG2", "CLC"),
  "iDC" = c("CD1A", "CD1B", "CD1C", "CD1E", "CLEC10A"),
  "Macrophages" = c("CD163", "CD68", "CD14", "ITGAM"),
  "Mast cells" = c("TPSAB1", "TPSB2", "CPA3", "MS4A2"),
  "Neutrophils" = c("FPR1", "SIGLEC5", "CSF3R", "FCGR3B"),
  "NK CD56bright" = c("IL2RB", "NCAM1", "KLRC1", "TNFRSF10A"),
  "NK CD56dim" = c("KIR2DL1", "KIR2DL3", "KIR3DL1"),
  "NK cells" = c("NCR1", "KLRB1", "CD160", "KLRK1"),
  "pDC" = c("LILRA4", "CLEC4C", "IL3RA"),
  "T helper cells" = c("CD4", "CD40LG"),
  "Tcm" = c("CCR7", "SELL", "TCF7"),
  "Tem" = c("CXCR3", "GZMK", "GZMA"),
  "TFH" = c("CXCR5", "BCL6", "IL21"),
  "Tgd" = c("TRGC1", "TRGC2", "TRDC", "TRDV2"),
  "Th1 cells" = c("TBX21", "IFNG", "TNF", "CXCR3"),
  "Th17 cells" = c("IL17A", "RORC", "IL23R", "CCR6"),
  "Th2 cells" = c("GATA3", "IL4", "IL5", "IL13"),
  "Treg" = c("FOXP3", "IL2RA", "CTLA4", "IKZF2"),
  "Monocytes" = c("CD14", "CD300E", "LYZ")
)

# 6.2 运行 ssGSEA
message("运行 ssGSEA 计算免疫评分...")
ssgsea_cache_file <- "processed_data/ssgsea_imm_res.rds"

if (file.exists(ssgsea_cache_file)) {
  message(paste("读取缓存的 ssGSEA 结果:", ssgsea_cache_file))
  imm_res <- readRDS(ssgsea_cache_file)
} else {
  if (exists("ssgseaParam")) {
    # 新版 GSVA (>= 1.52) 使用 ssgseaParam
    gspar_imm <- ssgseaParam(exprData = expr_mat, geneSets = immune_genes_list)
    imm_res <- gsva(gspar_imm, verbose = TRUE, BPPARAM = bp_param)
  } else {
    # 旧版
    imm_res <- gsva(expr_mat, immune_genes_list, method = "ssgsea", 
                    kcdf = "Gaussian", verbose = TRUE, parallel.sz = 4)
  }
  saveRDS(imm_res, ssgsea_cache_file)
}

# 6.3 准备绘图数据 (全部 24 个，或者选 Top 20)
# 这里直接用全部 24 个，或者为了满足"20个"的要求，去掉几个方差最小的
message("准备免疫雷达图数据...")
imm_mat <- imm_res[, common_samps]

# 选 Top 20 变化最大的细胞类型 (如果需要严格限制为20个)
vars <- apply(imm_mat, 1, var)
top_imm_cells <- names(sort(vars, decreasing = TRUE)[1:min(20, length(immune_genes_list))])
imm_mat <- imm_mat[top_imm_cells, ]

# 计算均值
radar_data_i <- data.frame(t(imm_mat)) %>%
  mutate(Group = group_info) %>%
  group_by(Group) %>%
  summarise(across(everything(), mean)) %>%
  arrange(Group) # High, Low

radar_df_i <- as.data.frame(radar_data_i[, -1])
rownames(radar_df_i) <- radar_data_i$Group

# 6.4 绘图
draw_radar_plot(radar_df_i, "Immune Cell Infiltration (High vs Low Risk)", "Fig9I_Immune_Radar")

message("全部完成！")
