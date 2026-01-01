# ==============================================================================
# 脚本名称: 04_fig9G_gsva.R
# 功能: Fig 9G: High/Low MFRS 组间 GSVA 富集分析热图
# 策略: 
#   1. 读取 TPM 矩阵 (Log2转换后) 用于 Risk Score 计算和 GSVA 分析
#   2. 计算 Risk Score 并分组
#   3. 获取 MSigDB C2 (Curated) 基因集
#   4. 运行 GSVA 分析
#   5. 计算 GSVA 分数与 Risk Score 的相关性
#   6. 选取与 Risk Score 显著正相关的通路绘制热图
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

packages <- c("survival", "tidyverse", "data.table", "survminer", "pheatmap", "msigdbr")
for (pkg in packages) {
  ensure_package(pkg)
}

# 安装 Bioconductor 包
bioc_packages <- c("GSVA", "clusterProfiler")
for (pkg in bioc_packages) {
  ensure_bioc_package(pkg)
}

setwd("d:/Works/R_learning/project")

# ------------------------------------------------------------------------------
# 2. 数据读取与预处理
# ------------------------------------------------------------------------------
message("正在读取数据...")

# 读取 TPM 数据 (用于 GSVA 和 Risk Score)
if (!file.exists("processed_data/expression_matrix.csv")) {
  stop("找不到 processed_data/expression_matrix.csv")
}
message("读取 TPM 矩阵...")
tpm_df <- fread("processed_data/expression_matrix.csv", data.table = FALSE)
rownames(tpm_df) <- tpm_df[, 1]
tpm_df <- tpm_df[, -1]

# 读取临床数据
clin_df <- read.csv("processed_data/clinical_cleaned.csv")

# 确保样本交集
common_samples <- intersect(colnames(tpm_df), clin_df$case_submitter_id)
message(paste("共有样本数:", length(common_samples)))

tpm_data <- tpm_df[, common_samples]
clin_data <- clin_df[clin_df$case_submitter_id %in% common_samples, ]
rownames(clin_data) <- clin_data$case_submitter_id
clin_data <- clin_data[common_samples, ]

# ------------------------------------------------------------------------------
# 3. 计算 Risk Score 和分组
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

# ------------------------------------------------------------------------------
# 4. 准备 GSVA 基因集 (MSigDB C2)
# ------------------------------------------------------------------------------
message("正在获取 MSigDB C2 基因集...")

# 定义本地 GMT 文件路径 (作为备选方案)
# 如果 msigdbr 失败，请手动下载 c2.all.v*.symbols.gmt 并放在 data 目录下
local_gmt_path <- "data/c2.all.v2023.2.Hs.symbols.gmt"

# 尝试获取基因集
gs_list <- tryCatch({
  # 优先尝试使用 msigdbr 在线获取
  message("尝试使用 msigdbr 获取基因集...")
  m_df <- msigdbr(species = "Homo sapiens", category = "C2")
  split(x = m_df$gene_symbol, f = m_df$gs_name)
  
}, error = function(e) {
  message("警告: msigdbr 连接失败 (可能是网络问题)。")
  
  # 如果失败，检查本地是否有 GMT 文件
  # 搜索 ../data (项目根目录) 或 data (当前目录) 下任何以 c2 开头并以 .gmt 结尾的文件
  gmt_files <- list.files("../data", pattern = "^c2.*\\.gmt$", full.names = TRUE)
  if (length(gmt_files) == 0) {
    gmt_files <- list.files("data", pattern = "^c2.*\\.gmt$", full.names = TRUE)
  }
  
  if (length(gmt_files) > 0) {
    message(paste("发现本地 GMT 文件，正在读取:", gmt_files[1]))
    c2_gmt <- clusterProfiler::read.gmt(gmt_files[1])
    split(c2_gmt$gene, c2_gmt$term)
  } else {
    message("错误: 未找到本地 GMT 文件，且无法连接 msigdbr。")
    message("请执行以下步骤解决:")
    message("1. 访问 GSEA 官网: https://www.gsea-msigdb.org/gsea/downloads.jsp")
    message("2. 下载 'c2.all.v*.Hs.symbols.gmt' 文件 (Human C2 Gene Sets)")
    message("3. 将文件放入项目的 'data' 文件夹中")
    stop("终止运行: 缺少基因集数据。")
  }
})

message(paste("获取到基因集数量:", length(gs_list)))

# ------------------------------------------------------------------------------
# 5. 运行 GSVA
# ------------------------------------------------------------------------------
message("正在运行 GSVA (这可能需要几分钟)...")

# 转换数据为矩阵
expr_mat <- as.matrix(tpm_data)

# 运行 GSVA
# kcdf="Gaussian" 适用于 log2 转换后的连续型数据 (如 Log2TPM)
# 注意: GSVA 1.52+ 版本 API 发生变化，需使用 gsvaParam
if (exists("gsvaParam")) {
  message("检测到新版 GSVA (>= 1.52)，使用 gsvaParam 方法...")
  
  # 1. 设置并行计算 (大幅提速)
  if (!require("BiocParallel", quietly = TRUE)) library(BiocParallel)
  
  # Windows 下建议使用 SnowParam，核心数设为 4
  n_cores <- 4
  if (.Platform$OS.type == "windows") {
    bp_param <- SnowParam(workers = n_cores, progressbar = TRUE)
  } else {
    bp_param <- MulticoreParam(workers = n_cores, progressbar = TRUE)
  }
  message(paste("已启用并行计算，核心数:", n_cores))

  # 2. 构造参数对象并过滤基因集
  # 关键优化: 仅分析基因数在 10 到 500 之间的基因集
  # 原因: 太小的集合统计效力低，太大的集合生物学意义模糊且计算极慢
  message("优化策略: 过滤基因集 (Min=10, Max=500)...")
  
  gspar <- gsvaParam(exprData = expr_mat, geneSets = gs_list, 
                     kcdf = "Gaussian",
                     minSize = 10, 
                     maxSize = 500)
  
  # 运行分析
  gsva_res <- gsva(gspar, verbose = TRUE, BPPARAM = bp_param)
  
} else {
  message("检测到旧版 GSVA，使用传统方法...")
  # 旧版也加上 min.sz 和 max.sz 限制
  gsva_res <- gsva(expr_mat, gs_list, kcdf = "Gaussian", verbose = TRUE, 
                   parallel.sz = 4, min.sz=10, max.sz=500)
}

# ------------------------------------------------------------------------------
# 6. 筛选与 Risk Score 相关的通路
# ------------------------------------------------------------------------------
message("正在计算相关性...")

# 提取 Risk Score
risk_scores <- analysis_data$risk_score
names(risk_scores) <- analysis_data$case_submitter_id

# 确保样本顺序一致
common_samps_gsva <- intersect(colnames(gsva_res), names(risk_scores))
gsva_mat <- gsva_res[, common_samps_gsva]
risk_vec <- risk_scores[common_samps_gsva]

# 计算每个通路与 Risk Score 的相关性
cor_res <- apply(gsva_mat, 1, function(x) {
  cor.test(x, risk_vec, method = "pearson")$estimate
})

# 排序
cor_res <- sort(cor_res, decreasing = TRUE)
  
# 选取 Top 20 正相关通路
top_pathways <- names(head(cor_res, 20))

message("Top 5 正相关通路:")
print(head(top_pathways, 5))

# 检查文中提到的特定通路是否在结果中 (可选)
target_pathways <- c("PEDERSEN_TARGETS_OF_611_CTF_ISOFORM_OF_ERBB2",
                     "REACTOME_CELL_EXTRACELLULAR_MATRIX_INTERACTIONS",
                     "BRUECKNER_TARGETS_OF_MIRLET7A3_DN",
                     "ELVIDGE_HYPOXIA_UP")

# 如果这些通路存在于 GSVA 结果中，优先展示它们
available_targets <- intersect(target_pathways, rownames(gsva_mat))
if (length(available_targets) > 0) {
  message("发现文中提到的目标通路:")
  print(available_targets)
  # 将目标通路加入绘图列表 (去重)
  plot_pathways <- unique(c(available_targets, top_pathways))
  # 截取前 25 个以防太多
  if(length(plot_pathways) > 25) plot_pathways <- plot_pathways[1:25]
} else {
  plot_pathways <- top_pathways
}

# ------------------------------------------------------------------------------
# 7. 绘制相关性热图 (修正为下三角相关性矩阵)
# ------------------------------------------------------------------------------
message("正在绘制相关性热图...")

ensure_package("corrplot")

# 准备数据: 样本为行，变量(通路+RiskScore)为列
# 使用之前筛选好的 plot_pathways
gsva_plot_data <- t(gsva_mat[plot_pathways, ])

# 合并 Risk Score
plot_data <- data.frame(gsva_plot_data)
# 确保行名一致
plot_data <- plot_data[names(risk_vec), ]
plot_data$Risk_Score <- risk_vec

# 计算相关性矩阵
M <- cor(plot_data, method = "pearson")

# 颜色设置
col_corr <- colorRampPalette(c("navy", "white", "firebrick3"))(200)

if (!dir.exists("results")) dir.create("results")

# PDF
pdf("results/Fig9G_GSVA_Correlation.pdf", width = 12, height = 12)
corrplot(M, 
         method = "square",       # 方块形状
         type = "lower",          # 下三角
         tl.col = "black",        # 文本颜色
         tl.cex = 0.6,            # 文本大小
         col = col_corr,          # 颜色
         title = "Correlation between Risk Score and GSVA Pathways",
         mar = c(0,0,2,0),        # 调整边距
         diag = TRUE)             # 显示对角线
dev.off()

# PNG
png("results/Fig9G_GSVA_Correlation.png", width = 3000, height = 3000, res = 300)
corrplot(M, 
         method = "square", 
         type = "lower", 
         tl.col = "black", 
         tl.cex = 0.6, 
         col = col_corr,
         title = "Correlation between Risk Score and GSVA Pathways",
         mar = c(0,0,2,0),
         diag = TRUE)
dev.off()

message("Fig 9G 完成，已保存至 results/Fig9G_GSVA_Correlation.pdf")
