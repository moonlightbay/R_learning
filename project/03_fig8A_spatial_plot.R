# ==============================================================================
# 脚本名称: 03_fig8A_spatial_plot.R
# 功能: Fig 8A: 空间转录组特征图 (Spatial Feature Map)
#       展示 Tumor Border (Bdy), Malignant (Mal, Mal1), Normal (New) 的分布
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

# Seurat 是空间转录组分析的核心包
packages <- c("Seurat", "tidyverse", "patchwork", "hdf5r")
for (pkg in packages) {
  ensure_package(pkg)
}

setwd("d:/Works/R_learning/project")

# ------------------------------------------------------------------------------
# 2. 加载空间转录组数据
# ------------------------------------------------------------------------------
# 数据路径
data_dir <- "d:/Works/R_learning/data/fig2-fig6/spatial_transcriptomics"

message("正在加载 10x Visium 数据...")

# Seurat 的 Load10X_Spatial 默认期望的文件名是 "filtered_feature_bc_matrix.h5"
# 但您的目录下文件名带有前缀 "Visium_FFPE_Human_Cervical_Cancer_"
# 我们需要显式指定文件名

# 检查文件是否存在
h5_file <- file.path(data_dir, "Visium_FFPE_Human_Cervical_Cancer_filtered_feature_bc_matrix.h5")
if (!file.exists(h5_file)) {
  stop(paste("找不到 H5 文件:", h5_file))
}

# Load10X_Spatial 会自动寻找 spatial 文件夹，这在您的目录结构中是存在的
# 关键是 filename 参数要指向正确的 h5 文件
st_data <- Load10X_Spatial(
  data.dir = data_dir,
  filename = "Visium_FFPE_Human_Cervical_Cancer_filtered_feature_bc_matrix.h5",
  assay = "Spatial",
  slice = "slice1",
  filter.matrix = TRUE
)

message(paste("加载完成，细胞数:", ncol(st_data), "基因数:", nrow(st_data)))

# ------------------------------------------------------------------------------
# 3. 数据预处理 (SCTransform)
# ------------------------------------------------------------------------------
message("正在进行 SCTransform 标准化 (这可能需要几分钟)...")
# SCTransform 替代了传统的 NormalizeData, FindVariableFeatures, ScaleData
# 尝试使用 glmGamPoi 提速
if (requireNamespace("glmGamPoi", quietly = TRUE)) {
  message("检测到 glmGamPoi，将使用其加速 SCTransform...")
  st_data <- SCTransform(st_data, assay = "Spatial", method = "glmGamPoi", verbose = TRUE)
} else {
  message("未检测到 glmGamPoi，使用默认方法 (较慢)...")
  st_data <- SCTransform(st_data, assay = "Spatial", verbose = TRUE)
}

# ------------------------------------------------------------------------------
# 4. 降维与聚类
# ------------------------------------------------------------------------------
message("正在进行降维与聚类...")

# PCA
st_data <- RunPCA(st_data, assay = "SCT", verbose = FALSE)

# UMAP (用于降维可视化，虽然空间图主要看 SpatialDimPlot)
st_data <- RunUMAP(st_data, reduction = "pca", dims = 1:30)

# FindNeighbors
st_data <- FindNeighbors(st_data, reduction = "pca", dims = 1:30)

# FindClusters
# resolution 决定了聚类的细致程度。
# 原文有 4 类 (New, Bdy, Mal, Mal1)，我们尝试 0.3 - 0.5 左右的分辨率
st_data <- FindClusters(st_data, verbose = FALSE, resolution = 0.4)

# ------------------------------------------------------------------------------
# 5. 初步可视化 (查看聚类结果)
# ------------------------------------------------------------------------------
# 保存一个带有数字标签的图，供用户对照
p_raw <- SpatialDimPlot(st_data, label = TRUE, label.size = 3, pt.size.factor = 1.6) + 
  ggtitle("Unsupervised Clustering (Check IDs)")

if (!dir.exists("results")) dir.create("results")
ggsave("results/Fig8A_Spatial_Clustering_Raw.pdf", p_raw, width = 8, height = 8)
# 修正: ggsave 的 width/height 默认单位是英寸(in)，800英寸太大了
# 如果要指定像素，需要配合 dpi 参数，或者直接用 png() 函数
# 这里我们改回合理的英寸尺寸
ggsave("results/Fig8A_Spatial_Clustering_Raw.png", p_raw, width = 10, height = 10, dpi = 300)

message("初步聚类图已保存至 results/Fig8A_Spatial_Clustering_Raw.png")
message("请打开该图片，对比论文原图，确定数字 ID 与 'New', 'Bdy', 'Mal', 'Mal1' 的对应关系。")

# ------------------------------------------------------------------------------
# 5.1 辅助注释：绘制标记基因 (Marker Genes)
# ------------------------------------------------------------------------------
# 为了帮助确定哪个聚类对应 Mal, Normal, Bdy，我们绘制一些经典标记物
# 宫颈癌常见标记:
# 上皮/肿瘤: EPCAM, KRT18, KRT19, CDKN2A (p16)
# 间质/成纤维: COL1A1, ACTA2, VIM
# 免疫: PTPRC (CD45), CD3D, CD68
# 正常鳞状上皮: KRT5, KRT14

markers <- c("EPCAM", "KRT18", "COL1A1", "PTPRC", "KRT5", "CDKN2A")
# 检查这些基因是否存在于数据中
valid_markers <- markers[markers %in% rownames(st_data)]

if (length(valid_markers) > 0) {
  message("正在绘制标记基因图以辅助注释...")
  p_markers <- SpatialFeaturePlot(st_data, features = valid_markers, ncol = 3)
  ggsave("results/Fig8A_Marker_Expression.png", p_markers, width = 12, height = 8, dpi = 300)
} else {
  message("未找到指定的标记基因，跳过 Marker 绘图。")
}

# ------------------------------------------------------------------------------
# 6. 注释与重命名 (需要用户修改此处!)
# ------------------------------------------------------------------------------
# 假设的对应关系 (请根据实际生成的 Raw 图进行修改)
# 论文中的标签为: Normal (原New), Bdy (Border), Mal (Malignant), Mal1
# 原文提到使用了 Cottrazm 进行边界界定，这里我们使用 Seurat 聚类近似
# 例如: 0 -> Mal, 1 -> Normal, 2 -> Bdy, 3 -> Mal1

# 获取当前的聚类 ID
current_ids <- levels(st_data)
message(paste("当前聚类 ID:", paste(current_ids, collapse = ", ")))

# --- [用户修改区域开始] ---
# 请根据 Raw 图的分布修改下面的对应关系
# 这是一个占位符逻辑，实际运行后需要您来指定
# new_labels <- c("Mal", "Normal", "Bdy", "Mal1", "Other") 
# names(new_labels) <- levels(st_data)
# st_data <- RenameIdents(st_data, new_labels)
# --- [用户修改区域结束] ---

# 为了演示，我们暂时直接使用数字 ID 绘图，
# 但将图例标题改为 "Regions"
st_data$region_label <- Idents(st_data)

# ------------------------------------------------------------------------------
# 7. 最终绘图 (模拟 Fig 8A 风格)
# ------------------------------------------------------------------------------
# 自定义颜色 (参考常见的病理染色或论文配色)
# Mal/Mal1 通常用红色/橙色系，Normal 用绿色/蓝色系，Bdy 用中间色
# 这里定义一个调色板，如果类别名改了，这里也会自动适配
# 扩展颜色以支持更多聚类 (防止 "Insufficient values" 报错)
my_cols <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628", "#F781BF", "#999999")

p_final <- SpatialDimPlot(st_data, label = FALSE, pt.size.factor = 1.6, alpha = 0.8) +
  scale_fill_manual(values = my_cols) +
  labs(fill = "Regions") +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold")
  ) +
  ggtitle("Spatial Feature Map")

ggsave("results/Fig8A_Spatial_Feature_Map.pdf", p_final, width = 8, height = 8)
ggsave("results/Fig8A_Spatial_Feature_Map.png", p_final, width = 10, height = 10, dpi = 300)

message("Fig 8A 完成，已保存至 results/Fig8A_Spatial_Feature_Map.pdf")
