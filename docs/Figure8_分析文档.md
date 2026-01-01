# Figure 8 分析文档：肿瘤边界特征与预后模型验证

---

## 图8A（1）：空间转录组特征图

### （1）代码文件名
**03_fig8A_scRNA_process.R** 和 **03_fig8A_spatial_plot.R**

**说明**：图8A对应两个代码文件，分别处理scRNA-seq参考数据和空间转录组数据。

### （2）依赖的包

**scRNA_process.R**:
- `Seurat` - 单细胞数据分析核心包
- `tidyverse` - 数据处理
- `patchwork` - 图片组合
- `Matrix` - 稀疏矩阵处理

**spatial_plot.R**:
- `Seurat` - 空间转录组分析
- `tidyverse` - 数据处理
- `hdf5r` - 读取H5文件

### （3）方法分析

#### **文件1：03_fig8A_scRNA_process.R**

**目标**：处理单细胞RNA测序参考数据集（E-MTAB-12305），为空间转录组反卷积做准备。

**步骤**：

1. **加载10x数据**
   ```r
   counts <- Read10X(data.dir = sample_path)
   sc_obj <- CreateSeuratObject(counts = counts, project = sample, 
                                min.cells = 3, min.features = 200)
   ```
   读取10x格式数据（barcodes.tsv.gz, features.tsv.gz, matrix.mtx.gz）

2. **质控（QC）与过滤**
   - 基础过滤：保留至少在3个细胞中表达的基因
   - 保留至少有200个基因表达的细胞
   - （可能包含）线粒体基因过滤、双联体去除等

3. **降维聚类**
   - 标准化
   - 寻找高变基因
   - PCA降维
   - 聚类和UMAP可视化
   - 用于识别成纤维细胞亚群（如C0 MYH11+）

**输出**：预处理后的单细胞参考数据，供空间反卷积使用

#### **文件2：03_fig8A_spatial_plot.R**

**目标**：使用空间转录组数据生成肿瘤边界特征图。

**步骤**：

1. **加载Visium空间数据**
   ```r
   h5_file <- file.path(data_dir, 
                        "Visium_FFPE_Human_Cervical_Cancer_filtered_feature_bc_matrix.h5")
   spatial_data <- Load10X_Spatial(data.dir = data_dir, filename = h5_file)
   ```

2. **标准化和降维**
   - SCTransform标准化
   - PCA + UMAP

3. **无监督聚类**
   ```r
   spatial_data <- FindNeighbors(spatial_data, reduction = "pca")
   spatial_data <- FindClusters(spatial_data, resolution = 0.5)
   ```
   识别不同组织区域

4. **注释组织区域**
   根据聚类结果和已知标记基因，手动注释：
   - **Mal/Mal1**：恶性区域
   - **Normal**：正常组织
   - **Bdy**：肿瘤边界

5. **绘制空间特征图**
   ```r
   SpatialDimPlot(spatial_data, label = TRUE, label.size = 3)
   ```
   在组织切片图像上叠加细胞类型标签

**输出**：
- `Fig8A_Spatial_Feature_Map.pdf/png`

**说明**：原文使用Cottrazm包进行肿瘤边界界定，由于缺少完整参考数据，此处采用无监督聚类近似识别。

### （4）结果比较
（留空）

---

## 图8B（2）：高低MFRS组生存曲线

### （1）代码文件名
**03_fig8B_km_mfrs.R**

### （2）依赖的包
- `survival` - 生存分析
- `survminer` - KM曲线绘制
- `tidyverse` - 数据处理

### （3）方法分析

**目标**：展示高MFRS组的总生存率显著低于低MFRS组。

**步骤**：

1. **计算MFRS（风险评分）**
   ```r
   risk_genes <- c("SEMA3C", "RHOB", "PPP1R14A", "PCP4", "FHL2", 
                   "DES", "CUX1", "CALD1", "CA12", "ACTN1")
   multi_formula <- as.formula(paste("Surv(OS_time, OS_status) ~", 
                                    paste(risk_genes, collapse = " + ")))
   multi_cox <- coxph(multi_formula, data = analysis_data)
   risk_score <- predict(multi_cox, type = "risk")
   ```

2. **按风险评分分组**
   ```r
   cutoff <- median(risk_score)
   analysis_data$MFRS_group <- ifelse(risk_score > cutoff, 
                                      "High MFRS", "Low MFRS")
   ```

3. **拟合KM模型**
   ```r
   fit <- survfit(Surv(OS_time/365, OS_status) ~ MFRS_group, 
                  data = analysis_data)
   ```

4. **绘制KM曲线**
   ```r
   ggsurvplot(fit,
              pval = TRUE,              # 显示P值
              pval.method = TRUE,       # 显示检验方法
              conf.int = TRUE,          # 显示置信区间
              risk.table = TRUE,        # 风险表
              palette = c("#E64B35", "#4DBBD5"))  # 高=红，低=蓝
   ```

**关键结果**：P < 0.0001，说明高低MFRS组生存差异极显著。

**输出**：
- `Fig8B_KM_Curve_MFRS.pdf/png`

### （4）结果比较
（留空）

---

## 图8C（3）：10个预后基因表达雷达图

### （1）代码文件名
**03_fig8C_radar_plot.R**

### （2）依赖的包
- `fmsb` - 雷达图绘制
- `survival` - 风险评分
- `tidyverse` - 数据处理

### （3）方法分析

**目标**：用雷达图（蜘蛛网图）展示10个预后基因在高低MFRS组的表达模式差异。

**步骤**：

1. **计算风险评分并分组**（与8B相同）

2. **计算每组的平均表达量**
   ```r
   mean_expr_high <- analysis_data %>%
     filter(risk_group == "High MFRS") %>%
     select(all_of(risk_genes)) %>%
     summarise(across(everything(), mean))
   
   mean_expr_low <- analysis_data %>%
     filter(risk_group == "Low MFRS") %>%
     select(all_of(risk_genes)) %>%
     summarise(across(everything(), mean))
   ```

3. **准备fmsb格式数据**
   ```r
   radar_data <- rbind(
     max_row,      # 第一行：最大值
     min_row,      # 第二行：最小值
     mean_expr_high,  # 第三行：High MFRS
     mean_expr_low    # 第四行：Low MFRS
   )
   ```
   fmsb包要求前两行为最大值和最小值

4. **绘制雷达图**
   ```r
   radarchart(radar_data,
              pcol = c("#E64B35", "#4DBBD5"),  # 线条颜色
              pfcol = c(rgb(0.9, 0.3, 0.2, 0.3), 
                       rgb(0.3, 0.7, 0.8, 0.3)),  # 填充颜色（半透明）
              plwd = 2,
              plty = 1,
              cglcol = "grey",
              axislabcol = "black",
              title = "Expression Pattern of 10 Risk Genes")
   ```

5. **添加图例**
   ```r
   legend("topright", 
          legend = c("High MFRS", "Low MFRS"),
          col = c("#E64B35", "#4DBBD5"),
          lwd = 2)
   ```

**关键观察**：论文提到只有DES和PCP4在低MFRS组表达更高，其他基因在高MFRS组表达更高。

**输出**：
- `Fig8C_Radar_Plot.pdf/png`

### （4）结果比较
（留空）

---

## 图8D（4）：临床特征分布饼图

### （1）代码文件名
**03_fig8D_clinical_pie.R**

### （2）依赖的包
- `tidyverse` - 数据处理和绘图
- `patchwork` / `cowplot` - 图片组合
- `survival` - 风险评分

### （3）方法分析

**目标**：分析预后基因在不同临床状态（年龄、种族、分期、分级）中的分布。

**步骤**：

1. **计算风险评分并分组**

2. **准备8个临床变量**
   - Status（生存状态：Alive/Dead）
   - Age（年龄分组）
   - Race（种族：White/Black/Asian等）
   - Stage（FIGO分期：I/II/III/IV）
   - Grade（肿瘤分级：G1/G2/G3）
   - T（T分期）
   - M（M分期）
   - N（N分期）

3. **为每个变量绘制环形饼图**
   ```r
   create_pie <- function(data, var_name, title) {
     plot_data <- data %>%
       group_by(Risk_group, !!sym(var_name)) %>%
       summarise(count = n()) %>%
       mutate(percentage = count / sum(count) * 100)
     
     ggplot(plot_data, aes(x = "", y = percentage, 
                           fill = !!sym(var_name))) +
       geom_bar(stat = "identity", width = 1) +
       coord_polar("y", start = 0) +
       facet_wrap(~ Risk_group) +
       theme_void() +
       labs(title = title)
   }
   ```

4. **布局排列**
   - 2行8列：每列代表一个临床变量
   - 每列两个饼图：上=High MFRS，下=Low MFRS
   - 颜色：每列使用不同色系，深浅表示不同类别

**输出**：
- `Fig8D_Clinical_Pie.pdf/png`

### （4）结果比较
（留空）

---

## 图8E（5）：免疫药物IC50预测箱线图

### （1）代码文件名
**03_fig8E_drug_sensitivity.R**

### （2）依赖的包
- `pRRophetic` - 药物敏感性预测（需从GitHub安装）
- `survival` - 风险评分
- `tidyverse` - 数据处理
- `ggpubr` - 箱线图

### （3）方法分析

**目标**：预测高低MFRS组对4种免疫相关药物的敏感性（IC50值）。

**步骤**：

1. **计算风险评分并分组**

2. **准备表达矩阵**
   ```r
   expr_mat <- as.matrix(t(tpm_df))  # 转置：行=样本，列=基因
   ```

3. **使用pRRophetic预测IC50**
   ```r
   drugs <- c("Parthenolide", "Obatoclax Mesylate", "FTI.277", "AZD8055")
   
   ic50_results <- list()
   for (drug in drugs) {
     ic50 <- pRRopheticPredict(
       testMatrix = expr_mat,
       drug = drug,
       tissueType = "all",
       batchCorrect = "eb"
     )
     ic50_results[[drug]] <- ic50
   }
   ```
   pRRophetic基于GDSC数据库训练的模型预测药物反应

4. **合并结果并可视化**
   ```r
   plot_data <- data.frame(
     Sample = names(ic50),
     IC50 = ic50,
     Drug = drug,
     Risk_group = risk_group
   )
   
   ggplot(plot_data, aes(x = Risk_group, y = IC50, fill = Risk_group)) +
     geom_boxplot() +
     geom_jitter(width = 0.2, alpha = 0.5) +
     facet_wrap(~ Drug, scales = "free_y") +
     stat_compare_means() +
     scale_fill_manual(values = c("High MFRS" = "#E64B35", 
                                  "Low MFRS" = "#4DBBD5"))
   ```

**关键结果**：论文指出所有4种药物的IC50值在高MFRS组中较低（更敏感），低MFRS组较高（耐药）。

**输出**：
- `Fig8E_Drug_Sensitivity.pdf/png`

**注意**：如果未安装pRRophetic包，脚本会使用模拟数据展示绘图逻辑。

### （4）结果比较
（留空）

---

## 总结

Figure 8通过多模态数据验证和扩展了预后模型：

1. **空间验证**（8A）：空间转录组展示肿瘤边界特征
2. **生存验证**（8B）：独立数据集验证高低MFRS组生存差异
3. **表达模式**（8C）：雷达图直观展示10个基因的表达差异
4. **临床关联**（8D）：分析预后基因与临床特征的分布关系
5. **治疗意义**（8E）：预测药物敏感性，提示潜在治疗策略

特别地，图8A是唯一需要多个代码文件的子图，体现了空间转录组分析的复杂性。
