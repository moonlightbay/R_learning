# Figure 9 分析文档：高低MFRS组间差异基因富集分析

---

## 图9A（1）：差异基因热图

### （1）代码文件名
**04_fig9A_heatmap.R**

### （2）依赖的包
- `DESeq2` - 差异表达分析（Bioconductor）
- `pheatmap` - 热图绘制
- `survival` - 风险评分
- `tidyverse` - 数据处理

### （3）方法分析

**目标**：识别高低MFRS组间的差异表达基因（DEGs），并绘制热图展示表达谱。

**步骤**：

1. **混合数据策略**
   - **TPM矩阵**：用于计算风险评分和样本分组
   - **Raw Counts矩阵**：用于DESeq2差异分析

   ```r
   tpm_df <- fread("processed_data/expression_matrix.csv")  # Log2(TPM+1)
   counts_df <- fread("processed_data/expression_matrix_counts.csv")  # Raw counts
   ```

2. **计算风险评分并分组**
   ```r
   multi_cox <- coxph(Surv(OS_time, OS_status) ~ Gene1 + ... + Gene10)
   risk_score <- predict(multi_cox, type = "risk")
   risk_group <- ifelse(risk_score > median(risk_score), "High", "Low")
   ```

3. **构建DESeq2数据对象**
   ```r
   dds <- DESeqDataSetFromMatrix(
     countData = counts_mat,  # 行=基因，列=样本
     colData = col_data,      # 样本信息（含risk_group）
     design = ~ risk_group
   )
   ```

4. **运行DESeq2分析**
   ```r
   dds <- DESeq(dds)
   res <- results(dds, contrast = c("risk_group", "High", "Low"))
   ```
   计算High vs Low的差异表达

5. **筛选显著差异基因**
   ```r
   sig_genes <- res %>%
     filter(padj < 0.05 & abs(log2FoldChange) > 1)
   ```
   标准：校正后P值<0.05 且 |Log2FC| > 1

6. **选取Top基因绘制热图**
   ```r
   top_genes <- sig_genes %>%
     arrange(padj) %>%
     head(30)  # 选取最显著的30个基因
   
   plot_mat <- assay(vst(dds))[top_genes$gene, ]  # VST标准化
   
   pheatmap(plot_mat,
            scale = "row",
            cluster_rows = TRUE,
            cluster_cols = TRUE,
            annotation_col = annotation_col,  # 标注High/Low组
            color = colorRampPalette(c("blue", "white", "red"))(100))
   ```

**输出**：
- `Fig9A_Heatmap.pdf/png`
- `DEG_results.csv`：所有差异基因列表

### （4）结果比较
（留空）

---

## 图9B（2）：差异基因火山图

### （1）代码文件名
**04_fig9B_volcano.R**

### （2）依赖的包
- `DESeq2` - 差异表达分析
- `ggplot2` - 火山图绘制
- `ggrepel` - 基因标签防重叠
- `tidyverse` - 数据处理

### （3）方法分析

**目标**：用火山图展示所有基因的差异表达情况，区分上调和下调基因。

**步骤**：

1. **运行DESeq2分析**（与9A相同流程）

2. **准备火山图数据**
   ```r
   volcano_data <- as.data.frame(res) %>%
     mutate(
       log2FC = log2FoldChange,
       negLogP = -log10(pvalue),
       Significant = case_when(
         padj < 0.05 & log2FC > 1 ~ "Up",
         padj < 0.05 & log2FC < -1 ~ "Down",
         TRUE ~ "NS"
       )
     )
   ```

3. **标记关键基因**
   ```r
   top_up <- volcano_data %>% filter(Significant == "Up") %>% 
             arrange(padj) %>% head(10)
   top_down <- volcano_data %>% filter(Significant == "Down") %>% 
               arrange(padj) %>% head(10)
   label_genes <- rbind(top_up, top_down)
   ```

4. **绘制火山图**
   ```r
   ggplot(volcano_data, aes(x = log2FC, y = negLogP)) +
     geom_point(aes(color = Significant), alpha = 0.6, size = 1.5) +
     scale_color_manual(values = c("Up" = "red", "Down" = "blue", "NS" = "grey")) +
     geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey30") +
     geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey30") +
     geom_text_repel(data = label_genes, aes(label = gene), 
                     size = 3, max.overlaps = 20) +
     labs(x = "Log2 Fold Change (High vs Low)",
          y = "-Log10(P value)",
          title = "Volcano Plot: DEGs between High and Low MFRS")
   ```

**输出**：
- `Fig9B_Volcano.pdf/png`

### （4）结果比较
（留空）

---

## 图9C/D/E（3/4/5）：GO和KEGG富集分析

### （1）代码文件名
**04_fig9C_E_enrichment.R**

### （2）依赖的包
- `DESeq2` - 差异分析
- `clusterProfiler` - 富集分析核心包（Bioconductor）
- `org.Hs.eg.db` - 人类基因注释数据库
- `enrichplot` - 富集结果可视化
- `tidyverse` - 数据处理

### （3）方法分析

**目标**：对差异基因进行功能富集分析，揭示高低MFRS组的生物学差异。

**步骤**：

1. **获取差异基因列表**（使用9A/9B的DESeq2结果）

2. **基因ID转换**
   ```r
   gene_entrez <- bitr(sig_genes$gene,
                       fromType = "SYMBOL",
                       toType = "ENTREZID",
                       OrgDb = org.Hs.eg.db)
   ```
   将基因Symbol转换为Entrez ID（富集分析所需）

3. **GO-BP富集分析**（图9C）
   ```r
   ego_bp <- enrichGO(
     gene = gene_entrez$ENTREZID,
     OrgDb = org.Hs.eg.db,
     ont = "BP",  # Biological Process
     pAdjustMethod = "BH",
     pvalueCutoff = 0.05,
     qvalueCutoff = 0.05
   )
   
   barplot(ego_bp, showCategory = 15, title = "GO-BP Enrichment")
   ```
   
   **论文结果**：主要富集于消化过程、中间丝组织等

4. **GO-CC富集分析**（图9D）
   ```r
   ego_cc <- enrichGO(
     gene = gene_entrez$ENTREZID,
     OrgDb = org.Hs.eg.db,
     ont = "CC",  # Cellular Component
     pAdjustMethod = "BH",
     pvalueCutoff = 0.05
   )
   
   barplot(ego_cc, showCategory = 15, title = "GO-CC Enrichment")
   ```
   
   **论文结果**：质膜外泌成分、刷状缘、肌动蛋白细胞突触聚集体等

5. **KEGG通路富集**（图9E）
   ```r
   kegg <- enrichKEGG(
     gene = gene_entrez$ENTREZID,
     organism = "hsa",  # Homo sapiens
     pAdjustMethod = "BH",
     pvalueCutoff = 0.05
   )
   
   barplot(kegg, showCategory = 15, title = "KEGG Pathway Enrichment")
   ```
   
   **论文结果**：
   - 年轻型糖尿病发病机制
   - Ras信号传导
   - 金黄色葡萄球菌感染
   - 碳水化合物消化和吸收
   - 癌症中的转录失调

**绘图参数**：
- X轴：GeneRatio（富集基因数/通路总基因数）
- Y轴：通路名称
- 颜色：P值
- 点大小：Count（富集基因数）

**输出**：
- `Fig9C_GO_BP.pdf/png`
- `Fig9D_GO_CC.pdf/png`
- `Fig9E_KEGG.pdf/png`

### （4）结果比较
（留空）

---

## 图9F（6）：GSEA富集分析

### （1）代码文件名
**04_fig9F_gsea.R**

### （2）依赖的包
- `DESeq2` - 差异分析
- `clusterProfiler` - GSEA分析
- `org.Hs.eg.db` - 基因注释
- `enrichplot` - 结果可视化

### （3）方法分析

**目标**：使用GSEA（Gene Set Enrichment Analysis）识别在高低MFRS组间系统性差异的基因集。

**步骤**：

1. **构建基因排序列表**
   ```r
   # 使用所有基因（不仅是显著差异基因）
   gene_list <- res$log2FoldChange
   names(gene_list) <- res$gene
   gene_list <- sort(gene_list, decreasing = TRUE)  # 按Log2FC降序排列
   ```

2. **基因ID转换**
   ```r
   gene_list_entrez <- bitr(names(gene_list),
                            fromType = "SYMBOL",
                            toType = "ENTREZID",
                            OrgDb = org.Hs.eg.db)
   # 更新gene_list的名称为Entrez ID
   ```

3. **运行GSEA（GO-BP）**
   ```r
   gsea_go <- gseGO(
     geneList = gene_list,
     OrgDb = org.Hs.eg.db,
     ont = "BP",
     pvalueCutoff = 0.05,
     pAdjustMethod = "BH"
   )
   ```

4. **提取Top通路**
   - **上调通路**（NES > 0）：
     - 胶原分解代谢过程
     - 伤口愈合
     - 表皮细胞扩散
     - 内胚层形成
     - 细胞外基质分解
   
   - **下调通路**（NES < 0）：
     - 钾离子稳态
     - 杀死其他生物体的细胞
     - 轴丝组装

5. **绘制GSEA图**
   ```r
   # 选取Top 5 上调和Top 5 下调通路
   top_pathways <- rbind(
     head(gsea_go[gsea_go$NES > 0, ], 5),
     head(gsea_go[gsea_go$NES < 0, ], 5)
   )
   
   gseaplot2(gsea_go, geneSetID = top_pathways$ID, 
             title = "GSEA Enrichment Plot")
   ```

**输出**：
- `Fig9F_GSEA.pdf/png`
- GSEA详细结果表

### （4）结果比较
（留空）

---

## 图9G（7）：GSVA富集分析热图

### （1）代码文件名
**04_fig9G_gsva.R**

### （2）依赖的包
- `GSVA` - 基因集变异分析（Bioconductor）
- `msigdbr` - MSigDB基因集数据库
- `pheatmap` - 热图绘制
- `survival` - 风险评分

### （3）方法分析

**目标**：计算每个样本在多个基因集上的富集分数（GSVA），并与风险评分相关性分析。

**步骤**：

1. **准备表达矩阵**
   ```r
   expr_mat <- as.matrix(tpm_df)  # Log2(TPM+1)，行=基因，列=样本
   ```

2. **获取MSigDB基因集**
   ```r
   gene_sets <- msigdbr(species = "Homo sapiens", category = "C2")  # Curated
   gene_sets_list <- split(gene_sets$gene_symbol, gene_sets$gs_name)
   ```

3. **运行GSVA**
   ```r
   gsva_res <- gsva(
     expr = expr_mat,
     gset.idx.list = gene_sets_list,
     method = "gsva",
     parallel.sz = 4  # 并行计算
   )
   ```
   输出矩阵：行=基因集，列=样本，值=富集分数

4. **计算与Risk Score的相关性**
   ```r
   cor_results <- apply(gsva_res, 1, function(pathway_score) {
     cor.test(pathway_score, risk_score, method = "pearson")
   })
   
   # 筛选显著正相关通路（cor > 0.3, P < 0.05）
   sig_pathways <- cor_results %>%
     filter(cor > 0.3 & pvalue < 0.05) %>%
     arrange(desc(cor)) %>%
     head(20)
   ```

5. **绘制热图**
   ```r
   plot_mat <- gsva_res[sig_pathways$pathway, ]
   
   # 样本按Risk Score排序
   sample_order <- order(risk_score)
   plot_mat <- plot_mat[, sample_order]
   
   pheatmap(plot_mat,
            scale = "row",
            cluster_rows = TRUE,
            cluster_cols = FALSE,
            annotation_col = data.frame(Risk_Score = risk_score[sample_order]),
            color = colorRampPalette(c("blue", "white", "red"))(100),
            show_colnames = FALSE)
   ```

**论文结果**：以下通路与风险评分显著正相关：
- `PEDERSEN_TARGETS_OF_611_CTF_ISOFORM_OF_ERBB2`
- `REACTOME_CELL_EXTRACELLULAR_MATRIX_INTERACTIONS`
- `BRUECKNER_TARGETS_OF_MIRLET7A3_DN`
- `ELVIDGE_HYPOXIA_UP`

**输出**：
- `Fig9G_GSVA_Heatmap.pdf/png`

### （4）结果比较
（留空）

---

## 图9H（8）：差异富集通路雷达图

### （1）代码文件名
**04_fig9H_I_radar.R**

### （2）依赖的包
- `GSVA` - 基因集分析
- `fmsb` - 雷达图绘制
- `msigdbr` - 基因集数据库
- `tidyverse` - 数据处理

### （3）方法分析

**目标**：用雷达图展示高低MFRS组间显著富集的Top 20通路差异。

**步骤**：

1. **运行GSVA分析**（与9G相同流程）

2. **识别组间差异显著的通路**
   ```r
   pathway_diff <- apply(gsva_res, 1, function(pathway_score) {
     t_test <- t.test(pathway_score[high_samples], 
                     pathway_score[low_samples])
     data.frame(
       mean_high = mean(pathway_score[high_samples]),
       mean_low = mean(pathway_score[low_samples]),
       pvalue = t_test$p.value
     )
   })
   
   # 选取Top 20显著通路
   top_pathways <- pathway_diff %>%
     filter(pvalue < 0.05) %>%
     mutate(diff = abs(mean_high - mean_low)) %>%
     arrange(desc(diff)) %>%
     head(20)
   ```

3. **准备雷达图数据**
   ```r
   radar_data <- rbind(
     max_values,  # 第一行：最大值
     min_values,  # 第二行：最小值
     high_scores, # 第三行：High MFRS平均分数
     low_scores   # 第四行：Low MFRS平均分数
   )
   colnames(radar_data) <- top_pathways$pathway
   ```

4. **绘制雷达图**
   ```r
   radarchart(radar_data,
              pcol = c("#E64B35", "#4DBBD5"),
              pfcol = c(rgb(0.9, 0.3, 0.2, 0.3), rgb(0.3, 0.7, 0.8, 0.3)),
              plwd = 2,
              cglcol = "grey",
              title = "Top 20 Differential Pathways")
   ```

**论文结果**：高MFRS组显著富集的通路包括：
- `GOBP_BRANCHED_CHAIN_AMINO_ACID_TRANSPORT`（支链氨基酸转运）
- `GOBP_REGULATION_OF_VENTRICULAR_CARDIAC_MUSCLE_CELL_ACTION_POTENTIAL`
- `GOBP_NEGATIVE_REGULATION_OF_CYCLIC_NUCLEOTIDE_PHOSPHODIESTERASE_ACTIVITY`
- `GOBP_PROTEIN_LOCALIZATION_TO_PLASMA_MEMBRANE_RAFT`
- `GOBP_INHIBITORY_SYNAPSE_ASSEMBLY`
- `GOBP_NEGATIVE_REGULATION_OF_GLUTAMATE_SECRETION`

**输出**：
- `Fig9H_Pathway_Radar.pdf/png`

### （4）结果比较
（留空）

---

## 图9I（9）：免疫细胞浸润雷达图

### （1）代码文件名
**04_fig9H_I_radar.R**

### （2）依赖的包
- `GSVA` - 免疫细胞评分
- `fmsb` - 雷达图
- `msigdbr` / 自定义基因集 - 免疫细胞标记基因

### （3）方法分析

**目标**：比较高低MFRS组的免疫细胞浸润景观。

**步骤**：

1. **定义免疫细胞标记基因集**
   ```r
   # 使用预定义的22种免疫细胞类型标记基因（LM22或类似）
   immune_gene_sets <- list(
     "Activated_DC" = c("gene1", "gene2", ...),
     "M0_Macrophages" = c(...),
     "Activated_NK_cells" = c(...),
     ...
   )
   ```

2. **使用GSVA计算免疫细胞评分**
   ```r
   immune_scores <- gsva(
     expr = expr_mat,
     gset.idx.list = immune_gene_sets,
     method = "gsva"
   )
   ```

3. **计算每组的平均免疫细胞评分**
   ```r
   mean_scores_high <- rowMeans(immune_scores[, high_samples])
   mean_scores_low <- rowMeans(immune_scores[, low_samples])
   ```

4. **选择Top 20差异免疫细胞**
   ```r
   immune_diff <- abs(mean_scores_high - mean_scores_low)
   top_immune <- names(sort(immune_diff, decreasing = TRUE))[1:20]
   ```

5. **绘制雷达图**
   ```r
   radar_data <- rbind(
     max_values,
     min_values,
     mean_scores_high[top_immune],
     mean_scores_low[top_immune]
   )
   
   radarchart(radar_data, 
              pcol = c("#E64B35", "#4DBBD5"),
              title = "Immune Cell Infiltration Landscape")
   ```

**论文结果**：高MFRS组中以下免疫细胞更为丰富：
- 活化的树突状细胞（Activated dendritic cells）
- M0巨噬细胞（M0 macrophages）
- 活化的自然杀伤细胞（Activated NK cells）
- 静息自然杀伤细胞（Resting NK cells）

**输出**：
- `Fig9I_Immune_Radar.pdf/png`

### （4）结果比较
（留空）

---

## 总结

Figure 9通过多层次富集分析全面揭示了高低MFRS组的生物学差异：

1. **差异基因识别**（9A-9B）：热图+火山图展示DEGs
2. **功能注释**（9C-9E）：GO-BP/CC + KEGG揭示功能差异
3. **通路层面**（9F）：GSEA识别系统性变化的生物学过程
4. **样本层面**（9G）：GSVA计算每个样本的通路活性
5. **综合比较**（9H-9I）：雷达图直观展示通路和免疫细胞景观差异

这些分析共同提示高MFRS组与细胞外基质重塑、免疫激活和代谢重编程相关，为理解预后差异的分子机制提供了重要线索。
