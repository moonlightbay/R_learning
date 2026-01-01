# Figure 7 分析文档：宫颈癌成纤维细胞预后模型

---

## 图7A（1）：单因素Cox回归分析森林图

### （1）代码文件名
**02_fig7A_univariate_cox.R**

### （2）依赖的包
- `survival` - Cox回归分析核心包
- `survminer` - 生存曲线可视化
- `tidyverse` - 数据处理（包含ggplot2）
- `data.table` - 快速数据读取

### （3）方法分析

**目标**：对13个候选基因进行单因素Cox回归分析，识别与预后显著相关的基因。

**步骤**：

1. **数据读取与合并**
   - 读取表达矩阵（Log2转换后的TPM数据）
   - 读取清洗后的临床数据
   - 转置表达矩阵并与临床数据按`case_submitter_id`合并
   - 过滤掉`OS_time`为NA或≤0的样本

2. **基因筛选策略**
   ```r
   custom_genes <- c("SMTN", "SEMA3C", "RHOB", "PPP1R14A", "PDLIM7", "PCP4",
                     "FHL2", "DSTN", "DES", "CUX1", "CALD1", "CA12", "ACTN1")
   ```
   使用论文中识别的13个C0 MYH11+ Fibroblasts相关基因。

3. **批量单因素Cox回归**
   ```r
   run_cox <- function(gene) {
     formula <- as.formula(paste0("Surv(OS_time, OS_status) ~ `", gene, "`"))
     cox_fit <- coxph(formula, data = analysis_data)
     # 提取HR、置信区间和P值
   }
   ```
   对每个基因单独构建模型：`Surv(time, status) ~ Gene`

4. **森林图绘制**
   - 使用ggplot2绘制
   - X轴：Hazard Ratio（HR）
   - 误差棒：95%置信区间（绿色）
   - 点：菱形，颜色映射P值（黑色→蓝色→绿色渐变）
   - 虚线：HR=1参考线，区分保护因子（HR<1）和风险因子（HR>1）
   - 基因按用户指定顺序逆序排列（从上到下）

**输出**：
- `Fig7A_ForestPlot.pdf/png`：森林图
- `univariate_cox_sig_genes.csv`：显著基因列表

### （4）结果比较
（留空）

---

## 图7B（2）：Lasso回归筛选基因

### （1）代码文件名
**02_fig7B_lasso.R**

### （2）依赖的包
- `glmnet` - Lasso回归核心包
- `survival` - 生存分析
- `tidyverse` - 数据处理
- `ggpubr` - 图片组合

### （3）方法分析

**目标**：使用Lasso-Cox回归从13个候选基因中进一步筛选，消除多重共线性。

**步骤**：

1. **准备Lasso输入矩阵**
   ```r
   x <- as.matrix(analysis_data[, exist_genes])  # 基因表达矩阵
   y <- Surv(analysis_data$OS_time, analysis_data$OS_status)  # 生存对象
   ```

2. **拟合Lasso模型**
   ```r
   fit <- glmnet(x, y, family = "cox")  # 系数路径
   cv_fit <- cv.glmnet(x, y, family = "cox", nfolds = 10)  # 交叉验证
   ```
   - `fit`：用于绘制系数路径图
   - `cv_fit`：通过10折交叉验证寻找最佳Lambda

3. **绘制两张子图**
   
   **子图A：Lasso系数路径**
   - X轴：Log Lambda
   - Y轴：系数值
   - 不同颜色线条代表不同基因的系数变化轨迹

   **子图B：交叉验证误差**
   - X轴：Log Lambda
   - Y轴：Partial Likelihood Deviance（偏离度）
   - 误差棒：交叉验证的标准误
   - 虚线：标注`lambda.min`（最小误差对应的Lambda）

4. **提取最终筛选的基因**
   ```r
   coef_min <- coef(cv_fit, s = "lambda.min")
   active_genes <- rownames(coef_min)[which(coef_min != 0)]
   ```
   使用`lambda.min`提取非零系数的基因。

**输出**：
- `Fig7B_Lasso.pdf/png`：组合图（包含A和B两个子图）
- `lasso_selected_genes.csv`：筛选出的基因及其系数

### （4）结果比较
（留空）

---

## 图7C（3）：多因素Cox回归分析

### （1）代码文件名
**02_fig7C_multivariate_cox.R.R**

### （2）依赖的包
- `survival` - Cox回归
- `survminer` - 生存分析可视化
- `tidyverse` - 数据处理

### （3）方法分析

**目标**：对最终选定的10个基因进行多因素Cox回归，评估在调整其他基因后每个基因的独立预后价值。

**步骤**：

1. **定义最终10个基因**
   ```r
   final_genes <- c("SEMA3C", "RHOB", "PPP1R14A", "PCP4", "FHL2", 
                    "DES", "CUX1", "CALD1", "CA12", "ACTN1")
   ```

2. **构建多因素Cox模型**
   ```r
   multi_formula <- as.formula(paste("Surv(OS_time, OS_status) ~", 
                                    paste(final_genes, collapse = " + ")))
   multi_cox <- coxph(multi_formula, data = analysis_data)
   ```
   公式形式：`Surv(time, status) ~ Gene1 + Gene2 + ... + Gene10`

3. **提取结果**
   - HR：风险比（exp(coef)）
   - 95%置信区间：lower和upper
   - P值：多因素模型中的显著性

4. **绘制森林图**
   - 风格与图7A一致
   - 标题改为"Multivariate Cox Regression Analysis"
   - 展示调整后的独立效应

**输出**：
- `Fig7C_Multivariate_ForestPlot.pdf/png`
- `multivariate_cox_results.csv`

### （4）结果比较
（留空）

---

## 图7D（4）：风险基因系数柱状图

### （1）代码文件名
**02_fig7D_7E_risk_score.R**

### （2）依赖的包
- `survival` - Cox回归
- `tidyverse` - 数据处理和绘图
- `ggpubr` - 图片组合

### （3）方法分析

**目标**：可视化10个风险基因在多因素模型中的系数大小。

**步骤**：

1. **构建多因素Cox模型**（与7C相同）

2. **提取系数和P值**
   ```r
   coef_data <- data.frame(
     Gene = rownames(multi_summary$coefficients),
     Coef = multi_summary$coefficients[, "coef"],  # log(HR)
     P_val = multi_summary$coefficients[, "Pr(>|z|)"]
   )
   ```

3. **计算-log10(P值)**
   ```r
   coef_data$logP <- -log10(coef_data$P_val)
   ```

4. **绘制柱状图**
   - X轴：基因名称
   - Y轴：系数值（Coefficient，即log Hazard Ratio）
   - 填充色：按-log10(P值)映射颜色（黑→蓝→绿渐变）
   - `scale_fill_gradientn`设置limits为c(0, 2.5)复刻论文

**输出**：
- `Fig7D_Coefficients.pdf/png`

### （4）结果比较
（留空）

---

## 图7E（5）：风险评分曲线与生存状态图

### （1）代码文件名
**02_fig7D_7E_risk_score.R**

### （2）依赖的包
- `survival` - 风险评分计算
- `tidyverse` - 数据处理
- `ggpubr` - 图片拼接

### （3）方法分析

**目标**：展示患者风险评分分布和生存状态，验证高风险组死亡率更高。

**步骤**：

1. **计算风险评分**
   ```r
   risk_score <- predict(multi_cox, type = "risk", newdata = analysis_data)
   ```
   `type="risk"`返回相对风险（Hazard Ratio）

2. **确定分组阈值**
   ```r
   n_survivors <- sum(analysis_data$OS_status == 0)  # 存活人数
   cutoff_score <- temp_sorted[n_survivors]  # 以存活人数为界
   ```
   特殊设计：使用存活患者数作为分界点，而非传统的中位数

3. **划分高低风险组**
   ```r
   analysis_data$risk_group <- ifelse(risk_score > cutoff_score, 
                                      "High MFRS", "Low MFRS")
   ```

4. **绘制两个子图**

   **子图1：风险评分曲线**
   - X轴：患者排名（按风险评分从小到大排序）
   - Y轴：风险评分
   - 区域填充：高风险区（红色）和低风险区（蓝色）
   - 竖虚线：标注存活人数位置

   **子图2：生存状态散点图**
   - X轴：患者排名
   - Y轴：生存时间（转换为年）
   - 点颜色：红色=死亡，绿色=存活
   - 竖虚线：与子图1对齐

**输出**：
- `Fig7E_RiskScore_Survival.pdf/png`（上下拼接）

### （4）结果比较
（留空）

---

## 图7F（6）：差异表达热图

### （1）代码文件名
**02_fig7F_7G_heatmap_roc.R**

### （2）依赖的包
- `pheatmap` - 热图绘制
- `survival` - 风险评分
- `tidyverse` - 数据处理

### （3）方法分析

**目标**：展示10个风险基因在高低风险组间的表达差异模式。

**步骤**：

1. **计算风险评分并分组**（与7E相同逻辑）

2. **准备热图矩阵**
   ```r
   heatmap_data <- analysis_data %>% arrange(risk_score)
   plot_mat <- t(heatmap_data[, final_genes])  # 转置：行=基因，列=样本
   ```

3. **准备注释条**
   ```r
   annotation_col <- data.frame(
     Group = factor(heatmap_data$risk_group, levels = c("Low MFRS", "High MFRS"))
   )
   ann_colors <- list(Group = c("Low MFRS" = "blue", "High MFRS" = "#E64B35"))
   ```

4. **绘制热图**
   ```r
   pheatmap(plot_mat,
            scale = "row",           # 按行标准化（Z-score）
            cluster_rows = TRUE,     # 聚类基因
            cluster_cols = FALSE,    # 不聚类样本（保持风险评分排序）
            show_colnames = FALSE,   # 隐藏样本名
            annotation_col = annotation_col,
            color = colorRampPalette(c("blue", "white", "#E64B35"))(100))
   ```

**输出**：
- `Fig7F_Heatmap.pdf/png`

### （4）结果比较
（留空）

---

## 图7G（7）：时间依赖性ROC曲线

### （1）代码文件名
**02_fig7F_7G_heatmap_roc.R**

### （2）依赖的包
- `timeROC` - 时间依赖性ROC分析
- `survival` - 生存分析
- `tidyverse` - 数据处理

### （3）方法分析

**目标**：评估预后模型在1年、3年、5年生存预测的准确性。

**步骤**：

1. **计算时间依赖性ROC**
   ```r
   times_of_interest <- c(365, 365*3, 365*5)  # 1/3/5年（天数）
   roc_res <- timeROC(
     T = analysis_data$OS_time,
     delta = analysis_data$OS_status,
     marker = analysis_data$risk_score,
     cause = 1,  # 1=死亡
     times = times_of_interest,
     iid = TRUE  # 计算置信区间
   )
   ```

2. **提取TPR和FPR**
   - TPR：真阳性率（灵敏度）
   - FPR：假阳性率（1-特异度）
   - 对每个时间点提取ROC曲线坐标

3. **提取AUC值**
   ```r
   auc_1y <- round(roc_res$AUC[1], 2)
   auc_3y <- round(roc_res$AUC[2], 2)
   auc_5y <- round(roc_res$AUC[3], 2)
   ```

4. **绘制ROC曲线**
   - 三条曲线：1年（红色）、3年（蓝色）、5年（绿色）
   - 对角虚线：随机分类器参考线
   - 图例：显示AUC值

**输出**：
- `Fig7G_ROC.pdf/png`

### （4）结果比较
（留空）

---

## 图7H（8）：PCA散点图

### （1）代码文件名
**02_fig7H_pca.R**

### （2）依赖的包
- `survival` - 风险评分
- `tidyverse` - 数据处理和绘图

### （3）方法分析

**目标**：通过主成分分析（PCA）展示高低风险组在全基因组表达谱上的差异。

**步骤**：

1. **准备PCA输入数据**
   ```r
   gene_cols <- setdiff(colnames(analysis_data), non_gene_cols)
   pca_input <- analysis_data[, gene_cols]  # 全基因组表达矩阵
   ```

2. **运行PCA**
   ```r
   pca_res <- prcomp(pca_input, scale. = FALSE)
   ```
   `scale.=FALSE`：不标准化，保持论文中的大坐标范围

3. **提取PC坐标**
   ```r
   pca_data <- as.data.frame(pca_res$x)
   pca_data$risk_group <- analysis_data$risk_group
   ```

4. **绘制散点图**
   - X轴：PC1（标注方差解释比例）
   - Y轴：PC2（标注方差解释比例）
   - 点颜色：高风险（红色）vs 低风险（蓝色）
   - 椭圆：95%置信椭圆，展示组间分离程度

**输出**：
- `Fig7H_PCA.pdf/png`

### （4）结果比较
（留空）

---

## 图7I（9）：SDC1基因生存曲线

### （1）代码文件名
**02_fig7I_km.R**

### （2）依赖的包
- `survival` - 生存分析
- `survminer` - Kaplan-Meier曲线绘制
- `tidyverse` - 数据处理

### （3）方法分析

**目标**：展示关键基因SDC1高低表达组的生存差异。

**步骤**：

1. **根据SDC1表达量分组**
   ```r
   cutoff_sdc1 <- median(analysis_data$SDC1)
   analysis_data$SDC1_group <- ifelse(analysis_data$SDC1 > cutoff_sdc1, 
                                      "High SDC1", "Low SDC1")
   ```
   使用中位数作为分组阈值

2. **拟合KM模型**
   ```r
   fit <- survfit(Surv(OS_time/365, OS_status) ~ SDC1_group, data = analysis_data)
   ```
   生存时间转换为年

3. **绘制KM曲线**
   ```r
   ggsurvplot(fit,
              pval = TRUE,              # 显示P值（Log-rank检验）
              risk.table = TRUE,        # 显示风险表
              surv.median.line = "hv",  # 标注中位生存时间
              palette = c("#E64B35", "#4DBBD5"))  # 高=红，低=蓝
   ```

**输出**：
- `Fig7I_KM_Curve_SDC1.pdf/png`

### （4）结果比较
（留空）

---

## 图7J（10）：临床因素多因素Cox回归

### （1）代码文件名
**02_fig7J_multivariate_clinical.R**

### （2）依赖的包
- `survival` - Cox回归
- `tidyverse` - 数据处理

### （3）方法分析

**目标**：证明风险评分（MFRS）是独立于临床因素的预后因子。

**步骤**：

1. **计算风险评分并分组**

2. **处理临床变量**
   - 年龄：按中位数分为High Age和Low Age
   - 种族：因子化（如White, Black, Asian等）
   - 分期：使用清洗后的FIGO Stage（I, II, III, IV）
   - T/M/N分期：使用清洗后的TNM分期
   - Grade：肿瘤分级

3. **构建多因素Cox模型**
   ```r
   formula <- Surv(OS_time, OS_status) ~ Risk_group + Age_group + 
                                        Race + Stage + Grade + T + M + N
   multi_cox <- coxph(formula, data = analysis_data)
   ```
   包含风险组和所有可用临床变量

4. **绘制森林图**
   - 风格与图7A/7C一致
   - 重点关注Risk_group的HR和P值
   - 如果Risk_group显著（P<0.05），说明MFRS是独立预后因子

**输出**：
- `Fig7J_Multivariate_Clinical_ForestPlot.pdf/png`

### （4）结果比较
（留空）

---

## 图7K（11）：列线图（Nomogram）

### （1）代码文件名
**02_fig7K_nomogram.R**

### （2）依赖的包
- `rms` - 列线图绘制
- `survival` - 生存分析
- `tidyverse` - 数据处理

### （3）方法分析

**目标**：构建基于MFRS和临床因素的预测列线图，用于临床应用。

**步骤**：

1. **数据准备**（与7J相同）

2. **使用rms包构建Cox模型**
   ```r
   dd <- datadist(nomogram_data)
   options(datadist = 'dd')
   
   cox_nom <- cph(Surv(OS_time, OS_status) ~ risk_score + Age_group + 
                                              Race + Grade, 
                  data = nomogram_data, x = TRUE, y = TRUE)
   ```
   `cph`是rms包的Cox模型函数，支持列线图

3. **绘制列线图**
   ```r
   nom <- nomogram(cox_nom,
                   fun = list(function(x) survest(cox_nom, times=365, lp=x)$surv,
                             function(x) survest(cox_nom, times=365*3, lp=x)$surv,
                             function(x) survest(cox_nom, times=365*5, lp=x)$surv),
                   fun.at = c(0.9, 0.8, 0.7, 0.5, 0.3),
                   funlabel = c("1-Year Survival", "3-Year Survival", "5-Year Survival"))
   ```

4. **解读列线图**
   - 顶部：Points刻度尺
   - 中间：各变量刻度（风险评分、年龄、种族、分级）
   - 底部：1/3/5年生存概率预测

**输出**：
- `Fig7K_Nomogram.pdf/png`

### （4）结果比较
（留空）

---

## 图7L（12）：相关性矩阵图

### （1）代码文件名
**02_fig7L_correlation.R**

### （2）依赖的包
- `GGally` - 相关性矩阵绘图
- `survival` - 风险评分
- `tidyverse` - 数据处理

### （3）方法分析

**目标**：展示10个预后基因、OS（生存时间）和Risk（风险评分）之间的相关性。

**步骤**：

1. **准备绘图数据**
   ```r
   plot_data <- analysis_data %>%
     select(OS_time, all_of(final_genes), risk_score) %>%
     rename(OS = OS_time, Risk = risk_score)
   ```
   共12个变量：OS + 10个基因 + Risk

2. **定义自定义绘图函数**
   - **右上三角**：热图+数值（相关系数）
     - 颜色：蓝-白-红渐变
     - 文本：显示相关系数值
   - **左下三角**：散点图
     - 展示两两变量的线性关系
   - **对角线**：密度图或变量名

3. **使用GGally::ggpairs绘制**
   ```r
   ggpairs(plot_data,
           upper = list(continuous = upper_fn),  # 自定义右上
           lower = list(continuous = "points"),   # 左下散点图
           diag = list(continuous = "densityDiag"))  # 对角线密度图
   ```

**输出**：
- `Fig7L_Correlation.pdf/png`

### （4）结果比较
（留空）

---

## 图7M（13）：基因表达差异箱线图

### （1）代码文件名
**02_fig7M_expression_diff.R**

### （2）依赖的包
- `survival` - 风险评分
- `tidyverse` - 数据处理
- `ggpubr` - 箱线图和统计检验

### （3）方法分析

**目标**：展示8个风险基因在高低MFRS组间的表达差异。

**步骤**：

1. **计算风险评分并分组**

2. **选择目标8个基因**
   ```r
   target_genes <- c("ACTN1", "FHL2", "CA12", "PPP1R14A", 
                     "CALD1", "RHOB", "CUX1", "SEMA3C")
   ```

3. **准备长格式数据**
   ```r
   plot_data_long <- analysis_data %>%
     select(Risk_group, all_of(target_genes)) %>%
     pivot_longer(cols = all_of(target_genes), 
                  names_to = "Gene", values_to = "Expression")
   ```

4. **绘制箱线图组合**
   - 布局：2行4列
   - 每个基因一个子图
   - X轴：High MFRS（红）vs Low MFRS（蓝）
   - Y轴：基因表达量（Log2 TPM）
   - 添加散点：`geom_jitter`
   - 显著性检验：`stat_compare_means`（Wilcoxon检验）

**输出**：
- `Fig7M_Expression_Diff.pdf/png`

### （4）结果比较
（留空）

---

## 总结

Figure 7完整地展示了基于TCGA宫颈癌数据构建预后模型的全流程：

1. **基因筛选**（7A-7C）：单因素Cox → Lasso回归 → 多因素Cox
2. **模型展示**（7D-7E）：系数可视化 → 风险评分分布
3. **模型验证**（7F-7H）：热图 → ROC曲线 → PCA
4. **独立验证**（7I-7J）：关键基因生存曲线 → 临床因素校正
5. **临床应用**（7K）：列线图预测工具
6. **机制探索**（7L-7M）：基因相关性 → 表达差异

所有子图共同构成了一个完整的、可重现的预后模型研究。
