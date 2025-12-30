# Fig 7 复现计划：宫颈癌成纤维细胞预后模型 (Prognostic Model)

## 1. 目标描述
**核心任务**：基于 TCGA 宫颈癌 (CC) 批量转录组数据，构建预后模型，并复现论文 **Figure 7** 中的关键图表。
**复现目标 (建议选 3-5 张)**：
1.  **Fig 7E**: 高低风险组的生存曲线 (Kaplan-Meier Curve) —— *最核心，必须复现*
2.  **Fig 7G/H**: 模型的预测准确性 (ROC Curve, 1/3/5年) —— *验证模型好坏的标准图*
3.  **Fig 7M**: 风险基因在不同分组的表达差异 (Boxplot) —— *展示基因层面差异*
4.  **Fig 7J**: 多因素 Cox 回归森林图 (Multivariate Cox) —— *证明独立预后价值*
5.  **Fig 7A/B**: Lasso 回归筛选变量的过程图 —— *如果有余力可做*

## 2. 数据集解析
所有数据位于 `data/fig7-9/prognostic/` 目录下。

| 文件/文件夹 | 用途 | 关键操作 |
|---|---|---|
| `data/` (文件夹) | 存放原始表达量数据。包含大量 UUID 命名的子文件夹，每个文件夹内有一个 `.tsv` 文件。 | 需要遍历读取所有 `.tsv`，提取 `unstranded` (未链特异性) 或 `tpm_unstranded` 列，合并为一个大的表达矩阵。 |
| `gdc_sample_sheet.tsv` | 样本映射表。连接 **文件名** (File Name) 和 **样本ID** (Sample ID)。 | 用于将 `data/` 里的乱码文件名转换为可读的 TCGA 样本 ID。 |
| `clinical.tsv` | 临床信息表。包含生存时间 (days_to_death/last_follow_up)、生存状态 (vital_status) 等。 | 需要清洗生存时间，将 "Alive/Dead" 转换为 0/1 状态，并与表达矩阵对齐。 |

## 3. 复现方法论 (Methodology)

### 步骤一：数据清洗与合并 (Data Cleaning)
由于原始数据分散在数百个文件夹中，首先需要构建一个**整合矩阵**。
*   **输入**: `data/` 下的所有 `.tsv` 文件, `gdc_sample_sheet.tsv`, `clinical.tsv`
*   **输出**: 
    1.  `expression_matrix.csv` (行=基因, 列=样本)
    2.  `survival_data.csv` (样本ID, 生存时间, 生存状态, 临床协变量)

### 步骤二：模型构建 (Model Construction)
论文使用了 **Lasso-Cox 回归** 方法筛选基因并构建模型。
1.  **差异分析/预筛选**: (可选) 先筛选出在 Tumor vs Normal 中差异表达的基因，或者直接使用文中提到的 "C0 subtype marker genes" (如果文中未提供具体列表，则直接对全基因组或高变基因进行 Lasso)。
2.  **单因素 Cox**: 找出与生存显著相关的基因。
3.  **Lasso 回归**: 消除多重共线性，筛选出最关键的 ~10 个基因 (Fig 7A-D)。
4.  **计算风险评分 (Risk Score)**: $Score = \sum (Expression \times Coefficient)$。
5.  **分组**: 根据 Risk Score 的中位数或最佳截断值 (cutoff)，将病人分为 High Risk 和 Low Risk 组。

### 步骤三：可视化 (Visualization)
使用 R 语言包进行绘图：
*   **生存分析**: `survival`, `survminer` (用于 Fig 7E)
*   **ROC曲线**: `timeROC` 或 `survivalROC` (用于 Fig 7G/H)
*   **差异比较**: `ggplot2`, `ggpubr` (用于 Fig 7M)
*   **森林图**: `survival`, `survminer` (用于 Fig 7J)

## 4. 关键 R 包
```r
library(tidyverse)  # 数据处理
library(survival)   # 生存分析核心
library(survminer)  # 生存曲线可视化
library(glmnet)     # Lasso 回归
library(timeROC)    # 时间依赖性 ROC
```

## 5. 论文原文摘录 (Reference)
> "We adopted a strategy of utilizing marker genes of key fibroblast subtypes as predictor genes."
> "Univariate Cox and Lasso regression analyses were performed using the 'survival' R package."
> "Risk score = (Expressed Gene1 × Coefficient1) + ... + (Expressed GeneN × CoefficientN)"
> "Key prognostic significance of C0 MYH11 +fibroblasts in cervical cancer patients. A A forest plot showed the results of univariate Cox analysis, identifying 13 C0 MYH11 +Fibroblasts-related genes associated with cervical cancer prognosis, with P≤ 0.05. The reference line (HR =1) distinguished protective factors (HR <1) from risk factors (HR >1). B Lasso regression was applied to select 10 genes contributing to the risk score, and the results were displayed in a lambda plot (lambda.min =0.017). C A forest plot depicted the final 10 genes associated with cervical cancer prognosis. D A bar plot showed the coefficient values of the 10 risk genes. E A curve plot illustrated the risk scores in high MFRS and low MFRS groups, with scatter plots displaying survival/death events over time for the high MFRS and low MFRS groups. F A heatmap displayed the differential expression of the 10 risk genes between the high MFRS and low MFRS groups. G ROC curves showed the AUC values for 1-year (AUC =0.80), 3-year (AUC =0.89), and 5-year (AUC =0.85) survival predictions. H A scatter plot displayed the distribution of genes along PC1 and PC2 for the high MFRS and low MFRS groups. I Kaplan–Meier survival curves showed survival differences between the high and low SDC1 expression groups. J A forest plot presented the results of multivariate Cox analysis for clinical factors and risk scores in the training cohort. K The nomogram model based on the MFRS was constructed, including race, age, and grade. L Scatter and heatmaps showed the correlation between the 10 prognostic genes and OS and MFRS. M Box and scatter plots illustrated the expression differences of 8 genes from the 10 risk genes between the high MFRS and low MFRS groups"