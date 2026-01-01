# 代码分析文档目录

本目录包含project目录下所有代码文件的详细分析文档。

## 文档列表

### 1. 数据准备分析_01_data_prep.md
**对应代码**: `01_data_prep.py`

数据准备是整个分析流程的第一步，负责：
- 读取TCGA宫颈癌原始数据
- 清洗和标准化临床信息
- 整合表达矩阵（TPM和Raw Counts）
- 输出后续分析所需的标准化数据

**关键输出**：
- `clinical_cleaned.csv`
- `expression_matrix.csv` (Log2转换的TPM)
- `expression_matrix_counts.csv` (Raw Counts)

---

### 2. Figure7_分析文档.md
**对应代码**: `02_fig7A_univariate_cox.R` ~ `02_fig7M_expression_diff.R` (共13个文件)

**主题**: 宫颈癌成纤维细胞预后模型构建

包含13个子图分析：

| 子图 | 序号 | 代码文件 | 主要内容 |
|------|------|----------|----------|
| 7A | 1 | 02_fig7A_univariate_cox.R | 单因素Cox回归森林图 |
| 7B | 2 | 02_fig7B_lasso.R | Lasso回归筛选基因 |
| 7C | 3 | 02_fig7C_multivariate_cox.R.R | 多因素Cox回归 |
| 7D | 4 | 02_fig7D_7E_risk_score.R | 风险基因系数柱状图 |
| 7E | 5 | 02_fig7D_7E_risk_score.R | 风险评分曲线与生存状态 |
| 7F | 6 | 02_fig7F_7G_heatmap_roc.R | 差异表达热图 |
| 7G | 7 | 02_fig7F_7G_heatmap_roc.R | 时间依赖性ROC曲线 |
| 7H | 8 | 02_fig7H_pca.R | PCA散点图 |
| 7I | 9 | 02_fig7I_km.R | SDC1基因生存曲线 |
| 7J | 10 | 02_fig7J_multivariate_clinical.R | 临床因素多因素Cox |
| 7K | 11 | 02_fig7K_nomogram.R | 列线图（Nomogram） |
| 7L | 12 | 02_fig7L_correlation.R | 相关性矩阵图 |
| 7M | 13 | 02_fig7M_expression_diff.R | 基因表达差异箱线图 |

**分析流程**: 基因筛选 → 模型构建 → 模型验证 → 临床应用 → 机制探索

---

### 3. Figure8_分析文档.md
**对应代码**: `03_fig8A_scRNA_process.R` ~ `03_fig8E_drug_sensitivity.R` (共6个文件)

**主题**: 肿瘤边界特征与预后模型验证

包含5个子图分析：

| 子图 | 序号 | 代码文件 | 主要内容 |
|------|------|----------|----------|
| 8A | 1 | 03_fig8A_scRNA_process.R<br>03_fig8A_spatial_plot.R | 空间转录组特征图<br>（注：唯一对应多个代码文件的子图） |
| 8B | 2 | 03_fig8B_km_mfrs.R | 高低MFRS组生存曲线 |
| 8C | 3 | 03_fig8C_radar_plot.R | 10个预后基因表达雷达图 |
| 8D | 4 | 03_fig8D_clinical_pie.R | 临床特征分布饼图 |
| 8E | 5 | 03_fig8E_drug_sensitivity.R | 免疫药物IC50预测 |

**分析重点**: 多模态数据验证（空间转录组 + 生存分析 + 药物敏感性）

---

### 4. Figure9_分析文档.md
**对应代码**: `04_fig9A_heatmap.R` ~ `04_fig9H_I_radar.R` (共7个文件)

**主题**: 高低MFRS组间差异基因富集分析

包含9个子图分析：

| 子图 | 序号 | 代码文件 | 主要内容 |
|------|------|----------|----------|
| 9A | 1 | 04_fig9A_heatmap.R | 差异基因热图 |
| 9B | 2 | 04_fig9B_volcano.R | 差异基因火山图 |
| 9C | 3 | 04_fig9C_E_enrichment.R | GO-BP富集分析 |
| 9D | 4 | 04_fig9C_E_enrichment.R | GO-CC富集分析 |
| 9E | 5 | 04_fig9C_E_enrichment.R | KEGG通路富集 |
| 9F | 6 | 04_fig9F_gsea.R | GSEA富集分析 |
| 9G | 7 | 04_fig9G_gsva.R | GSVA富集分析热图 |
| 9H | 8 | 04_fig9H_I_radar.R | 差异富集通路雷达图 |
| 9I | 9 | 04_fig9H_I_radar.R | 免疫细胞浸润雷达图 |

**分析重点**: 差异基因识别 → 功能注释 → 通路分析 → 免疫景观

---

## 文档结构说明

每个子图的分析遵循统一的格式：

### （1）代码文件名
列出对应的R或Python文件

### （2）依赖的包
- 列出所有使用的R包或Python库
- 注明是CRAN包还是Bioconductor包
- 说明特殊安装要求（如GitHub安装）

### （3）方法分析
**完整、清晰、简洁**的分步说明：
- 整体流程概述
- 详细步骤分解
- 关键代码片段
- 参数说明
- 设计要点
- 输出文件

### （4）结果比较
预留空白部分，供用户填写与论文结果的对比

---

## 特殊说明

### 多文件对应关系

1. **一个代码文件对应多个子图**:
   - `02_fig7D_7E_risk_score.R` → 图7D和7E
   - `02_fig7F_7G_heatmap_roc.R` → 图7F和7G
   - `04_fig9C_E_enrichment.R` → 图9C、9D和9E
   - `04_fig9H_I_radar.R` → 图9H和9I

2. **多个代码文件对应一个子图**:
   - 图8A ← `03_fig8A_scRNA_process.R` + `03_fig8A_spatial_plot.R`

### 数据流向

```
01_data_prep.py (数据准备)
    ↓
    ├─→ expression_matrix.csv (TPM) ──→ 用于生存分析和风险评分
    ├─→ expression_matrix_counts.csv ──→ 用于DESeq2差异分析
    └─→ clinical_cleaned.csv ──────────→ 用于临床信息整合
         ↓
    Fig 7 系列 (预后模型构建)
         ↓
    Fig 8 系列 (模型验证)
         ↓
    Fig 9 系列 (富集分析)
```

---

## 使用建议

1. **初学者**: 按顺序阅读（数据准备 → Fig7 → Fig8 → Fig9）
2. **关注特定分析**: 直接跳转到对应的Figure文档
3. **代码复现**: 结合文档中的方法分析和关键代码片段
4. **结果验证**: 在"结果比较"部分记录自己的复现结果

---

## 统计信息

- **总文档数**: 4个
- **总代码文件数**: 24个 (1个.py + 23个.R)
- **总子图数**: 27个 (数据准备 + 13个Fig7 + 5个Fig8 + 9个Fig9)
- **总文档行数**: ~1800行

---

## 更新日志

- 2026-01-01: 完成所有文档的初始版本
