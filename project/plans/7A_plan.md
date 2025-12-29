# Figure 7A 复现方案：单因素 Cox 回归森林图

## 1. 目标描述
**图表类型**：森林图 (Forest Plot)
**内容**：展示单因素 Cox 回归分析的结果。
**目的**：识别与宫颈癌患者预后显著相关的基因（特别是 C0 亚型相关的成纤维细胞标记基因）。
**统计指标**：
- **P value**: 显著性水平 (通常筛选 P < 0.05)
- **Hazard Ratio (HR)**: 风险比。HR > 1 为风险因子 (Risk Factor)，HR < 1 为保护因子 (Protective Factor)。
- **95% CI**: 置信区间。

## 2. 数据准备
需要使用之前预处理生成的两个文件：
1.  `project/processed_data/expression_matrix.csv`: 基因表达矩阵 (行=基因, 列=样本, 已做 Log2 转换)。
2.  `project/processed_data/clinical_cleaned.csv`: 临床数据 (包含 `OS_time`, `OS_status`)。

## 3. 分析流程
1.  **数据加载与合并**：
    - 读取表达矩阵，转置为 (行=样本, 列=基因)。
    - 读取临床数据。
    - 通过样本 ID (`case_submitter_id`) 将两者合并。
    - 过滤掉生存时间或状态缺失的样本。

2.  **基因筛选 (关键步骤)**：
    - 论文提到使用 "C0 MYH11+ Fibroblasts-related genes"。
    - **策略**：用户已从论文图片中识别出 13 个关键基因。
    - **基因列表**: `SMTN`, `SEMA3C`, `RHOB`, `PPP1R14A`, `PDLIM7`, `PCP4`, `FHL2`, `DSTN`, `DES`, `CUX1`, `CALD1`, `CA12`, `ACTN1`。
    - **验证**: 代码将自动检查这些基因是否存在于我们的表达矩阵中。

3.  **单因素 Cox 回归**：
    - 对每个候选基因，构建 Cox 模型：`Surv(time, status) ~ gene_expression`。
    - 提取 HR, 95% CI, P值。

4.  **绘图**：
    - 使用 `survminer` 包的 `ggforest` 函数（最简单）或 `ggplot2` 手绘（更灵活）。
    - 按照论文样式：左侧显示基因名和 P值，中间是森林图（点代表 HR，线代表 CI），右侧显示 HR 数值。

## 4. R 代码实现
脚本将保存为 `project/02_fig7A_univariate_cox.R`。

注意：已完成！大致和论文图一致。对于较为微小的不同，解释如下：
可能的原因：
1. 生存数据更新 (Follow-up Data)
2. 样本量差异（不同的数据筛选标准或缺失值处理）
3. 数据处理流程：我使用的是 TPM 数据并取了 Log2 转换，论文可能使用 FPKM 或其他标准化方法。