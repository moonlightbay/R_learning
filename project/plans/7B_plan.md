7B 应用 Lasso 回归选择 10 个对风险评分有贡献的基因，结果显示在 lambda 图中 (lambda.min = 0.017)。
## 2. 数据准备
需要使用之前预处理生成的两个文件：
1.  `project/processed_data/expression_matrix.csv`: 基因表达矩阵 (行=基因, 列=样本, 已做 Log2 转换)。
2.  `project/processed_data/clinical_cleaned.csv`: 临床数据 (包含 `OS_time`, `OS_status`)。