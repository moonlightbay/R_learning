import os
import pandas as pd
import numpy as np
import glob

# --- 配置路径 ---
BASE_DIR = r"d:\Works\R_learning\data\fig7-9\prognostic"
DATA_DIR = os.path.join(BASE_DIR, "data")
SAMPLE_SHEET_PATH = os.path.join(BASE_DIR, "gdc_sample_sheet.tsv")
CLINICAL_PATH = os.path.join(BASE_DIR, "clinical.tsv")
OUTPUT_DIR = r"d:\Works\R_learning\project\processed_data"

if not os.path.exists(OUTPUT_DIR):
    os.makedirs(OUTPUT_DIR)

print("正在读取 Sample Sheet...")
sample_sheet = pd.read_csv(SAMPLE_SHEET_PATH, sep="\t")
# 建立 File ID -> Case ID 的映射 (根据用户指示修改)
# 注意：这将直接使用病人ID (如 TCGA-XX-XXXX) 作为列名，方便与临床数据对齐
file_id_to_case = dict(zip(sample_sheet['File ID'], sample_sheet['Case ID']))
file_name_to_case = dict(zip(sample_sheet['File Name'], sample_sheet['Case ID']))

print("正在读取 Clinical Data...")
# 增加 na_values 参数处理 TCGA 常见的缺失值标记
# 用户特别指出存在 "'--" 这种格式的缺失值
clinical = pd.read_csv(CLINICAL_PATH, sep="\t", na_values=["--", "not reported", "[Not Applicable]", "[Not Available]", "[Discrepancy]", "'--"])

# --- 列名适配 (根据用户反馈修正) ---
# 定义可能的列名映射：原始列名 -> 标准化列名
# TCGA GDC 下载的数据通常带有 'category.field' 格式的前缀
possible_mappings = {
    'cases.submitter_id': 'case_submitter_id',
    'submitter_id': 'case_submitter_id',
    
    'demographic.vital_status': 'vital_status',
    'vital_status': 'vital_status',
    
    'demographic.days_to_death': 'days_to_death',
    'days_to_death': 'days_to_death',
    
    'diagnoses.days_to_last_follow_up': 'days_to_last_follow_up',
    'days_to_last_follow_up': 'days_to_last_follow_up',
    
    'demographic.age_at_index': 'age_at_index',
    'age_at_index': 'age_at_index',
    
    'demographic.race': 'race',
    'race': 'race',
    
    'diagnoses.ajcc_pathologic_stage': 'ajcc_pathologic_stage',
    'ajcc_pathologic_stage': 'ajcc_pathologic_stage'
}

# 筛选出实际存在的列并重命名
found_cols = [c for c in possible_mappings.keys() if c in clinical.columns]
clinical_clean = clinical[found_cols].copy()
clinical_clean = clinical_clean.rename(columns=possible_mappings)

# 去除重复列 (以防映射后有重复，例如同时存在 race 和 demographic.race)
clinical_clean = clinical_clean.loc[:, ~clinical_clean.columns.duplicated()]

# --- 新增：去除重复病人 (关键修正) ---
# TCGA 临床数据中同一个病人可能有多个条目（对应不同治疗阶段），但生存数据通常是一样的
# 我们只需要保留每个病人的一条记录，否则会造成统计膨胀
if 'case_submitter_id' in clinical_clean.columns:
    before_dedup = len(clinical_clean)
    clinical_clean = clinical_clean.drop_duplicates(subset=['case_submitter_id'], keep='first')
    print(f"已去除重复病人记录: {before_dedup} -> {len(clinical_clean)}")
# ------------------------------------

print(f"已提取并重命名临床列: {clinical_clean.columns.tolist()}")

# 处理生存时间：优先用 days_to_death，如果是 Alive 则用 days_to_last_follow_up
def get_survival_time(row):
    # 尝试获取 days_to_death
    d_death = row.get('days_to_death')
    if pd.notna(d_death):
        return float(d_death)
    
    # 尝试获取 days_to_last_follow_up
    d_last = row.get('days_to_last_follow_up')
    if pd.notna(d_last):
        return float(d_last)
        
    return None # 无法获取生存时间，后续在R中过滤

def get_status(row):
    status = row.get('vital_status')
    if pd.isna(status): return None
    if status == 'Dead': return 1
    if status == 'Alive': return 0
    return None

clinical_clean['OS_time'] = clinical_clean.apply(get_survival_time, axis=1)
clinical_clean['OS_status'] = clinical_clean.apply(get_status, axis=1)
# 保存清洗后的临床数据
clinical_clean.to_csv(os.path.join(OUTPUT_DIR, "clinical_cleaned.csv"), index=False)

print("正在遍历数据文件夹并合并表达矩阵 (这可能需要几分钟)...")
# 查找所有 .tsv 文件
all_files = glob.glob(os.path.join(DATA_DIR, "**", "*.tsv"), recursive=True)

merged_data = {}

for i, file_path in enumerate(all_files):
    file_name = os.path.basename(file_path)
    
    # 尝试从文件名匹配 Case ID
    case_id = file_name_to_case.get(file_name)
    
    # 如果文件名匹配不到，尝试用父文件夹名 (File ID)
    if not case_id:
        parent_folder = os.path.basename(os.path.dirname(file_path))
        case_id = file_id_to_case.get(parent_folder)
        
    if not case_id:
        continue # 无法识别样本，跳过

    # 读取表达量文件 (跳过注释行)
    # TCGA STAR Counts 通常格式: gene_id, gene_name, gene_type, unstranded, ...
    try:
        df = pd.read_csv(file_path, sep="\t", comment='#', header=0)
        
        # --- 数据清洗关键步骤 ---
        # 1. 过滤掉 TCGA 统计行 (N_unmapped, N_multimapping 等)
        # 这些行的 gene_id 通常以 'N_' 开头
        if 'gene_id' in df.columns:
            df = df[~df['gene_id'].str.startswith('N_')]
            
        # 2. 确保 gene_name 存在且不为空
        if 'gene_name' in df.columns:
            df = df.dropna(subset=['gene_name'])
        # -----------------------

        # 假设第4列是 unstranded counts (根据TCGA标准)
        # 如果列名不同，请检查文件内容。通常是 'unstranded' 或 'tpm_unstranded'
        if 'tpm_unstranded' in df.columns:
             # 优先使用 TPM (Transcripts Per Million) 进行生存分析，因为它已经归一化了测序深度和基因长度
            counts = df.set_index('gene_name')['tpm_unstranded']
        elif 'unstranded' in df.columns:
            counts = df.set_index('gene_name')['unstranded']
        else:
            # 尝试取第4列
            counts = df.set_index(df.columns[1]).iloc[:, 2] 
            
        # 处理重复基因名
        # ---------------------------------------------------------
        # 方案 A: 仅保留第一个 (原先的处理方式 - 已注释)
        # 缺点：可能会随机丢弃该基因的其他转录本信息
        # counts = counts[~counts.index.duplicated(keep='first')]
        
        # 方案 B: 分组取平均 (推荐 - 新增)
        # 优点：综合了该基因所有条目的表达量，更稳健
        counts = counts.groupby(level=0).mean()
        # ---------------------------------------------------------
        
        merged_data[case_id] = counts
    except Exception as e:
        print(f"Error reading {file_name}: {e}")

    if i % 50 == 0:
        print(f"已处理 {i}/{len(all_files)} 个文件...")

print("正在构建最终矩阵...")
expression_matrix = pd.DataFrame(merged_data)
# 填充缺失值为0 (针对某个样本完全缺失某个基因的情况)
expression_matrix = expression_matrix.fillna(0)

# --- 数据转换 (Log Transformation) ---
# ---------------------------------------------------------
# 方案 A: 原始 TPM (不推荐 - 已注释)
# 缺点：数据分布高度偏态，极值会影响 Lasso/Cox 模型性能
# pass 

# 方案 B: Log2(TPM + 1) 转换 (推荐 - 新增)
# 优点：稳定方差，使数据接近正态分布，符合统计模型假设
print("正在进行 Log2(TPM + 1) 转换...")
expression_matrix = np.log2(expression_matrix + 1)
# ---------------------------------------------------------

# 保存
expression_matrix.to_csv(os.path.join(OUTPUT_DIR, "expression_matrix.csv"))
print(f"完成！数据已保存至 {OUTPUT_DIR}")