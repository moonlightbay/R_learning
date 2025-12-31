# Figure 9: 高/低 MFRS 组间差异基因的富集分析

## 图表说明

**图 9** 展示了高/低 MFRS 组间差异基因的富集分析结果，包含以下子图：

- **图 A**: 热图显示了高/低 MFRS 组间差异基因的表达谱
- **图 B**: 火山图展示了高/低 MFRS 组间差异上调和差异下调的基因
- **图 C-E**: 条形图展示了高/低 MFRS 组间差异基因的 GO-BP (生物过程)、GO-CC (细胞组分) 和 KEGG 通路分析结果
- **图 F**: 基于差异上调和下调基因的详细 GSEA 富集结果
- **图 G**: 热图显示了高/低 MFRS 组间差异基因集的 GSVA 分析结果
- **图 H**: 雷达图比较了高/低 MFRS 组间通过不同富集方法获得的条目差异
- **图 I**: 雷达图显示了高/低 MFRS 组间不同状态下各种免疫细胞类型含量的差异

---

## 详细分析

### 高/低 MFRS 组间的富集分析

为进一步探索高/低 MFRS 组间的功能差异，我们分析了两组间的差异表达基因 (DEGs) 并研究了它们所富集的通路。

#### 1. 差异表达基因识别

我们首先识别了高/低 MFRS 组间的差异表达基因 (图 9A, B)。

#### 2. 功能富集分析

随后进行了 GO-BP、GO-CC 和 KEGG 通路分析 (图 9C-E)：

**GO-BP (生物过程) 分析**:
- 主要富集在消化过程、中间丝组织等相关通路

**GO-CC (细胞组分) 分析**:
- 富集于质膜外泌成分、刷状缘、肌动蛋白细胞突触聚集体、连环蛋白复合体、中间丝细胞骨架和膜外泌成分

**KEGG 通路分析**:
- 突出显示了与以下通路相关的基因：
  - 年轻型糖尿病发病机制 (maturity onset diabetes of the young)
  - Ras 信号传导
  - 金黄色葡萄球菌感染
  - 碳水化合物消化和吸收
  - 癌症中的转录失调

#### 3. 基因集富集分析 (GSEA)

我们进行了 GSEA 分析，结果显示：

**上调基因** 主要富集在以下过程：
- 胶原分解代谢过程
- 伤口愈合
- 表皮细胞扩散
- 内胚层形成
- 细胞外基质分解

**下调基因** 主要富集在：
- 钾离子稳态
- 杀死其他生物体的细胞
- 轴丝组装
- 轴丝组装 (图 9F)

#### 4. 基因集变异分析 (GSVA)

我们将几个富集通路与我们构建的预后模型的风险评分进行了相关性分析 (图 9G)。值得注意的是，以下术语与风险评分呈显著正相关：
- `PEDERSEN_TARGETS_OF_611_CTF_ISOFORM_OF_ERBB2`
- `REACTOME_CELL_EXTRACELLULAR_MATRIX_INTERACTIONS`
- `BRUECKNER_TARGETS_OF_MIRLET7A3_DN`
- `ELVIDGE_HYPOXIA_UP`

这表明这些通路具有潜在的统计学意义。

#### 5. 雷达图分析

为进一步研究高/低 MFRS 组间显著富集的通路，我们使用雷达图进行分析 (图 9H)。结果显示，以下术语在高 MFRS 组中显著富集：

- `GOBP_BRANCHED_CHAIN_AMINO_ACID_TRANSPORT` (支链氨基酸转运)
- `GOBP_REGULATION_OF_VENTRICULAR_CARDIAC_MUSCLE_CELL_ACTION_POTENTIAL` (心室心肌细胞动作电位调节)
- `GOBP_NEGATIVE_REGULATION_OF_CYCLIC_NUCLEOTIDE_PHOSPHODIESTERASE_ACTIVITY` (环核苷酸磷酸二酯酶活性的负向调节)
- `GOBP_PROTEIN_LOCALIZATION_TO_PLASMA_MEMBRANE_RAFT` (蛋白向质膜筏的定位)
- `GOBP_INHIBITORY_SYNAPSE_ASSEMBLY` (抑制性突触组装)
- `GOBP_NEGATIVE_REGULATION_OF_GLUTAMATE_SECRETION` (谷氨酸分泌的负向调节)

#### 6. 免疫细胞景观分析

最后，为了探索两组间的免疫细胞景观，我们发现高 MFRS 组中以下免疫细胞更为丰富 (图 9I)：

- 活化的树突状细胞 (activated dendritic cells)
- M0 巨噬细胞 (M0 macrophages)
- 活化的自然杀伤细胞 (activated NK cells)
- 静息自然杀伤细胞 (resting NK cells)

---

## 总结

本分析通过多种富集分析方法，全面揭示了高/低 MFRS 组间在基因表达、通路富集和免疫微环境方面的显著差异。这些发现为进一步理解 MFRS 相关的肿瘤生物学机制和潜在的免疫治疗靶点提供了重要线索。
