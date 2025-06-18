# DHM DepMap Analysis

本项目基于 DepMap 数据库分析二氢杨梅素（DHM）的药物敏感性及其与基因表达的关联，揭示 DHM 作用的潜在分子机制。

## 项目背景

DHM 是一种天然黄酮类化合物，具有抗炎、抗氧化和代谢调节作用。本项目旨在：

- 评估不同组织来源细胞系对 DHM 的敏感性；
- 探索关键代谢/炎症相关基因（SIRT1、AMPK、NF-κB）在肝细胞中的依赖性与表达特征；
- 结合 CRISPR 敲除与表达数据，分析潜在的功能关联。

## 数据来源

- [DepMap Public ](https://depmap.org/portal/)  
  - `sample_info.csv` 
  - `CCLE_expression.csv` 
  - `CRISPR_gene_effect.csv` 
  - `Drug_sensitivity_(PRISM_Repurposing_Primary_Screen)_19Q4.csv`  
## 分析内容与图表输出

1. **小提琴图**：展示 "liver"、"lung"、"kidney"、"skin" 细胞系对 DHM 的敏感性差异；
2. **柱状图**：肝细胞系中 SIRT1、PRKAA1、NFKB1 的 CRISPR 依赖性平均评分；
3. **热图**：目标基因在肝细胞系中的表达情况；
4. **散点图**：SIRT1 依赖性与 DHM 敏感性的 Pearson 相关性；
5. **表达分布图**：SIRT1 在不同组织来源细胞系中的表达比较。

## 目录结构

```plaintext
depmap_dhm_analysis/
├── depmap_dhm_analysis.py      # 主程序文件
├── data/                      # 数据文件夹
│   ├── sample_info.csv        # 细胞系样本元数据
│   ├── CCLE_expression.csv    # 基因表达数据
│   ├── CRISPR_gene_effect.csv # CRISPR 依赖性数据
│   └── Drug_sensitivity.csv   # 药物敏感性数据
├── results/                   # 结果文件夹
│   ├── 1_violinplot_dhm_by_tissue.png  # 组织敏感性小提琴图
│   ├── 2_bar_crispr_dependency.png       # CRISPR 依赖性条形图
│   ├── 3_heatmap_target_expression.png   # 基因表达热图
│   ├── 4_scatter_sirt1_vs_dhm.png        # SIRT1 相关性散点图
│   └── 5_boxplot_sirt1_by_tissue.png     # SIRT1 表达箱线图
└── README.md                   # 项目说明文档
```

## 使用方法

```bash
# 1. 安装依赖（需 Python 3.8+）
pip install pandas numpy matplotlib seaborn scipy

# 2. 将 DepMap 文件放入 ./depmap/ 文件夹下

# 3. 运行主程序
python depmap_dhm_analysis.py
# 4. 结果将保存在 ./文件夹中，分析结果将输出到控制台
```

## 分析流程

1. **数据加载与预处理**：
   - 加载样本信息、基因表达、CRISPR 依赖性和药物敏感性数据
   - 筛选包含 DHM 的药物敏感性列
   - 合并数据并筛选肝脏来源细胞系
2. **组织特异性敏感性分析**：
   - 绘制不同组织细胞系的 DHM 敏感性小提琴图
   - 计算各组织平均敏感性及组间差异 p 值
   - 识别各组织中对 DHM 最敏感的细胞系
3. **CRISPR 依赖性分析**：
   - 提取 SIRT1、PRKAA1、NFKB1 基因的 CRISPR 依赖性数据
   - 计算肝脏细胞系中这些基因的平均依赖性分数
   - 绘制条形图展示依赖性结果
4. **基因表达与敏感性关联分析**：
   - 分析 SIRT1 基因 CRISPR 依赖性与 DHM 敏感性的相关性
   - 计算全基因组基因表达与 DHM 敏感性的皮尔逊相关系数
   - 筛选相关性最高的前5个基因

## 输出结果

### 控制台输出

```plaintext
dihydromyricetin (BRD:BRD-K01614093-001-02-6)

Data loaded. 568 samples with DHM data.
Total liver-derived cell lines in sample_info: 54
Liver-derived cell lines with DHM data: 33

Mean ± SD and p-values comparing DHM sensitivity across tissues:
Liver: 0.06 ± 0.54 (p = 0.7420)
Lung: 0.15 ± 0.52 (p = 0.2369)
Kidney: 0.09 ± 0.41 (p = 0.9855)
Skin: 0.01 ± 0.38 (p = 0.2218)

Most DHM-sensitive cell line per tissue:
    sample_collection_site cell_line_name  logfold_change
220                 kidney          A-704       -0.556862
74                   liver       SU.86.86       -1.660440
28                    lung       STM91-01       -1.275652
296                   skin          EBC-1       -0.667397

Average CRISPR dependency scores:
Gene
NFKB1     0.043778
PRKAA1    0.002072
SIRT1     0.011826

SIRT1 CRISPR dependency vs DHM sensitivity: r=-0.03, p=0.5220

Top 5 genes most correlated with DHM sensitivity:
OR5B21 (219968)     0.765277
RGS7BP (401190)     0.691619
MYCBPAP (84073)     0.664929
BPIFA1 (51297)      0.657830
PLA2G2D (26279)     0.655268
```

### 图表输出

- `1_violinplot_dhm_by_tissue.png`：不同组织 DHM 敏感性小提琴图
- `2_bar_crispr_dependency.png`：靶点基因 CRISPR 依赖性条形图
- `3_heatmap_target_expression.png`：靶点基因表达热图
- `4_scatter_sirt1_vs_dhm.png`：SIRT1 依赖性与 DHM 敏感性散点图
- `5_boxplot_sirt1_by_tissue.png`：SIRT1 表达组织分布箱线图

## 项目贡献

本项目由 Jiajun LEI 完成，用于Med5018 课程结课作业。





