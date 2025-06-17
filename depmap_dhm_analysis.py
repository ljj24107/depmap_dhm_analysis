import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import os

plt.rcParams['font.size'] = 15

class DHMDepMapAnalyzer:
    """DepMap数据库分析类"""

    def __init__(self, data_dir):
        """初始化数据目录与数据存储变量"""
        self.data_dir = data_dir
        self.sample_info = None
        self.ccle_expression = None
        self.crispr_gene_effect = None
        self.drug_sensitivity = None
        self.meta = None
        self.dhm_col = None
        self.liver_lines = None

    def load_data(self):
        """加载DepMap数据文件"""
        try:
            self.sample_info = pd.read_csv(os.path.join(self.data_dir, "sample_info.csv"))
            self.ccle_expression = pd.read_csv(os.path.join(self.data_dir, "CCLE_expression.csv"), index_col=0)
            self.crispr_gene_effect = pd.read_csv(os.path.join(self.data_dir, "CRISPR_gene_effect.csv"), index_col=0)
            drug = pd.read_csv(os.path.join(self.data_dir, "Drug_sensitivity.csv"))
        except Exception as e:
            raise RuntimeError("数据文件读取失败") from e

        # 标准化Drug字段
        drug = drug.rename(columns={drug.columns[0]: "DepMap_ID"}) #重命名第一列，方便后面计算

        # 查找包含 "dihydromyricetin" 的列名（已知有且只有一个匹配的列）
        dhm_cols = next(col for col in drug.columns if "dihydromyricetin" in col.lower())
        self.dhm_col = "logfold_change"
        drug = drug.rename(columns={dhm_cols: self.dhm_col})
        print(dhm_cols)

        # 合并基本表
        self.drug_sensitivity = drug[["DepMap_ID", self.dhm_col]]
        self.meta = pd.merge(self.sample_info, self.drug_sensitivity, on="DepMap_ID", how="inner")
        print("\nData loaded. {} samples with DHM data.".format(len(self.meta)))
        #print(self.meta)

    def count_liver_cell_lines(self):
        """筛选肝脏特异性细胞系"""
        all_liver = self.sample_info[self.sample_info['sample_collection_site'] == 'liver']
        print(f"Total liver-derived cell lines in sample_info: {len(all_liver)}")
        self.liver_lines = self.meta[self.meta['sample_collection_site'] == 'liver']
        print(f"Liver-derived cell lines with DHM data: {len(self.liver_lines)}")

    def plot_violin_by_tissue(self):
        """生成不同组织来源依据药物敏感性的小提琴图"""

        organs = ["liver", "lung", "kidney", "skin"]
        df = self.meta.copy()
        df = df[df["sample_collection_site"].isin(organs)]

        plt.figure(figsize=(8, 6))
        sns.violinplot(x="sample_collection_site", y=self.dhm_col, data=df,
                       palette="Set2", hue="sample_collection_site", legend=False)
        plt.title("DHM sensitivity by tissue origin")
        plt.ylabel("DHM sensitivity (logFC)")
        plt.xlabel(None)
        plt.tight_layout()
        plt.savefig("1_violinplot_dhm_by_tissue.png")
        plt.close()

        print("\nMean ± SD and p-values comparing DHM sensitivity across tissues:")
        for organ in organs:
            group = df[df['sample_collection_site'] == organ][self.dhm_col]
            others = df[df['sample_collection_site'] != organ][self.dhm_col]
            mean = group.mean()
            std = group.std()
            _, p = stats.ttest_ind(group, others, equal_var=False, nan_policy='omit')
            print(f"{organ.capitalize()}: {mean:.2f} ± {std:.2f} (p = {p:.4f})")

        # 找出每种组织来源最敏感的细胞系
        print("\nMost DHM-sensitive cell line per tissue:")
        top1 = df.loc[df.groupby("sample_collection_site")[self.dhm_col].idxmin()]
        print(top1[["sample_collection_site", "cell_line_name", self.dhm_col]])

    def plot_crispr_dependency_bar(self):
        """画目标基因CRISPR依赖性柱状图"""

        gene_symbols = ["SIRT1", "PRKAA1", "NFKB1"]
        matched_cols = []
        for g in gene_symbols:
            matches = [col for col in self.crispr_gene_effect.columns if col.startswith(f"{g} ")]
            if matches:
                matched_cols.append(matches[0])
        if not matched_cols:
            return

        crispr_subset = self.crispr_gene_effect[matched_cols].copy()
        crispr_subset["DepMap_ID"] = self.crispr_gene_effect.index
        data = crispr_subset[crispr_subset["DepMap_ID"].isin(self.liver_lines["DepMap_ID"])]
        rename_dict = {col: col.split(" ")[0] for col in matched_cols}
        data = data.rename(columns=rename_dict)
        #print(data)

        melted = data.melt(id_vars="DepMap_ID", value_vars=gene_symbols, var_name="Gene", value_name="Score")
        #print(melted)
        plt.figure(figsize=(6, 6))
        sns.barplot(data=melted, x="Gene", y="Score", errorbar="sd", capsize=0.2,
                    hue="Gene", palette="Pastel2_r", legend=False)
        plt.title("CRISPR dependency scores in liver cell lines",fontsize=12)
        plt.tight_layout()
        plt.xlabel(None)
        plt.ylabel("Score", fontsize=15)
        plt.savefig("2_bar_crispr_dependency.png")
        plt.close()

        means = melted.groupby("Gene")["Score"].mean()
        print("\nAverage CRISPR dependency scores:")
        print(means)

    def plot_expression_heatmap(self):
        """肝细胞系目标基因表达热图"""

        genes = ["SIRT1", "PRKAA1", "NFKB1"]
        expr_subset = self.ccle_expression.loc[self.ccle_expression.index.isin(self.liver_lines["DepMap_ID"])]
        matched_cols = [col for g in genes for col in self.ccle_expression.columns if col.startswith(g)]
        expr_data = expr_subset[matched_cols]
        # 对每个基因的表达量做标准化
        expr_scaled = (expr_data - expr_data.mean()) / expr_data.std()

        xticklabels = [col.split(" ")[0] for col in matched_cols]
        plt.figure(figsize=(5, 8))
        sns.heatmap(expr_scaled, cmap="RdBu_r", xticklabels=xticklabels, yticklabels=True)
        plt.title("Expression heatmap of target genes in liver cell lines", fontsize=12)
        plt.ylabel("Cell lines")
        plt.yticks(fontsize=11)
        plt.xticks(fontsize=12)
        plt.tight_layout()
        plt.savefig("3_heatmap_target_expression.png")
        plt.close()

    def plot_sirt1_dependency_vs_dhm(self):
        """SIRT1依赖性与DHM敏感性相关性散点图"""

        crispr = self.crispr_gene_effect.reset_index().rename(columns={"index": "DepMap_ID"})
        sirt1_col = [col for col in crispr.columns if "SIRT1" in col][0]
        gene_name = sirt1_col.split(" ")[0]
        sirt1_score = crispr[["DepMap_ID", sirt1_col]]
        merged = pd.merge(self.meta, sirt1_score, on="DepMap_ID")
        # 数据清洗
        merged = merged.replace([np.inf, -np.inf], np.nan).dropna(subset=[self.dhm_col, sirt1_col])

        if merged.shape[0] < 3:
            return

        r, p = stats.pearsonr(merged[sirt1_col], merged[self.dhm_col])
        plt.figure(figsize=(6, 5))
        sns.regplot(x=sirt1_col, y=self.dhm_col, data=merged)
        plt.title(f"{gene_name} CRISPR dependency vs DHM sensitivity\n(r={r:.2f}, p={p:.3g})", fontsize=12)
        plt.xlabel("CRISPR dependency")
        plt.ylabel("DHM sensitivity")
        plt.tight_layout()
        plt.savefig("4_catter_sirt1_vs_dhm.png")
        plt.close()

        print(f"\nSIRT1 CRISPR dependency vs DHM sensitivity: r={r:.2f}, p={p:.4f}")

    def plot_sirt1_expression_by_tissue(self):
        """不同组织SIRT1表达箱型图"""

        organs = ["liver", "lung", "kidney", "skin"]
        sirt1_col = [col for col in self.ccle_expression.columns if col.startswith("SIRT1")][0]
        expr_subset = self.ccle_expression[[sirt1_col]].copy()
        expr_subset["DepMap_ID"] = expr_subset.index
        merged = pd.merge(self.sample_info, expr_subset, on="DepMap_ID")
        merged = merged[merged["sample_collection_site"].isin(organs)]

        plt.figure(figsize=(8, 8))
        sns.boxplot(x="sample_collection_site", y=sirt1_col, data=merged, palette="RdYlBu",
                    hue="sample_collection_site", legend=False)
        #plt.xticks(rotation=45)
        plt.title("SIRT1 expression across selected tissues")
        plt.xlabel(None)
        plt.ylabel("Relative Expression")
        plt.tight_layout()
        plt.savefig("5_boxplot_sirt1_by_tissue.png")
        plt.close()

    def correlate_expression_with_dhm(self):
        """筛选肝细胞系中与DHM敏感性最相关的基因"""

        expr_subset = self.ccle_expression.loc[self.ccle_expression.index.isin(self.liver_lines["DepMap_ID"])].copy()
        dhm_scores = self.liver_lines.set_index("DepMap_ID")[self.dhm_col]
        common_ids = expr_subset.index.intersection(dhm_scores.index)
        expr_subset = expr_subset.loc[common_ids]
        dhm_scores = dhm_scores.loc[common_ids]

        # 数据清洗，去掉缺失值
        mask = ~(expr_subset.isna().any(axis=1) | dhm_scores.isna())
        expr_subset = expr_subset[mask]
        dhm_scores = dhm_scores[mask]

        expr_subset = expr_subset.loc[:, expr_subset.std() > 0]

        def safe_pearsonr(col):
            if col.isnull().any() or dhm_scores.isnull().any():
                return np.nan
            try:
                return stats.pearsonr(col, dhm_scores)[0]
            except Exception:
                return np.nan

        correlations = expr_subset.apply(safe_pearsonr)
        top_corr = correlations.abs().sort_values(ascending=False).head(5)
        print("\nTop 5 genes most correlated with DHM sensitivity:")
        print(top_corr)

    def run_all(self):
        self.load_data()
        self.count_liver_cell_lines()
        self.plot_violin_by_tissue()
        self.plot_crispr_dependency_bar()
        self.plot_expression_heatmap()
        self.plot_sirt1_dependency_vs_dhm()
        self.plot_sirt1_expression_by_tissue()
        self.correlate_expression_with_dhm()

if __name__ == "__main__":
    analyzer = DHMDepMapAnalyzer(data_dir="./depmap")
    analyzer.run_all()
