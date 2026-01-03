#!/usr/bin/env python3

# ============================================================
# 1. Import libraries
# ============================================================
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.patches as mpatches
from sklearn.metrics import pairwise_distances
from skbio import DistanceMatrix
from skbio.stats.ordination import pcoa
from skbio.stats.distance import permanova

# ============================================================
# 2. Confidence ellipse
# ============================================================
def plot_confidence_ellipse(x, y, ax, n_std=1.96, **kwargs):
    """Draw a 95% confidence ellipse."""
    if len(x) < 2:
        return
    cov = np.cov(x, y)
    mean_x, mean_y = np.mean(x), np.mean(y)
    lambda_, v = np.linalg.eig(cov)
    lambda_ = np.sqrt(lambda_)
    ell = mpatches.Ellipse(
        xy=(mean_x, mean_y),
        width=lambda_[0] * n_std * 2,
        height=lambda_[1] * n_std * 2,
        angle=np.rad2deg(np.arccos(v[0, 0])),
        **kwargs
    )
    ax.add_patch(ell)

# ============================================================
# 3. User settings – your HPC paths
# ============================================================
tpm_file = "/ourdisk/hpc/rnafold/dywang/dont_archive/On_process/V49/GeneExpression/filtered_normalized_Sidbers_TPM.tsv"
psm_file = "/ourdisk/hpc/rnafold/dywang/dont_archive/On_process/V45/Sidbers/Proteome/Sindbis_Pro/filtered_normalized_Human_Sidbers_PSM.tsv"
meta_file = "/ourdisk/hpc/rnafold/dywang/dont_archive/On_process/V49/GeneExpression/MetaData.txt"

figure_file = "/ourdisk/hpc/rnafold/dywang/dont_archive/On_process/V49/GeneExpression/Bray-Curtis/PCoA_BrayCurtis.tif"
permanova_file = "/ourdisk/hpc/rnafold/dywang/dont_archive/On_process/V49/GeneExpression/Bray-Curtis/PERMANOVA_results.tsv"

dpi = 300
figsize = (6, 5)

# Color themes
consistent_colors = {"Control": "#56B4E0", "Treated": "indianred"}
panel_backgrounds = {"TPM": "#fff8dc", "PSM": "#e8f5e9"}

plt.rcParams.update({
    "font.size": 10,
    "axes.titlesize": 14,
    "axes.labelsize": 11,
    "legend.fontsize": 10,
})

# ============================================================
# 4. Read input data
# ============================================================
tpm_df = pd.read_csv(tpm_file, sep="\t", index_col=0)
psm_df = pd.read_csv(psm_file, sep="\t", index_col=0)
meta_df = pd.read_csv(meta_file, sep="\t")

meta_df.columns = meta_df.columns.str.strip()
for col in meta_df.columns:
    if meta_df[col].dtype == "object":
        meta_df[col] = meta_df[col].astype(str).str.strip()

tpm_df.columns = tpm_df.columns.str.strip()
psm_df.columns = psm_df.columns.str.strip()

print("TPM columns found:", list(tpm_df.columns))
print("PSM columns found:", list(psm_df.columns))
print("MetaData columns:", list(meta_df.columns))

# ============================================================
# 4a. Define grouping for samples
# ============================================================
def infer_group_from_name(col_name: str) -> str:
    """Fallback rule: use column name to infer Control vs Treated."""
    if "ctrl" in col_name.lower():
        return "Control"
    else:
        return "Treated"

has_explicit = {"Sample_TPM", "Sample_PSM", "dex"}.issubset(meta_df.columns)

if has_explicit:
    print("Using explicit Sample_TPM / Sample_PSM / dex mapping from MetaData.txt")
    meta_tpm = meta_df.dropna(subset=["Sample_TPM"])
    meta_psm = meta_df.dropna(subset=["Sample_PSM"])

    sample_group_map_tpm = dict(zip(meta_tpm["Sample_TPM"], meta_tpm["dex"]))
    sample_group_map_psm = dict(zip(meta_psm["Sample_PSM"], meta_psm["dex"]))

    tpm_samples = [c for c in tpm_df.columns if c in sample_group_map_tpm]
    psm_samples = [c for c in psm_df.columns if c in sample_group_map_psm]

    print("Matched TPM samples:", tpm_samples)
    print("Matched PSM samples:", psm_samples)

    if not tpm_samples:
        raise ValueError("No TPM columns match Sample_TPM in MetaData.txt.")
    if not psm_samples:
        raise ValueError("No PSM columns match Sample_PSM in MetaData.txt.")

    tpm_df = tpm_df[tpm_samples]
    psm_df = psm_df[psm_samples]
else:
    print(
        "MetaData.txt does not have Sample_TPM / Sample_PSM columns; "
        "inferring groups directly from column names (contains 'ctrl' = Control)."
    )
    sample_group_map_tpm = {col: infer_group_from_name(col) for col in tpm_df.columns}
    sample_group_map_psm = {col: infer_group_from_name(col) for col in psm_df.columns}

print("Sample → group map (TPM):", sample_group_map_tpm)
print("Sample → group map (PSM):", sample_group_map_psm)

# ============================================================
# 5. Compute Bray–Curtis distance & PCoA
# ============================================================
def compute_pcoa_with_variance(df, sample_group_map):
    """Compute Bray–Curtis distance, run PCoA, and return coords + variance."""
    dist = pairwise_distances(df.T, metric="braycurtis")
    dm = DistanceMatrix(dist, ids=df.columns)
    ord_res = pcoa(dm)

    coords = ord_res.samples.iloc[:, :2].copy()
    coords.columns = ["PCoA1", "PCoA2"]
    coords["Group"] = [sample_group_map[s] for s in coords.index]

    var1 = ord_res.proportion_explained.iloc[0] * 100
    var2 = ord_res.proportion_explained.iloc[1] * 100
    return coords, dist, var1, var2

tpm_plot_data, tpm_dist, tpm_var1, tpm_var2 = compute_pcoa_with_variance(
    tpm_df, sample_group_map_tpm
)
psm_plot_data, psm_dist, psm_var1, psm_var2 = compute_pcoa_with_variance(
    psm_df, sample_group_map_psm
)

# ============================================================
# 6. PERMANOVA
# ============================================================
results = []
for name, dist_matrix, meta_data in [
    ("TPM", tpm_dist, tpm_plot_data),
    ("PSM", psm_dist, psm_plot_data),
]:
    dm = DistanceMatrix(dist_matrix, ids=meta_data.index)
    grouping = meta_data["Group"]
    permanova_result = permanova(dm, grouping, permutations=999)

    try:
        test_stat = float(permanova_result["test statistic"])
        p_val = float(permanova_result["p-value"])
        n_perm = int(permanova_result["number of permutations"])
    except Exception as e:
        print("Warning: falling back to positional indexing for PERMANOVA:", e)
        test_stat = float(permanova_result.iloc[1])
        p_val = float(permanova_result.iloc[2])
        n_perm = int(permanova_result.iloc[3])

    results.append({
        "Dataset": name,
        "Test_statistic_F": test_stat,
        "R2": np.nan,
        "p_value": p_val,
        "Permutations": n_perm
    })

permanova_df = pd.DataFrame(results)
permanova_df.to_csv(permanova_file, sep="\t", index=False)
print("PERMANOVA results saved as:", permanova_file)
print(permanova_df)

# ============================================================
# 7. Plot PCoA figures
# ============================================================
fig, axes = plt.subplots(1, 2, figsize=figsize, sharex=True, sharey=True)

# make the two axes occupy more of the figure area
fig.subplots_adjust(
    left=0.10,   # space from left edge
    right=0.80,  # leave room on right for legend
    bottom=0.15,
    top=0.80,    # leave some space for suptitle
    wspace=0.25  # horizontal space between the two plots
)

# --- TPM panel ---
ax0 = axes[0]
ax0.set_facecolor(panel_backgrounds["TPM"])
sns.scatterplot(
    data=tpm_plot_data,
    x="PCoA1",
    y="PCoA2",
    hue="Group",
    s=60,
    ax=ax0,
    palette=consistent_colors,
    marker="o",
    legend=False
)
for group, subset in tpm_plot_data.groupby("Group"):
    plot_confidence_ellipse(
        subset["PCoA1"], subset["PCoA2"], ax0,
        edgecolor=consistent_colors[group],
        facecolor=consistent_colors[group],
        alpha=0.25,
        linewidth=2
    )
ax0.set_title("TPM (Transcriptomic Level)")
ax0.set_xlabel(f"PCoA 1 ({tpm_var1:.2f}% variance)")
ax0.set_ylabel(f"PCoA 2 ({tpm_var2:.2f}% variance)")
ax0.grid(True)

# --- PSM panel ---
ax1 = axes[1]
ax1.set_facecolor(panel_backgrounds["PSM"])
sns.scatterplot(
    data=psm_plot_data,
    x="PCoA1",
    y="PCoA2",
    hue="Group",
    s=60,
    ax=ax1,
    palette=consistent_colors,
    marker="o"
)
for group, subset in psm_plot_data.groupby("Group"):
    plot_confidence_ellipse(
        subset["PCoA1"], subset["PCoA2"], ax1,
        edgecolor=consistent_colors[group],
        facecolor=consistent_colors[group],
        alpha=0.25,
        linewidth=2
    )
ax1.set_title("PSM (Proteomic Level)")
ax1.set_xlabel(f"PCoA 1 ({psm_var1:.2f}% variance)")
ax1.set_ylabel("")  # avoid duplicate y-label
ax1.grid(True)

# --- Make axes auto-scale to include ellipses, then pad & equalize ---
for ax in axes:
    ax.relim()              # recompute limits using *all* artists
    ax.autoscale_view()     # expand limits to fit them
    ax.margins(0.15)        # add 15% padding around data
    ax.set_aspect("equal", adjustable="box")

# --- Legend ---
handles, labels = ax1.get_legend_handles_labels()
new_handles = [
    plt.Line2D([], [], color=consistent_colors[label], marker='o', linestyle='', markersize=10)
    for label in labels
]
ax1.legend(
    new_handles,
    labels,
    title="Group",
    loc="center left",
    bbox_to_anchor=(1.02, 0.5),
    borderaxespad=0,
    frameon=False
)

plt.suptitle("PCoA (Bray–Curtis) — Control vs Treated at Two Biological Levels", fontsize=14)

# removed tight_layout so our manual subplots_adjust controls sizing
# plt.tight_layout(rect=[0, 0, 0.88, 0.96])

plt.savefig(figure_file, dpi=dpi, bbox_inches="tight", format="tif")
plt.show()

print("Figure saved as:", figure_file)
