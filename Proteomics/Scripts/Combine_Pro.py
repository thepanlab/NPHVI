import os
import pandas as pd
from functools import reduce

# Directory with protein files
directory = "/ourdisk/hpc/rnafold/dywang/dont_archive/On_process/V45/Sidbers/Proteome/Sindbis_Pro/Protein"

# Standard mapping: filename keyword → final column name
sample_mapping = {
    "ctrl1": "PSMs_ctrl1",
    "ctrl2": "PSMs_ctrl2",
    "ctrl3": "PSMs_ctrl3",
    "S1": "PSMs_S1",
    "S2": "PSMs_S2",
    "S3": "PSMs_S3"
}

# List to hold processed DataFrames
dfs = []

# Process each .txt file in directory
for filename in os.listdir(directory):
    if filename.endswith(".txt"):
        filepath = os.path.join(directory, filename)

        # Load only needed columns
        df = pd.read_csv(filepath, sep='\t', quotechar='"',
                         usecols=["Description", "Number of PSMs"])

        # Match filename to mapping
        sample_name = os.path.splitext(filename)[0]
        new_colname = None
        for key, mapped in sample_mapping.items():
            if key in sample_name:
                new_colname = mapped
                break

        if new_colname is None:
            raise ValueError(f"Filename {filename} not recognized (needs ctrl1/2/3 or S1/2/3).")

        # Rename
        df.rename(columns={"Number of PSMs": new_colname}, inplace=True)
        dfs.append(df)

# Merge on Description
combined_df = reduce(lambda l, r: pd.merge(l, r, on="Description", how="outer"), dfs)

# Fill missing with 0
combined_df.fillna(0, inplace=True)

# Ensure all standard columns exist (even if missing in input)
for col in sample_mapping.values():
    if col not in combined_df.columns:
        combined_df[col] = 0

# Save combined
combined_output = os.path.join(directory, "/ourdisk/hpc/rnafold/dywang/dont_archive/On_process/V45/Sidbers/Proteome/Sindbis_Pro/combined_protein_psms.tsv")
combined_df.to_csv(combined_output, sep='\t', index=False)
print(f"Combined file saved to: {combined_output}")

# --- Human vs Virus split ---
def is_human(desc):
    return isinstance(desc, str) and "ENSG" in desc and "ENST" in desc

human_df = combined_df[combined_df["Description"].apply(is_human)].copy()
virus_df = combined_df[~combined_df["Description"].apply(is_human)].copy()

# Columns
control_cols = ["PSMs_ctrl1", "PSMs_ctrl2", "PSMs_ctrl3"]
treatment_cols = ["PSMs_S1", "PSMs_S2", "PSMs_S3"]

# Filter: keep if ≥2 non-zero in control or treatment
def keep_row(row):
    control_nonzero = sum(row[c] != 0 for c in control_cols if c in row)
    treatment_nonzero = sum(row[c] != 0 for c in treatment_cols if c in row)
    return control_nonzero >= 2 or treatment_nonzero >= 2

human_df = human_df[human_df.apply(keep_row, axis=1)].copy()

# Normalize columns so each sums to mean of sums
psm_cols = control_cols + treatment_cols
col_sums = human_df[psm_cols].sum()
mean_total = col_sums.mean()

for col in psm_cols:
    if col_sums[col] != 0:
        human_df[col] = human_df[col] * (mean_total / col_sums[col])
    else:
        human_df[col] = 0

# Round to integer
human_df[psm_cols] = human_df[psm_cols].round().astype(int)

# --- Post-adjustment to fix rounding drift ---
final_sums = human_df[psm_cols].sum()
target_sum = int(round(mean_total))  # common target for all columns

for col in psm_cols:
    diff = target_sum - final_sums[col]
    if diff > 0:
        # add +1 to top 'diff' rows with largest values
        idx = human_df[col].nlargest(diff).index
        human_df.loc[idx, col] += 1
    elif diff < 0:
        # subtract -1 from top 'abs(diff)' rows with largest values
        idx = human_df[col].nlargest(abs(diff)).index
        human_df.loc[idx, col] -= 1

# Reorder: Description, controls, treatments
ordered_cols = ["Description"] + control_cols + treatment_cols
human_df = human_df[ordered_cols]

# Save outputs
human_output = os.path.join(directory, "/ourdisk/hpc/rnafold/dywang/dont_archive/On_process/V45/Sidbers/Proteome/Sindbis_Pro/combined_human_protein_psms.tsv")
human_df.to_csv(human_output, sep='\t', index=False)
print(f"Filtered & normalized human protein file saved to: {human_output}")

virus_output = os.path.join(directory, "/ourdisk/hpc/rnafold/dywang/dont_archive/On_process/V45/Sidbers/Proteome/Sindbis_Pro/combined_virus_protein_psms.tsv")
virus_df.to_csv(virus_output, sep='\t', index=False)
print(f"Virus protein file saved to: {virus_output}")
