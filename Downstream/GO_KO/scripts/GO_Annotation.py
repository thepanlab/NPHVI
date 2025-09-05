import pandas as pd
import requests
from goatools.obo_parser import GODag

# === STEP 1: Load your input file ===
input_file = "/ourdisk/hpc/rnafold/dywang/dont_archive/On_process/V45/Sidbers/Nanopore/GeneExpression/Deseq_new/Normalized_filtered_Sidbers_TPM.tsv"
print(f" Loading input file: {input_file}")
df = pd.read_csv(input_file, sep="\t")

# Extract Ensembl Gene IDs (ENSG) from the ID column
df['ENSG'] = df['ID'].str.extract(r'(ENSG\d+)')

# === STEP 2: Map Ensembl Gene IDs to NCBI Gene IDs using MyGene.info ===
def map_ensembl_to_ncbi(ensg_ids):
    url = "https://mygene.info/v3/query"
    mapping = {}
    batch_size = 1000
    for i in range(0, len(ensg_ids), batch_size):
        batch = ensg_ids[i:i + batch_size]
        params = {
            "q": ",".join(batch),
            "scopes": "ensembl.gene",
            "fields": "entrezgene",
            "species": "human"
        }
        res = requests.post(url, data=params)
        if res.ok:
            data = res.json()
            for item in data:
                if "entrezgene" in item:
                    mapping[item["query"]] = item["entrezgene"]
    return mapping

print(" Mapping Ensembl Gene IDs to NCBI Gene IDs using MyGene.info...")
ensg_unique = df['ENSG'].dropna().unique().tolist()
ensembl_to_ncbi = map_ensembl_to_ncbi(ensg_unique)
df['NCBI_temp'] = df['ENSG'].map(ensembl_to_ncbi)

# === STEP 3: Load gene2go and filter for human ===
gene2go_file = "/ourdisk/hpc/rnafold/dywang/dont_archive/On_process/V45/Sidbers/GO_File/gene2go"
print(f" Loading gene2go from: {gene2go_file}")
gene2go = pd.read_csv(gene2go_file, sep="\t", comment="#", header=None,
                      names=["tax_id", "GeneID", "GO_ID", "Evidence", "Qualifier", "GO_term", "PubMed", "Category"])
gene2go = gene2go[gene2go["tax_id"] == 9606]

# === STEP 4: Merge with expression data ===
print(" Merging GO annotations with expression data...")
df['NCBI_temp'] = df['NCBI_temp'].astype(str)
gene2go['GeneID'] = gene2go['GeneID'].astype(str)
merged = df.merge(gene2go, left_on="NCBI_temp", right_on="GeneID", how="left")

# === STEP 5: Fill in missing Category using go-basic.obo ===
obo_file = "/ourdisk/hpc/rnafold/dywang/dont_archive/On_process/V45/Sidbers/GO_File/go-basic.obo"
print(f" Loading ontology from: {obo_file}")
go_dag = GODag(obo_file, optional_attrs=['namespace'])

# Map GO_ID → Category
go_to_cat = {}
for go_id, term in go_dag.items():
    ns = term.namespace
    if ns == "biological_process":
        go_to_cat[go_id] = "BP"
    elif ns == "molecular_function":
        go_to_cat[go_id] = "MF"
    elif ns == "cellular_component":
        go_to_cat[go_id] = "CC"

# Use original category if available, else assign from ontology
print("🛠 Filling in missing GO categories...")
merged['Category'] = merged['Category'].map({
    "P": "BP", "F": "MF", "C": "CC",
    "p": "BP", "f": "MF", "c": "CC"
})
merged['Category'] = merged.apply(
    lambda row: go_to_cat.get(row['GO_ID'], row['Category']),
    axis=1
)

# === DEBUG: Check mapped categories ===
print(" Final Category counts:")
print(merged["Category"].value_counts(dropna=False))

# === STEP 6: Reorder columns ===
expression_cols = [col for col in df.columns if col not in ['ID', 'ENSG', 'NCBI_temp']]
final_cols = ['ID', 'ENSG', 'GO_ID', 'GO_term', 'Category'] + expression_cols
final_df = merged[final_cols]

# === STEP 7: Save full result ===
output_base = "/ourdisk/hpc/rnafold/dywang/dont_archive/On_process/V45/Sidbers/Nanopore/GeneExpression/Deseq_new/"
full_out = output_base + "Normalized_filtered_Sidbers_TPM_GO.tsv"
final_df.to_csv(full_out, sep="\t", index=False)
print(f" Full GO annotation saved to: {full_out}")

# === STEP 8: Split by Category ===
bp_df = final_df[final_df['Category'] == 'BP']
mf_df = final_df[final_df['Category'] == 'MF']
cc_df = final_df[final_df['Category'] == 'CC']

bp_df.to_csv(output_base + "Normalized_filtered_Sidbers_TPM_GO_BP.tsv", sep="\t", index=False)
mf_df.to_csv(output_base + "Normalized_filtered_Sidbers_TPM_GO_MF.tsv", sep="\t", index=False)
cc_df.to_csv(output_base + "Normalized_filtered_Sidbers_TPM_GO_CC.tsv", sep="\t", index=False)

print(" Separated GO annotations into:")
print(f"   - BP: {output_base}Normalized_filtered_Sidbers_TPM_GO_BP.tsv")
print(f"   - MF: {output_base}Normalized_filtered_Sidbers_TPM_GO_MF.tsv")
print(f"   - CC: {output_base}Normalized_filtered_Sidbers_TPM_GO_CC.tsv")
