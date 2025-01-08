import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import pearsonr

# Step 1: Load your datasets
transcript_data_path = '/ourdisk/hpc/rnafold/dywang/dont_archive/On_process/V45/Sidbers/Visualization/Figures/Correlation/Nanopore_Sidbers_TPM.csv'
proteome_data_path = '/ourdisk/hpc/rnafold/dywang/dont_archive/On_process/V45/Sidbers/Visualization/Figures/Correlation/m6A_consolidate_filtered.csv'

transcript_df = pd.read_csv(transcript_data_path)
proteome_df = pd.read_csv(proteome_data_path)

# Step 2: Preprocess the ID columns to extract the common ENSG ID
def extract_ensg_ids(id_column, delimiter="|"):
    return id_column.str.extract(r'(ENSG\d+)')[0]

transcript_df['Processed_ID'] = extract_ensg_ids(transcript_df['ID'])
proteome_df['Processed_ID'] = extract_ensg_ids(proteome_df['ID'])

# Step 3: Merge datasets based on ENSG IDs
merged_data = pd.merge(transcript_df, proteome_df, on='Processed_ID', suffixes=('_transcript', '_proteome'))

# Step 4: Select the relevant transcript and proteome columns for analysis
transcript_cols = ['ctrl1_transcript', 'ctrl2_transcript', 'ctrl3_transcript', 'S1_transcript', 'S2_transcript', 'S3_transcript']
proteome_cols = ['ctrl1_proteome', 'ctrl2_proteome', 'ctrl3_proteome', 'S1_proteome', 'S2_proteome', 'S3_proteome']

# Step 5: Log-transform the data (using log1p to handle zeros)
log_transcript_data = np.log1p(merged_data[transcript_cols].values.flatten())
log_proteome_data = np.log1p(merged_data[proteome_cols].values.flatten())

# Step 6: Calculate Pearson correlation coefficient and save to .txt file
def calculate_correlation_and_save(log_transcript_data, log_proteome_data, output_file):
    # Pearson correlation
    correlation_coefficient, p_value = pearsonr(log_transcript_data, log_proteome_data)
    
    # Write correlation results to a text file
    with open(output_file, 'w') as file:
        file.write(f"Pearson Correlation Coefficient: {correlation_coefficient}\n")
        file.write(f"P-value: {p_value}\n")
    print(f"Results saved to {output_file}")

# Step 7: Create the hexbin plot without the fitting line, and save as .tif with 300 dpi
def create_hexbin_plot_without_fit_and_save(log_transcript_data, log_proteome_data, figsize=(10, 8), linewidth=1.5, vmin=None, vmax=None, output_file='/ourdisk/hpc/rnafold/dywang/dont_archive/On_process/V45/Sidbers/Visualization/Figures/Correlation/TPM_m6A_log_transformed_no_fit.tif'):
    plt.figure(figsize=figsize)
    
    # Customizing hexbin plot with different color and line width, setting color bar range with vmin and vmax
    hb = plt.hexbin(log_transcript_data, log_proteome_data, gridsize=50, cmap='Blues', mincnt=1, edgecolors='white', linewidths=linewidth, vmin=vmin, vmax=vmax)
    
    # Customizing the color bar
    cbar = plt.colorbar(hb)
    cbar.set_label('Counts', fontsize=12)
    cbar.ax.tick_params(labelsize=10)
    
    # Add labels and title
    plt.title('Hexbin Plot (Log-Transformed Data)', fontsize=14)
    plt.xlabel('Log-Transformed TPM Data', fontsize=12)
    plt.ylabel('Log-Transformed m6A data', fontsize=12)

    # Save the figure as a .tif file with 300 dpi resolution
    plt.savefig(output_file, format='tif', dpi=300)
    plt.show()

# Define output file paths
correlation_output_file = '/ourdisk/hpc/rnafold/dywang/dont_archive/On_process/V45/Sidbers/Visualization/Figures/Correlation/TPM_m6A_log_transformed_correlation_coefficient.txt'
hexbin_plot_output_file = '/ourdisk/hpc/rnafold/dywang/dont_archive/On_process/V45/Sidbers/Visualization/Figures/Correlation/TPM_m6A_log_transformed_no_fit.tif'

# Step 8: Run the correlation calculation and save the result to the .txt file
calculate_correlation_and_save(log_proteome_data, log_transcript_data, correlation_output_file)

# Step 9: Generate and save the hexbin plot without the fitting line
create_hexbin_plot_without_fit_and_save(log_transcript_data, log_proteome_data, figsize=(5, 4), linewidth=2, vmin=0, vmax=100, output_file=hexbin_plot_output_file)
