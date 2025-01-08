import glob
import pandas as pd

# Specify the path to your directory
path = '/ourdisk/hpc/rnafold/dywang/dont_archive/On_process/V45/Sidbers/Proteome/Abundance/'  # Updated for the provided example file location

# Get a list of all files in the directory
files = glob.glob(path + '*.txt')
files.sort()

# Function to split records with multiple IDs and assign average totalPeakHeight
def split_and_average(df):
    rows = []
    for _, row in df.iterrows():
        ids = row['locus']
        if isinstance(ids, str) and ids.startswith('{') and ids.endswith('}'):
            ids = ids[1:-1].split(',')
            num_ids = len(ids)
            avg_peak_height = row['totalPeakHeight'] / num_ids
            for single_id in ids:
                new_row = row.copy()
                new_row['locus'] = single_id.strip()
                new_row['totalPeakHeight'] = avg_peak_height
                rows.append(new_row)
        else:
            rows.append(row)
    return pd.DataFrame(rows)

# Initialize final dataframe with columns for the desired output
final_df = pd.DataFrame(columns=['locus', 'description'])

# Process each file and combine
for file in files:
    # Read the file into a DataFrame
    df = pd.read_csv(file, sep='\t', comment='#')
    
    # Process to split records with multiple IDs and assign average totalPeakHeight
    processed_df = split_and_average(df)
    
    # Summarize totalPeakHeight
    summed_df = processed_df.groupby('locus').agg({
        'totalPeakHeight': 'sum',
        'description': 'first'  # Assume 'description' is the same for the same 'locus'
    }).reset_index()

    # Get identifier from filename
    parts = file.split("/")[-1].split(".")[0].split('_')
    identifier = parts[1] + '_' + parts[2]

    # Rename the 'totalPeakHeight' column
    summed_df.rename(columns={'totalPeakHeight': f'total_Peak_Height_{identifier}'}, inplace=True)

    # Merge the dataframes on locus and description
    if final_df.empty:
        final_df = summed_df
    else:
        final_df = pd.merge(final_df, summed_df, on=['locus', 'description'], how='outer')

# Replace NA values with 0
final_df = final_df.fillna(0)

# Ensure 'description' is the second column
columns = final_df.columns.tolist()
columns.remove('description')
columns.insert(1, 'description')
final_df = final_df[columns]

# Save the final DataFrame to a new text file with tab-separated values
output_file_path_final = '/ourdisk/hpc/rnafold/dywang/dont_archive/On_process/V45/Sidbers/Proteome/combined_proteome.csv'
final_df.to_csv(output_file_path_final, sep='\t', index=False)

print(f"Combined data saved to: {output_file_path_final}")
