import pandas as pd

# Load the CSV file
file_path = '/ourdisk/hpc/rnafold/dywang/dont_archive/On_process/V45/Sidbers/Nanopore/GeneExpression/consolidated_TPM.csv'
data = pd.read_csv(file_path)

# Define the control and treated columns
control_columns = ['TPM_s1_a549', 'TPM_c1a549', 'TPM_mocka59']
treated_columns = ['TPM_c2_a549', 'TPM_sva549', 'TPM_s2_a549']

# Create a function that will convert a row's values to 0 if less than 2 are non-zero
def apply_zero_condition(row):
    if (row != 0).sum() < 2:
        return [0] * len(row)
    else:
        return row

# Apply this function to both control and treated groups
data[control_columns] = data[control_columns].apply(apply_zero_condition, axis=1)
data[treated_columns] = data[treated_columns].apply(apply_zero_condition, axis=1)

# Now, filter out rows where all values are zero in both control and treated groups
filtered_data = data[(data[control_columns].sum(axis=1) > 0) | (data[treated_columns].sum(axis=1) > 0)]

# Convert all the numeric values to integers
filtered_data_int = filtered_data.astype({col: 'int' for col in filtered_data.columns if col != 'ID'})

# Save the filtered data with integer values to a new CSV file
final_output_file_path = '/ourdisk/hpc/rnafold/dywang/dont_archive/On_process/V45/Sidbers/Nanopore/GeneExpression/filtered_Sidbers_TPM.csv'
filtered_data_int.to_csv(final_output_file_path, index=False)
