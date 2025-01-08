import pandas as pd
import matplotlib.pyplot as plt
import os

# Load the data
file_path = '/ourdisk/hpc/rnafold/dywang/dont_archive/On_process/V45/Sidbers/Nanopore/GeneExpression/Deseq/Top_10_Unique_Pathways_1.csv'  # Update with your actual file path
data = pd.read_csv(file_path)

# Prepare data for plotting
grouped_data = data.groupby('Pathway')['log2FoldChange'].mean().sort_values()

# Separate data into increased and decreased fold changes
# Sort increased values in ascending order so they appear at the top in the plot
increased_values = grouped_data[grouped_data > 0].sort_values()
decreased_values = grouped_data[grouped_data <= 0].sort_values(ascending=False)

# Set up the figure size and other plot parameters
fig, ax = plt.subplots(figsize=(4, 4))  # Adjust the figure size here

# Plot negative values in blue with larger markers and thicker lines (decreased at the bottom)
ax.hlines(y=decreased_values.index, xmin=0, xmax=decreased_values, color='#6699FF', linewidth=3)
ax.plot(decreased_values, decreased_values.index, 'o', color='#6699FF', markersize=8)

# Plot positive values in red with larger markers and thicker lines (increased at the top)
ax.hlines(y=increased_values.index, xmin=0, xmax=increased_values, color='#FF6666', linewidth=3)
ax.plot(increased_values, increased_values.index, 'o', color='#FF6666', markersize=8)

# Labels and grid
ax.set_xlabel("log2FoldChange")
ax.set_ylabel("Pathway")
plt.title("Log2 Fold Change by GO Description (Increased at Top, Decreased at Bottom)")
plt.grid(axis='x', linestyle='--', alpha=0.7)

# Specify the output directory and filename
output_dir = "/ourdisk/hpc/rnafold/dywang/dont_archive/On_process/V45/Sidbers/Nanopore/GeneExpression/Deseq/"  # Update with your desired output directory
output_filename = "Top10_KO.tif"

# Ensure the output directory exists
os.makedirs(output_dir, exist_ok=True)

# Save the figure as a .tif file with 300 dpi
output_path = os.path.join(output_dir, output_filename)
plt.savefig(output_path, format="tif", dpi=300)

plt.show()
