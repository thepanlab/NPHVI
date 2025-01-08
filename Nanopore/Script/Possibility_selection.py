import csv

input_filename = '/scratch/dywang/Nanopore/m6A/c1a549/m6Anet/c1a549.csv'
output_filename = '/scratch/dywang/Nanopore/m6A/selection/c1a549.csv'
threshold = 0.95

with open(input_filename, 'r', newline='') as input_file, \
        open(output_filename, 'w', newline='') as output_file:

    reader = csv.reader(input_file, delimiter=',')  # Updated delimiter to comma for CSV
    writer = csv.writer(output_file, delimiter=',')  # Updated delimiter to comma for CSV

    # Find the index of the 'probability_modified' column
    header = next(reader)
    probability_index = header.index('probability_modified')

    # Write the header to the output file
    writer.writerow(header)

    for row in reader:
        try:
            probability = float(row[probability_index])
            if probability > threshold:
                writer.writerow(row)
        except (ValueError, IndexError):
            # This handles any conversion errors or missing data, skipping problematic rows
            pass
