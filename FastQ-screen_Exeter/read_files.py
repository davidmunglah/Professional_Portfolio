import pandas as pd
import os

# Define the path where the FastQ-Screen files are stored
directory_path = '/Users/davidmunglah/Documents/Exeter_Bioinformatician/Bioinf_OCT24'

def parse_fastq_screen(file_path):
    """
    Parses a FastQ-Screen output file and extracts relevant mapping statistics.

    Args:
    - file_path (str): The path to the FastQ-Screen output file.

    Returns:
    - dict: A dictionary containing genome mapping statistics with genome names as keys
            and a tuple of (mapped reads, percentage mapped) as values.
    """
    # Initialize an empty dictionary to store results
    summary = {}

    # Read the file and skip comment lines (denoted by '#')
    with open(file_path, 'r') as file:
        lines = file.readlines()
        for idx, line in enumerate(lines[2:], start=3):  # Start from line 3 to skip the header
            columns = line.strip().split('\t')

            # Ensure there are enough columns in the line to process
            if len(columns) < 4:
                print(f"Warning: Skipping line {idx} due to unexpected format in file {file_path}: {line}")
                continue  # Skip this line if it doesn't have the expected number of columns

            try:
                genome = columns[0]
                reads_processed = int(columns[1])  # Total number of reads processed
                unmapped = int(columns[2])  # Number of unmapped reads
                percent_unmapped = float(columns[3])  # Percentage of unmapped reads

                # Calculate number of mapped reads and percentage mapped
                mapped_reads = reads_processed - unmapped
                percent_mapped = 100 - percent_unmapped

                # Store the result as a string in the format "mapped_reads (percent_mapped%)"
                summary[genome] = f"{mapped_reads} ({percent_mapped:.2f}%)"
            except ValueError as e:
                print(f"Warning: Skipping line {idx} due to data conversion error in file {file_path}: {line} - Error: {e}")
                continue  # Skip lines where data conversion fails
    
    return summary

def generate_summary_table(directory_path):
    """
    Generates a summary table for all FastQ-Screen files in a given directory.

    Args:
    - directory_path (str): Path to the directory containing FastQ-Screen files.

    Returns:
    - DataFrame: A Pandas DataFrame with samples as rows, genomes as columns, and 
                "reads (percentage mapped)" as the cell values.
    """
    # List all files in the directory that match the naming pattern
    files = [f for f in os.listdir(directory_path) if f.endswith('_fastp_screen.txt')]

    # Initialize an empty DataFrame to store the final summary
    all_data = {}

    for file in files:
        # Get the sample name from the filename (strip off the file extension)
        sample_name = file.split('_R1_001_fastp_screen.txt')[0]

        # Parse the file and get the mapping summary
        file_path = os.path.join(directory_path, file)
        summary = parse_fastq_screen(file_path)

        # Add the summary to the final table using the sample name as the key
        all_data[sample_name] = summary

    # Convert the dictionary of results into a Pandas DataFrame
    df = pd.DataFrame(all_data).T  # Transpose to get samples as rows, genomes as columns

    # Fill missing values with '0 (0%)' if any genome is missing for a sample
    df.fillna('0 (0%)', inplace=True)

    return df

# Generate the summary table
summary_df = generate_summary_table(directory_path)

# Save the table to a CSV file
output_file = '/Users/davidmunglah/Professional_Portfolio/results/contamination_summary.csv'
summary_df.to_csv(output_file)

# Output the summary to view
print(summary_df)