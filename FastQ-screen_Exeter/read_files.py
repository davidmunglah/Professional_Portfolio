import pandas as pd
import os
import seaborn as sns
import matplotlib.pyplot as plt

# Define the path where the FastQ-Screen files are stored
directory_path = '/Users/davidmunglah/Documents/Exeter_Bioinformatician/Bioinf_OCT24'
output_dir = '/Users/davidmunglah/Professional_Portfolio/results'

def parse_fastq_screen(file_path):
    summary = {}
    with open(file_path, 'r') as file:
        lines = file.readlines()
        for idx, line in enumerate(lines[2:], start=3):  # Start from line 3 to skip the header
            columns = line.strip().split('\t')
            if len(columns) < 4:
                print(f"Warning: Skipping line {idx} due to unexpected format: {line}")
                continue
            try:
                genome = columns[0]
                reads_processed = int(columns[1])
                unmapped = int(columns[2])
                percent_unmapped = float(columns[3])

                mapped_reads = reads_processed - unmapped
                percent_mapped = 100 - percent_unmapped
                summary[genome] = f"{mapped_reads} ({percent_mapped:.2f}%)"
            except ValueError as e:
                print(f"Warning: Skipping line {idx} due to error: {e}")
                continue
    return summary


def generate_summary_table(directory_path):
    files = [f for f in os.listdir(directory_path) if f.endswith('_fastp_screen.txt')]
    all_data = {}
    for file in files:
        sample_name = file.split('_R1_001_fastp_screen.txt')[0]
        file_path = os.path.join(directory_path, file)
        summary = parse_fastq_screen(file_path)
        all_data[sample_name] = summary
    df = pd.DataFrame(all_data).T
    df.fillna('0 (0%)', inplace=True)
    return df

# Save summary table to CSV
summary_df = generate_summary_table(directory_path)
summary_df.to_csv(os.path.join(output_dir, 'contamination_summary.csv'))

def extract_percentages(df):
    return df.applymap(lambda x: float(x.split('(')[-1].replace('%', '').replace(')', '')))

def plot_genome_bars(df, output_dir):
    percentage_df = extract_percentages(df)
    for genome in percentage_df.columns:
        plt.figure(figsize=(8, 6))
        sns.barplot(x=percentage_df.index, y=percentage_df[genome], palette='Blues_d')
        plt.xticks(rotation=45, ha='right')
        plt.title(f'Percentage of Reads Mapped to {genome}')
        plt.ylabel('Percentage of Mapped Reads')
        plt.xlabel('Sample Name')
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, f'barplot_{genome}.png'), dpi=300)
        plt.close()

# Generate and save bar plots for each genome
plot_genome_bars(summary_df, output_dir)

# Output the summary to view
print(summary_df)

# Save the percentage values for long-term trend tracking
percentage_df = extract_percentages(summary_df)
percentage_df.to_csv(os.path.join(output_dir, 'contamination_percentages.csv'))

# Plot all genomes together in a single plot (clear with legend)
def plot_combined_bars(df, output_dir):
    percentage_df = extract_percentages(df)
    plt.figure(figsize=(12, 8))
    bar_width = 0.1  # Adjust bar width
    indices = range(len(percentage_df))

    for i, genome in enumerate(percentage_df.columns):
        plt.bar([index + bar_width * i for index in indices], percentage_df[genome],
                width=bar_width, label=genome)
    
    plt.xticks([r + bar_width * len(percentage_df.columns) / 2 for r in indices],
               percentage_df.index, rotation=45, ha='right')
    
    plt.ylabel('Percentage of Mapped Reads')
    plt.title('Genome Mapping Percentage per Sample')
    plt.legend(loc='best', bbox_to_anchor=(1.05, 1), borderaxespad=0.)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'combined_genome_mapping.png'), dpi=300)
    plt.show()

# Generate combined plot with legend
plot_combined_bars(summary_df, output_dir)

