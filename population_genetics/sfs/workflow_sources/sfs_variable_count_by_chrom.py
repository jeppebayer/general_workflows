import gzip
import csv
import sys
from collections import defaultdict

def parse_vcf(vcf_file, population_id, output1, output2, max_count=100):
    allele_counts = []
    summary = defaultdict(lambda: defaultdict(int))  # Nested dictionary to hold summary counts per chromosome
    with gzip.open(vcf_file, 'rt') if vcf_file.endswith('.gz') else open(vcf_file, 'r') as f:
        with open(output1, 'w', newline='') as out1:
            writer1 = csv.writer(out1, delimiter='\t')
            writer1.writerow(["population", "chromosome", "position", "reference_count", "alternative_count", "invariable", "fold_count"])

            for line in f:
                if line.startswith('#'):
                    if line.startswith('#CHROM'):
                        headers = line.strip().split('\t')
                        pop_idx = headers.index(population_id)  # Find index of selected population
                    continue

                cols = line.strip().split('\t')
                chromosome = cols[0].replace(".1", "")  # Remove the redundant ".1" suffix
                position = cols[1]
                genotype_info = cols[pop_idx]

                # Split genotype data and extract the allele counts (assuming format like 0/0/.../1:...:272,2)
                alleles_part = genotype_info.split(':')[0]
                alleles = alleles_part.split('/')

                # Count the number of references (0) and alternatives (1)
                ref_count = alleles.count('0')
                alt_count = alleles.count('1')

                # Check if the site is invariable (all alleles are 0 or all alleles are 1)
                invariable = ref_count == 0 or alt_count == 0

                # Calculate fold count (smaller value between ref_count and alt_count)
                fold_count = min(ref_count, alt_count)

                # Write to output1
                writer1.writerow([population_id, chromosome, position, ref_count, alt_count, str(invariable), fold_count])

                # Update summary counts per chromosome
                summary[chromosome][f"alt_count_{alt_count}"] += 1
                summary[chromosome][f"fold_count_{fold_count}"] += 1

    # Write detailed summary to output2
    with open(output2, 'w', newline='') as out2:
        writer2 = csv.writer(out2, delimiter='\t')
        writer2.writerow(["chromosome", "count_type", "count_number", "number_of_sites"])

        # Iterate through chromosomes and all combinations of count_type and count_number
        for chromosome in summary:
            for count_type in ["alt_count", "fold_count"]:
                for count_number in range(max_count + 1):
                    # Generate the key based on count_type and count_number
                    key = f"{count_type}_{count_number}"
                    count_value = summary[chromosome].get(key, 0)  # Get count value or 0 if not present
                    writer2.writerow([chromosome, count_type, count_number, count_value])

        # Optional: Add invariable site counts if needed (not included in the main logic based on requirements)

# Example usage:
vcf_file = sys.argv[1]  # "your_file.vcf"  # Provide the path to your VCF file
population_id = sys.argv[2]  # "pop1"  # Specify the population ID
output1 = sys.argv[3]  # f"{population_id}_variants_count.tsv"  # Path to output file 1
output2 = sys.argv[4]  # f"output2_summary.tsv"  # Path to output file 2

parse_vcf(vcf_file, population_id, output1, output2)

print(f"Output written to {output1} and {output2}")

