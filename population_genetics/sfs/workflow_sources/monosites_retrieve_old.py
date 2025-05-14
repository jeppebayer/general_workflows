import sys
import pysam
import pandas as pd
import numpy as np
from collections import defaultdict

def compress_and_index_vcf(vcf_file, compressed_vcf_file):
    """
    Compress and index a VCF file using pysam.tabix_compress and pysam.tabix_index.
    """
    # Compress the VCF file with tabix_compress
    pysam.tabix_compress(vcf_file, compressed_vcf_file, force=True)
    print(f"Compressed VCF saved to {compressed_vcf_file}")
    
    # Index the compressed VCF file with tabix_index
    pysam.tabix_index(compressed_vcf_file, preset="vcf", force=True)
    print(f"Index created for {compressed_vcf_file}")

def query_vcf_by_region(vcf, chrom, start, end):
    """
    Query a bgzipped and indexed VCF file for a specific region.
    """
    return list(vcf.fetch(chrom, start, end))

def generate_summary_and_subtracted_bed(bed_file, vcf_file, output_summary_file, output_bed_file):
    """
    Process a BED file and VCF file to generate summary and subtracted BED regions.
    """
    # Read BED file and group by chromosome
    bed_df = pd.read_csv(bed_file, sep="\t", header=None, names=["chrom", "start", "end"])
    bed_df["region_size"] = bed_df["end"] - bed_df["start"]

    # Preload BED regions into a dictionary keyed by chromosome
    bed_regions = defaultdict(list)
    for _, row in bed_df.iterrows():
        bed_regions[row["chrom"]].append((row["start"], row["end"]))

    # Convert lists of regions to NumPy arrays for vectorized operations
    for chrom in bed_regions:
        bed_regions[chrom] = np.array(bed_regions[chrom])

    # Open the compressed VCF file with Tabix
    vcf = pysam.TabixFile(vcf_file)
    variant_counts = defaultdict(lambda: {"total_variants": 0, "variants_in_bed": 0})

    # Process each chromosome in the BED file
    new_regions = []
    for chrom, regions in bed_regions.items():
        if chrom not in vcf.contigs:
            continue

        # Process variants in the chromosome using Tabix fetch
        variant_positions = []
        for start, end in regions:
            for record in vcf.fetch(chrom, start, end):
                variant_positions.append(int(record.split("\t")[1]))  # Extract position

        # Update counts for summary
        variant_positions = sorted(set(variant_positions))
        for start, end in regions:
            in_region_count = sum(start <= pos < end for pos in variant_positions)
            variant_counts[chrom]["variants_in_bed"] += in_region_count
            variant_counts[chrom]["total_variants"] += len(variant_positions)

            # Subtract variants from BED regions
            current_start = start
            for pos in variant_positions:
                if start <= pos < end:
                    if current_start < pos:
                        new_regions.append([chrom, current_start, pos])
                    current_start = pos + 1
            if current_start < end:
                new_regions.append([chrom, current_start, end])

    # Prepare summary data
    summary_data = []
    for chrom in bed_regions.keys():
        total_region_size = bed_df[bed_df["chrom"] == chrom]["region_size"].sum()
        variants_in_bed = variant_counts[chrom]["variants_in_bed"]
        variants_outside_bed = variant_counts[chrom]["total_variants"] - variants_in_bed
        summary_data.append([chrom, total_region_size, variants_in_bed, variants_outside_bed])

    # Save the summary file
    summary_df = pd.DataFrame(summary_data, columns=["chromosome", "total_region_size", "variants_in_bed", "variants_outside_bed"])
    summary_df.to_csv(output_summary_file, sep="\t", index=False)
    print(f"Summary file saved to {output_summary_file}")

    # Save the new BED file
    new_bed_df = pd.DataFrame(new_regions, columns=["chrom", "start", "end"])
    new_bed_df.to_csv(output_bed_file, sep="\t", header=False, index=False)
    print(f"Subtracted BED file saved to {output_bed_file}")

# File paths (update these with your actual file paths)
original_vcf_file = sys.argv[2]# "your_vcf_file.vcf"  # Uncompressed VCF file
compressed_vcf_file = sys.argv[2]+".gz"#"your_vcf_file.vcf.gz"  # Compressed VCF file
bed_file_path = sys.argv[1]#"your_bed_file.bed"  # BED file
output_summary_path = sys.argv[3]#"output_summary.tsv"  # Output summary TSV file
output_bed_path = sys.argv[4]#"subtracted_bed_file.bed"  # Subtracted BED file

# Compress and index the VCF file
compress_and_index_vcf(original_vcf_file, compressed_vcf_file)

# Generate the summary and subtracted BED file
generate_summary_and_subtracted_bed(bed_file_path, compressed_vcf_file, output_summary_path, output_bed_path)

