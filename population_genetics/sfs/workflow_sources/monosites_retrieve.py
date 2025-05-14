import logging
import os
import bisect
import pysam
import pandas as pd
from collections import defaultdict

# Configure the logger
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
    handlers=[logging.StreamHandler()]
)

def check_and_prepare_vcf(vcf_file):
    """
    Check if the VCF file is gzipped and indexed. If not, compress and index it.
    """
    gzipped_vcf = f"{vcf_file}.gz"
    tbi_file = f"{gzipped_vcf}.tbi"

    if not os.path.exists(gzipped_vcf):
        logging.info(f"Compressing VCF file: {vcf_file}")
        pysam.tabix_compress(vcf_file, gzipped_vcf, force=True)
        logging.info(f"Compressed VCF saved to {gzipped_vcf}")

    if not os.path.exists(tbi_file):
        logging.info(f"Indexing gzipped VCF file: {gzipped_vcf}")
        pysam.tabix_index(gzipped_vcf, preset="vcf", force=True)
        logging.info(f"Index created for {gzipped_vcf}")

    return gzipped_vcf

def merge_bed_regions(bed_df):
    """
    Merge overlapping or adjacent BED regions to minimize fetch calls.
    """
    merged = []
    current_chrom, current_start, current_end = None, None, None

    for _, row in bed_df.iterrows():
        chrom, start, end = row["chrom"], row["start"], row["end"]

        if current_chrom is None or chrom != current_chrom or start > current_end:
            if current_chrom is not None:
                merged.append([current_chrom, current_start, current_end])
            current_chrom, current_start, current_end = chrom, start, end
        else:
            current_end = max(current_end, end)

    if current_chrom is not None:
        merged.append([current_chrom, current_start, current_end])

    return pd.DataFrame(merged, columns=["chrom", "start", "end"])

def get_variants_in_region(variant_positions, start, end):
    """
    Use binary search to retrieve the list of variants within a region [start, end).
    Returns a list of positions and the count of variants in the region.
    """
    left = bisect.bisect_left(variant_positions, start)
    right = bisect.bisect_left(variant_positions, end)
    return variant_positions[left:right], right - left

def subtract_variants_from_region(region_variants, chrom, start, end):
    """
    Subtract variants from a BED region [start, end) using the pre-retrieved list of variant positions.
    Returns the updated list of regions after removing variants.
    """
    new_regions = []
    current_start = start

    for pos in region_variants:
        if current_start < pos:
            new_regions.append(f"{chrom}\t{current_start}\t{pos}\n")
        current_start = pos + 1

    if current_start < end:
        new_regions.append(f"{chrom}\t{current_start}\t{end}\n")

    return new_regions

def generate_summary_and_subtracted_bed_optimized(bed_file, vcf_file, output_summary_file, output_bed_file):
    """
    Process a BED file and VCF file to generate summary and subtracted BED regions efficiently.
    """
    logging.info("Reading and merging BED file regions...")
    bed_df = pd.read_csv(bed_file, sep="\t", header=None, names=["chrom", "start", "end"])
    #bed_df = merge_bed_regions(bed_df)
    bed_df["region_size"] = bed_df["end"] - bed_df["start"]
    logging.info(f"Merged BED file contains {len(bed_df)} regions.")

    logging.info(f"Opening VCF file: {vcf_file}")
    vcf = pysam.TabixFile(vcf_file)

    with open(output_summary_file, "w") as summary_file, open(output_bed_file, "w") as bed_file:
        summary_file.write("chromosome\ttotal_region_size\tvariants_in_bed\tvariants_outside_bed\n")

        grouped_bed = bed_df.groupby("chrom")

        for chrom, regions in grouped_bed:
            logging.info(f"Processing chromosome: {chrom}")
            if chrom not in vcf.contigs:
                logging.warning(f"Chromosome {chrom} not found in VCF. Skipping.")
                continue

            logging.info(f"Fetching variants for chromosome: {chrom}")
            variant_positions = sorted(int(record.split("\t")[1]) for record in vcf.fetch(chrom))
            logging.info(f"Fetched {len(variant_positions)} unique variants for chromosome: {chrom}")

            total_region_size = 0
            variants_in_bed = 0
            new_regions = []

            for _, row in regions.iterrows():
                start, end = row["start"], row["end"]
                total_region_size += end - start

                # Use binary search once to get variants and their count in the region
                region_variants, in_region_count = get_variants_in_region(variant_positions, start, end)
                variants_in_bed += in_region_count

                # Generate subtracted BED regions using the retrieved variant list
                new_regions.extend(subtract_variants_from_region(region_variants, chrom, start, end))

            variants_outside_bed = len(variant_positions) - variants_in_bed
            summary_file.write(f"{chrom}\t{total_region_size}\t{variants_in_bed}\t{variants_outside_bed}\n")
            bed_file.writelines(new_regions)
            logging.info(f"Finished processing chromosome: {chrom}")

if __name__ == "__main__":
    import sys
    bed_file_path = sys.argv[1]
    vcf_file_path = sys.argv[2]
    output_summary_path = sys.argv[3]
    output_bed_path = sys.argv[4]
    gzipped_vcf_path = check_and_prepare_vcf(vcf_file_path)
    generate_summary_and_subtracted_bed_optimized(bed_file_path, gzipped_vcf_path, output_summary_path, output_bed_path)

