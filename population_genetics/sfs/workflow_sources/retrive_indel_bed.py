import gzip
import os
import sys

def open_vcf_file(file_path):
    """
    Open a VCF file, automatically detecting if it's gzipped or plain text.
    """
    try:
        with gzip.open(file_path, 'rt') as f:
            # Test if the file can be read as gzipped
            f.read(1)
        return gzip.open(file_path, 'rt')  # Return the file object if gzipped
    except (OSError, gzip.BadGzipFile):
        return open(file_path, 'r')  # Return the file object for plain text

def generate_bed_with_merged_regions(vcf_file, bed_file):
    with open_vcf_file(vcf_file) as vcf, open(bed_file, 'w') as bed:
        current_region = None  # To store the current merged region

        for line in vcf:
            # Skip header lines
            if line.startswith("#"):
                continue

            # Parse VCF line
            fields = line.strip().split('\t')
            chrom = fields[0]
            pos = int(fields[1])  # 1-based position
            info = fields[7]

            # Extract the TYPE field
            type_field = [entry.split('=')[1] for entry in info.split(';') if entry.startswith('TYPE=')]
            if not type_field:
                continue
            variant_types = type_field[0].split(',')

            # Check if "ins" or "del" is in the TYPE field
            if any(var_type in ['ins', 'del'] for var_type in variant_types):
                # Calculate region boundaries (0-based, inclusive start, exclusive end)
                start = max(0, pos - 6)  # 5bp upstream
                end = pos + 5           # 5bp downstream

                # Merge with current region if overlapping
                if current_region and current_region[0] == chrom and current_region[2] >= start:
                    # Extend the current region
                    current_region[2] = max(current_region[2], end)
                else:
                    # Write the previous region if it exists
                    if current_region:
                        bed.write(f"{current_region[0]}\t{current_region[1]}\t{current_region[2]}\n")
                    # Start a new region
                    current_region = [chrom, start, end]

        # Write the last region if it exists
        if current_region:
            bed.write(f"{current_region[0]}\t{current_region[1]}\t{current_region[2]}\n")

# Specify input VCF file and output BED file
vcf_file = sys.argv[1] #"variants.vcf.gz"  # Can be gzipped or plain text
bed_file = sys.argv[2] #"output_regions.bed"

# Generate the BED file with merged regions
generate_bed_with_merged_regions(vcf_file, bed_file)

