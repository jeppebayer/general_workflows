import sys
import pysam

def find_closest_variant_in_region(vcf_file, chrom, start, end):
    """
    Fetch all variants within a region and find the closest variant for each query position.
    
    Args:
    vcf_file (str): Path to the VCF file (should be bgzipped and indexed with tabix).
    chrom (str): Chromosome ID (e.g., 'chr1').
    start (int): Start position of the region (1-based after conversion).
    end (int): End position of the region (1-based).
    
    Returns:
    list of tuples: A list of (variant position, ref, alts, info) for all variants in the region.
    """
    vcf = pysam.VariantFile(vcf_file)
    variants = []

    for rec in vcf.fetch(chrom, start, end):
        var_pos = rec.pos
        ref = rec.ref
        alts = rec.alts
        info = rec.info
        info_str = ";".join(f"{key}={value}" for key, value in info.items())
        variants.append((var_pos, ref, alts, info_str))

    return variants

def process_queries(vcf_file, regions, output_file):
    """
    Process each position in the regions from the BED file and write the closest variants to an output file,
    skipping positions that directly hit a variant.
    
    Args:
    vcf_file (str): Path to the VCF file (bgzipped and indexed with tabix).
    regions (list of tuples): List of (chrom, start, end) regions to query.
    output_file (str): Path to the output file.
    """
    with open(output_file, 'w') as out_f:
        # Write header
        out_f.write("Query_Chrom\tQuery_Position\tNearest_Variant_Position\tNearest_Variant_Ref\tNearest_Variant_Alt\tNearest_Variant_Info\tDistance\n")
        
        # Process each region
        for chrom, start, end in regions:
            # Convert the 0-based start position from BED to 1-based for VCF
            start_1based = start + 1
            
            # Fetch all variants in the current region
            variants = find_closest_variant_in_region(vcf_file, chrom, start_1based, end)
            
            # Process each query position in the region
            for pos in range(start_1based, end + 1):  # Query each 1-based position
                closest_variant = None
                closest_distance = float('inf')
                
                # Check each variant and calculate the distance
                for var_pos, ref, alts, info in variants:
                    distance = abs(var_pos - pos)
                    
                    # Skip if the position hits a variant exactly
                    if var_pos == pos:
                        closest_variant = None
                        break
                    
                    # Find the closest variant
                    if distance < closest_distance:
                        closest_distance = distance
                        closest_variant = (var_pos, ref, alts, info)
                
                # Skip writing if the query hits a variant exactly
                if closest_variant is None:
                    continue
                
                # Write the closest variant to the output file
                var_pos, ref, alts, info = closest_variant
                alts_str = ",".join(alts)
                out_f.write(f"{chrom}\t{pos}\t{var_pos}\t{ref}\t{alts_str}\t{info}\t{closest_distance}\n")

def read_bed_file(bed_file):
    """
    Read a BED file and extract chromosome, start, and end positions for queries.
    
    Args:
    bed_file (str): Path to the BED file.
    
    Returns:
    list of tuples: List of (chrom, start, end) regions from the BED file.
    """
    regions = []
    with open(bed_file, 'r') as bed_f:
        for line in bed_f:
            if line.strip():
                chrom, start, end = line.strip().split()[:3]
                regions.append((chrom, int(start), int(end)))
    return regions

# Example usage
vcf_file = sys.argv[1]#"yourfile.vcf.gz"
bed_file = sys.argv[2]#"yourqueries.bed"
output_file = sys.argv[3]# "output.tsv"

# Read the BED file to get the regions
regions = read_bed_file(bed_file)

# Process the queries and write to the output file
process_queries(vcf_file, regions, output_file)

