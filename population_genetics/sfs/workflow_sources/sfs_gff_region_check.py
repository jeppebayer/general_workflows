import sys
import pandas as pd

def check_regions_length(gff_file, cds_bed_output, cds_report_output, exon_bed_output, exon_report_output):
    # Counters for CDS
    total_cds = 0
    failed_cds = 0
    
    # Counters for Exons
    total_exons = 0
    failed_exons = 0
    
    # Open output files
    with open(cds_bed_output, 'w') as cds_bed, open(cds_report_output, 'w') as cds_report, \
         open(exon_bed_output, 'w') as exon_bed, open(exon_report_output, 'w') as exon_report:
        
        # Open and read the GFF file
        with open(gff_file, 'r') as gff:
            for line in gff:
                if line.startswith("#"):
                    continue  # Skip comment lines
                
                fields = line.strip().split('\t')
                if len(fields) < 9:
                    continue  # Skip incomplete lines
                
                region_type = fields[2].lower()
                chrom = fields[0]
                start = int(fields[3])
                end = int(fields[4])
                length = end - start + 1
                
                # Process CDS regions
                if region_type == "cds":
                    total_cds += 1
                    if length % 3 == 0:
                        cds_bed.write(f"{chrom}\t{start - 1}\t{end}\n")  # BED format is 0-based
                    else:
                        failed_cds += 1
                        cds_report.write(f"Failed CDS: {chrom}\t{start}\t{end}\tLength: {length}\n")
                
                # Process Exon regions
                elif region_type == "exon":
                    total_exons += 1
                    if length % 3 == 0:
                        exon_bed.write(f"{chrom}\t{start - 1}\t{end}\n")  # BED format is 0-based
                    else:
                        failed_exons += 1
                        exon_report.write(f"Failed Exon: {chrom}\t{start}\t{end}\tLength: {length}\n")
        
        # Write summary information at the end of each report
        cds_report.write(f"\nTotal CDS Regions: {total_cds}\n")
        cds_report.write(f"Failed CDS Regions: {failed_cds}\n")
        
        exon_report.write(f"\nTotal Exons: {total_exons}\n")
        exon_report.write(f"Failed Exons: {failed_exons}\n")

# Example usage:
gff_file = sys.argv[1]# "path/to/your/annotation.gff"
cds_bed_output = sys.argv[2]# "cds_pass.bed"
cds_report_output = sys.argv[3]# "cds_fail_report.txt"
exon_bed_output = sys.argv[4]# "exon_pass.bed"
exon_report_output = sys.argv[5]# "exon_fail_report.txt"

check_regions_length(gff_file, cds_bed_output, cds_report_output, exon_bed_output, exon_report_output)

