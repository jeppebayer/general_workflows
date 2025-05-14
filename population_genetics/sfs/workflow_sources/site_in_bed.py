import sys
def filter_sites_within_bed_regions(bed_file, site_file, output_file):
    with open(bed_file, 'r') as bed, open(site_file, 'r') as sites, open(output_file, 'w') as output:
        # Read and write the header line from the site file
        header = sites.readline().strip()
        output.write(header + '\n')
        
        # Read the first region from the BED file
        bed_line = bed.readline()
        if not bed_line:
            return  # BED file is empty
        
        bed_fields = bed_line.strip().split()
        bed_chr, bed_start, bed_end = bed_fields[0], int(bed_fields[1]) + 1, int(bed_fields[2])  # Adjust start to 1-based

        for site_line in sites:
            site_fields = site_line.strip().split()
            site_chr, site_pos = site_fields[0], int(site_fields[1])
            
            # Move to the next region in the BED file until the current region is relevant
            while site_chr > bed_chr or (site_chr == bed_chr and site_pos > bed_end):
                bed_line = bed.readline()
                if not bed_line:
                    return  # Reached the end of the BED file
                bed_fields = bed_line.strip().split()
                bed_chr, bed_start, bed_end = bed_fields[0], int(bed_fields[1]) + 1, int(bed_fields[2])  # Adjust start to 1-based
            
            # Check if the site is within the current BED region
            if site_chr == bed_chr and bed_start <= site_pos <= bed_end:
                output.write(site_line)

# Usage
bed_file = sys.argv[1] #"regions.bed"  # Replace with the path to your BED file
site_file = sys.argv[2] #"sites.txt"  # Replace with the path to your site file
output_file = sys.argv[3] #"filtered_sites.txt"  # Replace with the desired output file path

filter_sites_within_bed_regions(bed_file, site_file, output_file)

