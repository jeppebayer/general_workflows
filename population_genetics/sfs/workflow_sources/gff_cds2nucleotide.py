import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Data import CodonTable
import csv

# Load the genome sequence from the FASTA file
def load_genome(fasta_file):
    genome_dict = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))
    return genome_dict

# Parse the GFF file and extract CDS regions by gene
def parse_gff(gff_file):
    cds_regions = {}
    with open(gff_file, 'r') as file:
        for line in file:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if fields[2] == 'CDS':
                chrom = fields[0]
                start = int(fields[3]) - 1  # Convert to 0-based index
                end = int(fields[4])
                strand = fields[6]
                attributes = fields[8]
                
                # Extract gene ID
                gene_id = None
                for attr in attributes.split(';'):
                    if attr.startswith("ID="):
                        gene_id = attr.split('=')[1]
                        break
                
                if gene_id:
                    if gene_id not in cds_regions:
                        cds_regions[gene_id] = {'chrom': chrom, 'strand': strand, 'regions': []}
                    cds_regions[gene_id]['regions'].append((start, end))
    return cds_regions

# Determine if mutations at a given nucleotide position would result in a synonymous or missense mutation
def calculate_mutation_probabilities(codon, codon_position):
    original_aa = codon.translate()
    synonymous_count = 0
    missense_count = 0

    # Mutation possibilities
    bases = ['A', 'T', 'C', 'G']
    original_base = codon[codon_position]
    
    for base in bases:
        if base != original_base:
            mutated_codon = codon[:codon_position] + base + codon[codon_position + 1:]
            mutated_aa = mutated_codon.translate()
            if mutated_aa == original_aa:
                synonymous_count += 1
            else:
                missense_count += 1

    total_mutations = synonymous_count + missense_count
    synonymous_prob = synonymous_count / total_mutations if total_mutations > 0 else 0
    missense_prob = missense_count / total_mutations if total_mutations > 0 else 0
    return missense_prob, synonymous_prob

# Retrieve and format the CDS sequences, along with codon and position information
def retrieve_cds_sequences(cds_regions, genome_dict, output_file):
    with open(output_file, 'w', newline='') as out_csv:
        writer = csv.writer(out_csv, delimiter='\t')
        writer.writerow(['Chromosome', 'Genome Position', 'Strand', 'Nucleotide', 'Codon', 'Gene', 'Codon Position', 'Missense Probability', 'Synonymous Probability'])
        
        for gene_id, info in cds_regions.items():
            chrom = info['chrom']
            strand = info['strand']
            regions = sorted(info['regions'], key=lambda x: x[0])  # Sort by start position
            full_cds_seq = ""

            # Build the full CDS sequence
            for start, end in regions:
                full_cds_seq += genome_dict[chrom].seq[start:end].upper()  # Ensure uppercase

            if strand == '-':
                full_cds_seq = full_cds_seq.reverse_complement().upper()

            # Divide the CDS sequence into codons and write each nucleotide with its codon and position
            for codon_start in range(0, len(full_cds_seq), 3):
                codon = full_cds_seq[codon_start:codon_start+3]
                if len(codon) == 3:
                    for i in range(3):
                        genome_pos = regions[0][0] + codon_start + i + 1 if strand == '+' else regions[-1][1] - codon_start - i
                        codon_position = i  # Position within the codon (0, 1, or 2)
                        missense_prob, synonymous_prob = calculate_mutation_probabilities(codon, i)
                        writer.writerow([chrom, genome_pos, strand, codon[i], codon, gene_id, codon_position, missense_prob, synonymous_prob])

# Define file paths
fasta_file = sys.argv[1] #"path/to/genome.fasta"
gff_file = sys.argv[2] #"path/to/annotations.gff"
output_file = sys.argv[3]#"path/to/output_file.tsv"

# Execute the function calls
genome_dict = load_genome(fasta_file)
cds_regions = parse_gff(gff_file)
retrieve_cds_sequences(cds_regions, genome_dict, output_file)

print("Output written to", output_file)

