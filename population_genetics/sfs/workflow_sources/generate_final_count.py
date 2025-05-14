import sys

def read_monosites_count(monosites_count):
    chrom_count_dict = {}
    chrom_counts = open(monosites_count)
    for line in chrom_counts:
        infos = line.strip("\n").split("\t")
        chrom_id = infos[0].split(".")[0]
        if chrom_id != "Chromosome":
            N_count = int(float(infos[1]))
            S_count = int(float(infos[2]))
        else:
            N_count = infos[1]
            S_count = infos[2]
        chrom_count_dict[chrom_id]=(N_count,S_count)
    return chrom_count_dict

def read_mono_count_intergenic(inter_mono):
    inter_mono_count = {}
    chrom_counts = open(inter_mono)
    for line in chrom_counts:
        infos = line.strip("\n").split("\t")
        chrom_id = infos[0].split(".")[0]
        if chrom_id != "chromosome":
            total_region_size = int(infos[1])
            variant_in_region = int(infos[2])
            mono_count = total_region_size - variant_in_region
        else:
            mono_count = "header"
        inter_mono_count[chrom_id] = mono_count
    return inter_mono_count

monosites_file_gene = sys.argv[1]
monosites_file_inter = sys.argv[2]

varaint_N_summary = sys.argv[3]
variant_S_summary = sys.argv[4]
variant_inter_summary = sys.argv[5]

chrom_mono_count = read_monosites_count(monosites_file_gene)
inter_mono_count = read_mono_count_intergenic(monosites_file_inter)

N_sum = open(varaint_N_summary)
N_sum_out = open(sys.argv[6],"w")
S_sum = open(variant_S_summary)
S_sum_out = open(sys.argv[7],"w")
inter_sum = open(variant_inter_summary)
inter_sum_out = open(sys.argv[8],"w")

for line in N_sum:
    infos = line.strip("\n").split("\t")
    chrom_id = infos[0]
    count_type = infos[1]
    if chrom_id != "chromosome":
        count_number = int(infos[2])
        number_of_sites = int(infos[3])
    else:
        count_numer = infos[2]
        number_of_sites = infos[3]
    if (chrom_id in chrom_mono_count) and count_number == 0:
        sites_update = number_of_sites + chrom_mono_count[chrom_id][0]
        info_update = [chrom_id,count_type,str(count_number),str(sites_update)]
        N_sum_out.write('\t'.join(info_update) + '\n')
    else:
        N_sum_out.write(line)
for line in S_sum:
    infos = line.strip("\n").split("\t")
    chrom_id = infos[0]
    count_type = infos[1]
    if chrom_id != "chromosome":
        count_number = int(infos[2])
        number_of_sites = int(infos[3])
    else:
        count_numer = infos[2]
        number_of_sites = infos[3]
    if (chrom_id in chrom_mono_count) and count_number == 0:
        sites_update = number_of_sites + chrom_mono_count[chrom_id][1]
        info_update = [chrom_id,count_type,str(count_number),str(sites_update)]
        S_sum_out.write('\t'.join(info_update) + '\n')
    else:
        S_sum_out.write(line)

for line in inter_sum:
    infos = line.strip("\n").split("\t")
    chrom_id = infos[0]
    count_type = infos[1]
    if chrom_id != "chromosome":
        count_number = int(infos[2])
        number_of_sites = int(infos[3])
    else:
        count_numer = infos[2]
        number_of_sites = infos[3]
    if (chrom_id in chrom_mono_count) and count_number == 0:
        sites_update = number_of_sites + inter_mono_count[chrom_id]
        info_update = [chrom_id,count_type,str(count_number),str(sites_update)]
        inter_sum_out.write('\t'.join(info_update) + '\n')
    else:
        inter_sum_out.write(line)
