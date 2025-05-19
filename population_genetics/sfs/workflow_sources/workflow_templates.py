#!/bin/env python3
from gwf import AnonymousTarget
import os, glob


def annotate_N_S_sites(work_path: str, script_path: str, spid: str, genome: str, gff: str, log_path: str):
    inputs = {"genome":genome,
            "gff":gff}
    outputs = {
            "CDS_summary":f"{work_path}/{spid}_species_data/{spid}_cds_nucleotide.tsv",
            "log":f"{log_path}/{spid}_cds_per_nucleotide_annotation.DONE"}
    options = {
            'cores':1,
            'memory':'8g',
            'walltime':'12:00:00',
            'account':'EcoGenetics'
        }
    spec = f"""
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate biopython
    echo "START: $(date)"
    echo "JobID: $SLURM_JOBID"
    mkdir -p {work_path}/{spid}_species_data
    mkdir -p {log_path}
    cd {work_path}/{spid}_species_data
    python {script_path}/gff_cds2nucleotide.py {genome} {gff} {spid}_cds_nucleotide.tsv
    echo "FINISH: $(date)"
    jobinfo $SLURM_JOBID
    echo done > {outputs["log"]}
    """
    return AnonymousTarget(inputs=inputs,outputs=outputs,options=options,spec=spec)

def split_annotated_vcfs(work_path: str, ann_vcf: str, spid: str, pop: str, name:str, log_path: str):
    """
    ann_vcf: snpEff annotated vcf file
    """
    inputs = {"ann": ann_vcf}
    outputs = {"log": f"{log_path}/{name}.log",
            "N_vcf":f"{work_path}/{pop}/vcf_by_type/{pop}_N.vcf",
            "S_vcf":f"{work_path}/{pop}/vcf_by_type/{pop}_S.vcf",
            "inter_vcf":f"{work_path}/{pop}/vcf_by_type/{pop}_intergenic.vcf",
            "inter_trim_vcf":f"{work_path}/{pop}/vcf_by_type/{pop}_intergenic_trimmed.vcf",
            "inter_close_vcf":f"{work_path}/{pop}/vcf_by_type/{pop}_intergenic_close_to_gene.vcf"}
    options = {
              'cores':6,
              'memory':'12g',
              'walltime':'12:00:00',
              'account':"EcoGenetics"}
    spec = f"""
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate bcftools
    echo "START: $(date)"
    echo "JobID: $SLURM_JOBID"
    mkdir -p {work_path}/{pop}/vcf_by_type
    mkdir -p {log_path}
    cd {work_path}/{pop}/vcf_by_type
    ### Retrive missense variants ###
    echo {pop} > pop_name.txt
    bcftools view -i 'INFO/ANN[*] ~ "missense_variant" || INFO/ANN[*] ~ "nonsense_variant"' -S pop_name.txt --threads {options['cores']} {inputs["ann"]} > {pop}_N.vcf
    ### Retreive synonymous variants ###
    bcftools view -i 'INFO/ANN[*] ~ "synonymous_variant"' -S pop_name.txt --threads {options['cores']} {inputs["ann"]} > {pop}_S.vcf
    ### Retrive intergenic regions ###
    bcftools view -i 'INFO/ANN[*] ~ "intergenic_region"' -S pop_name.txt --threads {options['cores']} {inputs["ann"]} > {pop}_intergenic.vcf
    ### Retrive intergenic regions, trimmed for downstream and upstream 5000bp defaulat###
    bcftools view -i 'INFO/ANN[*] ~ "intergenic_region" && INFO/ANN[*] !~ "upstream_gene_variant" && INFO/ANN[*] !~ "downstream_gene_variant"' {inputs["ann"]} -S pop_name.txt --threads {options['cores']} > {pop}_intergenic_trimmed.vcf
    ### Retrive intergenic regions, but only for downstream and upstream 5000bp defaulat###
    bcftools view -i 'INFO/ANN[*] ~ "intergenic_region" && INFO/ANN[*] ~ "upstream_gene_variant" && INFO/ANN[*] ~ "downstream_gene_variant"' {inputs["ann"]} -S pop_name.txt --threads {options['cores']} > {pop}_intergenic_close_to_gene.vcf

    ### Check variants within population, remove invariable sites ###
    #########
    echo "FINISH: $(date)"
    jobinfo $SLURM_JOBID
    echo done > {outputs["log"]}
    """
    return AnonymousTarget(inputs=inputs,outputs=outputs,options=options,spec=spec)
def sfs_variants_single_pop(work_path: str, type_code: str, spid: str, pop: str, name:str, log_path: str, count_script_path: str):
    inputs = [f"{work_path}/{pop}/vcf_by_type/{pop}_{type_code}.vcf"]
    outputs = [f"{log_path}/{name}.log",
            f"{work_path}/{pop}/count/{pop}_{type_code}_count_single.tsv",
            f"{work_path}/{pop}/count/{pop}_{type_code}_count_summary.tsv"] 
    options = {
              'cores':1,
              'memory':'2g',
              'walltime':'6:00:00',
              'account':"EcoGenetics"}
    spec = f"""
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate biopython
    echo "START: $(date)"
    echo "JobID: $SLURM_JOBID"
    mkdir -p {work_path}/{pop}/count
    mkdir -p {log_path}
    cd {work_path}/{pop}/count
    ### Count and summary ###
    echo "python {count_script_path}/sfs_variable_count_by_chrom.py {inputs[0]} {pop} {outputs[1]} {outputs[2]}"
    python {count_script_path}/sfs_variable_count_by_chrom.py {inputs[0]} {pop} {outputs[1]} {outputs[2]}
    ########
    echo "FINISH: $(date)"
    jobinfo $SLURM_JOBID
    echo done > {outputs[0]}
    """
    return AnonymousTarget(inputs=inputs,outputs=outputs,options=options,spec=spec)



def pop_bam2bed(work_path: str, bam_file: str, mindp: int, maxdp: int, spid: str, pop: str, name: str, log_path: str):
    inputs = {"bam":bam_file}
    outputs = {
            "log":f"{log_path}/{name}.log",
            "bed":f"{work_path}/{pop}/monosites/beds/{spid}_{pop}_bam_cov_pass.bed"
            }
    options = {
              'cores':4,
              'memory':'32g',
              'walltime':'12:00:00',
              'account':"EcoGenetics"}
    spec = f"""
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate samtools117
    echo "START: $(date)"
    echo "JobID: $SLURM_JOBID"
    mkdir -p {work_path}/{pop}/monosites/beds
    cd {work_path}/{pop}/monosites/beds
    samtools depth -@ 4 -a {bam_file} > {spid}_{pop}_bam_cov.txt
    # Extract min and max coverage thresholds from the TSV file
    min_coverage={mindp} 
    max_coverage={maxdp}
    awk -v min_coverage="$min_coverage" -v max_coverage="$max_coverage" \
    '$3 >= min_coverage && $3 <= max_coverage' {spid}_{pop}_bam_cov.txt > {spid}_{pop}_bam_cov_pass.txt
    awk '{{if (NR==1) {{chrom=$1; start=$2; end=$2}} else if ($2==end+1 && $1==chrom) {{end=$2}} else {{print chrom"\t"start-1"\t"end; chrom=$1; start=$2; end=$2}}}} END {{print chrom"\t"start-1"\t"end}}' {spid}_{pop}_bam_cov_pass.txt > {spid}_{pop}_bam_cov_pass.bed
    echo done > {outputs["log"]}
    echo "START: $(date)"
    jobinfo $SLURM_JOBID
    """
    return AnonymousTarget(inputs=inputs,outputs=outputs,options=options,spec=spec)

def pop_indel2bed(work_path: str, raw_vcf: str, spid: str, pop: str, name: str, log_path: str, script_path: str):
    inputs = {
            "vcf":raw_vcf
            }
    outputs = {
            "log":f"{log_path}/{name}.log",
            "bed":f"{work_path}/{pop}/monosites/beds/{spid}_{pop}_indel.bed"
            }
    options = {
              'cores':1,
              'memory':'2g',
              'walltime':'4:00:00',
              'account':"EcoGenetics"}
    spec = f"""
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate biopython
    echo "START: $(date)"
    echo "JobID: $SLURM_JOBID"
    mkdir -p {work_path}/{pop}/monosites/beds
    cd {work_path}/{pop}/monosites/beds
    python {script_path}/retrive_indel_bed.py {raw_vcf} {spid}_{pop}_indel.bed
    echo done > {outputs["log"]}
    echo "START: $(date)"
    jobinfo $SLURM_JOBID
    """
    return AnonymousTarget(inputs=inputs,outputs=outputs,options=options,spec=spec)

def pop_region_retreive(work_path:str, spid: str, pop: str, name: str, log_path: str, gene_bed: str, inter_bed:str):
    inputs = {
        "gene_bed":f"{gene_bed}",
        "intergenic_bed":f"{inter_bed}",
        "cov_bed":f"{work_path}/{pop}/monosites/beds/{spid}_{pop}_bam_cov_pass.bed",
        "indel_bed":f"{work_path}/{pop}/monosites/beds/{spid}_{pop}_indel.bed"
        }
    outputs = {
        "gene_pass_bed":f"{work_path}/{pop}/monosites/beds/{spid}_{pop}_gene_pass_merged.bed",
        "intergenic_pass_bed":f"{work_path}/{pop}/monosites/beds/{spid}_{pop}_intergenic_pass_merged.bed",
        "log":f"{log_path}/{name}.log"
        }
    options = {
              'cores':1,
              'memory':'2g',
              'walltime':'4:00:00',
              'account':"EcoGenetics"}
    spec = f"""
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate bedtools
    echo "START: $(date)"
    echo "JobID: $SLURM_JOBID"
    mkdir -p {work_path}/{pop}/monosites/beds
    cd {work_path}/{pop}/monosites/beds
    sort -k1,1 -k2,2n {inputs["gene_bed"]} > {spid}_{pop}_gene_sort.bed
    sort -k1,1 -k2,2n {inputs["intergenic_bed"]} > {spid}_{pop}_intergenic_sort.bed
    sort -k1,1 -k2,2n {inputs["cov_bed"]} > {spid}_{pop}_cov_sort.bed
    bedtools intersect -a {spid}_{pop}_gene_sort.bed -b {spid}_{pop}_cov_sort.bed > {spid}_{pop}_gene_cov.bed
    bedtools intersect -a {spid}_{pop}_intergenic_sort.bed -b {spid}_{pop}_cov_sort.bed > {spid}_{pop}_intergenic_cov.bed
    sort -k1,1 -k2,2n {spid}_{pop}_gene_cov.bed > {spid}_{pop}_gene_cov_sort.bed
    sort -k1,1 -k2,2n {spid}_{pop}_intergenic_cov.bed > {spid}_{pop}_intergenic_cov_sort.bed
    sort -k1,1 -k2,2n {inputs["indel_bed"]} > {spid}_{pop}_indel_sort.bed
    bedtools subtract -a {spid}_{pop}_gene_cov_sort.bed -b {spid}_{pop}_indel_sort.bed > {spid}_{pop}_gene_pass.bed
    bedtools subtract -a {spid}_{pop}_intergenic_cov_sort.bed -b {spid}_{pop}_indel_sort.bed > {spid}_{pop}_intergenic_pass.bed
    sort -k1,1 -k2,2n {spid}_{pop}_gene_pass.bed > {spid}_{pop}_gene_pass_sort.bed
    sort -k1,1 -k2,2n {spid}_{pop}_intergenic_pass.bed > {spid}_{pop}_intergenic_pass_sort.bed
    bedtools merge -i {spid}_{pop}_gene_pass_sort.bed > {spid}_{pop}_gene_pass_merged.bed
    bedtools merge -i {spid}_{pop}_intergenic_pass_sort.bed > {spid}_{pop}_intergenic_pass_merged.bed
    #sort -u {spid}_{pop}_gene_pass_merged.bed > {spid}_{pop}_gene_pass_final.bed
    #sort -u {spid}_{pop}_intergenic_pass_merged.bed >  {spid}_{pop}_intergenic_pass_final.bed
    rm {spid}_{pop}_gene_sort.bed
    rm {spid}_{pop}_intergenic_sort.bed
    rm {spid}_{pop}_gene_cov.bed
    rm {spid}_{pop}_intergenic_cov.bed
    rm {spid}_{pop}_gene_cov_sort.bed
    rm {spid}_{pop}_intergenic_cov_sort.bed
    rm {spid}_{pop}_indel_sort.bed
    echo done > {outputs["log"]}
    echo "START: $(date)"
    jobinfo $SLURM_JOBID
    """
    return AnonymousTarget(inputs=inputs,outputs=outputs,options=options,spec=spec)

def monosites_retrival(work_path:str, spid: str, pop: str, name: str,script_path: str, log_path: str):
    inputs = {
        "gene_bed":f"{work_path}/{pop}/monosites/beds/{spid}_{pop}_gene_pass_merged.bed",
        "intergenic_bed":f"{work_path}/{pop}/monosites/beds/{spid}_{pop}_intergenic_pass_merged.bed",
        "N_vcf":f"{work_path}/{pop}/vcf_by_type/{pop}_N.vcf",
        "S_vcf":f"{work_path}/{pop}/vcf_by_type/{pop}_S.vcf",
        "inter_vcf":f"{work_path}/{pop}/vcf_by_type/{pop}_intergenic.vcf"
    }
    outputs = {
        "gene_bed_invariant":f"{work_path}/{pop}/monosites/summary/{spid}_{pop}_gene_monosites.bed",
        "intergenic_bed_invariant":f"{work_path}/{pop}/monosites/summary/{spid}_{pop}_intergenic_monosites.bed",
        "gene_summary_tsv":f"{work_path}/{pop}/monosites/summary/{spid}_{pop}_gene_callable_summary.tsv",
        "intergenic_summary_tsv":f"{work_path}/{pop}/monosites/summary/{spid}_{pop}_intergenic_callable_summary.tsv",
        "log":f"{log_path}/{name}.log"
    }
    options = {
              'cores':1,
              'memory':'16g',
              'walltime':'12:00:00',
              'account':"EcoGenetics"}
    spec = f"""
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate pysam
    echo "START: $(date)"
    echo "JobID: $SLURM_JOBID"
    mkdir -p {work_path}/{pop}/monosites/summary
    cd {work_path}/{pop}/monosites/summary
    conda activate bcftools
    bgzip -c {inputs["N_vcf"]} > N_vcf.gz
    bgzip -c {inputs["S_vcf"]} > S_vcf.gz
    tabix -p vcf N_vcf.gz
    tabix -p vcf S_vcf.gz
    bcftools concat -a N_vcf.gz S_vcf.gz -Oz -o gene.vcf.gz
    bcftools sort gene.vcf.gz -Oz -o sorted_gene.vcf.gz
    bgzip -d sorted_gene.vcf.gz
    conda activate pysam
    python {script_path}/monosites_retrieve.py {inputs["gene_bed"]} sorted_gene.vcf {spid}_{pop}_gene_callable_summary.tsv {spid}_{pop}_gene_monosites.bed
    python {script_path}/monosites_retrieve.py {inputs["intergenic_bed"]} {inputs["inter_vcf"]} {spid}_{pop}_intergenic_callable_summary.tsv {spid}_{pop}_intergenic_monosites.bed
    echo done > {outputs["log"]}
    echo "START: $(date)"
    jobinfo $SLURM_JOBID
    """
    return AnonymousTarget(inputs=inputs,outputs=outputs,options=options,spec=spec)

def monosites_gene_sites_summary(work_path:str, spid: str, pop: str, name: str,script_path: str, log_path: str):
    inputs = {
            "site_ann":f"{work_path}/{spid}_species_data/{spid}_cds_nucleotide.tsv",
            "gene_bed":f"{work_path}/{pop}/monosites/summary/{spid}_{pop}_gene_monosites.bed"
            }
    outputs = {
            "log":f"{log_path}/{name}.log",
            "NS_count":f"{work_path}/{pop}/monosites/summary/{spid}_{pop}_cds_monosites_NS_count.tsv"
            }
    options = {
              'cores':1,
              'memory':'16g',
              'walltime':'4:00:00',
              'account':"EcoGenetics"}
    spec = f"""
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate biopython
    echo "START: $(date)"
    echo "JobID: $SLURM_JOBID"
    mkdir -p {work_path}/{pop}/monosites/summary
    cd {work_path}/{pop}/monosites/summary
    sed '1d' {inputs["site_ann"]} > headless_site_ann.tsv
    head -n 1 {inputs["site_ann"]} > header_site_ann.tsv
    sort -k1,1 -k2,2n headless_site_ann.tsv > headless_{spid}_site_ann_sorted.tsv
    cat header_site_ann.tsv headless_{spid}_site_ann_sorted.tsv > {spid}_site_ann_sorted.tsv
    rm headless_site_ann.tsv
    rm header_site_ann.tsv
    rm headless_{spid}_site_ann_sorted.tsv
    python {script_path}/site_in_bed.py {inputs["gene_bed"]} {spid}_site_ann_sorted.tsv {spid}_{pop}_callable_region_cds_monosites.tsv
    python {script_path}/count_monosites_NS.py {spid}_{pop}_callable_region_cds_monosites.tsv {spid}_{pop}_cds_monosites_NS_count.tsv
    echo done > {outputs["log"]}
    echo "START: $(date)"
    jobinfo $SLURM_JOBID
    """
    return AnonymousTarget(inputs=inputs,outputs=outputs,options=options,spec=spec)

def generate_final_count(work_path:str, spid: str, pop: str, name: str,script_path: str, log_path: str):
    inputs = {
            "NS_count":f"{work_path}/{pop}/monosites/summary/{spid}_{pop}_cds_monosites_NS_count.tsv",
            "inter_genic_count":f"{work_path}/{pop}/monosites/summary/{spid}_{pop}_intergenic_callable_summary.tsv",
            "N_variant":f"{work_path}/{pop}/count/{pop}_N_count_summary.tsv",
            "S_variant":f"{work_path}/{pop}/count/{pop}_S_count_summary.tsv",
            "inter_variant":f"{work_path}/{pop}/count/{pop}_intergenic_count_summary.tsv"
            }
    outputs = {
            "N_results":f"{work_path}/{pop}/results/N_site_sfs.tsv",
            "S_results":f"{work_path}/{pop}/results/S_site_sfs.tsv",
            "inter_results":f"{work_path}/{pop}/results/inter_site_sfs.tsv",
            "log":f"{log_path}/{name}.log"
            }
    options = {
              'cores':1,
              'memory':'2g',
              'walltime':'1:00:00',
              'account':"EcoGenetics"}
    spec = f"""
    CONDA_BASE=$(conda info --base)
    source $CONDA_BASE/etc/profile.d/conda.sh
    conda activate biopython
    echo "START: $(date)"
    echo "JobID: $SLURM_JOBID"
    mkdir -p {work_path}/{pop}/results
    cd {work_path}/{pop}/results
    python {script_path}/generate_final_count.py {inputs["NS_count"]} {inputs["inter_genic_count"]} {inputs["N_variant"]} {inputs["S_variant"]} {inputs["inter_variant"]} {outputs["N_results"]} {outputs["S_results"]} {outputs["inter_results"]}
    echo done > {outputs["log"]}
    echo "START: $(date)"
    jobinfo $SLURM_JOBID
    """
    return AnonymousTarget(inputs=inputs,outputs=outputs,options=options,spec=spec)

