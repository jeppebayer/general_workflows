from gwf import Workflow
from gwf.workflow import collect
import os, yaml, glob, sys
from workflow_templates import *

def purge_dup_for_assembly_workflow(config_file: str):
    """
    """
    # --------------------------------------------------
    #                  Configuration
    # --------------------------------------------------
    CONFIG = yaml.safe_load(open(config_file))
    ACCOUNT: str = CONFIG['account']
    ANN_FILE: str = CONFIG['ann_vcf']
    RAW_VCF: str = CONFIG['raw_vcf']
    GFF: str = CONFIG['annotation_gff']
    GENOME: str = CONFIG['genome_fasta']
    POP_BAM: str = CONFIG['population_bam_alignment']
    POP_DEPTH: str = CONFIG['population_varaint_depth_summary']
    TASK_ID: str = CONFIG['task_id']
    SPECIES_ID: str = CONFIG['species_id']
    POP_ID: str = CONFIG['population']
    WORK_DIR: str = CONFIG['working_directory_path']
    OUTPUT_DIR: str = CONFIG['output_directory_path']
    LOG_DIR: str = CONFIG['log_directory_path']
    DATA_DIR: str = CONFIG['data_directory_path']
    HELP_SCRIPTS_PATH: str = CONFIG['help_scripts_path']
    GENE_BED: str = CONFIG['gene_bed']
    INTER_BED: str = CONFIG['intergenic_bed']
    # --------------------------------------------------
    #                  Workflow
    # --------------------------------------------------

    gwf = Workflow(
        defaults={'account': ACCOUNT}
    )
    # summarise N and S information per nucleotide in the CDS regions of species reference genome
    gwf.target_from_template(
        name = f'annotate_CDS_nucleotide_{SPECIES_ID}',
        template = annotate_N_S_sites(
            work_path = WORK_DIR,
            script_path = HELP_SCRIPTS_PATH,
            spid = SPECIES_ID,
            genome = GENOME,
            gff = GFF,
            log_path = LOG_DIR
            )
        )
    # Split a annoated vcf file from snpEff by types
    gwf.target_from_template(
        name = f'sfs_{TASK_ID}_{SPECIES_ID}_{POP_ID}',
        template = split_annotated_vcfs(
            work_path = WORK_DIR,
            ann_vcf = ANN_FILE,
            spid = SPECIES_ID,
            pop = POP_ID,
            name = f'sfs_{TASK_ID}_{SPECIES_ID}_{POP_ID}',
            log_path = LOG_DIR
        )
    )

    # Count number of variable sites per type 
    for type_code in ["N","S","intergenic","intergenic_trimmed","intergenic_close_to_gene"]:
        gwf.target_from_template(
            name = f'sfs_{SPECIES_ID}_{POP_ID}_count_{type_code}',
            template = sfs_variants_single_pop(
                work_path = WORK_DIR,
                type_code = type_code,
                spid = SPECIES_ID,
                pop = POP_ID,
                name = f'sfs_{SPECIES_ID}_{POP_ID}_count_{type_code}',
                log_path = LOG_DIR,
                count_script_path = HELP_SCRIPTS_PATH
            )
        )
    # Generate bed file based on Coverage of bam files for callable regions
    gwf.target_from_template(
        name = f'sfs_{SPECIES_ID}_{POP_ID}_bam_cov_pass_bed',
        template = pop_bam2bed(
            work_path = WORK_DIR,
            bam_file = POP_BAM,
            depth_threshold = POP_DEPTH,
            spid = SPECIES_ID,
            pop = POP_ID,
            name = f'sfs_{SPECIES_ID}_{POP_ID}_bam_cov_pass_bed',
            log_path = LOG_DIR
            )
        )
    # Generate bed file based on indel variants from the raw vcf file
    gwf.target_from_template(
        name = f'sfs_{SPECIES_ID}_{POP_ID}_indel_bed',
        template = pop_indel2bed(
            work_path = WORK_DIR,
            raw_vcf = RAW_VCF,
            spid = SPECIES_ID,
            pop = POP_ID,
            name = f'sfs_{SPECIES_ID}_{POP_ID}_indel_bed',
            log_path = LOG_DIR,
            script_path = HELP_SCRIPTS_PATH
            )
        )
    # combine bed files to account all callable regions
    gwf.target_from_template(
        name = f'sfs_{SPECIES_ID}_{POP_ID}_callable_bed',
        template = pop_region_retreive(
            work_path = WORK_DIR,
            spid = SPECIES_ID,
            pop = POP_ID,
            name =  f'sfs_{SPECIES_ID}_{POP_ID}_callable_bed',
            log_path = LOG_DIR,
            gene_bed = GENE_BED,
            inter_bed = INTER_BED
            )
        )
    # retrival callable bed summary and invariant beds
    gwf.target_from_template(
        name = f'sfs_{SPECIES_ID}_{POP_ID}_callable_summary',
        template = monosites_retrival(
            work_path = WORK_DIR,
            spid = SPECIES_ID,
            pop = POP_ID,
            name = f'sfs_{SPECIES_ID}_{POP_ID}_callable_summary',
            script_path = HELP_SCRIPTS_PATH,
            log_path = LOG_DIR
        )
    )
    # count NS site number of callable gene regions
    gwf.target_from_template(
        name = f'sfs_{SPECIES_ID}_{POP_ID}_cds_monosites_NS_count',
        template = monosites_gene_sites_summary(
            work_path = WORK_DIR,
            spid = SPECIES_ID,
            pop = POP_ID,
            name = f'sfs_{SPECIES_ID}_{POP_ID}_cds_monosites_NS_count',
            script_path = HELP_SCRIPTS_PATH,
            log_path = LOG_DIR
            )
        )
    # generate final count tsv for SFS (including invariant sites in the callable regions)
    gwf.target_from_template(
        name = f'sfs_{SPECIES_ID}_{POP_ID}_final_tsv_generate',
        template = generate_final_count(
            work_path = WORK_DIR,
            spid = SPECIES_ID,
            pop = POP_ID,
            name = f'sfs_{SPECIES_ID}_{POP_ID}_final_tsv_generate',
            script_path = HELP_SCRIPTS_PATH,
            log_path = LOG_DIR)
        )
    return gwf
