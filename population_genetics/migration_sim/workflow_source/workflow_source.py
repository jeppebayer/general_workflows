#!/bin/env python3
from gwf import Workflow
from gwf.workflow import collect
import os, yaml, glob, sys, re
import subprocess as sp
from workflow_templates import *

def migration_simulation(config_file: str = glob.glob('*config.y*ml')[0]):
    """
    Workflow: Makes 2dSFS for each population pair, and simulates present and past migration rates at specified intervals, using FastSimCoal.
    
    :param str config_file:
        Configuration file containing pre-defined set of variables
    """
    # --------------------------------------------------
    #                  Configuration
    # --------------------------------------------------
    
    # EDIT CONFIGS
    CONFIG = yaml.safe_load(open(config_file)) #opens yml file for reading
    ACCOUNT: str = CONFIG['account'] # gets the setting after account as string
    SPECIES_NAME: str = CONFIG['species_name']
    OUTPUT_DIR: str = CONFIG['output_directory_path']
    WORK_DIR: str = CONFIG['working_directory_path']
    #VCF_BASE_FOLDER: str = CONFIG['vcf_base_folder']
    TAXONOMY: str = CONFIG['taxonomic_group']
    #BED_PATH: str = CONFIG['bed_files_path']
    #GENOME_PATH: str = CONFIG['reference_genome_path']
    #ANNOTATION_GTF: str = CONFIG['annotation_gtf']
    #VCF_FILES: list = CONFIG['vcf_lists']
    #COLLECTION_SITES_FILE: str = CONFIG['collection_sites_file']
    #neutral_position_count_file: str = CONFIG['count_of_positions_file']
    MAX_MIG_DIVIDE: int = CONFIG['max_generations_migration_divide']
    MIN_MIG_DIVIDE: int = CONFIG['min_generations_migration_divide']
    VCF_FILE: str = CONFIG['vcf_file']
    


    # --------------------------------------------------
    #                  Workflow
    # --------------------------------------------------
    
    gwf = Workflow(
        defaults={'account': ACCOUNT}
    ) # what does this do?
    

    ### MAKE DIRECTORIES
    #print("making directories")
    # make directories
    new_wd=f'{WORK_DIR}/migration_sim/{TAXONOMY}/{SPECIES_NAME.replace(" ","_")}/fsc/'
    new_out_mig=f'{OUTPUT_DIR}/migration_sim/{TAXONOMY}/{SPECIES_NAME.replace(" ","_")}/'
    if not os.path.isdir(new_wd):
        os.makedirs(new_wd)
    if not os.path.isdir(new_out_mig):
        os.makedirs(new_out_mig)

    # defines new files/directories base on config
    #bed_genes = glob.glob(f'{BED_PATH}/*_genomic.genes.bed')[0]
    #bed_intergenes = glob.glob(f'{BED_PATH}/*_genomic.intergenic.bed')[0]
    #bed_repeats = glob.glob(f'{BED_PATH}/*_genomic.repeats.bed')[0]
    #bed_neutral = f'{new_wd}/{os.path.basename(bed_genes).replace("genes","neutral")}'
    # f'{bed_repeats.replace("repeats","neutral")}'
    

    #####################################
    ### MAKE data ready for 2dSFS BED ###
    #####################################

    
	# Parse individual IDs and their counters
	
    #awk '/^#CHROM/ {{for(i=10; i<=NF; i++) print $i, "\t", i-9; exit}}' {inputs['vcf_file']} > {outputs['parsed_pops_file']}
    #args = ["awk", r'/^#CHROM/ {for(i=10; i<=NF; i++) print $i, "\t", i-9; exit}', VCF_FILE]
    #cmd = """awk '/^#CHROM/ {for(i=10; i<=NF; i++) print $i, "\t", i-9; exit}' VCF_FILE"""
    #parsed_pops = sp.Popen(args, stdin = sp.PIPE, stdout = sp.PIPE, stderr = sp.PIPE, shell = True)
    #parsed_pops.wait()
    awk_parse_pops_from_vcf(VCF_FILE, f'{new_wd}/2dSFS/{species_abbreviation(SPECIES_NAME)}_parsed_pops.tsv')
    
    #args = ["awk", r'{OFS="\t"; print $2,$4,$5,$6}', "B3LYPD.txt"]
    #p = sp.Popen(args, stdin = sp.PIPE, stdout = sp.PIPE, stderr = sp.PIPE )
    #print(p.stdout.readline()) # will give you the first line of the awk output



	# Generate all possible pairs where counter2 >= counter1
	#awk -v OFS='\t' 'NR==FNR {{id[NR] = $1; counter[NR] = $2; next}}
	 #	 {{for (i=1; i<FNR; i++) if (counter[i] <= $2) print id[i], $1, counter[i], $2}}' {outputs['parsed_pops_file']} {outputs['parsed_pops_file']} > {outputs['pop_pairs_file']}
		 # $ids_output_file $ids_output_file > $pairs_output_file


    # Define target making genome.fna.fai
    make_ready_for_2dSFS_pairs = gwf.target_from_template(
            name='prepare_data_sfs',
            template=prepare_data(
                vcf_file=VCF_FILE,
                output_prefix=species_abbreviation(SPECIES_NAME),
                working_dir=f'{new_wd}/2dSFS/')
        )

    #################
    ### RUN 2dSFS ###
    #################
    # create input dict based on output from above
    # make function that returns dicst of pop pairs
    # read pop{outputs['pop_pairs_file']}
# neq idea:# make
    # make dictionary directly from parsing the names of directories that will be made in previous template

#    input_dict_list = create_input_dict_2dSFS(pair_file = make_ready_for_2dSFS_pairs.outputs['pop_pairs_file'])
    #names = [d for d in input_dict_list['name']]
    #print(input_dict_list)
    
    #neutral_vcf_files_runtemplate_map = gwf.map(
     #   #name=create_input_dict_2dSFS_target(pair_file = make_ready_for_2dSFS_pairs.outputs['pop_pairs_file']),
      #  template_func = pair_2DSFS_map_target,
       # inputs = create_input_dict_2dSFS_target(pair_file = make_ready_for_2dSFS_pairs.outputs['pop_pairs_file']),
        #extra = {'out2_file': make_ready_for_2dSFS_pairs.outputs['outfile_dat'], 
         #        'working_directory': f'{new_wd}/2dSFS/'})
    

    return gwf
