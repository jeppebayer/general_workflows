#!/bin/env python3
from gwf import Workflow
from gwf.workflow import collect
import os, yaml, glob, sys, re
from workflow_templates import *

def fst_and_pi_wf(config_file: str = glob.glob('*config.y*ml')[0]):
    """
    Workflow: Estimates pi and Fst from filtered and neutral :format:`VCF` using equations desribed in :script: popoolation :script: fst-sliding.
    
    :param str config_file:
        Configuration file containing pre-defined set of variables
    """
    # --------------------------------------------------
    #                  Configuration
    # --------------------------------------------------
    
    CONFIG = yaml.safe_load(open(config_file)) #opens yml file for reading
    ACCOUNT: str = CONFIG['account'] # gets the setting after account as string
    SPECIES_NAME: str = CONFIG['species_name']
    OUTPUT_DIR: str = CONFIG['output_directory_path']
    WORK_DIR: str = CONFIG['working_directory_path']
    #VCF_BASE_FOLDER: str = CONFIG['vcf_base_folder']
    TAXONOMY: str = CONFIG['taxonomic_group']
    BED_PATH: str = CONFIG['bed_files_path']
    GENOME_PATH: str = CONFIG['reference_genome_path']
    #ANNOTATION_GTF: str = CONFIG['annotation_gtf']
    VCF_FILE: str = CONFIG['vcf_file']
    #COLLECTION_SITES_FILE: str = CONFIG['collection_sites_file']
    COUNT_POSITIONS_FILE: str = CONFIG['count_of_positions_file']


    # --------------------------------------------------
    #                  Workflow
    # --------------------------------------------------
    
    gwf = Workflow(
        defaults={'account': ACCOUNT}
    ) # what does this do?
    

    ### MAKE DIRECTORIES
    #print("making directories")
    # make directories
    new_wd_fst_pi=f'{WORK_DIR}/fst_pi/{TAXONOMY}/{SPECIES_NAME.replace(" ","_")}/intermediary_files'
    new_out_pi=f'{OUTPUT_DIR}/pi/{TAXONOMY}/{SPECIES_NAME.replace(" ","_")}/'
    new_out_fst=f'{OUTPUT_DIR}/fst/{TAXONOMY}/{SPECIES_NAME.replace(" ","_")}/'
    #if not os.path.isdir(new_wd):
    #    os.makedirs(new_wd)

    # defines new files/directories base on config
    bed_genes = glob.glob(f'{BED_PATH}/*_genomic.genes.bed')[0]
    ##bed_intergenes = glob.glob(f'{BED_PATH}/*_genomic.intergenic.bed')[0]
    bed_repeats = glob.glob(f'{BED_PATH}/*_genomic.repeats.bed')[0]
    bed_neutral = f'{os.path.dirname(bed_genes)}/{os.path.basename(bed_genes).replace("genes","neutral")}'
    print(bed_neutral)
    ## f'{bed_repeats.replace("repeats","neutral")}'
    


    # Define target making genome.fna.fai, if not already present
    #make_fna_fai_target = gwf.target_from_template(
    #        name='make_genome_fai',
     #       template=make_genome_fai(
      #          ref_genome_file=GENOME_PATH,
       #         fasta_fai_output=f'{new_wd}/{GENOME_PATH.split('/')[-1].replace(".fna", ".fna.fai")}')
        #)

    ########################
    ### MAKE NEUTRAL BED ###
    ########################

    # Define target with input genes and repeats bed
        # output neutral bed
    make_neutral_bed_target = gwf.target_from_template(
        name='bed_make_neutral_bed',
        template=make_neutral_bed_imporved(
            genes_bed_file = bed_genes,
            repeats_bed_file = bed_repeats,
            reference_genome_file = GENOME_PATH,
            neutral_bed_out = bed_neutral, 
            working_dir = new_wd_fst_pi)
    )

    # make varying intervals of intergenic region beds
    #interbed_wd = f'{WORK_DIR}/fst_pi/{TAXONOMY}/{SPECIES_NAME.replace(" ","_")}/intergene_neutrality'
    #percentage_list = list(range(0,10+1))
    #bp_list = list(range(0,5000+1,500))
    #output_list = [ os.path.join(interbed_wd,'intergenic_' + str(n) + '_percent_rem.bed') for n in percentage_list]

    #make_neutral_bed_inter = gwf.target_from_template(
     #   name='make_various_intergene_beds',
      #  template=make_various_intergene_beds(
       #     intergenes_bed_file = bed_intergenes,
        #    repeats_bed_file = bed_repeats,
         #   working_directory = interbed_wd,
          #  genome_path = make_fna_fai_target.outputs['genome_fai'])
    #)

    # The several bed outputs should go in here, and we map over them to create netural beds
    # make a list of dictionaries with input intergenic beds and output neutral bed names
    #bed_neutral = f'{new_wd}/{os.path.basename(bed_intergenes).replace("genes","neutral")}'
    #make_neutral_bed_target_fromIntergene = gwf.target_from_template(
     #   name='make_neutral_bed_from_intergene',
      #  template=make_neutral_bed_inter(
       #     intergenes_bed_file = bed_intergenes,
        #    repeats_bed_file = bed_repeats,
         #   genome_bedstyle = make_fna_fai_target.outputs['genome_fai'],
          #  neutral_bed_out = bed_neutral)
    #)






# Target to deminish VCF to only contain neutral sites

    ###############################
    ###     Make neutral vcf    ###
    ###############################
    print("Making neutral version of VCF file")
    #input_dict=[{'vcf_file': f, 'working_directory': p} for f, p in zip(files_list, files_popdir)] # making combined dictionary of files and new wd paths using list comprehension and zip

    neutral_vcf_files_template = gwf.target_from_template(
        #name=make_neutral_vcf,
        name='vcf_make_neutral',
        template=make_neutral_vcf_improved(
            vcf_file = VCF_FILE,
            working_directory = f'{new_wd_fst_pi}/neutral_vcf_temp/',
            neutral_bed = make_neutral_bed_target.outputs['neutral_bed']
        )
    )


    ###############################
    ###     CAlC Allele Freq    ###
    ###############################
        
    # get allele frequencies
    print("Re-calculating allele frequencies for each population")

    # add template with single run: ( can be paralellized if need be - requires modifications )
    #recalculate_AF
    recalculate_AF_template = gwf.target_from_template(
        name = 'recalculate_AF',
        template = recalculate_AF_improved(
            input_vcf = neutral_vcf_files_template.outputs['neutral_vcf'],
            working_dir = new_wd_fst_pi, 
            species_short = species_abbreviation(SPECIES_NAME))
    )



    # get the output with allele freqs:
    


    ###########################
    ###     CAlC Fst        ###
    ###########################
    print(recalculate_AF_template.outputs['allele_freq'])
    fst_from_AF_template = gwf.target_from_template(
        name = 'fst_calc_from_AF',
        template = fst_calc_from_AF_improved(
            allele_freq_file = recalculate_AF_template.outputs['allele_freq'],
            working_dir = new_wd_fst_pi, 
            species_short = species_abbreviation(SPECIES_NAME))
    )
    
    
    ##########################
    ###     CAlC pi        ###
    ##########################


    # need counting file fir this
    pi_calculation_all_pos = gwf.target_from_template(
        name = 'pi_calculation_all_positions',
        template = calculate_pi_template_improved(
            allele_freq_file = recalculate_AF_template.outputs['allele_freq'],
            working_directory = new_wd_fst_pi,
            count_file = COUNT_POSITIONS_FILE, 
            species_short = species_abbreviation(SPECIES_NAME))
    )


    return gwf
