#!/bin/env python3
from gwf import Workflow
from gwf.workflow import collect
import os, yaml, glob, sys, re
from workflow_templates import *


#def fst_and_pi_wf(config_file: str = glob.glob('*config.y*ml')[0]):
def fst_and_pi_wf(config_file: str, gwf):    
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
    AREA_TYPE: str = CONFIG['area_collected']
    BED_PATH: str = CONFIG['bed_files_path']
    GENOME_PATH: str = CONFIG['reference_genome_path']
    #ANNOTATION_GTF: str = CONFIG['annotation_gtf']
    VCF_FILE: str = CONFIG['vcf_file']
    VCF_AN: int = CONFIG['vcf_AN']
    INDELS_INFO_FILE: str = CONFIG['bcf_tools_stats_indels_info']
    BAM_PATH: str = CONFIG['bam_general_path']
    #COLLECTION_SITES_FILE: str = CONFIG['collection_sites_file']
    COUNT_POSITIONS_FILE: str = CONFIG['count_of_positions_file']
    PI_COLLECTION_FILE: str = CONFIG['output_collection_file_pi']
    RUN_PI: str = CONFIG['estimate_pi_yes_no']
    RUN_FST: str = CONFIG['estimate_fst_yes_no']
    FILT_MIN_COUNT: int = CONFIG['filter-sample-min-count']
    FILT_MIN_DEPTH: int = CONFIG['filter-sample-min-read-depth']
    FILT_MAX_DEPTH: int = CONFIG['filter-sample-max-read-depth']
    FILT_WINDOW: int = CONFIG['window-interval-width']
    FILT_POOL_SIZE: int = CONFIG['pool-sizes']
    EXCLUDE_POPS: list = CONFIG['exclude_list']
    # --------------------------------------------------
    #                  Workflow
    # --------------------------------------------------
    
    ##gwf = Workflow(
    ##    defaults={'account': ACCOUNT}
    ##) # what does this do?
    

    ### MAKE DIRECTORIES
    #print("making directories")
    # make directories
    #new_wd_fst_pi=f'{WORK_DIR}/fst_pi/{TAXONOMY}/{AREA_TYPE.replace(" ","_")}/{SPECIES_NAME.replace(" ","_")}/intermediary_files'
    new_wd_fst_pi=f'{WORK_DIR}/fst_pi/{TAXONOMY}/{SPECIES_NAME.replace(" ","_")}/{AREA_TYPE.replace(" ","_")}/intermediary_files'
    if not os.path.isdir(new_wd_fst_pi):
        os.makedirs(new_wd_fst_pi)
    if RUN_PI:
        #new_out_pi=f'{OUTPUT_DIR}/pi/{TAXONOMY}/{AREA_TYPE.replace(" ","_")}/{SPECIES_NAME.replace(" ","_")}/'
        new_out_pi=f'{OUTPUT_DIR}/pi/{TAXONOMY}/{SPECIES_NAME.replace(" ","_")}/{AREA_TYPE.replace(" ","_")}/'
        if not os.path.isdir(new_out_pi):
            os.makedirs(new_out_pi)
    else:
        print("Pi is not being estimated. Activate pi calculation in config file, if you want pi estimated.")
    if RUN_FST:
        #new_out_fst=f'{OUTPUT_DIR}/fst/{TAXONOMY}/{AREA_TYPE.replace(" ","_")}/{SPECIES_NAME.replace(" ","_")}/'
        new_out_fst=f'{OUTPUT_DIR}/fst/{TAXONOMY}/{SPECIES_NAME.replace(" ","_")}/{AREA_TYPE.replace(" ","_")}/'
        if not os.path.isdir(new_out_fst):
            os.makedirs(new_out_fst)
    else:
        print("Fst is not being estimated. Activate Fst calculation in config file, if you want Fst estimated.")
    
    #if not os.path.isdir(new_wd):
    #    os.makedirs(new_wd)

    # defines new files/directories base on config
    #print(f'{BED_PATH}')
    #print(f'{BED_PATH}/EG_OrcVil_23072024_genomic.genes.bed')
    #print(glob.glob(f'{BED_PATH}/*_genomic.genes.bed'))
    if RUN_PI or RUN_FST:
        print(f'Running: {TAXONOMY} - {SPECIES_NAME.replace(" ","_")} - {AREA_TYPE.replace(" ","_")}')
        bed_genes = glob.glob(f'{BED_PATH}/*_genomic.genes.bed')[0]
        ##bed_intergenes = glob.glob(f'{BED_PATH}/*_genomic.intergenic.bed')[0]
        bed_repeats = glob.glob(f'{BED_PATH}/*_genomic.repeats.bed')[0]
        bed_neutral = f'{os.path.dirname(bed_genes)}/{os.path.basename(bed_genes).replace("_genomic.genes.bed", "_" + AREA_TYPE.replace(" ","_") + "_genomic.genes.bed").replace("genes","neutral")}'
        #print(bed_neutral)
        ## f'{bed_repeats.replace("repeats","neutral")}'
    


    ########################
    ### MAKE NEUTRAL BED ###
    ########################

    # Make neutral bed. For both pi and fst.     

    if RUN_PI or RUN_FST:
        make_neutral_bed_target = gwf.target_from_template(
            name=f'bed_make_neutral_bed_{species_abbreviation(SPECIES_NAME)}_{AREA_TYPE.replace(" ","_")}',
            template=make_neutral_bed_improved(
                genes_bed_file = bed_genes,
                repeats_bed_file = bed_repeats,
                reference_genome_file = GENOME_PATH,
                neutral_bed_out = bed_neutral, 
                working_dir = new_wd_fst_pi)
        )



# Target to deminish VCF to only contain neutral sites

    ###############################
    ###     Make neutral vcf    ###
    ###############################
    if RUN_PI or RUN_FST:
        #print("Making neutral version of VCF file")
        #input_dict=[{'vcf_file': f, 'working_directory': p} for f, p in zip(files_list, files_popdir)] # making combined dictionary of files and new wd paths using list comprehension and zip
        #print(VCF_FILE)
        neutral_vcf_files_template = gwf.target_from_template(
            name=f'vcf_make_neutral_{species_abbreviation(SPECIES_NAME)}_{AREA_TYPE.replace(" ","_")}',
            template=make_neutral_vcf_improved(
                vcf_file = VCF_FILE,
                working_directory = f'{new_wd_fst_pi}/neutral_vcf_temp/',
                neutral_bed = make_neutral_bed_target.outputs['neutral_bed']
            )
        )
        
    
    ###############################
    ###     Make neutral BAM    ###
    ###############################
    if RUN_PI:
        #print("Making neutral version of BAM file")
        #input_dict=[{'vcf_file': f, 'working_directory': p} for f, p in zip(files_list, files_popdir)] # making combined dictionary of files and new wd paths using list comprehension and zip
        
        #import json
        #from pathlib import Path
        #if not os.path.isdir(f'{WORK_DIR}/fst_pi/{TAXONOMY}/bam_list_check'):
        #    os.makedirs(f'{WORK_DIR}/fst_pi/{TAXONOMY}/bam_list_check')
        #STATE_FILE = Path(f'{WORK_DIR}/fst_pi/{TAXONOMY}/bam_list_check/bam_dict_state.json')

        ## Load previous state
        #if STATE_FILE.exists():
        #    with STATE_FILE.open("r") as f:
        #        all_bam_dicts = json.load(f)
        #else:
        #    all_bam_dicts = []

        # Your new bam_dict_list for this run
        bam_dict_list = find_bam_files(BAM_PATH)  # however it's generated

        ## Check for duplicates
        #for new_dict in bam_dict_list:
        #    if new_dict in all_bam_dicts:
        #        print("Duplicate found:", new_dict)
        #        print("Not running duplicate.")
        #    else:
        #        all_bam_dicts.append(new_dict)
                

        ## Save updated state
        #with STATE_FILE.open("w") as f:
        #    json.dump(all_bam_dicts, f, indent=2)

                    
        #bam_dict_list = find_bam_files(BAM_PATH)
        #print((bam_dict_list))
        neutral_bam_files_map = gwf.map(
            name=create_bam_run_name,
            #name=f'bam_make_neutral_{species_abbreviation(SPECIES_NAME)}_{AREA_TYPE.replace(" ","_")}',
            template_func = make_neutral_bam_improved,
            inputs = bam_dict_list,
            extra = {'neutral_bed': make_neutral_bed_target.outputs['neutral_bed'],
                    'working_directory': f'{new_wd_fst_pi}/neutral_bam_temp/'})



    ###############################
    ###     CAlC Allele Freq    ###
    ###############################
    if RUN_PI or RUN_FST:
        # get allele frequencies
        #print("Re-calculating allele frequencies for each population")

        # add template with single run: ( can be paralellized if need be - requires modifications )
        #recalculate_AF
        recalculate_AF_template = gwf.target_from_template(
            name = f'recalculate_AF_{species_abbreviation(SPECIES_NAME)}_{AREA_TYPE.replace(" ","_")}',
            template = recalculate_AF_improved(
                input_vcf = neutral_vcf_files_template.outputs['neutral_vcf'],
                working_dir = new_wd_fst_pi, 
                species_short = species_abbreviation(SPECIES_NAME), 
                VCF_AN = VCF_AN)
        )



    ###########################
    ###     CAlC Fst        ###
    ###########################
    if RUN_FST:
        #print(recalculate_AF_template.outputs['allele_freq'])
            # only run pairwise fst calculations, if more than one pop:
        #print(vcf_column_count(VCF_FILE))
    
        if vcf_column_count(VCF_FILE) > 10: 
            # Checking if there is more than one population
            # #print("enter if")
            fst_from_AF_template = gwf.target_from_template(
                name = f'fst_popoolation_calc_from_AF_{species_abbreviation(SPECIES_NAME)}_{AREA_TYPE.replace(" ","_")}',
                template = fst_calc_from_AF_improved(
                    allele_freq_file = recalculate_AF_template.outputs['allele_freq'],
                    working_dir = new_wd_fst_pi,
                    output_directory = new_out_fst,
                    species_short = species_abbreviation(SPECIES_NAME),
                    landcover_type = AREA_TYPE.replace(" ","_") )
            )
            
            ##################################
            ###     CAlC HUDSON Fst        ###
            ##################################
            #print(recalculate_AF_template.outputs['allele_freq'])
            fst_HUD_from_AF_template = gwf.target_from_template(
                name = f'fst_HUD_calc_from_AF_{species_abbreviation(SPECIES_NAME)}_{AREA_TYPE.replace(" ","_")}',
                template = fst_HUDSON_calc_from_AF_improved(
                    allele_freq_file = recalculate_AF_template.outputs['allele_freq'],
                    working_dir = new_wd_fst_pi, 
                    output_directory = new_out_fst,
                    species_short = species_abbreviation(SPECIES_NAME),
                    landcover_type = AREA_TYPE.replace(" ","_"))
            )


            ############################################
            ###     CAlC Hivert/Poolfstat Fst        ###
            ############################################
            
            fst_HIVERT_from_AF_template = gwf.target_from_template(
                name = f'fst_HIVERT_calc_from_ACount_{species_abbreviation(SPECIES_NAME)}_{AREA_TYPE.replace(" ","_")}',
                template = fst_HIVERT_calc_from_AF_improved(
                    allele_count_file = recalculate_AF_template.outputs['allele_count'],
                    positions_file = recalculate_AF_template.outputs['allele_positions'],
                    working_dir = new_wd_fst_pi, 
                    output_directory = new_out_fst,
                    species_short = species_abbreviation(SPECIES_NAME),
                    landcover_type = AREA_TYPE.replace(" ","_"))
            )
        else:
            print("Only one population. Not estimating Fst.")


    ##########################
    ###     CAlC pi        ###
    ##########################

    if RUN_PI:
        # need counting file for this
        pi_calculation_all_pos = gwf.target_from_template(
            name = f'pi_calculation_all_positions_{species_abbreviation(SPECIES_NAME)}_{AREA_TYPE.replace(" ","_")}',
            template = calculate_pi_template_improved(
                allele_freq_file = recalculate_AF_template.outputs['allele_freq'],
                allele_count_file = recalculate_AF_template.outputs['allele_count'],
                indels_file = INDELS_INFO_FILE,
                working_directory = new_wd_fst_pi,
                output_directory = new_out_pi,
                count_file = COUNT_POSITIONS_FILE, 
                species_short = species_abbreviation(SPECIES_NAME),
                landcover_type = AREA_TYPE.replace(" ","_"))
        )

        # Pi from greenedalf
        #neutral_bam_files_map.output['neutral_bam']
        #bam_dict_list = find_bam_files(BAM_PATH)
        # no I should assemble them as in gwf
        #print(collect(neutral_bam_files_map.outputs, ['neutral_bam']))
        
        pi_grenedalf_target = gwf.target_from_template(
            #name=create_grened_run_name,
            name=f'pi_greened_{species_abbreviation(SPECIES_NAME)}_{AREA_TYPE.replace(" ","_")}',
            template = pi_grenedalf_template(
                    neutral_bams=collect(neutral_bam_files_map.outputs, ['neutral_bam']),
                    working_directory = f'{new_wd_fst_pi}/neutral_bam_temp/',
                    ouput_dir = new_out_pi,
                    FILT_MIN_COUNT = FILT_MIN_COUNT,
                    FILT_MIN_DEPTH = FILT_MIN_DEPTH,
                    FILT_MAX_DEPTH = FILT_MAX_DEPTH,
                    FILT_WINDOW = FILT_WINDOW,
                    FILT_POOL_SIZE = FILT_POOL_SIZE,
                    species_short = species_abbreviation(SPECIES_NAME),
                    landcover_type = AREA_TYPE.replace(" ","_"))
        )
        
        

        ####################################
        ###     Add pi to collection     ###
        ####################################

        #print(classify_land_use(AREA_TYPE.replace(" ","_")))
        # need counting file for this
        #print("Collecting estimates:" + PI_COLLECTION_FILE)
        pi_collection_greenedalf = gwf.target_from_template(
            name = f'pi_nei_collection_{species_abbreviation(SPECIES_NAME)}_{AREA_TYPE.replace(" ","_")}',
            template = add_pi_to_collection_file(
                mean_file = pi_calculation_all_pos.outputs['pi_mean'],
                working_directory = new_wd_fst_pi,
                collection_output_file = PI_COLLECTION_FILE,
                taxonomy = TAXONOMY.replace(" ", "_"),
                species_short = species_abbreviation(SPECIES_NAME),
                landcover_type = AREA_TYPE.replace(" ","_"),
                python_ignore_list = EXCLUDE_POPS)
        )
        
        #print("Collecting estimates:" + PI_COLLECTION_FILE.replace(".txt", "_greenedalf.txt"))
        pi_collection_greenedalf = gwf.target_from_template(
            name = f'pi_greened_collection_{species_abbreviation(SPECIES_NAME)}_{AREA_TYPE.replace(" ","_")}',
            template = add_pi_to_collection_file(
                mean_file = pi_grenedalf_target.outputs['neutral_pi_greenedalf_mean'],
                working_directory = new_wd_fst_pi,
                collection_output_file = PI_COLLECTION_FILE.replace(".txt", "_greenedalf.txt"),
                taxonomy = TAXONOMY.replace(" ", "_"),
                species_short = species_abbreviation(SPECIES_NAME),
                landcover_type = AREA_TYPE.replace(" ","_"),
                python_ignore_list = EXCLUDE_POPS)
        )


        #print("Collecting estimates:" + PI_COLLECTION_FILE.replace(".txt", "_greenedalf.txt").replace("all_pi","all_tajD"))
        pi_collection_greenedalf = gwf.target_from_template(
            name = f'TajD_greened_collection_{species_abbreviation(SPECIES_NAME)}_{AREA_TYPE.replace(" ","_")}',
            template = add_pi_to_collection_file(
                mean_file = pi_grenedalf_target.outputs['neutral_tajD_greenedalf_mean'],
                working_directory = new_wd_fst_pi,
                collection_output_file = PI_COLLECTION_FILE.replace(".txt", "_greenedalf.txt").replace("all_pi","all_tajD"),
                taxonomy = TAXONOMY.replace(" ", "_"),
                species_short = species_abbreviation(SPECIES_NAME),
                landcover_type = AREA_TYPE.replace(" ","_"),
                python_ignore_list = EXCLUDE_POPS)
        )

        #print("Collecting estimates:" + PI_COLLECTION_FILE.replace(".txt", "_greenedalf.txt").replace("all_pi","all_wattersonTheta"))
        pi_collection_greenedalf = gwf.target_from_template(
            name = f'Watterson_greened_collection_{species_abbreviation(SPECIES_NAME)}_{AREA_TYPE.replace(" ","_")}',
            template = add_pi_to_collection_file(
                mean_file = pi_grenedalf_target.outputs['neutral_watTet_greenedalf_mean'],
                working_directory = new_wd_fst_pi,
                collection_output_file = PI_COLLECTION_FILE.replace(".txt", "_greenedalf.txt").replace("all_pi","all_wattersonTheta"),
                taxonomy = TAXONOMY.replace(" ", "_"),
                species_short = species_abbreviation(SPECIES_NAME),
                landcover_type = AREA_TYPE.replace(" ","_"),
                python_ignore_list = EXCLUDE_POPS)
        )
        
    print(" ")
    return gwf
