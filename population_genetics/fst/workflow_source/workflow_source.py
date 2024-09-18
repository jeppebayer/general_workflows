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
    ANNOTATION_GTF: str = CONFIG['annotation_gtf']
    VCF_FILES: list = CONFIG['vcf_lists']
    COLLECTION_SITES_FILE: str = CONFIG['collection_sites_file']
    neutral_position_count_file: str = CONFIG['count_of_positions_file']


    # --------------------------------------------------
    #                  Workflow
    # --------------------------------------------------
    
    gwf = Workflow(
        defaults={'account': ACCOUNT}
    ) # what does this do?
    

    ### MAKE DIRECTORIES
    #print("making directories")
    # make directories
    new_wd=f'{WORK_DIR}/fst_pi/{TAXONOMY}/{SPECIES_NAME.replace(" ","_")}/genome_or_annot'
    new_out_pi=f'{OUTPUT_DIR}/pi/{TAXONOMY}/{SPECIES_NAME.replace(" ","_")}/'
    new_out_fst=f'{OUTPUT_DIR}/fst/{TAXONOMY}/{SPECIES_NAME.replace(" ","_")}/'
    if not os.path.isdir(new_wd):
        os.makedirs(new_wd)

    # defines new files/directories base on config
    bed_genes = glob.glob(f'{BED_PATH}/*_genomic.genes.bed')[0]
    bed_intergenes = glob.glob(f'{BED_PATH}/*_genomic.intergenic.bed')[0]
    bed_repeats = glob.glob(f'{BED_PATH}/*_genomic.repeats.bed')[0]
    bed_neutral = f'{new_wd}/{os.path.basename(bed_genes).replace("genes","neutral")}'
    # f'{bed_repeats.replace("repeats","neutral")}'
    


    # Define target making genome.fna.fai
    make_fna_fai_target = gwf.target_from_template(
            name='make_genome_fai',
            template=make_genome_fai(
                ref_genome_file=GENOME_PATH,
                fasta_fai_output=f'{new_wd}/{GENOME_PATH.split('/')[-1].replace(".fna", ".fna.fai")}')
        )

    ########################
    ### MAKE NEUTRAL BED ###
    ########################

    # Define target with input genes and repeats bed
        # output neutral bed
    make_neutral_bed_target = gwf.target_from_template(
        name='make_neutral_bed',
        template=make_neutral_bed(
            genes_bed_file = bed_genes,
            repeats_bed_file = bed_repeats,
            genome_bedstyle = make_fna_fai_target.outputs['genome_fai'],
            neutral_bed_out = bed_neutral)
    )

    # make varying intervals of intergenic region beds
    interbed_wd = f'{WORK_DIR}/fst_pi/{TAXONOMY}/{SPECIES_NAME.replace(" ","_")}/intergene_neutrality'
    percentage_list = list(range(0,10+1))
    bp_list = list(range(0,5000+1,500))
    output_list = [ os.path.join(interbed_wd,'intergenic_' + str(n) + '_percent_rem.bed') for n in percentage_list]

    make_neutral_bed_inter = gwf.target_from_template(
        name='make_various_intergene_beds',
        template=make_various_intergene_beds(
            intergenes_bed_file = bed_intergenes,
            repeats_bed_file = bed_repeats,
            working_directory = interbed_wd,
            percentage_list = percentage_list,
            bp_list = bp_list,
            output_list = output_list,
            genome_path = make_fna_fai_target.outputs['genome_fai'])
    )
    # The several bed outputs should go in here, and we map over them to create netural beds
    # make a list of dictionaries with input intergenic beds and output neutral bed names
    bed_neutral = f'{new_wd}/{os.path.basename(bed_intergenes).replace("genes","neutral")}'
    make_neutral_bed_target_fromIntergene = gwf.target_from_template(
        name='make_neutral_bed_from_intergene',
        template=make_neutral_bed_inter(
            intergenes_bed_file = bed_intergenes,
            repeats_bed_file = bed_repeats,
            genome_bedstyle = make_fna_fai_target.outputs['genome_fai'],
            neutral_bed_out = bed_neutral)
    )

# Target to deminish VCF to only contain neutral sites
    # make list of dictionary with :new_wd + file directory name
    for GROUP in VCF_FILES:
        #if not GROUP['vcf_files_list']:
        #	continue
        # make this run on all groups, except if called "bad_samples"
        if not GROUP['vcf_files_list']:
            continue
        if GROUP['group_name'] == 'bad_samples':
            continue
        if GROUP['group_name'] == 'grassland':
            #print('yes')
            print(f'Group being submitted: {GROUP['group_name']}')
            ########################
            ### MAKE NEUTRAL VCF ###
            ########################

            files_list = GROUP['vcf_files_list'] 
            #print(files_list)
            # files_popdir=[ f'{new_wd}/{os.path.dirname(el).split("/")[-2]}' for el in files_list ]  # via list comprehension extract pop directory name and append it to the wkdirectory path
            files_popdir=[ f'{os.path.dirname(el)}/neutral_vcf/' for el in files_list ]  # via list comprehension extract pop directory name and append it to the wkdirectory path
            
            input_dict=[{'vcf_file': f, 'working_directory': p} for f, p in zip(files_list, files_popdir)] # making combined dictionary of files and new wd paths using list comprehension and zip

            neutral_vcf_files_runtemplate_map = gwf.map(
                    #name=make_neutral_vcf,
                    template_func = make_neutral_vcf,
                    inputs = input_dict,
                    extra = {'neutral_bed': make_neutral_bed_target.outputs['neutral_bed']})
        #else:
            #print("no")
            
            ########################
            ### GET ALLELE FREQ  ###
            ########################

        # Making ready for next map:
            neutral_vcf_list = neutral_vcf_files_runtemplate_map.outputs, ['neutral_vcf'] # get list of lists with output from former map
            neutral_vcf_list = [x['neutral_vcf'] for x in neutral_vcf_list[0]]	# make list of output vcf files, via accesing output. only file paths, not as dictionary or something else

        # MAKING lists of new working directories, to output i personal folder
            new_wd_alfq = f'{new_wd.replace("/fst_pi/","/allele_frequencies/").replace("genome_or_annot","")}/'
            wd_list = [ f'{new_wd_alfq}{re.split("//|/",el)[-5]}/{re.split("//|/",el)[-4]}' for el in neutral_vcf_list ]  
                # via list comprehension extract pop directory name (and name of dir before) and append it to the wkdirectory path
            
            input_dict = [{'vcf_file': f, 'working_directory': p} for f, p in zip(neutral_vcf_list, wd_list)] # making combined dictionary of files and new wd paths using list comprehension and zip

        # Map over neutral vcfs to get allele frequency
            get_allele_freq_runtemplate_map = gwf.map(
                    #name=make_neutral_vcf,
                    template_func = extract_allele_frq,
                    inputs = input_dict)
            #print("Allele freqs run")


            #################################################################################
            ### Calculate pi on allele frequencies per population on all possible sites   ###
            #################################################################################
                # per neutral site

            freq_files = collect(get_allele_freq_runtemplate_map.outputs, ['allele_frq_file']) # dict with allele_frq_files: [list of files]
            freq_files = freq_files['allele_frq_files']		# changing into list instead of dict
            new_wd_fst_pi = f'{WORK_DIR}/fst_pi/{TAXONOMY}/{SPECIES_NAME.replace(" ","_")}/'

            pi_calculation_all_pos = gwf.target_from_template(
                name = 'pi_calculation_all_positions',
                template = calculate_pi_template(
                    allele_freq_files = freq_files,
                    working_directory = new_wd_fst_pi)
            )
                # outputs two files: 
                # 		one with: scaf 0pos 0pos pi1 pi2 pi3 pi4 ...
                # 		one with: scaf 0pos 0pos 'name' and 'score'
                #				name consists of a comma separated list of population names
                #				score consists of a comma separated list of pi scores

            pi_remodelling = gwf.target_from_template(
                name = 'pi_remodelling_file',
                template = modify_pi_file_template(
                    sorted_pi_file = pi_calculation_all_pos.outputs['sorted_pi_file'],
                    working_directory = new_wd_fst_pi)
            )
            #print("Remodelled pi")

            #  change to bed file combination thingy
            print(make_neutral_bed_target.outputs['neutral_bed'])
            print(glob.glob(os.path.join(BED_PATH,'*genomic.genes.bed'))[0])
            print(glob.glob(os.path.join(BED_PATH,'*genomic.repeats.bed'))[0])
            print(pi_remodelling.outputs['pi_all_pops_bed'])
            print(ANNOTATION_GTF)
            print(GENOME_PATH)
            pi_rearrangement_all_pos = gwf.target_from_template(
                name = 'pi_add_context',
                template = add_context_info_pi(
                    pi_bedfile = pi_remodelling.outputs['pi_all_pops_bed'],
                    working_directory = new_wd_fst_pi,
                    output_directory = new_out_pi, 
                    species_gtf = ANNOTATION_GTF, 
                    neutral_bed = make_neutral_bed_target.outputs['neutral_bed'],
                    genes_bed = glob.glob(os.path.join(BED_PATH,'*genomic.genes.bed'))[0],
                    repeats_bed = glob.glob(os.path.join(BED_PATH,'*genomic.repeats.bed'))[0],
                    covered_sites_across_genome = f'/home/anneaa/EcoGenetics/people/Jeppe_Bayer/population_genetics/test_data/depthdist/multibam.test.merge.bed')
            )
                # modify

            # ADD PLOTTING SCRIPT?


            ################################################################
            ### Make ready for FST-calculation with a combined AF file   ###
            ################################################################
                # per neutral site
            # input : AF files
                # The output bed-style file from all files output from get_allele_freq_runtemplate_map
                    # Scaff 0-type_position 0-type_position variant_type AlleleFrequency 
            # wd:			new_wd_alfq
            paste_allele_freqs_all_pops = gwf.target_from_template(
                name = 'paste_af_files',
                template = paste_allele_freq(
                    allele_freq_files = freq_files,
                    working_directory = new_wd_alfq,
                    positions_type = 'all')
            )

            # OBS may be needed for fst calc.
            # But maybe keep it for potential construct?


            ################################################################
            ### Estimate FST for all sites based on combined AF file	 ###
            ################################################################
                # per neutral site
            # input : AF files
                # The output bed-style file from all files output from get_allele_freq_runtemplate_map
                    # Scaff 0-type_position 0-type_position variant_type AlleleFrequency 
            
            # define list of outputfiles:
            counter = 0
            outputlist_calculate_fst_pairs_target = []
            for number1 in range(len(freq_files)):
                #print(counter)
                for number2 in range(number1+1, len(freq_files)):
                    #print(number1)
                    #print(number2)
                    counter = counter + 1
                    iteration_number_pad = f"{counter:0>3}" 
                    #print(iteration_number_pad)
                    calculate_fst_pairs_target = gwf.target_from_template(
                        name = f'fst_estimation_{iteration_number_pad}',
                        template = calculate_fst_template(
                            allele_freq_file = paste_allele_freqs_all_pops.outputs['AF_all_pops'],
                            pop_index_1 = number1,
                            pop_index_2 = number2,
                            working_directory = new_wd_fst_pi,
                            output_file_name = f'fst_pair_{iteration_number_pad}_{number1}_{number2}_all_positions.fst')
                        )
                    outputlist_calculate_fst_pairs_target.append(calculate_fst_pairs_target.outputs)

            
            # paste fst. estimates together into one bigger file and calculate mean

            # collect makes it into dict, the last [] takes the dict element with a list in it, and returns the list
            fst_files_list = collect(outputlist_calculate_fst_pairs_target, ['pop_pair_fst'])['pop_pair_fsts']

            paste_fst_files = gwf.target_from_template(
                name = 'paste_fst_files',
                template = paste_fst_calc_mean(
                    fst_files = fst_files_list,
                    output_directory = new_out_fst,
                    species_short = species_abbreviation(SPECIES_NAME))
            )


            # Add plotting script?
                # For isolation by distance plot
                    #    add distance data - I already have this file somewhere, it should be placed in a logical way
                # Cladogram
                # needs mending:
            #plot_fst = gwf.target_from_template(
             #   name = 'fst_plots',
              #  template = fst_plots(
               #     fst_file = paste_fst_files.outputs['fst_allPos'],
                #    distance_file = COLLECTION_SITES_FILE, 
                 #   output_directory = new_out_fst,
                  #  species_short = species_abbreviation(SPECIES_NAME))
            #)



        #next 
        # Q: do we want to calculate pi/fst on scaffolds etc?
            # then I will need to add such awk calculations in a new or one of the existing template


    

    return gwf
