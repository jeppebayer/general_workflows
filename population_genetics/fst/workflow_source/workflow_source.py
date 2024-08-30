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
	VCF_FILES: list = CONFIG['vcf_lists']
	neutral_position_count_file: str = CONFIG['count_of_positions_file']


	# --------------------------------------------------
	#                  Workflow
	# --------------------------------------------------
	
	gwf = Workflow(
		defaults={'account': ACCOUNT}
	) # what does this do?
	

	### MAKE DIRECTORIES
	
	# make directories
	new_wd=f'{WORK_DIR}/fst_pi/{TAXONOMY}/{SPECIES_NAME.replace(" ","_")}/genome_or_annot'
	new_out_pi=f'{OUTPUT_DIR}/pi/{TAXONOMY}/{SPECIES_NAME.replace(" ","_")}/'
	new_out_fst=f'{OUTPUT_DIR}/fst/{TAXONOMY}/{SPECIES_NAME.replace(" ","_")}/'
	if not os.path.isdir(new_wd):
		os.makedirs(new_wd)

	# defines new files/directories base on config
	bed_genes = glob.glob(f'{BED_PATH}/*_genome.genes.bed')[0]
	bed_repeats = glob.glob(f'{BED_PATH}/*_genome.repeats.bed')[0]
	bed_neutral = f'{new_wd}/{os.path.basename(bed_genes).replace("genes","neutral")}'
	# f'{bed_repeats.replace("repeats","neutral")}'
	


	# Define target making genome.fna.fai
	make_fna_fai_target = gwf.target_from_template(
			name='make_genome_fai',
			template=make_genome_fai(
				ref_genome_file=GENOME_PATH,
				fasta_fai_output=f'{new_wd}/{GENOME_PATH.split('/')[-1].replace(".fna", ".fna.fai")}')
		)

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

# Target to deminish VCF to only contain neutral sites
# iterating files using .map

	# make list of dictionary with :new_wd + file directory name
	for GROUP in VCF_FILES:
		#if not GROUP['vcf_files_list']:
		#	continue

		if GROUP['group_name'] == 'good_samples':
			#print('yes')

			########################
			### MAKE NEUTRAL VCF ###
			########################

			files_list = GROUP['vcf_files_list'] # access the good samples GROUP and get the list of vcf files
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
			
			pi_rearrangement_all_pos = gwf.target_from_template(
				name = 'pi_rearrangement_all_pos',
				template = long_to_wide_pi(
					pi_sorted_file = pi_calculation_all_pos.outputs['pi_all_pops'],
					neutral_position_count_file = os.path.join(new_out_pi/os.path.basename(pi_calculation_all_pos.outputs['pi_all_pops']).replace("_long.pi", ".pi")))
			)


			#pi_rearrangement_all_pos = gwf.target_from_template(
			#	name = 'pi_rearrangement_all_pos',
			#	template = calculate_pi_template(
			#		pi_sorted_file = pi_calculation_all_pos.outputs['pi_all_pops'],
			#		neutral_position_count_file = new_out_pi/os.path.basename(pi_calculation_all_pos.outputs['pi_all_pops']).replace("_long.pi", ".pi"))
			#)
			#Make this template into a "add non-variable positions as well" 
			# potentially using a bed-file? cat and sort.
				# could maybe more successfully be added in the earlier step?
					# what about the na's in this file. should they indeed be 0?
					# if the variant does not exist in the population, then I guess it should count as a 0-pi-position, because it does not add to pi
					# if this is the case, which I think right now,then colMeans should be na.rm=FALSE. indicating that na is included in the calculation. or to be sure do: colSums/nrow




			# 
			# 
			# 
			# 
			# 	#pi_calculation_all_pos = gwf.target_from_template(
			#	name = 'pi_calculation_all_positions',
			#	template = calculate_pi_template(
			#		allele_freq_files = freq_files,
			#		working_directory = new_wd_fst_pi,
			#		output_directory = new_out_pi,
			#		neutral_position_count_file = neutral_position_count_file,
			#		positions_type = 'all')
			#)


			# OR: make this one without calculating the mean.
			# Calculate mean in new template
			# Done, but needs the inputfile of neutral_position_count (maybe make a test)


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
			#for file in freq_files:
			#	popul = os.path.basename( file ).split(".")[0]
			#	print(popul)
			#	counter1 = freq_files.index(file) # get file position in list
			#	print(counter1)
			#print(paste_allele_freqs_all_pops.outputs['AF_all_pops'])
			counter = 0
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

			# Done
			
			
		#next 
		# make template to add all fst_files together, and calculate mean fst
		# Q: do we want to calculate pi/fst on scaffolds etc?
			# then I will need to add such awk calculations in a new or one of the existing template
	

	return gwf
