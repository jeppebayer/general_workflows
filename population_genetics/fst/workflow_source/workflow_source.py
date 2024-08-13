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
	#OUTPUT_DIR: str = CONFIG['output_directory_path']
	WORK_DIR: str = CONFIG['working_directory_path']
	#VCF_BASE_FOLDER: str = CONFIG['vcf_base_folder']
	TAXONOMY: str = CONFIG['taxonomic_group']
	BED_PATH: str = CONFIG['bed_files_path']
	GENOME_PATH: str = CONFIG['reference_genome_path']
	VCF_FILES: list = CONFIG['vcf_lists']


	# --------------------------------------------------
	#                  Workflow
	# --------------------------------------------------
	
	gwf = Workflow(
		defaults={'account': ACCOUNT}
	) # what does this do?
	

	### MAKE DIRECTORIES
	
	# make directories
	new_wd=f'{WORK_DIR}/fst/{TAXONOMY}/{SPECIES_NAME.replace(" ","_")}/genome_or_annot'
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
			new_wd = f'{new_wd.replace("/fst/","/allele_frequencies/").replace("genome_or_annot","")}/'
			wd_list = [ f'{new_wd}{re.split("//|/",el)[-5]}/{re.split("//|/",el)[-4]}' for el in neutral_vcf_list ]  # via list comprehension extract pop directory name (and name of dir before) and append it to the wkdirectory path
			
			input_dict = [{'vcf_file': f, 'working_directory': p} for f, p in zip(neutral_vcf_list, wd_list)] # making combined dictionary of files and new wd paths using list comprehension and zip

		# Map over neutral vcfs to get allele frequency
			get_allele_freq_runtemplate_map = gwf.map(
					#name=make_neutral_vcf,
					template_func = extract_allele_frq,
					inputs = input_dict)

			##############################
			### GET COMMON ALLELE FRQ  ###
			##############################

			freq_files = collect(get_allele_freq_runtemplate_map.outputs, ['allele_frq_file']) # dict with allele_frq_files: [list of files]
			#print(freq_files)
			#print(list(freq_files.values())[0]) 	# to get first element of list within dictionary, get the values in from dict using values. then convert that to list. get first element of first element.[0]
			#print(len(list(freq_files.values())[0])) 	# to get first element of list within dictionary, get the values in from dict using values. then convert that to list. get first element of first element.[0]
			#print()
			new_wd = os.path.dirname(os.path.dirname(list(freq_files.values())[0][0]))
				# to get first element of list within dictionary, get the values in from dict using values. then convert that to list. get first element of first element.[0]
				# then take directory name twice, to get working directory
			count_files = len(list(freq_files.values())[0])
				# count number of files to add to outputfile
			freq_files = list(freq_files.values())[0] # made into a list

			make_neutral_bed_target = gwf.target_from_template(
				name='common_sites_allele_frq',
				template=common_sites_allele_frq(
					allele_freq_files = freq_files,
					working_directory = new_wd,
					files_count = count_files)
			)


			


			#print(neutral_vcf_list)
			#print(new_wd)
			#print(neutral_vcf_list)
			#print()
			#print(input_dict)
			#inputs_and_new_wd_dict = neutral_vcf_collect
			#result = [dict(item, **{'elem':'value'}) for item in neutral_vcf_list[0]]


			#get_allele_freqs = gwf.map(
			#	template_func = extract_allele_frq,
			#	inputs = inputs_and_new_wd_dict)


	# collect outputs from neutral_vcf_files_runtemplate_map
	# from those output should be defined as:
		# exported in : /home/anneaa/EcoGenetics/general_workflow_outputs/population_genetics/pi/
		# deciding where to put the intermediate


	
	# define target with input vcf and output neutral vcf
		# files needed:
			# bed file with genic positions, and repeats
			# vcf filtered : "/faststorage/project/EcoGenetics/general_workflow_outputs/population_genetics/vcf/collembola/Entomobrya_nicoleti/grassland/EntNic_aaRJ-C225/filtered_vcf/EntNic_aaRJ-C225.freebayes_n3_p100_minaltfrc0_minaltcnt2.bcftoolsfilter_SnpGap5_DP300_600_biallelic_AO1.vcf.gz"
				# add to config file
		# add to gwf.map() so it iterates over all the files

	#sample_site_type = f'a' #obs fix this ON MONDAY THIS IS WHERE I GOT TO
		# in this process I started listing VCF samples in config file, 
		# since this info can be derived from their path.
		# next is to access list, make new list or dict with sample_site_type information
	 #maybe it is not nessesaey, since I may just move the path into fst output

	

	return gwf
