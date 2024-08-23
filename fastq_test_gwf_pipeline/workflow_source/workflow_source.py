#!/bin/env python3
from gwf import Workflow
from gwf.workflow import collect
import os, yaml, glob, sys, re, datetime
from workflow_templates import *

def fastq_test_wf(config_file: str = glob.glob('*config.y*ml')[0]):
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
	#SPECIES_NAME: str = CONFIG['species_name']
	OUTPUT_DIR: str = CONFIG['output_directory_path']
	WORK_DIR: str = CONFIG['working_directory_path']
	#VCF_BASE_FOLDER: str = CONFIG['vcf_base_folder']
	#TAXONOMY: str = CONFIG['taxonomic_group']
	#BED_PATH: str = CONFIG['bed_files_path']
	#GENOME_PATH: str = CONFIG['reference_genome_path']
	fastq_directories_list: list = CONFIG['fastq_folder_list']


	# --------------------------------------------------
	#                  Workflow
	# --------------------------------------------------
	
	gwf = Workflow(
		defaults={'account': ACCOUNT}
	) # what does this do?
	

	### MAKE DIRECTORIES
	
	# wd
	if not os.path.isdir(WORK_DIR):
		os.makedirs(WORK_DIR)
	# od
	if not os.path.isdir(OUTPUT_DIR):
		os.makedirs(OUTPUT_DIR)



	
	### Find fastq files in directories
	#print(fastq_directories_list)
	print("Checking for fastq files (.fq.gz) recursively in these directories:")
	fastq_files_list = []
	for direct in fastq_directories_list:
		print(direct)
		fastq_files_list.extend([f for f in glob.glob(direct + '/**', recursive=True) if os.path.isfile(f) & f.endswith(".fq.gz") & (not bool(re.search("erroneous/", f )))]	)
		# somehow add if folder is not called "erroneous" maybe above works
	
	## modify the paths and file names for output
	pathlist = [os.path.dirname(file).split("BACKUP/")[-1] for file in fastq_files_list]
	out_pathlist = [ OUTPUT_DIR + '/' + pathis for pathis in pathlist]
	out_path_filelist = [os.path.basename(file).replace('.fq.gz', '.txt') for file in fastq_files_list]
	outdir_file_list = [os.path.join(dirname, filename) for dirname, filename in zip(out_pathlist, out_path_filelist)]
	#print(os.path.dirname(fastq_files_list[10]))
	#print(fastq_files_list[10])
	#print(outdir_file_list)

	## Make dictionary of input files, and modified file path for output folder
	input_dict = [{'fastq_file': f, 'filename_output': d} for f, d in zip(fastq_files_list, outdir_file_list)]
	#print(input_dict)
		# fastq_file
		# filename_output

	test_fastqfiles_map = gwf.map(
		#name=make_neutral_vcf,
		template_func = gzip_test,
		inputs = input_dict,
		extra = {})

# make template that lists all files put in erroneous folders, dependent on the other one to finish
# or make template that takes in outputs of the other template, and lists the names of the non-empty outputs, with fq.gz file ending. into file with date and time in name.
	date_time = datetime.datetime.now()
	date_time = date_time.strftime("%d-%m-%Y_%H-%M")

	# change keys of input dict
		# tested_files_output
		# filename_fastq
	input_dict_test = input_dict
	#print(input_dict_test[1])
	for elem in input_dict_test:
		elem['filename_fastq'] = elem.pop('fastq_file')
		elem['tested_files_output'] = elem.pop('filename_output')

	#print(input_dict_test[1])
	#input_dict_test['filename_fastq'] = input_dict_test.pop('fastq_file')
	#input_dict_test['tested_files_output'] = input_dict_test.pop('filename_output')
	
		# del dictionary[old_key]
	check_output_target = gwf.map(
		#name=make_neutral_vcf,
		template_func = check_output,
		inputs = input_dict_test,
		extra = {'filename_output': os.path.join(OUTPUT_DIR, f'erroneous_files_{date_time}.txt')})



	print("Preparing for second check for fastq integrity:")
	fastq_files_list_1 = [name for name in fastq_files_list if '_1.fq.gz' in name]
	fastq_files_list_2 = [name for name in fastq_files_list if '_2.fq.gz' in name]
	pathlist = [os.path.dirname(file).split("BACKUP/")[-1] for file in fastq_files_list_1]
	out_pathlist = [ OUTPUT_DIR + '/' + pathis for pathis in pathlist]
	out_new_name = [os.path.basename(file).replace('.fq.gz', '_fqInteg.txt') for file in fastq_files_list_1]
	outdir_file_list_test2 = [os.path.join(dirname, filename) for dirname, filename in zip(out_pathlist, out_new_name)]

	# obs. THese lists needs to be sorted. they are not pairing up correctly!

	input_dict_test2 = [{'forward': f, 'reverse': r, 'test_output': o} for f, r, o in zip(fastq_files_list_1.sort(), fastq_files_list_2.sort(), outdir_file_list_test2.sort())]
	
	 # fastq_info file_1.fastq.gz file_2.fastq.gz    
	check_output_target = gwf.map(
		#name=make_neutral_vcf,
		template_func = check_fq_integrity,
		inputs = input_dict_test2,
		extra = {'test_summary_file': os.path.join(OUTPUT_DIR, f'integrity_issues_file_{date_time}.txt')})


# additional templates to make:
	# possibly after first screening (being gzip screening)
		# continues only on good files.
	
	# check MD5 screening
	# Check for four line and same length of seq and quality score

# add template to check if all reverse and forward pairs are still present in main folders. if not, move loner to erroneous folder.

	# conda install bioconda::fastq_utils

	return gwf
