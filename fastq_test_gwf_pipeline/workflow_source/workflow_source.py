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
	fastq_files_list = []
	for direct in fastq_directories_list:
		print(direct)
		fastq_files_list.extend([f for f in glob.glob(direct + '/**', recursive=True) if os.path.isfile(f) & f.endswith(".fq.gz") & (not bool(re.search("erroneous/", f )))]	)
		# somehow add if folder is not called "erroneous"
	
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


	test_fastqfiles_map = gwf.map(
		#name=make_neutral_vcf,
		template_func = gzip_test,
		inputs = input_dict,
		extra = {'testresult_info_directory': OUTPUT_DIR})

# make template that lists all files put in erroneous folders, dependent on the other one to finish
# or make template that takes in outputs of the other template, and lists the names of the non-empty outputs, with fq.gz file ending. into file with date and time in name.
	date_time = datetime.datetime.now()
	date_time = date_time.strftime("%d-%m-%Y_%H-%M")
	#print(test_fastqfiles_map.outputs)
	# print(isinstance(collect(test_fastqfiles_map.outputs, ['output_file']), dict))
	#with open(os.path.join(OUTPUT_DIR, f'erroneous_files_{date_time}.txt'), 'w') as fp:
	#	pass

	#check_output_target = gwf.map(
	#	template_func = check_output,
	#	inputs = test_fastqfiles_map.outputs,
	#	extra = { 'filename_output': os.path.join(OUTPUT_DIR, f'erroneous_files_{date_time}.txt')}
	#)
	check_output_target = gwf.target_from_template(
		name= 'check_output_cat_list',
		template = check_output(
			test_files = collect(test_fastqfiles_map.outputs, ['output_file']),
			filename_output = os.path.join(OUTPUT_DIR, f'erroneous_files_{date_time}.txt'))
	)
	

# additional templates to make:
	# possibly after first screening (being gzip screening)
		# continues only on good files.
	
	# check MD5 screening
	# Check for four line and same length of seq and quality score

	return gwf
