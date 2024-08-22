#!/bin/env python3
from gwf import AnonymousTarget
import os, glob

########################## Functions and templates ##########################




def gzip_test(fastq_file: str, filename_output: str):
	"""
	Template: Test fastqfile for corruption with gzip -t
	
	Template I/O::
	
		inputs = {}
		outputs = {}
	
	:param
	"""
	inputs = {'fastq_to_test': fastq_file}
	outputs = {'output_file': filename_output}
	options = {
		'cores': 1,
		'memory': '5g',
		'walltime': '11:00:00'
	}
	spec = f"""
	# Sources environment 										OBS EDIT:
	source /home/"$USER"/.bashrc
	conda activate fastq_test_env

	if [ "$USER" == "jepe" ]; then
		source /home/"$USER"/.bashrc
		source activate popgen ########### OBS make dedicated env
	fi
	
	echo "START: $(date)"
	echo "JobID: $SLURM_JOBID"
	echo {inputs['fastq_to_test']}

	mkdir -p {os.path.dirname(outputs['output_file'])}
	gzip -t {inputs['fastq_to_test']} &>> {outputs['output_file']}
		# Outputs error messages and stdout to outputfile

	# last bit should be in another template, since for jobs that fail, they do not go to below part.


	echo "END: $(date)"
	echo "$(jobinfo "$SLURM_JOBID")"
	"""
	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)




def check_output(tested_files_output: dict, filename_fastq: str, filename_output: str):
	"""
	Template: Checks output of gzip template above, and adds filenames from temporary files to a new file.
	
	Template I/O::
	
		inputs = {}
		outputs = {}
	
	:param
	"""
	inputs = {	'test_output': tested_files_output,
		   'fastq_file_path': filename_fastq}
	outputs = {}
	options = {
		'cores': 1,
		'memory': '5g',
		'walltime': '11:00:00'
	}
	spec = f"""
	# Sources environment 										OBS EDIT:
	source /home/"$USER"/.bashrc
	conda activate fastq_test_env

	if [ "$USER" == "jepe" ]; then
		source /home/"$USER"/.bashrc
		source activate popgen ########### OBS make dedicated env
	fi

	echo "START: $(date)"
	echo "JobID: $SLURM_JOBID"
	
	#"_TEMPErr.txt"

	
	if [ -s {inputs['test_output']} ]; then
        # The file is not-empty.
		echo CORRUPTED -  {inputs['fastq_file_path']}
		## adding filename  to file with recognizable ending
		#echo {inputs['fastq_to_test']} >> {os.path.join(testresult_info_directory, os.path.basename(inputs['fastq_to_test'])).replace(".fq.gz", "_TEMPErr.txt")}
		
		# add line with filename to file
		echo {inputs['fastq_file_path']} >> {filename_output.replace(".txt",".txt.TEMP")}

		

		# Then make directory called "erroneous" within original file directory, and move it in there.
		echo Moving fastq to erroneous directory
		mkdir -p {os.path.dirname(inputs['fastq_file_path'])}/erroneous
		mv {inputs['fastq_file_path']} {os.path.dirname(inputs['fastq_file_path'])}/erroneous/
	else
        # The file is empty.
		echo All good -  {inputs['fastq_file_path']}
	fi

	#echo {inputs['test_output']}
	#cat {os.path.join(os.path.dirname(outputs['filename_output']), "*_TEMPErr.txt")} > {outputs['filename_output'].replace(".txt",".txt.TEMP")}

	# rm {os.path.join(os.path.dirname(outputs['filename_output']), "*_TEMPErr.txt")}

	#mv {outputs['filename_output'].replace(".txt",".txt.TEMP")} {outputs['filename_output']}




	mv {outputs['output_file']} {outputs['output_file']}
	



	echo "END: $(date)"
	echo "$(jobinfo "$SLURM_JOBID")"
	"""
	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)















