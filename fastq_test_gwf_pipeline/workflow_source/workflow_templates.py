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
	gzip -t {inputs['fastq_to_test']} &>> {outputs['output_file']} || exit 0
		# Outputs error messages and stdout to outputfile, and make sure it only exits with code 0, so dependencies will run anyway.

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
	outputs = {'check_done': tested_files_output.replace(".txt","_checkout.txt")}
	options = {
		'cores': 1,
		'memory': '1g',
		'walltime': '04:00:00'
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
		
		
		# add line with filename to file
		echo {inputs['fastq_file_path']} >> {filename_output}


		# Then make directory called "erroneous" within original file directory, and move it in there.
		echo Moving fastq to erroneous directory
		mkdir -p {os.path.dirname(inputs['fastq_file_path'])}/erroneous
		mv {inputs['fastq_file_path']} {os.path.dirname(inputs['fastq_file_path'])}/erroneous/

			# If input is forward read (_1.fq.gz), check if rw is present (_2.fq.gz), if that is the case, move it
		if [ "$(echo {inputs['fastq_file_path']} | grep -cim1 "_1.fq.gz")" -eq 1 ]; then
			echo file is forward read, checking if reverse read is present, and moving it.
			if [ -f {inputs['fastq_file_path'].replace("_1.fq.gz", "_2.fq.gz")} ]; then
				mv {inputs['fastq_file_path'].replace("_1.fq.gz", "_2.fq.gz")} {os.path.dirname(inputs['fastq_file_path'])}/erroneous/
			fi	
		elif [ "$(echo {inputs['fastq_file_path']} | grep -cim1 "_2.fq.gz")" -eq 1 ]; then
			echo file is reverse read, checking if forward read is present, and moving it.
			if [ -f {inputs['fastq_file_path'].replace("_2.fq.gz", "_1.fq.gz")} ]; then
				mv {inputs['fastq_file_path'].replace("_2.fq.gz", "_1.fq.gz")} {os.path.dirname(inputs['fastq_file_path'])}/erroneous/
			fi
		fi
	else
        # The file is empty.
		echo All good -  {inputs['fastq_file_path']}
	fi
	echo '' > {outputs['check_done']}

	echo "END: $(date)"
	echo "$(jobinfo "$SLURM_JOBID")"
	"""
	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)















