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
			# If input is reverse read (_1.fq.gz), check if fw is present (_2.fq.gz), if that is the case, move it
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
	echo.> {outputs['check_done']}

	echo "END: $(date)"
	echo "$(jobinfo "$SLURM_JOBID")"
	"""
	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def check_file_size_return_memor_function(file: str):
	size_is = os.path.getsize(file)
	if size_is <= 10000000000:
		memory_is = '16G'
	elif  size_is > 80000000000:
		memory_is = '90G'
	elif  size_is > 70000000000:
		memory_is = '80G'
	elif  size_is > 60000000000:
		memory_is = '70G'
	elif  size_is > 50000000000:
		memory_is = '60G'
	elif  size_is > 40000000000:
		memory_is = '50G'
	elif  size_is > 30000000000:
		memory_is = '40G'
	elif  size_is > 20000000000:
		memory_is = '30G'
	else:
		memory_is = '20G'
	return memory_is



def check_fq_integrity_pairedend(forward: str, reverse: str, test_output: str, test_summary_file: str):
	"""
	Template: Checks output of gzip template above, and adds filenames from temporary files to a new file.
	
	Template I/O::
	
		inputs = {}
		outputs = {}
	
	:param
	"""
	inputs = {	'forward_read': forward,
		   'reverse_read': reverse}
	outputs = {'check_done_file': test_output}
	options = {
		'cores': 1,
		'memory': check_file_size_return_memor_function(forward),
		'walltime': '11:00:00'
	}
	spec = f"""
	# Sources environment 										OBS EDIT:
	source /home/"$USER"/.bashrc
	conda activate fastq_test_env

		# some ran fine at 16G, others did ok with 20g. pairs with a combined size around 60G failed at 20
			# I could make a dynamic setting?
	if [ "$USER" == "jepe" ]; then
		source /home/"$USER"/.bashrc
		source activate popgen ########### OBS make dedicated env
	fi

	echo "START: $(date)"
	echo "JobID: $SLURM_JOBID"
	
	fastq_info {inputs['forward_read']} {inputs['reverse_read']} > {outputs['check_done_file']}

	if grep -q " " {outputs['check_done_file']}; then
		if [ -s {outputs['check_done_file']} ]; then
			# The file is not-empty.
			echo File or file pair not in regular fastq format - {inputs['forward_read']} {inputs['reverse_read']}
			#exit 
			
			# add line with filename to file
			echo "JobID: $SLURM_JOBID" {inputs['forward_read']} {inputs['reverse_read']} >> {test_summary_file}
		fi
	else
		# The file only contain line breaks or something, but not text. However it is not truly empty.
		# # make the file really be empty
		echo -n > {outputs['check_done_file']}
		echo All good -  {inputs['forward_read']} {inputs['reverse_read']} 
	fi

	echo "END: $(date)"
	echo "$(jobinfo "$SLURM_JOBID")"
	"""
	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)





def check_fq_integrity_singleend(forward: dict, test_output: str, test_summary_file: str):
	"""
	Template: Checks output of gzip template above, and adds filenames from temporary files to a new file.
	
	Template I/O::
	
		inputs = {}
		outputs = {}
	
	:param
	"""
	inputs = {	'forward_read': forward}
	outputs = {'check_done_file': test_output}
	options = {
		'cores': 1,
		'memory': '30g',
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
	
	fastq_info {inputs['forward_read']} > {outputs['check_done_file']}

	if grep -q " " {outputs['check_done_file']}; then
		if [ -s {outputs['check_done_file']} ]; then
			# The file is not-empty.
			echo File or file pair not in regular fastq format - {inputs['forward_read']}
			#exit 
			
			# add line with filename to file
			echo "JobID: $SLURM_JOBID" {inputs['forward_read']} >> {test_summary_file}
		fi
	else
		# The file only contain line breaks or something, but not text. However it is not truly empty.
		# # make the file really be empty
		echo -n > {outputs['check_done_file']}
		echo All good -  {inputs['forward_read']}
	fi

	echo "END: $(date)"
	echo "$(jobinfo "$SLURM_JOBID")"
	"""
	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)










