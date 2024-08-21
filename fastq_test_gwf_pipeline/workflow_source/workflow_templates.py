#!/bin/env python3
from gwf import AnonymousTarget
import os, glob

########################## Functions and templates ##########################




def gzip_test(fastq_file: str, filename_output: str, testresult_info_directory: str):
	"""
	Template: Test fastqfile for corruption with gzip -t
	
	Template I/O::
	
		inputs = {}
		outputs = {}
	
	:param
	"""
	inputs = {'fastq_to_test': fastq_file,
		   'information_outdir': testresult_info_directory}
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
	mkdir -p {os.path.dirname(outputs['output_file'])}
	gzip -t {inputs['fastq_to_test']} &>> {outputs['output_file']}
		# Outputs error messages and stdout to outputfile

	if [ -s {outputs['output_file']} ]; then
        # The file is not-empty.
		echo CORRUPTED - {outputs['output_file']} 
		# adding filename  to file with recognizable ending
		echo {inputs['fastq_to_test']} >> {os.path.join(inputs['information_outdir'], os.path.basename(inputs['fastq_to_test'])).replace(".fq.gz", "_TEMPErr.txt")}
			# Will then make dependent template to concatenate those
		
		# Then make directory called "erroneous" within original file directory, and move it in there.
		mkdir -p {os.path.dirname(inputs['fastq_to_test'])}/erroneous
		mv {inputs['fastq_to_test']} {os.path.dirname(inputs['fastq_to_test'])}/erroneous/
	else
        # The file is empty.
		echo All good - {outputs['output_file']}
	fi

	
	echo "END: $(date)"
	echo "$(jobinfo "$SLURM_JOBID")"
	"""
	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)




def check_output(test_files: str, filename_output: str):
	"""
	Template: Checks output of gzip template above, and adds filenames from temporary files to a new file.
	
	Template I/O::
	
		inputs = {}
		outputs = {}
	
	:param
	"""
	inputs = {	'test_outputs' : test_files}
	outputs = {'filename_output' : filename_output}
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

	echo {inputs['test_outputs']}
	cat {os.path.join(os.path.dirname(outputs['filename_output']), "*_TEMPErr.txt")} > {outputs['filename_output']}

	echo "END: $(date)"
	echo "$(jobinfo "$SLURM_JOBID")"
	"""
	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)















