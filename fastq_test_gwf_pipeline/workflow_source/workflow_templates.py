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
	mkdir -p {os.path.dirname(outputs['output_file'])}
	gzip -t {inputs['fastq_to_test']} > {outputs['output_file']}

	if [ -s {outputs['output_file']} ]; then
        # The file is not-empty.
		echo CORRUPTED - {outputs['output_file']} 
		# add to some list of files, where filename has date and time of gwf run in it.
			# maybe by making template where all output is collected and checked

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

















