#!/bin/env python3
from gwf import AnonymousTarget
import os, glob

########################## Functions ##########################

def species_abbreviation(species_name: str) -> str:
	"""Creates species abbreviation from species name.

	:param str species_name:
		Species name written as *genus* *species*"""
	genus, species = species_name.replace(' ', '_').split('_') # first changes space to underscore, then splits at underscore. could have split at space?
	genus = genus[0].upper() + genus[1:3] # returns first character ([0]) as upper case? apparently python 1:3 does not include the last number. it means "up until"
	species = species[0].upper() + species[1:3]
	return genus + species


def make_genome_fai(ref_genome_file: str, fasta_fai_output: str):
	"""
	Template: Create fai indexed genome file for bedtools.
	
	Template I/O::
	
		inputs = {}
		outputs = {}
	
	:param
	"""
	inputs = {'reference_genome_file': ref_genome_file}
	outputs = {'genome_fai': fasta_fai_output}
	options = {
		'cores': 1,
		'memory': '5g',
		'walltime': '11:00:00'
	}
	spec = f"""
	# Sources environment 										OBS EDIT:
	source /home/"$USER"/.bashrc
	conda activate ecogen_neutral_diversity_wf
	if [ "$USER" == "jepe" ]; then
		source /home/"$USER"/.bashrc
		source activate popgen ########### OBS make dedicated env
	fi

	echo "START: $(date)"
	echo "JobID: $SLURM_JOBID"

	samtools faidx {inputs['reference_genome_file']} -o {outputs['genome_fai']}
	
	echo "END: $(date)"
	echo "$(jobinfo "$SLURM_JOBID")"
	"""
	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)



def make_neutral_bed(genes_bed_file: str, repeats_bed_file: str, neutral_bed_out: str, genome_bedstyle: str):
	"""
	Template: Create netural bed file from genes and repeats bed for species.
	
	Template I/O::
	
		inputs = {}
		outputs = {}
	
	:param
	"""
	inputs = {'genes_bed': genes_bed_file,
		   	  'repeats_bed': repeats_bed_file,
			  'genome_bed_style': genome_bedstyle}
	outputs = {'neutral_bed': neutral_bed_out}
	options = {
		'cores': 1,
		'memory': '5g',
		'walltime': '11:00:00'
	}
	spec = f"""
	# Sources environment 										OBS EDIT:
	source /home/"$USER"/.bashrc
	conda activate ecogen_neutral_diversity_wf
	if [ "$USER" == "jepe" ]; then
		source /home/"$USER"/.bashrc
		source activate popgen ########### OBS make dedicated env
	fi

	echo "START: $(date)"
	echo "JobID: $SLURM_JOBID"

	cat {inputs['genes_bed']} > {outputs['neutral_bed'].replace("neutral.","neutralTEMP.")}
	cat {inputs['repeats_bed']} >> {outputs['neutral_bed'].replace("neutral.","neutralTEMP.")}
		
	bedtools sort -i {outputs['neutral_bed'].replace("neutral.","neutralTEMP.")} \
		> {outputs['neutral_bed'].replace("neutral.","neutralTEMP.").replace(".bed", "_sort.bed")}
		
	# Now get the reverse complement, ie. only neutral and non repetitive parts.
	bedtools complement \
		-i {outputs['neutral_bed'].replace("neutral.","neutralTEMP.").replace(".bed", "_sort.bed")} \
		-g {inputs['genome_bed_style']} \
		> {outputs['neutral_bed'].replace("neutral.","neutralTEMP.")}
			# when done test running, this should end up in same folder as genes_bed_file

	    # only keep neutral file
    rm {outputs['neutral_bed'].replace("neutral.","neutralTEMP.").replace(".bed", "_sort.bed")}
	mv {outputs['neutral_bed'].replace("neutral.","neutralTEMP.")} {outputs['neutral_bed']}
	
	
	echo "END: $(date)"
	echo "$(jobinfo "$SLURM_JOBID")"
	"""
	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)






def make_neutral_vcf(vcf_file: str, working_directory: str, neutral_bed: str):
	"""
	Template: Create netural vcf file from the neutral bedfile created by make_neutral_bed function/template.
	
	Template I/O::
	
		inputs = {}
		outputs = {}
	
	:param
	"""
	inputs = {'vcf_file': vcf_file,
		   	  'neutral_bed': neutral_bed}
	outputs = {'neutral_vcf': f'{working_directory}/{vcf_file.split("/")[-1].replace(".vcf","_neutral.vcf")}'}
	# change to just adjust vcf_file, if should be placed in same folder
	options = {
		'cores': 1,
		'memory': '5g',
		'walltime': '11:00:00'
	}
	spec = f"""
	# Sources environment 										OBS EDIT:
	source /home/"$USER"/.bashrc
	conda activate ecogen_neutral_diversity_wf
	if [ "$USER" == "jepe" ]; then
		source /home/"$USER"/.bashrc
		source activate popgen ########### OBS make dedicated env
	fi

	echo "START: $(date)"
	echo "JobID: $SLURM_JOBID"

	# making directories
	mkdir -p {working_directory}
	# mkdir -p {os.path.dirname(inputs['vcf_file'])}/neutral

	# copy vcf and make bgzip it
	cp {inputs['vcf_file']}	{working_directory}/{os.path.basename(inputs['vcf_file'])}

	bgzip -d {working_directory}/{os.path.basename(inputs['vcf_file'])}
	bgzip {working_directory}/{os.path.basename(inputs['vcf_file']).replace(".vcf.gz", ".vcf")}

	# bcftools index vcf
	bcftools index {working_directory}/{os.path.basename(inputs['vcf_file'])} 
		# produces .vcf.gz.csi
	
	# run bcftools command to get neutral file
	bcftools view --with-header --output-type z \
		-R {inputs['neutral_bed']} {working_directory}/{os.path.basename(inputs['vcf_file'])} \
		> {outputs['neutral_vcf']}

	rm {working_directory}/{os.path.basename(inputs['vcf_file'])}
	rm {working_directory}/{os.path.basename(inputs['vcf_file']).replace(".vcf.gz",".vcf.gz.csi")}
		

	# clean up in template: remove this below at some point:

	#vcftools --gzvcf {inputs['vcf_file']} --bed {inputs['neutral_bed']} --temp {os.path.dirname(outputs['neutral_vcf'])} --stdout | gzip -c > {os.path.dirname(outputs['neutral_vcf'])}/$folder_name/{os.path.basename(outputs['neutral_vcf'])}
	#vcftools --gzvcf {inputs['vcf_file']} --bed {inputs['neutral_bed']} --temp {os.path.dirname(outputs['neutral_vcf'])} --stdout --recode-INFO-all --recode > {os.path.dirname(outputs['neutral_vcf'])}/$folder_name/{os.path.basename(outputs['neutral_vcf'])}
	# finds polyploidy, despite filtering?
	# I cannot get it to work with vcftools? if I recode, then it only outputs info, nothing else
	# I will postpone making it work until Jeppe is here. Since I can then get access and do bcftools instead?
	# I will continue, and trying with the non-neutral vcf to get the frequencies
	
	echo "END: $(date)"
	echo "$(jobinfo "$SLURM_JOBID")"
	"""
	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)









