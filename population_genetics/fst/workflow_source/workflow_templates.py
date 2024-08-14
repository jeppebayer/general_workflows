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

	bcftools index {outputs['neutral_vcf']}

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


def awk_time_column_by_number(file, outfile, column, number):
    # Convert all elements to strings and concatenate them with a space as separator
	cmd = "awk -v column=" + column + " -v number=" + number + " '{ for (x=1; x<column; x++) { print $1, $column*number}}' " + file + " > " + outfile 
	os.system (cmd)


#def awk_runner(inputfile, outputfile):
 #   cmd = "awk 'FNR<=27{print;next} ++count%10==0{print;count}' " + inputfile + " > " + outputfile
  #  os.system(cmd)



def extract_allele_frq(vcf_file: str, working_directory: str):
	"""
	Template: Extract allele frequencies using bcftools on vcf files
	
	Template I/O::
	
		inputs = {}
		outputs = {}
	
	:param
	"""
	inputs = {'neutral_vcf_file': vcf_file}
	outputs = {'allele_frq_file': f'{working_directory}/{vcf_file.split("/")[-1].replace(".vcf.gz",".frq.bed")}'}
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

	mkdir -p {working_directory}
		# working directory is input as dictionary connected to the vcf file
		# should be in my local folder, not in wf_outputs
	
	# line example:
	# HiC_scaffold_1  2875    .       AACTA   AACTT   4156.93 PASS    AB=0.442539;ABP=19.7303;AC=48;AF=0.48;AN=100;AO=258;CIGAR=4M1X;DP=583;DPB=592.4;DPRA=0;EPP=9.60888;EPPR=14.9338;GTI=1;LEN=1;MEANALT=19;MQM=31.6705;MQMR=33.9386;NS=1;NUMALT=1;ODDS=0.119269;PAIRED=0.914729;PAIREDR=0.855596;PAO=6.16667;PQA=222.167;PQR=222.167;PRO=6.16667;QA=9347;QR=9649;RO=277;RPL=256;RPP=546.013;RPPR=338.914;RPR=2;RUN=1;SAF=143;SAP=9.60888;SAR=115;SRF=112;SRP=25.0308;SRR=165;TYPE=snp       GT:DP:AD:RO:QR:AO:QA    0/0/0/0/0/0/0/0/0/0/0/0/0/0/0/0/0/0/0/0/0/0/0/0/0/0/0/0/0/0/0/0/0/0/0/0/0/0/0/0/0/0/0/0/0/0/0/0/0/0/0/0/1/1/1/1/1/1/1/1/1/1/1/1/1/1/1/1/1/1/1/1/1/1/1/1/1/1/1/1/1/1/1/1/1/1/1/1/1/1/1/1/1/1/1/1/1/1/1/1:583:277,258:277:9649:258:9347
	##INFO=<ID=AC,Number=A,Type=Integer,Description="Total number of alternate alleles in called genotypes">
	##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">
	##INFO=<ID=AF,Number=A,Type=Float,Description="Estimated allele frequency in the range (0,1]">
	##INFO=<ID=TYPE,Number=A,Type=String,Description="The type of allele, either snp, mnp, ins, del, or complex.">

	bcftools query -f '%CHROM\t%POS0\t%POS0\t%TYPE\t%AF\n' {inputs['neutral_vcf_file']} > {outputs['allele_frq_file']}

	echo "END: $(date)"
	echo "$(jobinfo "$SLURM_JOBID")"
	"""
	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)



def concatenate_list_elements(input_list):
    # Convert all elements to strings and concatenate them with a space as separator
    result = ' '.join(str(element) for element in input_list)
    return result




def common_sites_allele_frq(allele_freq_files: list, working_directory: str, files_count = int):
	"""
	Template: Get common sites for all .frq.bed files outputted from extract_allele_frq()
	
	Template I/O::
	
		inputs = {}
		outputs = {}
	
	:param
	"""
	inputs = {'freq_files': allele_freq_files}
	outputs = {'common_allele_frq_file': f'{working_directory}/{allele_freq_files[0].split("/")[-1].split(".")[-2]}.{files_count}.frq'}
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

	mkdir -p {working_directory}
	
	echo {working_directory}/{allele_freq_files[0].split("/")[-1].split(".")[-2]}.{files_count}.frq

	
	bedtools intersect -sorted -wa -wb -f 1.0 -C -a {inputs['freq_files'][0]} -b {concatenate_list_elements(inputs['freq_files'][1:])} > {outputs['common_allele_frq_file']}

	echo "END: $(date)"
	echo "$(jobinfo "$SLURM_JOBID")"
	"""
	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

# target for combining vcf's keeping only common variants among all samples
	# 	potentiallt using bcftools "isec -n+2 -p output_dir *.vcf.gz"
		# "present in n files" variants with the -n parameter
		# Ensure that your VCF is decomposed and normalized (left aligned, parsimonious representation) before you do this though, multi-allelic variants and non-normalized entries can mess up comparison.
	# use BEDtools 'intersect' for the two original VCFs. bedtools intersect -a $vcffileA -b all_vcfs*
	# use VCFtools 'vcf-annotate' to add the 1000 Genomes rs numbers, then 'grep' to keep the variants that were annotated as such.


	# now I want to use an awk. How to do that within python? define function like above? then no f string thingy