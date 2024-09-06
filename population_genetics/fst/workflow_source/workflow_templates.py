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

	bcftools query -f '%CHROM\t%POS0\t%POS\t%TYPE\t%AF\n' {inputs['neutral_vcf_file']} > {outputs['allele_frq_file']}

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
	
	# Deprecated 

	Template I/O::
	
		inputs = {}
		outputs = {}
	
	:param
	"""
	inputs = {'freq_files': allele_freq_files}
	outputs = {'common_sites': f'{working_directory}/{allele_freq_files[0].split("/")[-1].split(".")[-2]}.{files_count}_count_overlap.frq'}
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
	
	echo {working_directory}/{inputs['freq_files'][0].split("/")[-1].split(".")[-2]}.{files_count}.frq

	
	number_of_files={len(inputs['freq_files'])}

	bedtools intersect -sorted -wa -f 1.0 -c -a {inputs['freq_files'][0]} -b {concatenate_list_elements(inputs['freq_files'][1:])} > {outputs['common_sites'].replace(".frq", ".frq.temp")}
	
	echo `wc -l {outputs['common_sites'].replace(".frq", ".frq.temp")}`

	number_of_files={len(inputs['freq_files'])}

	awk -v number_of_files=$number_of_files '$6==(number_of_files-1) {{print $0}}' {outputs['common_sites'].replace(".frq", ".frq.temp")} > {outputs['common_sites']}

	echo After keeping only common: `wc -l {outputs['common_sites']}`

	#rm {outputs['common_sites'].replace(".frq", ".frq.temp")}


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



	

# notes from monday:
#ok .
# estimate pi:
#   I need to estimate pi on all called variants, not just neutral, as I have done so far.
#   then I can make a file like this:

        # scaf  pos1    pos2    NsitesCov   comment pia pib pic pid pif ...
        # scaf1 5       2005    2000        CDS; gene;
        # scaf1 37920   37920   1        variant; neutral .1  .001    .002


def calculate_pi_template(allele_freq_files: list, working_directory: str):
	"""
	Template: Calculate pi from bed-style files with allele frequencies.
	They look like this: Scaff Start_0-type_position End_0-type_position Variant_type AlleleFrequency 
	positions_type: string to add to outputs filename. Is it all / common / another subset of positions?  eg. 'all' 'common'

	Outputs one file:
	pi per variant position: scaff pos pi_pop2	pi_pop3	pi_pop4 ...	


	Template I/O::
	
		inputs = {}
		outputs = {}
	
	:param
	"""
	inputs = {'freq_files_list': allele_freq_files }
	outputs = { 'pi_all_pops_pi': f'{working_directory}/pi_allPops_variant_positions.pi',
			'pi_all_pops_bed': f'{working_directory}/pi_allPops_variant_positions.bed'}
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
	mkdir -p {working_directory}/tmp


	# For every file with allelefrequencies
	echo -n > {working_directory}/tmp/pi.calc.tmp
	file_counter=0
	for file in {concatenate_list_elements(inputs['freq_files_list'])}; 
	do
	    file_counter=$((file_counter+1))

		popul=`basename $file|cut -d"." -f1`
		echo $popul

    	# add pi calculated at every line
		awk -v popul=$popul -v file_counter=$file_counter '{{ AF=$5;
			pi=1-(AF)^2-(1-AF)^2;
            print $1, $2, $3, popul, pi
		}}' $file >> {working_directory}/tmp/pi.calc.tmp
	done


	# sort file according to scaffold, position and population number
	sort -k 1,1 -k2,2n -k3,3n -k4,4 {working_directory}/tmp/pi.calc.tmp > {working_directory}/tmp/pi.calc.sort.tmp 
	rm {working_directory}/tmp/pi.calc.tmp

	# add header
	sed -i '1i chrom\tchromStart\tchromEnd\tpop_name\tpi' {working_directory}/tmp/pi.calc.sort.tmp 


	## section in python

	python

	import pandas as pd
	import numpy as np

	# Read data from a file (assuming 'data.txt' is the file name)
	data_file = {working_directory}/tmp/pi.calc.sort.tmp 
	df = pd.read_csv(data_file, sep='\s+')
	#df.info()

	# Pivot the data to wide format
	wide_df = df.pivot(index=['chrom', 'chromStart', 'chromEnd'], columns='pop_name', values='pi').reset_index()
	df = []	# reducing memory impact of dataframe

	# reindex to remove nan column
	wide_df = wide_df.reindex(index=wide_df.index.difference([np.nan]), columns=wide_df.columns.difference([np.nan]))
	#wide_df.info()

	# fill nas with 0
	wide_df_fill = wide_df.fillna('0')
	#wide_df_fill.info()

	wide_df = []	# reducing memory impact of dataframe

	# sort columns
	# separate dataset to sort on pop cols
	wide_df_fill_first = wide_df_fill[['chrom','chromStart','chromEnd']]
	wide_df_fill_second = wide_df_fill[wide_df_fill.columns.difference(['chrom','chromStart','chromEnd'])].sort_index(axis = 'columns', key=lambda col: col.str.lower())
	#wide_df_fill_second.info()

	# combine split sets
	wide_df_fill = pd.concat([ wide_df_fill_first, wide_df_fill_second ], ignore_index=False, axis=1)	#reusing df name

	wide_df_fill_second = []	# reducing memory impact of dataframe

	# Write output to a file or print it
	output_file = {outputs['pi_all_pops_pi']}
	# output_file = '/home/anneaa/EcoGenetics/people/anneaa/derived_dat_scripts/neutral_diversity_pipeline/fst_pi_gwf_intermediate_steps/fst_pi/Collembola/Entomobrya_nicoleti/pi_allPops_variant_positions.pi'
	wide_df_fill.to_csv(output_file, sep='\t', index=False)

	# make comma separated list for "name" and "score" columns in bed format
	wide_df_fill_concatPi = wide_df_fill.apply(lambda row: ", ".join(row.iloc[3:].astype(str)), axis = 'columns')
	wide_df_fill_concatPops = ", ".join(wide_df_fill.columns[3:].astype(str))

	wide_df_fill = []	# reducing memory impact of dataframe

	# make new pd dataframe of score and name columns
	new_df = pd.DataFrame({{
		'name': [wide_df_fill_concatPops] * len(wide_df_fill_concatPi),  # Repeat header for all rows
		'score': wide_df_fill_concatPi}})

	wide_df_fill_concatPops = []	# reducing memory impact of dataframe
	wide_df_fill_concatPi = []	# reducing memory impact of dataframe

		
	wide_df_bed = pd.concat([ wide_df_fill_first, new_df ], ignore_index=False, axis=1)

	# write to bedfile
	output_file = {outputs['pi_all_pops_bed']}
	# output_file = '/home/anneaa/EcoGenetics/people/anneaa/derived_dat_scripts/neutral_diversity_pipeline/fst_pi_gwf_intermediate_steps/fst_pi/Collembola/Entomobrya_nicoleti/pi_allPops_variant_positions.bed'
	wide_df_bed.to_csv(output_file, sep='\t', index = False, header = False)

	exit()

	## End of python

	
	echo "END: $(date)"
	echo "$(jobinfo "$SLURM_JOBID")"
	"""
	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)



def add_context_info_pi(pi_bedfile: str, working_directory: str, output_directory: str, species_gtf: str, bed_directory: str):
	"""
	Template: rearrange pi file, so we get this format:
	Scaff pos0 pos0 type pi_pop1 pi_pop2 pi_pop3 pi_pop4 ...	

	from this format:
	Scaff Start_0-type_position End_0-type_position Variant_type AlleleFrequency pop_nr pop_name
	
	Template I/O::
	
		inputs = {}
		outputs = {}
	
	:param
	"""
	inputs = {'pi_bed': pi_bedfile,
		   	'gtf_species': species_gtf, 
			'neutral_bed': f'{bed_directory}/*genomic.neutral.bed',
			'genes_bed': f'{bed_directory}/*genomic.genes.bed',
			'repeats_bed': f'{bed_directory}/*genomic.repeats.bed' }
	outputs = { 'pi_extended_bed': f'{output_directory}/' }
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

	mkdir -p {os.path.dirname(outputs['pi_extended_bed'])}

	# first make various bed-files, if they do not already exist.
	# hmm. no, This should be done in another template.
		# already made:
		# *genome.genes.bed
		# *genome.repeats.bed
		# *genome.neutral.bed
	# or I could stream it?
		# couldn't stream it in, for some reson
	mkdir -p {working_directory}/annotate_pi

	awk 'OFS="\t" {{if ($3 == "CDS") {{print $1, $4-1, $5, $3}}}}' {inputs['gtf_species']} > {working_directory}/annotate_pi/genomic_CDS.bed
	awk 'OFS="\t" {{if ($3 == "intron") {{print $1, $4-1, $5, $3}}}}' {inputs['gtf_species']} > {working_directory}/annotate_pi/genomic_intron.bed
	awk 'OFS="\t" {{if ($3 == "exon") {{print $1, $4-1, $5, $3}}}}' {inputs['gtf_species']} > {working_directory}/annotate_pi/genomic_exon.bed


	# First bed combination should be adding all covered sites
		# sites covered sufficiently in all pops (testfile): /home/anneaa/EcoGenetics/people/Jeppe_Bayer/population_genetics/test_data/depthdist/multibam.test.merge.bed
		# Now expand to have a file with per site info
		awk '{{ scaf = $1; start = $2; end = $3; OFS=FS="\t"; 
			 for (i = start; i < end; i++)
		        {{print scaf, i, i }}}}' /home/anneaa/EcoGenetics/people/Jeppe_Bayer/population_genetics/test_data/depthdist/multibam.test.merge.bed > /home/anneaa/EcoGenetics/people/anneaa/derived_dat_scripts/neutral_diversity_pipeline/fst_pi_gwf_intermediate_steps/fst_pi/Collembola/Entomobrya_nicoleti/expand_covSites.bed

					# consider using -sorted
	bedtools intersect -loj -wa -wb -sorted \
		-a /home/anneaa/EcoGenetics/people/anneaa/derived_dat_scripts/neutral_diversity_pipeline/fst_pi_gwf_intermediate_steps/fst_pi/Collembola/Entomobrya_nicoleti/expand_covSites.bed \
		-b /home/anneaa/EcoGenetics/people/anneaa/derived_dat_scripts/neutral_diversity_pipeline/fst_pi_gwf_intermediate_steps/fst_pi/Collembola/Entomobrya_nicoleti/pi_allPops_variant_positions_MOD.bed \
		> /home/anneaa/EcoGenetics/people/anneaa/derived_dat_scripts/neutral_diversity_pipeline/fst_pi_gwf_intermediate_steps/fst_pi/Collembola/Entomobrya_nicoleti/pi_allPops_all_positions_MOD.bed
	bedtools intersect -wao -sorted \
		-a /home/anneaa/EcoGenetics/people/anneaa/derived_dat_scripts/neutral_diversity_pipeline/fst_pi_gwf_intermediate_steps/fst_pi/Collembola/Entomobrya_nicoleti/expand_covSites.bed \
		-b /home/anneaa/EcoGenetics/people/anneaa/derived_dat_scripts/neutral_diversity_pipeline/fst_pi_gwf_intermediate_steps/fst_pi/Collembola/Entomobrya_nicoleti/pi_allPops_variant_positions_MOD.bed \
		> /home/anneaa/EcoGenetics/people/anneaa/derived_dat_scripts/neutral_diversity_pipeline/fst_pi_gwf_intermediate_steps/fst_pi/Collembola/Entomobrya_nicoleti/pi_allPops_all_positions_MOD.bed

	bedtools intersect -loj -wa -wb -sorted \
		-b /home/anneaa/EcoGenetics/people/anneaa/derived_dat_scripts/neutral_diversity_pipeline/fst_pi_gwf_intermediate_steps/fst_pi/Collembola/Entomobrya_nicoleti/expand_covSites.bed \
		-a /home/anneaa/EcoGenetics/people/anneaa/derived_dat_scripts/neutral_diversity_pipeline/fst_pi_gwf_intermediate_steps/fst_pi/Collembola/Entomobrya_nicoleti/pi_allPops_variant_positions_MOD.bed \
		> /home/anneaa/EcoGenetics/people/anneaa/derived_dat_scripts/neutral_diversity_pipeline/fst_pi_gwf_intermediate_steps/fst_pi/Collembola/Entomobrya_nicoleti/pi_allPops_all_positions_MOD_inversetest.bed

		# check if there is any missing, using the reverse (for each variant site, a site should be in the cov file, if not, something is up)
	if grep -e '-1' /home/anneaa/EcoGenetics/people/anneaa/derived_dat_scripts/neutral_diversity_pipeline/fst_pi_gwf_intermediate_steps/fst_pi/Collembola/Entomobrya_nicoleti/pi_allPops_all_positions_MOD_inversetest.bed |head -1
	then
		echo ERROR: Lines from the variants file, is not present in the file with sufficiently covered sites!
		echo Check this file for "-1", indicating missing match from bedtools intersect: /home/anneaa/EcoGenetics/people/anneaa/derived_dat_scripts/neutral_diversity_pipeline/fst_pi_gwf_intermediate_steps/fst_pi/Collembola/Entomobrya_nicoleti/pi_allPops_all_positions_MOD_inversetest.bed
	else
		rm /home/anneaa/EcoGenetics/people/anneaa/derived_dat_scripts/neutral_diversity_pipeline/fst_pi_gwf_intermediate_steps/fst_pi/Collembola/Entomobrya_nicoleti/pi_allPops_all_positions_MOD_inversetest.bed
	fi

	# modify intersect file to become standard bedfile again.

	awk 'NF == 8 {print }'
	


	# Use bedtools to add context to input pi bedfile
	tail -n 33339 /home/anneaa/EcoGenetics/people/anneaa/derived_dat_scripts/neutral_diversity_pipeline/fst_pi_gwf_intermediate_steps/fst_pi/Collembola/Entomobrya_nicoleti/pi_allPops_variant_positions.bed > /home/anneaa/EcoGenetics/people/anneaa/derived_dat_scripts/neutral_diversity_pipeline/fst_pi_gwf_intermediate_steps/fst_pi/Collembola/Entomobrya_nicoleti/pi_allPops_variant_positions_MOD.bed
	# Chekc if this bedfile is 0pos? it should be, but check.
	bedtools annotate -i /home/anneaa/EcoGenetics/people/anneaa/derived_dat_scripts/neutral_diversity_pipeline/fst_pi_gwf_intermediate_steps/fst_pi/Collembola/Entomobrya_nicoleti/pi_allPops_variant_positions_MOD.bed \
		-files {inputs['genes_bed']} \
			{inputs['repeats_bed']} \
			{working_directory}/annotate_pi/genomic_CDS.bed \
			{working_directory}/annotate_pi/genomic_intron.bed \
			{working_directory}/annotate_pi/genomic_exon.bed \
			{inputs['neutral_bed']}

	# column 6:11 will be: genes, repeats, CDS, intron, exon, neutral
	# if the number is above 0, the regions are overalpping
	bedtools annotate -i /home/anneaa/EcoGenetics/people/anneaa/derived_dat_scripts/neutral_diversity_pipeline/fst_pi_gwf_intermediate_steps/fst_pi/Collembola/Entomobrya_nicoleti/pi_allPops_variant_positions_MOD.bed \
		-files /home/anneaa/EcoGenetics/BACKUP/population_genetics/reference_genomes/collembola/Entomobrya_nicoleti/annotation/EG_EntNic_05092024_genomic.repeats.bed \
		/home/anneaa/EcoGenetics/BACKUP/population_genetics/reference_genomes/collembola/Entomobrya_nicoleti/annotation/EG_EntNic_05092024_genomic.genes.bed > /home/anneaa/EcoGenetics/people/anneaa/derived_dat_scripts/neutral_diversity_pipeline/fst_pi_gwf_intermediate_steps/fst_pi/Collembola/Entomobrya_nicoleti/pi_allPops_variant_positions_MOD_annot.bed
	
		
			# if the number is above 0, the regions are overalpping
	bedtools intersect -wa -wb -a /home/anneaa/EcoGenetics/people/anneaa/derived_dat_scripts/neutral_diversity_pipeline/fst_pi_gwf_intermediate_steps/fst_pi/Collembola/Entomobrya_nicoleti/pi_allPops_variant_positions_MOD.bed \
		-b /home/anneaa/EcoGenetics/BACKUP/population_genetics/reference_genomes/collembola/Entomobrya_nicoleti/annotation/EG_EntNic_05092024_genomic.repeats.bed |head
		\
		/home/anneaa/EcoGenetics/BACKUP/population_genetics/reference_genomes/collembola/Entomobrya_nicoleti/annotation/EG_EntNic_05092024_genomic.genes.bed > /home/anneaa/EcoGenetics/people/anneaa/derived_dat_scripts/neutral_diversity_pipeline/fst_pi_gwf_intermediate_steps/fst_pi/Collembola/Entomobrya_nicoleti/pi_allPops_variant_positions_MOD_annot.bed
	
		
	# run intersect to check if it is ok


	bedtools intersect {inputs['pi_bed']} ?

	























	python

	import pandas as pd

	# Read data from a file (assuming 'data.txt' is the file name)
	data_file = '/home/anneaa/EcoGenetics/people/anneaa/derived_dat_scripts/neutral_diversity_pipeline/fst_pi_gwf_intermediate_steps/fst_pi/Collembola/Entomobrya_nicoleti/tmp/pi.header_cat.test1'
	#data_file = '/home/anneaa/EcoGenetics/people/anneaa/derived_dat_scripts/neutral_diversity_pipeline/fst_pi_gwf_intermediate_steps/fst_pi/Collembola/Entomobrya_nicoleti/pi_allPops_variant_positions_long.pi'
	
	
	#df = pd.read_csv({inputs['sorted_pi']}, sep='\s+')
	df = pd.read_csv(data_file, sep='\s+')
	df
	# Pivot the data to wide format
	wide_df = df.pivot(index=['scaffold', '0pos', '0pos.1', 'type' ], columns='pop_name', values='pi').reset_index()
	wide_df.info()
	wide_df_fill = wide_df.fillna('na')

	# Write output to a file or print it
	#output_file = {outputs['pi_all_pops_positions']}
	output_file = '/home/anneaa/EcoGenetics/people/anneaa/derived_dat_scripts/neutral_diversity_pipeline/fst_pi_gwf_intermediate_steps/fst_pi/Collembola/Entomobrya_nicoleti/pi_allPops_variant_positions.pi'
	wide_df_fill.to_csv(output_file, sep='\t', index=False)

	exit()


	echo "END: $(date)"
	echo "$(jobinfo "$SLURM_JOBID")"
	"""
	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)






# needs work:

def calculate_mean_pi_template(pi_perpos_file: list, working_directory: str, output_directory: str, neutral_position_count_file: str, positions_type: str):
	"""
	Template: Calculate pi from bed-style files with allele frequencies.
	They look like this: Scaff Start_0-type_position End_0-type_position Variant_type AlleleFrequency 
	positions_type: string to add to outputs filename. Is it all / common / another subset of positions?  eg. 'all' 'common'

	Outputs two files:
	mean pi: pi_pop1	pi_pop2	pi_pop3	pi_pop4 ...	
	pi per variant position: scaff pos pi_pop2	pi_pop3	pi_pop4 ...	


	Template I/O::
	
		inputs = {}
		outputs = {}
	
	:param
	"""
	inputs = {'freq_files': allele_freq_files, 
		   		'number_of_positions_file': neutral_position_count_file}
	# neutral_position_count should probably be given as a file some way. talk to jeppe
	outputs = {'mean_pi_all_pops': f'{output_directory}/pi_mean_allPops_{positions_type}_positions.pi',
				'pi_all_pops': f'{output_directory}/pi_allPops_{positions_type}_positions.pi'}
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
	

	# For every file with allelefrequencies
	for file in {inputs['freq_files']}; 
	do
	    # add population name to file
		popul=`basename $file|cut -d"." -f1`
		echo $popul > {working_directory}/tmp/pi.$popul.tmp
		
    	# add pi calculated at every line
		awk '{{ first_pi=1-($5)^2-(1-$5)^2;
            print first_pi
		}}' $file >> {working_directory}/tmp/pi.$popul.tmp	
	done
	
	# add site estimates to one file with all pops
	paste -d'\t' {working_directory}/tmp/pi.*.tmp > {outputs['pi_all_pops']}

	rm {working_directory}/tmp/pi.*.tmp
	
	echo "END: $(date)"
	echo "$(jobinfo "$SLURM_JOBID")"
	"""
	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

"""

	neutral_position_count=`awk 'NR=2 {{print $4}}' {inputs['number_of_positions_file']}`
		# modify which column to take, according to Jeppes file. also if line is 1 or 2, depending on header?

	# calculate mean pi across all available positions, by deviding with given positions count.
		awk -v neutral_position_count=$neutral_position_count 'NR==1 {{print $0}};
			NR>1{{ sumpi+=$1 }};
			END {{print sumpi; print NR-1; print neutral_position_count; print sumpi/neutral_position_count}}
        ' {working_directory}/tmp/pi.$popul.tmp > {working_directory}/tmp/pi_mean.$popul.tmp
    	# prints 5 rows: Population, SUM_pi, N_variable_sites, N_variable_nonvariable_sites, MEAN_pi (SUM/N_variable_nonvariable_sites)
	# this should be modified to be run on each colum of file

	rm {working_directory}/tmp/pi_mean.*.tmp
	# add mean estimates to one file with all pops
	paste -d'\t' {working_directory}/tmp/pi_mean.*.tmp > {outputs['mean_pi_all_pops']}


"""



def paste_allele_freq(allele_freq_files: list, working_directory: str, positions_type: str):
	"""
	Template: make a combined file with all variant positions allele frequencies.
	Inputs look like this: Scaff Start_0-type_position End_0-type_position Variant_type AlleleFrequency 
	output like this: Scaff Start_0-type_position pop1_AF pop2_AF ...

	Positions_type: string to add to outputs filename. Is it all / common / another subset of positions?  eg. 'all' 'common'
	Template I/O::
	
		inputs = {}
		outputs = {}
	
	:param
	"""
	inputs = {'freq_files': allele_freq_files }
	# neutral_position_count should probably be given as a file some way. talk to jeppe
	outputs = { 'AF_all_pops': f'{working_directory}/allele_freq_allPops_{positions_type}_positions.txt'}
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
	mkdir -p {working_directory}/tmp
	
	# For every file with allelefrequencies
	echo {inputs['freq_files']}
	echo {concatenate_list_elements(inputs['freq_files'])}

	for file in {concatenate_list_elements(inputs['freq_files'])}; 
	do
	    # add population name to file
		popul=`basename $file|cut -d"." -f1`
		echo $popul > {working_directory}/tmp/allele_freq.$popul.tmp
		# add allelefrequencies
		awk '{{ print $5 }}' $file >> {working_directory}/tmp/allele_freq.$popul.tmp
	done

	# add site allele frequencies to one file with all pops
	paste -d'\t' {working_directory}/tmp/allele_freq.*.tmp > {outputs['AF_all_pops']}

	rm {working_directory}/tmp/allele_freq.*.tmp
	
	echo "END: $(date)"
	echo "$(jobinfo "$SLURM_JOBID")"
	"""
	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)



def calculate_fst_template(allele_freq_file: str, working_directory: str, output_file_name: str, pop_index_1: int, pop_index_2: int):
	"""
	Template: Calculate fst from bed-style files with allele frequencies.
	Inputs look like this: Scaff Start_0-type_position End_0-type_position Variant_type AlleleFrequency 
	
	Positions_type: string to add to outputs filename. Is it all / common / another subset of positions?  eg. 'all' 'common'
	Template I/O::
	
		inputs = {}
		outputs = {}
	
	:param
	"""
	inputs = {'af_file': allele_freq_file}
	outputs = {'pop_pair_fst': f'{working_directory}/{output_file_name}'}
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
	
	# For the input indices pairs:
	popul_1=`awk -v idx={pop_index_1} 'NR=1{{ print $idx; exit }}' {inputs['af_file']}`
	popul_2=`awk -v idx={pop_index_2} 'NR=1{{ print $idx; exit }}' {inputs['af_file']}`
	echo $popul_1 - $popul_2
	echo {pop_index_1} - {pop_index_2}

	# write to temporary outputfile
	echo "index: {pop_index_1}-{pop_index_2}" > {outputs['pop_pair_fst'].replace(".fst",".temp")}
	echo "colnr: {pop_index_1+1}-{pop_index_2+1}" > {outputs['pop_pair_fst'].replace(".fst",".temp")}
	#echo '$popul_1-$popul_2' >> {outputs['pop_pair_fst'].replace(".fst",".temp")}

	# calculate fst from index columns

	#awk -v firstpop={pop_index_1} -v secondpop={pop_index_2} 'NR==1{{print $firstpop "-" $secondpop; exit
	#	}}' {inputs['af_file']} >> {outputs['pop_pair_fst'].replace(".fst",".temp")}
	
	awk -v firstpop={pop_index_1+1} -v secondpop={pop_index_2+1} '
		NR==1{{print $firstpop "-" $secondpop }};
		NR>1{{
		first_pi=1-($firstpop)^2-(1-$firstpop)^2;
		second_pi=1-($secondpop)^2-(1-$secondpop)^2;
		
		pi_within=(first_pi + second_pi)/2;

		pi_total=1-(($firstpop + $secondpop)/2)^2-(((1-$firstpop) + (1-$secondpop))/2)^2;
		
		if (pi_total=="0")  # pi tot will be 0, thus zero devision
			{{ 	fst=NA; 
				#print $firstpop, $secondpop, first_pi, second_pi, pi_within, pi_total, fst;
				print fst;
				first_pi=0; second_pi=0; pi_within=0; pi_total=0; fst=NA	}}
		else 
			{{ 	fst=(pi_total-pi_within)/pi_total;
				#print $firstpop, $secondpop, first_pi, second_pi, pi_within, pi_total, fst;
				print fst;
				first_pi=0; second_pi=0; pi_within=0; pi_total=0; fst=NA  }}
		}}' {inputs['af_file']} >> {outputs['pop_pair_fst'].replace(".fst",".temp")}
	
	mv {outputs['pop_pair_fst'].replace(".fst",".temp")} {outputs['pop_pair_fst']}
	
	echo "END: $(date)"
	echo "$(jobinfo "$SLURM_JOBID")"
	"""
	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)