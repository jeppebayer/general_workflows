#!/bin/env python3
from gwf import AnonymousTarget
import os, glob


########################## Functions ##########################

def species_abbreviation(species_name: str) -> str:
	"""Function: Creates species abbreviation from species name.

	:param str species_name:
		Species name written as *genus* *species*"""
	genus, species = species_name.replace(' ', '_').split('_') # first changes space to underscore, then splits at underscore. could have split at space?
	genus = genus[0].upper() + genus[1:3] # returns first character ([0]) as upper case? apparently python 1:3 does not include the last number. it means "up until"
	species = species[0].upper() + species[1:3]
	return genus + species


def concatenate_list_elements(input_list):
	# Convert all elements to strings and concatenate them with a space as separator
	result = ' '.join(str(element) for element in input_list)
	return result


################### Templates #################

def recalculate_AF_improved(input_vcf: str, working_dir: str, species_short: str):
	"""
	Template: Recalculate Allele frequency for each population, and output it in a pop specific temporary file. Then paste them together
	
	Template I/O::
	
		inputs = {}
		outputs = {}
	
	:param
	"""
	inputs = {'vcf_file': input_vcf}
	outputs = {'allele_freq': f'{working_dir}/tmp/allele_freq/{species_short}_AF_allpops.frq'}
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

	# info on job
	echo "START: $(date)"
	echo "JobID: $SLURM_JOBID"

	########################
	### GET ALLELE FREQ  ###
	########################

	mkdir -p {working_dir}
	mkdir -p {working_dir}/tmp/
	mkdir -p {working_dir}/tmp/allele_freq/
	
	
	#CHROM


	if [[ {inputs['vcf_file']} =~ \\.gz$ ]]; then
		pop_line=`zcat {inputs['vcf_file']}| grep --max-count 1 "^\\#[A-Z]" `
	else
		pop_line=`grep --max-count 1 "^\\#[A-Z]" {inputs['vcf_file']}`
	fi

	echo $pop_line
	pop_line=`echo $pop_line | awk -F' ' '{{ for(i=10; i<=NF; i++) printf "%s%s", $i, (i<NF ? "\\t" : "\\n") }}'`
	echo $pop_line
	

	# loop over each pop:
	# while calculating allele frequencies for each population
	for popname in $pop_line
	do
		echo " "
		echo $popname
		bcftools view -Ou -s $popname {inputs['vcf_file']} | bcftools query -f '%CHROM\t%POS0\t%POS\t%TYPE\t%AC\t%AN\t%AF\n' | awk 'FS=OFS="\t" {{AC=$5; AN=$6; if(AN==0) print $0, "NA"; else if (AN > 0) print $0, AC/AN}}' > {working_dir}/tmp/allele_freq/{species_short}_AF_$popname.tmp
	done

	################################################
	### 		APPEND ALLELE FREQ TOGETHER  	 ###
	################################################

	# Append all population AF together

		# add one file to coming outfile
	# make header file
	echo -e "CHROM\tPOS0\tPOS\tTYPE\tAC\tAN\tAF" > {working_dir}/tmp/allele_freq/{species_short}_AF_header.tmp1
	for conseq_pops in `ls {working_dir}/tmp/allele_freq/{species_short}_AF_*.tmp`
	do 
		pop_is=`echo $conseq_pops | rev | cut -d"/" -f1|rev| cut -d"." -f 1`
		if [ $conseq_pops == `ls {working_dir}/tmp/allele_freq/{species_short}_AF_*.tmp|head -n1` ]
		then
			echo "appending allele frequency (AC/AN) of pop: $conseq_pops"
				# if first pop in list: make start files
				# header:
			awk -v pop_is=$pop_is 'FS=OFS="\t" {{print $0, pop_is}}' {working_dir}/tmp/allele_freq/{species_short}_AF_header.tmp1 > {working_dir}/tmp/allele_freq/{species_short}_AF_header.tmp2
				# data file:
			cat $conseq_pops > {working_dir}/tmp/allele_freq/{species_short}_AF_allpops.tmp1
			continue
		else
			echo "appending allele frequency (AC/AN) of pop: $conseq_pops"
				# make intermediary file with only AF
			awk '{{print $8}}' $conseq_pops > ${{conseq_pops/.tmp/.tmp2}}
				
				# paste pop into the file building up
			paste -d'\t' {working_dir}/tmp/allele_freq/{species_short}_AF_allpops.tmp1 ${{conseq_pops/.tmp/.tmp2}} > {working_dir}/tmp/allele_freq/{species_short}_AF_allpops.tmp2
				# change the file back to former naming for next pop
			mv {working_dir}/tmp/allele_freq/{species_short}_AF_allpops.tmp2 {working_dir}/tmp/allele_freq/{species_short}_AF_allpops.tmp1	

				# building header as well
			awk -v pop_is=$pop_is 'FS=OFS="\t" {{print $0, pop_is}}' {working_dir}/tmp/allele_freq/{species_short}_AF_header.tmp2 > {working_dir}/tmp/allele_freq/{species_short}_AF_header.tmp1
			mv {working_dir}/tmp/allele_freq/{species_short}_AF_header.tmp1 {working_dir}/tmp/allele_freq/{species_short}_AF_header.tmp2
				# add pop to header
		fi
	done
	mv {working_dir}/tmp/allele_freq/{species_short}_AF_header.tmp2 {working_dir}/tmp/allele_freq/{species_short}_AF_header.tmp1
	cat {working_dir}/tmp/allele_freq/{species_short}_AF_header.tmp1 {working_dir}/tmp/allele_freq/{species_short}_AF_allpops.tmp1 > {outputs['allele_freq']}




	echo "END: $(date)"
	echo "$(jobinfo "$SLURM_JOBID")"
	"""
	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)



def fst_calc_from_AF_improved(allele_freq_file: str, working_dir: str, species_short: str):
	"""
	Template: Recalculate Allele frequency for each population, and output it in a pop specific temporary file. Then paste them together
	
	Template I/O::
	
		inputs = {}
		outputs = {}
	
	:param
	"""
	inputs = {'allele_freq_file': allele_freq_file}
	outputs = {'fst_file': f'{working_dir}/tmp/fst/{species_short}_fst_allpairs.fst',
				'fst_file_mean': f'{working_dir}/tmp/fst/{species_short}_fst_allpairs_mean.fst'}
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

	# info on job
	echo "START: $(date)"
	echo "JobID: $SLURM_JOBID"


	mkdir -p {working_dir}/tmp/fst/

	#################
	### Calc Fst  ###
	#################

	NUMBER_MAINCOLS=7
	echo input al file: {inputs['allele_freq_file']}
	tot_cols=`awk 'NR == 1 {{print NF; exit}}' {inputs['allele_freq_file']}`
	echo $tot_cols
	for columns in `seq $((NUMBER_MAINCOLS+1)) $tot_cols`
	do
		echo $columns
		if [ $columns == $tot_cols ]
		then
			echo At last column - break
			break
		else
			for columns2 in `seq $((columns+1)) $tot_cols`
			do
				echo $columns2
				#popul_1=`awk -v idx=$columns '{{ print $idx; exit }}' {working_dir}/tmp/fst/AF_header.tmp1`
				#popul_2=`awk -v idx=$columns2 'NR=1{{ print $idx; exit }}' {working_dir}/tmp/fst/AF_header.tmp1`
				#echo $popul_1 - $popul_2

				# create output file name:
				fst_output_file_mod={outputs['fst_file']}
				#fst_output_file_mod=$fst_output_file.fst
				#fst_output_file_mod=${{fst_output_file_mod/.fst/"_pops_"$((columns-$NUMBER_MAINCOLS))"_"$((columns2-$NUMBER_MAINCOLS))".fst"}}
				fst_output_file_mod=${{fst_output_file_mod/.fst/"_pops_"$((columns-$NUMBER_MAINCOLS))"_"$((columns2-$NUMBER_MAINCOLS))".fst"}}
				###################
				### Calculate fst
				###################
				echo index: $columns-$columns2 > ${{fst_output_file_mod/.fst/.temp}}
				echo pops_nr: $((columns-$NUMBER_MAINCOLS))-$((columns2-$NUMBER_MAINCOLS)) >> ${{fst_output_file_mod/.fst/.temp}}
				#echo pops: $popul_1-$popul_2 >> ${{fst_output_file_mod/.fst/.temp}}
				
				awk -v firstpop=$columns -v secondpop=$columns2 '
					NR==1{{print $firstpop "-" $secondpop }};
					NR>1{{
					if ($firstpop == "" || $secondpop == "" || $firstpop == "NA" || $secondpop == "NA")  # if AF is na (aka empty field), put fst=NA
						{{ 	fst=NA; 
							print fst;
							first_pi=0; second_pi=0; pi_within=0; pi_total=0; fst=NA; next	}}

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
					}}' {inputs['allele_freq_file']} >> ${{fst_output_file_mod/.fst/.temp}}

					#mv ${{fst_output_file_mod/.fst/.temp}} $fst_output_file_mod
				
			done
		fi
	done

	#fst_output_file_new=allpairs_$fst_output_file.fst
	#paste -d'\t' {working_dir}/tmp/fst/$fst_output_file*.fst	> $working_directory/$fst_output_file_new
	paste -d'\t' {working_dir}/tmp/fst/{species_short}_fst_allpairs*.temp > {outputs['fst_file']}
	#{species_short}_fst_allpairs*.temp

	###############
	## Calculate mean
	##############

	awk 'BEGIN{{ FS = OFS = "\t"}};
		NR == 1 || NR==2 || NR==3 || NR==4 {{print $0}}
		NR > 2 {{
			for (i=1; i<=NF; i++) {{
				if ($i != "" || $i != "NA" ) {{
					count[i]++
					sum[i] += $i }};
			}}
		}}
		END {{
			for (i=1; i<=NF; i++) {{
				if (i == NF) {{
					printf "%f", sum[i] / count[i]  # No tab for the last column
				}} else {{
					printf "%f\\t", sum[i] / count[i]  # Tab between columns
				}}
			}}
		}}' {outputs['fst_file']} > {outputs['fst_file_mean']}
		
		# $working_directory/${{fst_output_file_new/.fst/_mean.fst}}
	
	
	echo "END: $(date)"
	echo "$(jobinfo "$SLURM_JOBID")"
	"""
	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)





def make_neutral_bed_improved(genes_bed_file: str, repeats_bed_file: str, neutral_bed_out: str, reference_genome_file: str, working_dir: str):
	"""
	Template: Create netural bed file from genes and repeats bed for species.
	
	Template I/O::
	
		inputs = {}
		outputs = {}
	
	:param
	"""
	inputs = {'genes_bed': genes_bed_file,
		   	  'repeats_bed': repeats_bed_file,
			  'reference_genome_file': reference_genome_file}
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

	# make fasta.fai index for bedtools
	mkdir -p {working_dir}/tmp/genome_stuff
	samtools faidx {inputs['reference_genome_file']} -o {working_dir}/tmp/genome_stuff/{os.path.basename(inputs['reference_genome_file']).replace(".fna", ".fna.fai")}
	

	cat {inputs['genes_bed']} {inputs['repeats_bed']} > {outputs['neutral_bed'].replace("neutral.","neutralTEMP1.")}
	##cat {inputs['repeats_bed']} >> {outputs['neutral_bed'].replace("neutral.","neutralTEMP.")}
	
	bedtools sort -i {outputs['neutral_bed'].replace("neutral.","neutralTEMP1.")} > {outputs['neutral_bed'].replace("neutral.","neutralTEMP1.")}

	bedtools merge -i {outputs['neutral_bed'].replace("neutral.","neutralTEMP1.")} > {outputs['neutral_bed'].replace("neutral.","neutralTEMP.")}
	rm {outputs['neutral_bed'].replace("neutral.","neutralTEMP1.")}
	
	bedtools sort -i {outputs['neutral_bed'].replace("neutral.","neutralTEMP.")} \
		> {outputs['neutral_bed'].replace("neutral.","neutralTEMP.").replace(".bed", "_sort.bed")}
		
	# Now get the reverse complement, ie. only neutral and non repetitive parts.
	bedtools complement \
		-i {outputs['neutral_bed'].replace("neutral.","neutralTEMP.").replace(".bed", "_sort.bed")} \
		-g {working_dir}/tmp/genome_stuff/{os.path.basename(inputs['reference_genome_file']).replace(".fna", ".fna.fai")} \
		> {outputs['neutral_bed'].replace("neutral.","neutralTEMP.")}
			# when done test running, this should end up in same folder as genes_bed_file

		# only keep neutral file
	rm {outputs['neutral_bed'].replace("neutral.","neutralTEMP.").replace(".bed", "_sort.bed")}
	mv {outputs['neutral_bed'].replace("neutral.","neutralTEMP.")} {outputs['neutral_bed']}
	
	
	echo "END: $(date)"
	echo "$(jobinfo "$SLURM_JOBID")"
	"""
	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)




def make_neutral_vcf_improved(vcf_file: str, working_directory: str, neutral_bed: str):
	"""
	Template: Create netural vcf file from the neutral bedfile created by make_neutral_bed function/template.
		Assuming bgzipped vcf
	Template I/O::
	
		inputs = {}
		outputs = {}
	
	:param
	"""
	inputs = {'vcf_file': vcf_file,
		   	  'neutral_bed': neutral_bed}
	outputs = {'neutral_vcf': f'{os.path.dirname(vcf_file)}/neutral_vcf/{os.path.basename(vcf_file).replace(".vcf", "_neutral.vcf")}',
				'neutral_vcf_index': f'{os.path.dirname(vcf_file)}/neutral_vcf/{os.path.basename(vcf_file).replace(".vcf.gz", "_neutral.vcf.gz.csi")}'}
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
	mkdir -p {os.path.dirname(outputs['neutral_vcf'])}
	# mkdir -p {os.path.dirname(inputs['vcf_file'])}/neutral

	# copy vcf into /neutral_vcf/ and make bgzip it
	cp {inputs['vcf_file']}	{os.path.dirname(outputs['neutral_vcf'])}
	#cp {inputs['vcf_file']} {working_directory}/{os.path.basename(inputs['vcf_file'])}

	#bgzip -d {working_directory}/{os.path.basename(inputs['vcf_file'])}
	#bgzip {working_directory}/{os.path.basename(inputs['vcf_file']).replace(".vcf.gz", ".vcf")}

	# bcftools index vcf
	#bcftools index {working_directory}/{os.path.basename(inputs['vcf_file'])} 
	bcftools index {os.path.dirname(outputs['neutral_vcf'])}/{os.path.basename(inputs['vcf_file'])} 
		# produces .vcf.gz.csi
	
	# run bcftools command to get neutral file
	bcftools view --with-header --output-type z \
		-R {inputs['neutral_bed']} {os.path.dirname(outputs['neutral_vcf'])}/{os.path.basename(inputs['vcf_file'])} \
		> {outputs['neutral_vcf']}

	bcftools index {outputs['neutral_vcf']} -o {outputs['neutral_vcf_index']}
	
	echo May I remove: {os.path.dirname(outputs['neutral_vcf'])}/{os.path.basename(inputs['vcf_file'])}
	echo May I remove: {os.path.dirname(outputs['neutral_vcf'])}/{os.path.basename(inputs['vcf_file']).replace(".gz", ".gz.csi")}

	echo "END: $(date)"
	echo "$(jobinfo "$SLURM_JOBID")"
	"""
	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)





# notes from monday:
#ok .
# estimate pi:
#   I need to estimate pi on all called variants, not just neutral, as I have done so far.
#   then I can make a file like this:

		# scaf  pos1    pos2    NsitesCov   comment pia pib pic pid pif ...
		# scaf1 5       2005    2000        CDS; gene;
		# scaf1 37920   37920   1        variant; neutral .1  .001    .002


def calculate_pi_template_improved(allele_freq_file: list, working_directory: str, count_file: str, species_short: str):
	"""
	Template: Calculate pi from bed-style files with allele frequencies.
	They look like this: Scaff Start_0-type_position End_0-type_position Variant_type AlleleFrequency 
	positions_type: string to add to outputs filename. Is it all / common / another subset of positions?  eg. 'all' 'common'

	Outputs one file with sorted variant positions for all populations :
	pi per variant position: scaff pos pi


	Template I/O::
	
		inputs = {}
		outputs = {}
	
	:param
	"""
	inputs = {'al_freq_file': allele_freq_file, 
				'count_file': count_file }
	outputs = { 'pi': f'{working_directory}/tmp/pi/{species_short}_pi.pi',
				'pi_mean': f'{working_directory}/tmp/pi/{species_short}_pi_mean.pi'}
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
	mkdir -p {working_directory}/tmp/pi/
	

	# get count from file
	echo get line from count file
	line_is=`grep "within_threshold" {inputs['count_file']} | grep "all" | grep "whole_genome" | grep "intergenic_excl_repeats"`
	echo $line_is
	total_division_nr=`echo $line_is| rev | awk '{{print $1}}' | rev`
	echo $total_division_nr


	# initiate pi files
	echo initiate pi files
	#echo -n > {working_directory}/tmp/pi/{species_short}_pi_calc.tmp
	#echo -n > {working_directory}/tmp/pi/{species_short}_pi_calc_01aa.tmp
	awk '{{print $1, $2, $3}}' {inputs['al_freq_file']} > {working_directory}/tmp/pi/{species_short}_pi_calc_01aa.tmp

	# built for loop, each iteration making a pop specific pi calc.
	#file_counter=0
	count_fields=`awk 'NR==1 {{print NF; exit}}' {inputs['al_freq_file']}`
	echo number of fields $count_fields
	echo run for loop for pi
	for number in `seq 8 $count_fields`
	do
		echo $number
		popul=`awk -v number=$number 'NR==1 {{print $number}}' {inputs['al_freq_file']}`
		echo $popul
		echo $popul > {working_directory}/tmp/pi/{species_short}_pi_calc_$popul.tmp

		# add pi calculated at every line
		awk -v popul=$popul -v number=$number '
			NR == 1 {{print $number}}
			NR > 1 {{ AF=$number;
				pi=1-(AF)^2-(1-AF)^2;
				print pi
			}}' {inputs['al_freq_file']} >> {working_directory}/tmp/pi/{species_short}_pi_calc_$popul.tmp
		
		# calc mean pi
		awk -v total_division_nr=$total_division_nr '
				NR==1{{print $1}} 
				NR>1{{ sum += $1; n++ }} 
				END {{ if (total_division_nr > 0) print sum / total_division_nr; }}
				' {working_directory}/tmp/pi/{species_short}_pi_calc_$popul.tmp > {working_directory}/tmp/pi/{species_short}_pi_calc_mean_$popul.tmp 
	done

	paste -d'\\t' {working_directory}/tmp/pi/{species_short}_pi_calc_* > {outputs['pi']}
	paste -d'\\t' {working_directory}/tmp/pi/{species_short}_pi_calc_01aa.tmp {working_directory}/tmp/pi/{species_short}_pi_calc_mean* > {outputs['pi_mean']}


	echo "END: $(date)"
	echo "$(jobinfo "$SLURM_JOBID")"
	"""
	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


