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
				
				# create output file name:
				fst_output_file_mod={outputs['fst_file']}
				fst_output_file_mod=${{fst_output_file_mod/.fst/"_pops_"$((columns-$NUMBER_MAINCOLS))"_"$((columns2-$NUMBER_MAINCOLS))".fst"}}
				###################
				### Calculate fst
				###################
				echo index: $columns-$columns2 > ${{fst_output_file_mod/.fst/.temp}}
				echo pops_nr: $((columns-$NUMBER_MAINCOLS))-$((columns2-$NUMBER_MAINCOLS)) >> ${{fst_output_file_mod/.fst/.temp}}
				
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
	#paste -d'\t' {working_dir}/tmp/fst/$fst_output_file*.fst > $working_directory/$fst_output_file_new
	paste -d'\t' {working_dir}/tmp/fst/{species_short}_fst_allpairs*.temp > {outputs['fst_file']}
	#{species_short}_fst_allpairs*.temp

	###############
	## Calculate mean
	##############

	awk 'BEGIN{{ FS = OFS = "\t"}};
		NR == 1 || NR==2 || NR==3 {{print $0}}
		NR > 3 {{
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


def make_various_intergene_beds(intergenes_bed_file: str, repeats_bed_file: str, working_directory: str, genome_path: str):
	"""
	Template: Create netural bed file from genes and repeats bed for species.
	
	Template I/O::
	
		inputs = {}
		outputs = {}
	
	:param
	"""
	inputs = {'intergenes_bed': intergenes_bed_file, 
		   'bed_stype_genome':genome_path,
		   'repeats_bed': repeats_bed_file}
	outputs = {'percent_graph': f'{working_directory}/snps_neutral_region_percentage_rem.pdf',
			'nbases_graph': f'{working_directory}/snps_neutral_region_Nbases_rem.pdf'}
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

	# for each interval: 
	# 	make new intergene bed: 
	# 	make new neutral bed : 
	# 	intersect with variants called - for each pop?
	# 	count lines

	mkdir -p {working_directory}
	echo Starting R script
	../../../scripts/intergenic_length_investigations.r {inputs['intergenes_bed']} {working_directory}
	echo Finihsed R script

	#rename '.temp' '' {working_directory}/intergenic_*_percent_rem.bed.temp 	
	#rename '.temp' '' {working_directory}/intergenic_*_bases_rem.bed.temp 

	echo Making new neutral bed file and counting population variants
	# 	take the complement of repeats
	bedtools complement -i {repeats_bed_file} -g {genome_path} > {working_directory}/repeat_complement_temp.bed

	for pop_filt_vcf in ...
	do
		pop_vcf_base=`basename $pop_filt_vcf`
		echo -n {working_directory}/${{pop_vcf_base/.vcf/.varcount}}
		for percent_file in {working_directory}/intergenic_*_percent_rem.bed.temp
		do
			#	Overlaps between percentage intergene file and repeat complement
			percent_file_name=`basename $percent_file`
			bedtools intersect -a $percent_file -b {working_directory}/repeat_complement_temp.bed > {working_directory}/${{percent_file_name/.bed.temp/_neutral.bed}}
				# this is new neutral file
		
			# get overlap with population filtered VCF file, and count lines (only output lines where vcf overlaps)
			count=`bedtools intersect -a $pop_filt_vcf -b {working_directory}/${{percent_file_name/.bed.temp/_neutral.bed}} |wc -l | awk '{{print $1}}'`
			echo -e $percent_file_name'\\t'$count >> {working_directory}/${{pop_vcf_base/.vcf/.varcount}}
		done

		for bases_file in {working_directory}/intergenic_*_bases_rem.bed.temp
		do
			#	Overlaps between percentage intergene file and repeat complement
			bases_file_name=`basename $bases_file`
			bedtools intersect -a $bases_file -b {working_directory}/repeat_complement_temp.bed > {working_directory}/${{bases_file_name/.bed.temp/_neutral.bed}}
				# this is new neutral file
		
			# get overlap with population filtered VCF file, and count lines (only output lines where vcf overlaps)
			count=`bedtools intersect -a $pop_filt_vcf -b {working_directory}/${{bases_file_name/.bed.temp/_neutral.bed}} |wc -l| awk '{{print $1}}'`
			echo -e $bases_file_name'\\t'$count >> {working_directory}/${{pop_vcf_base/.vcf/.varcount}}
		done
	done

	# maybe the output should be a graph? with all pops on same graph, using alpha coloring. lines combining y values. plot(type = 'l')
	# two graphs, one with percentage, one with n bases

	../../../scripts/plot_neutral_regions_saturation.r varcount {working_directory}
		# makes the output files in temp version
	
	mv {working_directory}/snps_neutral_region_percentage_rem.pdf {outputs['percent_graph']}
	mv {working_directory}/snps_neutral_region_Nbases_rem.pdf {outputs['nbases_graph']}

	echo "END: $(date)"
	echo "$(jobinfo "$SLURM_JOBID")"
	"""
	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def make_neutral_bed_inters(intergenes_bed_file: str, repeats_bed_file: str, neutral_bed_out: str, genome_bedstyle: str):
	"""
	Template: Create netural bed file from genes and repeats bed for species.
	
	Template I/O::
	
		inputs = {}
		outputs = {}
	
	:param
	"""
	inputs = {'intergenes_bed': intergenes_bed_file,
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

	cat {inputs['intergenes_bed']} > {outputs['neutral_bed'].replace("neutral.","neutralTEMP.")}
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



def make_neutral_bed_imporved(genes_bed_file: str, repeats_bed_file: str, neutral_bed_out: str, reference_genome_file: str, working_dir: str):
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

	bcftools query -f '%CHROM\\t%POS0\\t%POS\\t%TYPE\\t%AF\\n' {inputs['neutral_vcf_file']} > {outputs['allele_frq_file']}

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
	paste -d'\\t' {working_directory}/tmp/pi/{species_short}_pi_calc_mean* > {outputs['pi_mean']}


	echo "END: $(date)"
	echo "$(jobinfo "$SLURM_JOBID")"
	"""
	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)



def modify_pi_file_template(sorted_pi_file: str, working_directory: str):
	"""
	Template: Remodel pi file
	They look like this: Scaff Start_0-type_position End_0-type_position Variant_type AlleleFrequency 
	positions_type: string to add to outputs filename. Is it all / common / another subset of positions?  eg. 'all' 'common'

	Outputs one file:
	pi per variant position: scaff pos pi_pop2	pi_pop3	pi_pop4 ...	


	Template I/O::
	
		inputs = {}
		outputs = {}
	
	:param
	"""
	inputs = {'sorted_pi_file': sorted_pi_file }
	outputs = { 'pi_all_pops_pi': f'{working_directory}/pi/pi_allPops_variant_positions.pi',
			'pi_all_pops_bed': f'{working_directory}/pi/pi_allPops_variant_positions.bed'}
	options = {
		'cores': 1,
		'memory': '120g',
		'walltime': '11:59:00'
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
	mkdir -p {working_directory}/pi
	
	echo Starting to change the format from long to wide, and output both wide and a bed-format, using awk.

	sorted_pi_file={inputs['sorted_pi_file']}
	output_wide_pi={outputs['pi_all_pops_pi'].replace(".pi", ".pi.temp")}
	output_bed_pi={outputs['pi_all_pops_bed'].replace(".bed", ".bed.temp")}

	
	awk -v working_directory={working_directory} 'NR==1 {{print $0}}; NR>1 {{print > working_directory"/pi/" $1 "_sorted.pi.tmp"}}' $sorted_pi_file


	for file in {working_directory}/pi/*scaff*_sorted.pi.tmp
	do
		echo $file
		awk 'BEGIN {{
				OFS="\\t";
				IGNORECASE = 1;  # Make comparisons case-insensitive for sorting
			}}
			# For each line of input, store values
			NR>1{{
					# Create a unique key from chrom, chromStart, chromEnd
				key = $1 OFS $2 OFS $3;
					# Store the value of pi for each pop_name under that key
				data[key][$4] = $5;
					# Track the unique pop_names to use as column headers later
				if (!seen[$4]++) {{
					pop_names[count++] = $4	}}
				# Track the unique keys (chrom, chromStart, chromEnd) for sorting rows
				if (!seen_key[key]++) {{
					keys[++key_count] = key	}}
			}}
		END {{
				# Sort pop_name headers alphabetically
			asort(pop_names);
				# Sort rows based on the keys (which are based on chrom, chromStart, chromEnd)
			asort(keys);
				# Print the header
			header = "chrom\\tchromStart\\tchromEnd";
			for (i = 1; i <= count; i++) {{
				header = header OFS pop_names[i];
			}}
			print header;
				# Print the rows in sorted order
			for (i = 1; i <= key_count; i++) {{
				key = keys[i];
				row = key;
				for (j = 1; j <= count; j++) {{
					pop = pop_names[j];
					row = row OFS (data[key][pop] ? data[key][pop] : "0");
				}}
				print row;
			}}
		}}' $file |sort -k1,1 -k2,2n -k3,3n > ${{file/_sorted.pi.tmp/_sorted_wide.pi.tmp}}
	done
	
	# cat files
	file=`ls {working_directory}/pi/*scaff*_sorted_wide.pi.tmp |head -n1`
	awk 'NR==1{{print $0 }}' $file > $output_wide_pi 
	cat {working_directory}/pi/*scaff*_sorted_wide.pi.tmp >> $output_wide_pi

	
	# make bed_file
	awk '{{OFS=FS="\\t"; 
	   	{{printf "%s\\t%d\\t%d\\tpi\\t", $1, $2, $3}}; 
		for (i = 4; i <= NF; i++){{ if(i==NF)
			{{printf "%s", $i}} else{{printf "%s,", $i}}}}; 
		printf "\\n"}}' $output_wide_pi > $output_bed_pi
	
	
	mv $output_wide_pi {outputs['pi_all_pops_pi']}
	mv $output_bed_pi {outputs['pi_all_pops_bed']}
	
	echo "END: $(date)"
	echo "$(jobinfo "$SLURM_JOBID")"
	"""
	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)




def add_context_info_pi(pi_bedfile: str, working_directory: str, output_directory: str, species_gtf: str, covered_sites_across_genome: str, neutral_bed: str, genes_bed: str, repeats_bed: str):
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
			'neutral_bed': neutral_bed,
			'genes_bed': genes_bed,
			'repeats_bed': repeats_bed, 
			'covered_sites': covered_sites_across_genome }
	outputs = { 'pi_extended_bed': f'{output_directory}/{"_".join(os.path.basename(species_gtf).split("_")[:2])}_pi_annot_allSites.pi',
			'pi_mean_extended_bed': f'{output_directory}/{"_".join(os.path.basename(species_gtf).split("_")[:2])}_pi_mean.pi' }
	options = {
		'cores': 1,
		'memory': '30g',
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
	mkdir -p {working_directory}/annotate_pi

	awk 'OFS="\\t" {{if ($3 == "CDS") {{print $1, $4-1, $5, $3}}}}' {inputs['gtf_species']} > {working_directory}/annotate_pi/genomic_CDS.bed
	awk 'OFS="\\t" {{if ($3 == "intron") {{print $1, $4-1, $5, $3}}}}' {inputs['gtf_species']} > {working_directory}/annotate_pi/genomic_intron.bed
	awk 'OFS="\\t" {{if ($3 == "exon") {{print $1, $4-1, $5, $3}}}}' {inputs['gtf_species']} > {working_directory}/annotate_pi/genomic_exon.bed


	# First bed combination should be adding all covered sites
		# sites covered sufficiently in all pops (testfile): /home/anneaa/EcoGenetics/people/Jeppe_Bayer/population_genetics/test_data/depthdist/multibam.test.merge.bed
		# Now expand to have a file with per site info
		
	awk '{{ scaf = $1; start = $2; end = $3; OFS=FS="\\t"; 
			for (i = start; i < end; i++)
			{{print scaf, i, i+1 }}}}' {inputs['covered_sites']} \
				> {working_directory}/annotate_pi/expand_covered_sites.bed.temp
	#/home/anneaa/EcoGenetics/people/anneaa/derived_dat_scripts/neutral_diversity_pipeline/fst_pi_gwf_intermediate_steps/fst_pi/Collembola/Entomobrya_nicoleti/expand_covSites.bed

		
	bedtools intersect -loj -wa -wb -sorted \
		-a {working_directory}/annotate_pi/expand_covered_sites.bed.temp \
		-b {inputs['pi_bed']} \
		> {working_directory}/annotate_pi/pi_allPops_all_positions.bed.temp

		
	# check reverse result to see if we only have variant sites that are within coverage thresholds
	bedtools intersect -loj -wa -wb -sorted \
		-b {working_directory}/annotate_pi/expand_covered_sites.bed.temp \
		-a {inputs['pi_bed']} \
		> {working_directory}/annotate_pi/pi_allPops_all_positions_inversetest.bed.temp

		
	# check if there is any missing, using the reverse (for each variant site, a site should be in the cov file, if not, something is up)
	if grep -e '-1' {working_directory}/annotate_pi/pi_allPops_all_positions_inversetest.bed.temp |head -1
	then
		echo ERROR: Lines from the variants file, is not present in the file with sufficiently covered sites!
		echo Check this file for "-1", indicating missing match from bedtools intersect: {working_directory}/annotate_pi/pi_allPops_all_positions_inversetest.bed.temp
	else
		rm {working_directory}/annotate_pi/pi_allPops_all_positions_inversetest.bed.temp
	fi

	
	# modify intersect file to become  bedfile again.
		# do this while adding 0 pi values to all positions that were within coverage thresholds - ( only variants have pi estimates here )

	pops=`awk 'BEGIN {{ FS = "\\t" }}; NF == 8 && $8 != "." {{print $7; exit}}' {working_directory}/annotate_pi/pi_allPops_all_positions.bed.temp`
	zeros=`awk 'BEGIN {{ FS = "\\t" }}; NF == 8 && $8 != "." {{print $8; exit}}' {working_directory}/annotate_pi/pi_allPops_all_positions.bed.temp| awk '{{print NF}}'`
	echo Number of populations: $zeros
		# counts the number of populations
	zeros=`awk -v number=$zeros 'BEGIN{{ for(c=1; c<=number ;c++) printf "0%s", (c < number ? ", " : ""); printf "\\n"}}'`
		# make the zero-string for the new bed file
	
	
	awk -v pops="$pops" -v zeros="$zeros" 'BEGIN {{OFS = FS = "\\t" }};
		NF == 8 && $8 != "."  {{print $4, $5, $6, $7, $8 }};
		NF == 8 && $8 == "." {{print $1, $2, $3, pops, zeros}}' \
		{working_directory}/annotate_pi/pi_allPops_all_positions.bed.temp \
		> {working_directory}/annotate_pi/pi_allPops_all_positions_addPiZeroPos.bed.temp

		

	# Now, Use bedtools to add context to input pi bedfile
	#tail -n 33339 /home/anneaa/EcoGenetics/people/anneaa/derived_dat_scripts/neutral_diversity_pipeline/fst_pi_gwf_intermediate_steps/fst_pi/Collembola/Entomobrya_nicoleti/pi_allPops_variant_positions.bed > /home/anneaa/EcoGenetics/people/anneaa/derived_dat_scripts/neutral_diversity_pipeline/fst_pi_gwf_intermediate_steps/fst_pi/Collembola/Entomobrya_nicoleti/pi_allPops_variant_positions_MOD.bed
	# Chekc if this bedfile is 0pos? it should be, but check.
	#/home/anneaa/EcoGenetics/people/anneaa/derived_dat_scripts/neutral_diversity_pipeline/fst_pi_gwf_intermediate_steps/fst_pi/Collembola/Entomobrya_nicoleti/pi_allPops_all_positions_zeroadd.bed

	bedtools annotate -i {working_directory}/annotate_pi/pi_allPops_all_positions_addPiZeroPos.bed.temp \
		-files {inputs['genes_bed'][0]} \
		{inputs['repeats_bed'][0]} \
		{working_directory}/annotate_pi/genomic_CDS.bed \
		{working_directory}/annotate_pi/genomic_intron.bed \
		{working_directory}/annotate_pi/genomic_exon.bed \
		{inputs['neutral_bed'][0]} \
		> {working_directory}/annotate_pi/pi_allPops_all_positions_addPiZeroPos_annot.bed.temp
	# column 6:11 will be: genes, repeats, CDS, intron, exon, neutral
	# if the number is above 0, the regions are overalpping	( there should only be 0.000 or 1.000)

	
	# NEXT
	# The name column should be annotation column instead:
	# 	Remove the pops from the name column they are similar in all fields after all.O
	#		potentially just store one file where the order of populations are clear. (export one field of one line)		
	#	Add column 6:11 to name column, as tab separated, adding names for annotation: gene, rep, cds, int, ex, neu
	
	awk 'BEGIN {{FS="\\t"; OFS="\\t"}} {{
		string = "";    # Initialize the string variable

		# Check if field 6 > 0 and add "gene"
		if ($6 > 0) string = "gene";
		
		# Check if field 7 > 0 and add "repeat"
		if ($7 > 0) {{
			if (string != "") string = string ", repeat";
			else string = "repeat"	}}

		# Check if field 8 > 0 and add "exon"
		if ($8 > 0) {{
			if (string != "") string = string ", exon";
			else string = "exon"	}}

		# Check if field 9 > 0 and add "intron"
		if ($9 > 0) {{
			if (string != "") string = string ", intron";
			else string = "intron"	}}

		# Check if field 10 > 0 and add "neutral"
		if ($10 > 0) {{
			if (string != "") string = string ", neutral";
			else string = "neutral"	}}

		# Check if field 11 > 0 and add "upstream"
		if ($11 > 0) {{
			if (string != "") string = string ", upstream";
			else string = "upstream"}}

		print $1, $2, $3, string, $5    # Output fields $1, $2, $3, the string, and $5
	}}' {working_directory}/annotate_pi/pi_allPops_all_positions_addPiZeroPos_annot.bed.temp \
		> {working_directory}/annotate_pi/pi_allPops_all_positions_addPiZeroPos_annot_condenced.bed.temp

	
	# calculates average pi:

	awk -v pops=$pops 'BEGIN {{FS = OFS = "\\t"}} 
		{{	split($5, values, ", ");	# Split field 5 (comma-separated values) into an array
			if (NR == 1) {{				# Initialize the sum array on the first run
				for (i = 1; i <= length(values); i++) {{
					sum[i] = 0	}}		}}
			
			# Add the values of field 5 to the corresponding positions in the sum array
			for (i = 1; i <= length(values); i++) {{
				sum[i] += values[i]		}}
			
			# Increment the count of lines
			count++		}}
		END {{
			# Calculate the average for each element in the sum array
			for (i = 1; i <= length(sum); i++) {{
				avg[i] = sum[i] / count;
				#avg[i] = sprintf("%.2f", avg[i]);  # Format the average to 2 decimal places
			}}
			
			# Join the averages into a comma-separated string
			avg_values = avg[1];
			for (i = 2; i <= length(avg); i++) {{
				avg_values = avg_values ", " avg[i];
			}}
			
			# Output the last line's fields $1, $2, $3, the calculated averages for field 5, $6 and the rest of the fields
			print pops \\n avg_values;
		}}' {working_directory}/annotate_pi/pi_allPops_all_positions_addPiZeroPos_annot_condenced.bed.temp > {outputs['pi_mean_extended_bed']}


	mv {working_directory}/annotate_pi/pi_allPops_all_positions_addPiZeroPos_annot_condenced.bed.temp {outputs['pi_extended_bed']}
	
	rm {working_directory}/annotate_pi/*.temp



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
	paste -d'\\t' {working_directory}/tmp/pi.*.tmp > {outputs['pi_all_pops']}

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
	paste -d'\\t' {working_directory}/tmp/pi_mean.*.tmp > {outputs['mean_pi_all_pops']}


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
	paste -d'\\t' {working_directory}/tmp/allele_freq.*.tmp > {outputs['AF_all_pops']}

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
	outputs = {'pop_pair_fst': f'{working_directory}/fst_pair_estimates/{output_file_name}'}
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
	mkdir -p {working_directory}/fst_pair_estimates/
	
	# For the input indices pairs:
	popul_1=`awk -v idx={pop_index_1} 'NR=1{{ print $idx; exit }}' {inputs['af_file']}`
	popul_2=`awk -v idx={pop_index_2} 'NR=1{{ print $idx; exit }}' {inputs['af_file']}`
	echo $popul_1 - $popul_2
	echo {pop_index_1} - {pop_index_2}

	# write to temporary outputfile
	echo "index: {pop_index_1}-{pop_index_2}" > {outputs['pop_pair_fst'].replace(".fst",".temp")}
	echo "colnr: {pop_index_1+1}-{pop_index_2+1}" > {outputs['pop_pair_fst'].replace(".fst",".temp")}
	

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


def paste_fst_calc_mean(fst_files: list, output_directory: str, species_short: str):
	"""
	Template: make a combined file with all variant positions fst calculations.
	
	Positions_type: string to add to outputs filename. Is it all / common / another subset of positions?  eg. 'all' 'common'
	Template I/O::
	
		inputs = {}
		outputs = {}
	
	:param
	"""
	inputs = {'fst_files': fst_files }
	# neutral_position_count should probably be given as a file some way. talk to jeppe
	outputs = { 'fst_allPos': f'{os.path.join(output_directory, f'{species_short}_fst_allPos.fst')}',
			'fst_mean': f'{os.path.join(output_directory, f'{species_short}_fst_mean.fst')}'}
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

	mkdir -p {output_directory}
	
	# add fst estimantes to one file with all pops
	paste -d'\\t' {concatenate_list_elements(inputs['fst_files'])} > {outputs['fst_allPos']}

	# calculate mean:
	awk 'BEGIN{{ FS = OFS = "\\t"}};
		NR == 1 || NR==2 {{print $0}}
		NR > 2 {{ 
			for (i=1; i<=NF; i++) {{
				if ($i != "") {{
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
	}}' {outputs['fst_allPos']} > {outputs['fst_mean']}

	
	echo "END: $(date)"
	echo "$(jobinfo "$SLURM_JOBID")"
	"""
	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


# When making plotting scripts, how to call R scrips with inputs from python cod eblocks?
	# # subprocess.call ("/pathto/MyrScript.r")
# # No, I just do it as I would from commandline (within the f sentence in gwf)
	# Rscript rscript.r input1 input2
		# with bin/bash/R	(#! /usr/bin/Rscript)
	

def fst_plots(fst_file: str, distance_file: str, output_directory: str, species_short: str):
	"""
	Template: make a combined file with all variant positions fst calculations.
	
	Positions_type: string to add to outputs filename. Is it all / common / another subset of positions?  eg. 'all' 'common'
	Template I/O::
	
		inputs = {}
		outputs = {}
	
	:param
	"""
	inputs = {'fst_file': fst_file,
		   		'distance': distance_file  }
	# neutral_position_count should probably be given as a file some way. talk to jeppe
	outputs = { 'ibd': f'{os.path.join(output_directory, f'{species_short}_fst_ibd.png')}',
			'cladogram_neigbor': f'{os.path.join(output_directory, f'{species_short}_fst_cladogram_neighbor.png')}',
			'cladogram_UPGMA': f'{os.path.join(output_directory, f'{species_short}_fst_cladogram_UPGMA.png')}'}
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

	Rscript ../../../workflow_source/scripts/fst_plots.r {inputs['fst_file']} {inputs['distance']} {outputs['ibd']} {outputs['cladogram_neigbor']} {outputs['cladogram_UPGMA']}
		# now make script
	
	echo "END: $(date)"
	echo "$(jobinfo "$SLURM_JOBID")"
	"""
	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


