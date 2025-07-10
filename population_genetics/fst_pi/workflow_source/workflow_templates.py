#!/bin/env python3
from gwf import AnonymousTarget
import os, glob
import pandas as pd
import gzip
#from scripts import collect_pi_estimates
#os.environ['OPENBLAS_NUM_THREADS'] = '18'

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

def find_bam_files(directory):
    bam_file_list = []
    for root, dirs, files in os.walk(directory):
        for file in files:
            if file.endswith('.bam'):
                if "filtered" in file:
                    full_path = os.path.join(root, file)
                    bam_file_list.append({'bam_file': full_path})
    return bam_file_list

def get_dict_key_elements(d: dict, key: str) -> str:
    """Returns all elements of d[key] as a space-separated string."""
    try:
        values = d[key]
        if isinstance(values, (list, tuple, set)):
            return " ".join(str(v) for v in values)
        else:
            return str(values)
    except KeyError:
        return f"Key '{key}' not found in dictionary."
	

def pi_mean_greened(input_file_greened: str, species_short: str, landscape_type: str, outplace:str) -> str:
    """Returns sum and mean of pi, wattersons theta and tajimas D, and exports them to each file with said name."""
    
    import pandas as pd
    # Load the input TSV
    df = pd.read_csv(input_file_greened, sep="\t")
    
    # Clean column headers: remove specific string from column names
    clean_str = "filtered_"  # ← replace this with the string you want to remove
    df.columns = [col.replace(clean_str, "") for col in df.columns]
    
    # Define the patterns to search for
    patterns = ["theta_pi", "theta_watterson", "tajimas_d"]
    
    # Process each pattern
    for pattern in patterns:
        # Select matching columns
        matching_cols = [col for col in df.columns if pattern in col]
        
        if not matching_cols:
            print(f"No columns matched for pattern: {pattern}")
            continue
        
        # Subset the DataFrame
        subset = df[matching_cols]
        
        # Compute sum and mean
        summed = subset.sum()
        averaged = subset.mean()
        counted = subset.count()
        
        # Combine into a new DataFrame
        result_df = pd.DataFrame([counted, summed, averaged], index=["count", "sum", "average"])
        
        # Write to output
        result_df.to_csv(f"{outplace}/{species_short}_{landscape_type}_grenedalf_{pattern}_summary.tsv", sep="\t")
        print(f"Summary written to {outplace}/{species_short}_{landscape_type}_grenedalf_{pattern}_summary.tsv, later changed to : {outplace}/{species_short}_{landscape_type}_grenedalf_{pattern}_mean.tsv")


def classify_land_use(input_str):
    # Normalize: lowercase and remove underscores/spaces
    normalized = input_str.lower().replace("_", "").replace(" ", "")
    
    if normalized == "conservationagriculture":
        return "CA"
    elif normalized == "conventionalagriculture":
        return "K"
    elif normalized == "grassland":
        return "GR"
    else:
        return None  # or raise an error if unknown    

import pandas as pd



def parse_indel_distribution(stats_file):
    """
    Parses indelLength lines from bcftools stats output and computes average indel length.

    Returns:
        df (pd.DataFrame): Table of indel sizes and counts
        avg_length (float): Weighted average indel length
    """
    indel_data = []
    with open(stats_file) as f:
        for line in f:
            if line.startswith("IDD"):
                #print(line)
                parts = line.strip().split()
                length = abs(int(parts[2]))
                count = int(parts[3])
                #print(count)
                indel_data.append((length, count))
                #print(indel_data)
    
    if not indel_data:
        raise ValueError("No indelLength data found in file.")
    
    # Create DataFrame
    df = pd.DataFrame(indel_data, columns=["indel_length", "count"])
    # Compute weighted average
    total_count = df["count"].sum()
    weighted_sum = (df["indel_length"] * df["count"]).sum()
    avg_length = weighted_sum / total_count if total_count else 0
    avg_length = float(avg_length)
    return avg_length

            
def vcf_column_count(filename):
    with gzip.open(filename, 'rt') as f:
        for line in f:
            if line.startswith('#CHROM'):
                columns = line.strip().split('\t')
                return len(columns)
    return 0  # if no #CHROM line found


############### Naming functions #########################

def create_bam_run_name(idx: str, target: AnonymousTarget) -> str:
	#pair_name=target.inputs['name_pops']
	return f'bam_neutral_{os.path.dirname(target.outputs["neutral_bam"]).replace("-","_").split("/")[-5]}_{os.path.basename(target.outputs["neutral_bam"]).replace("-","_")}'


def create_grened_run_name(idx: str, target: AnonymousTarget) -> str:
	#pair_name=target.inputs['name_pops']
	return f'pi_grened_{os.path.basename(target.outputs["parameter_value_likelihood_file"]).replace("-","_")}'

################### Templates #################

def recalculate_AF_improved(input_vcf: str, working_dir: str, species_short: str, VCF_AN: int):
	"""
	Template: Recalculate Allele frequency for each population, and output it in a pop specific temporary file. Then paste them together
	
	Template I/O::
	
		inputs = {}
		outputs = {}
	
	:param
	"""
	inputs = {'vcf_file': input_vcf}
	outputs = {'allele_freq': f'{working_dir}/tmp/allele_freq/{species_short}_AF_allpops.frq',
			'allele_count': f'{working_dir}/tmp/allele_freq/{species_short}_Acount_allpops.count',
			'allele_positions': f'{working_dir}/tmp/allele_freq/{species_short}_positions.pos'}
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
	vcf_var={VCF_AN}
	re='^[0-9]+$' 
		# regular expression saying it should start with a number between 0-9 and match one or more (+) of those, then end the string
	counter=0
	for popname in $pop_line
	do
		echo " "
		echo $popname
		if [[ -n "$vcf_var" && ( $vcf_var =~ $re && "$vcf_var" -gt 0 )]]
		then
			# if variable is set and is a number and greater than 0
			echo AN variable is set. Using manual AN instead of letting BCFtools calculate it.
			bcftools view -Ou -s $popname {inputs['vcf_file']} | bcftools query -f '%CHROM\t%POS0\t%POS\t%TYPE\t%AC\t%AN\n' | awk -F"\t" -v VCF_AN=$vcf_var '{{AC=$5; AN=VCF_AN; if(AN==0) print "NA"; else if (AN > 0) print AC/AN}}' > {working_dir}/tmp/allele_freq/{species_short}_AF_$popname.tmp
			bcftools view -Ou -s $popname {inputs['vcf_file']} | bcftools query -f '%CHROM\t%POS0\t%POS\t%TYPE\t%AC\t%AN\n' | awk -F"\t" -v VCF_AN=$vcf_var 'BEGIN {{OFS = ","}} {{AC=$5; AN=VCF_AN; if(AN==0) print 0,0 ; else if (AN > 0) print AN-AC, AC}}' > {working_dir}/tmp/allele_freq/{species_short}_Acount_$popname.tmp
			if [[ couter -eq 0 ]]
			then
				# print out positions file
				bcftools view -Ou -s $popname {inputs['vcf_file']} | bcftools query -f '%CHROM\t%POS\n' > {outputs['allele_positions'].replace(".gz", "")}
			fi
		else
			echo AN variable not set. Letting BCFtools calcuate it.
			bcftools view -Ou -s $popname {inputs['vcf_file']} | bcftools query -f '%CHROM\t%POS0\t%POS\t%TYPE\t%AC\t%AN\n' | awk -F"\t" '{{AC=$5; AN=$6; if(AN==0) print "NA"; else if (AN > 0) print AC/AN}}' > {working_dir}/tmp/allele_freq/{species_short}_AF_$popname.tmp
			bcftools view -Ou -s $popname {inputs['vcf_file']} | bcftools query -f '%CHROM\t%POS0\t%POS\t%TYPE\t%AC\t%AN\n' | awk -F"\t" 'BEGIN {{OFS=","}} {{AC=$5; AN=$6; if(AN==0) print 0,0; else if (AN > 0) print AN-AC, AC}}' > {working_dir}/tmp/allele_freq/{species_short}_Acount_$popname.tmp
			if [[ couter -eq 0 ]]
			then
				# print out positions file
				bcftools view -Ou -s $popname {inputs['vcf_file']} | bcftools query -f '%CHROM\t%POS\n' > {outputs['allele_positions'].replace(".gz", "")}
			fi
		fi
	done

	################################################
	### 		APPEND ALLELE FREQ TOGETHER  	 ###
	################################################

	# Append all population AF together

		# add one file to coming outfile
	# make header file
	# echo -e "CHROM\tPOS0\tPOS\tTYPE\tAC\tAN" > {working_dir}/tmp/allele_freq/{species_short}_AF_header.tmp1

 	#paste 
	paste -d'\t' {working_dir}/tmp/allele_freq/{species_short}_AF_*.tmp > {working_dir}/tmp/allele_freq/{species_short}_AF_allpops.tmp1
	paste -d'\t' {working_dir}/tmp/allele_freq/{species_short}_Acount_*.tmp > {working_dir}/tmp/allele_freq/{species_short}_Acount_allpops.tmp1
	
	# get header (pop names)
	ls {working_dir}/tmp/allele_freq/{species_short}_AF_*.tmp| rev |cut -d _ -f 1|cut -d . -f 2|rev|tr "\n" "\t" | sed 's/$/\\n/' > {working_dir}/tmp/allele_freq/{species_short}_header.tmp1
	

	# cat them together
	cat {working_dir}/tmp/allele_freq/{species_short}_header.tmp1 {working_dir}/tmp/allele_freq/{species_short}_AF_allpops.tmp1 > {outputs['allele_freq']}
	cat {working_dir}/tmp/allele_freq/{species_short}_header.tmp1 {working_dir}/tmp/allele_freq/{species_short}_Acount_allpops.tmp1 > {outputs['allele_count']}


	# Clean up
	rm -f {working_dir}/tmp/allele_freq/{species_short}_AF_*.tmp
	rm -f {working_dir}/tmp/allele_freq/{species_short}_AF_*.tmp1
	rm -f {working_dir}/tmp/allele_freq/{species_short}_AF_*.tmp2
	rm -f {working_dir}/tmp/allele_freq/{species_short}_Acount_*.tmp
	rm -f {working_dir}/tmp/allele_freq/{species_short}_Acount_*.tmp1
	#gzip -q {outputs['allele_positions'].replace(".gz", "")}

	echo "END: $(date)"
	echo "$(jobinfo "$SLURM_JOBID")"
	"""
	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)



def fst_calc_from_AF_improved(allele_freq_file: str, working_dir: str, species_short: str, output_directory: str, landcover_type: str, allele_count_file: str):
	"""
	Template: Recalculate Allele frequency for each population, and output it in a pop specific temporary file. Then paste them together
	
	Template I/O::
	
		inputs = {}
		outputs = {}
	
	:param
	"""
	inputs = {'allele_freq_file': allele_freq_file,
           'allele_count_file': allele_count_file}
	outputs = {'fst_file': f'{output_directory}/{species_short}_{landcover_type}_fst_allpairs_popoolation.fst',
				'fst_file_mean': f'{output_directory}/{species_short}_{landcover_type}_fst_allpairs_popoolation_reg_mean.fst',
                'fst_file_PImean_fst_mean': f'{output_directory}/{species_short}_{landcover_type}_fst_allpairs_popoolation_mean.fst'}
	options = {
		'cores': 1,
		'memory': '5g',
		'walltime': '23:00:00'
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
	mkdir -p `dirname {outputs['fst_file']}`

	#################
	### Calc Fst  ###
	#################

	#NUMBER_MAINCOLS=7
	NUMBER_MAINCOLS=1	# new AF format
	echo input al file: {inputs['allele_freq_file']}
	tot_cols=`awk 'NR == 1 {{print NF; exit}}' {inputs['allele_freq_file']}`
	echo $tot_cols
	for columns in `seq $((NUMBER_MAINCOLS)) $tot_cols`
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
				# fst_output_file_mod={working_dir}/tmp/fst/{species_short}_fst_allpairs.fst
				fst_output_file_mod={working_dir}/tmp/fst/{os.path.basename(outputs['fst_file'])}				
				fst_output_file_mod=${{fst_output_file_mod/.fst/"_pops_"$((columns-$NUMBER_MAINCOLS))"_"$((columns2-$NUMBER_MAINCOLS))".fst"}}
				rm -f ${{fst_output_file_mod/.fst/.temp1}}
				###################
				### Calculate fst
				###################
                 	# get pop AN
				AN_bcf_firstpop=`awk -v columns=$columns 'NR==2 {{
						split($columns, a, ",")
						sum = a[1] + a[2]
						print sum
					}}' {inputs['allele_count_file']}`
				AN_bcf_secondpop=`awk -v columns2=$columns2 'NR==2 {{
						split($columns2, a, ",")
						sum = a[1] + a[2]
						print sum
					}}' {inputs['allele_count_file']}`
				echo "AN (Chr. number) is calculated as (pop1 and pop2): $AN_bcf_firstpop and $AN_bcf_secondpop"
								
				awk -v firstpop=$columns -v secondpop=$columns2 -v AN_bcf_firstpop=$AN_bcf_firstpop -v AN_bcf_secondpop=$AN_bcf_secondpop '
					NR==1{{print $firstpop ":" $secondpop }};
					NR>1{{
					if ($firstpop == "" || $secondpop == "" || $firstpop == "NA" || $secondpop == "NA")  # if AF is na (aka empty field), put fst=NA
						{{ 	fst="NA"; 
							print fst, fst, fst;
							first_pi=0; second_pi=0; pi_within=0; pi_total=0; fst=NA; next	}}

					#first_pi=1-($firstpop)^2-(1-$firstpop)^2;
					#second_pi=1-($secondpop)^2-(1-$secondpop)^2;
                    
                    first_pi=2*($firstpop)*(1-$firstpop)*(AN_bcf_firstpop/(AN_bcf_firstpop-1));
					second_pi=2*($secondpop)*(1-$secondpop)*(AN_bcf_secondpop/(AN_bcf_secondpop-1));
					
                    # pi=2*(AF)*(1-AF)*(AN_bcf/(AN_bcf-1))
                         
					pi_within=(first_pi + second_pi)/2;
					AN_avg=(AN_bcf_firstpop+AN_bcf_secondpop)/2
					pi_total=(1-(($firstpop + $secondpop)/2)^2-(((1-$firstpop) + (1-$secondpop))/2)^2)*(AN_avg/(AN_avg-1));
					
					if (pi_total=="0")  # pi tot will be 0, thus zero devision
						{{ 	fst="NA"; 
                            print fst, fst, fst; 
                                #which means pi within and total will also be printed empty
							first_pi=0; second_pi=0; pi_within=0; pi_total=0; fst=NA	}}
					else 
						{{ 	fst=(pi_total-pi_within)/pi_total;
							#print $firstpop, $secondpop, first_pi, second_pi, pi_within, pi_total, fst;
							#print fst;
                            print fst, pi_within, pi_total;
							first_pi=0; second_pi=0; pi_within=0; pi_total=0; fst=NA  }}
					}}' {inputs['allele_freq_file']} >> ${{fst_output_file_mod/.fst/.temp1}}

                         
                         # create fst calc based on mean total and witin pi calcs
                    awk 'NR == 1{{print $1}};
                        NR > 1{{print $2, $3}}' ${{fst_output_file_mod/.fst/.temp1}} > ${{fst_output_file_mod/.fst/.tempPI}}
						
                        # output regular fst
                    awk '{{print $1}}' ${{fst_output_file_mod/.fst/.temp1}} > ${{fst_output_file_mod/.fst/.temp}}

                    
                    
                    	# calculate mean total pi and mean wintin pi
                    #BEGIN{{ FS = OFS = "\t"}};
                	awk ' NR == 1 {{print $1}}
						NR > 1 {{
							for (i=1; i<=NF; i++) {{
								if ($i != "" && $i != "NA" ) {{
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
						}}' ${{fst_output_file_mod/.fst/.tempPI}} > ${{fst_output_file_mod/.fst/.tempPIavg}}
						
                              
                            # calculate mean fst based on mean within pi and mean total pi
                        awk 'NR == 1 {{ print $1 }};
                            NR > 1 {{pi_within=$1; pi_total=$2;
								if (pi_total=="0")  # pi tot will be 0, thus zero devision
									{{ 	fst=NA; 
										#print $firstpop, $secondpop, first_pi, second_pi, pi_within, pi_total, fst;
										print fst;
										first_pi=0; second_pi=0; pi_within=0; pi_total=0; fst=NA	}}
								else 
									{{ 	fst=(pi_total-pi_within)/pi_total;
										#print $firstpop, $secondpop, first_pi, second_pi, pi_within, pi_total, fst;
										print fst;
										#print fst, pi_within, pi_total;
										first_pi=0; second_pi=0; pi_within=0; pi_total=0; fst=NA  }}
								}}' ${{fst_output_file_mod/.fst/.tempPIavg}} > ${{fst_output_file_mod/.fst/.tempPIavg1}}
                         
                               
			done
		fi
	done

	paste -d'\t' {working_dir}/tmp/fst/{os.path.basename(outputs['fst_file']).replace(".fst", "*.temp")} > {outputs['fst_file']}
	paste -d'\t' {working_dir}/tmp/fst/{os.path.basename(outputs['fst_file']).replace(".fst", "*.tempPIavg1")} > {outputs['fst_file_PImean_fst_mean']}
     
	###############
	## Calculate mean
	##############
	awk 'BEGIN{{ FS = OFS = "\t"}};
		NR == 1 {{print $0}}
		NR > 2 {{
			for (i=1; i<=NF; i++) {{
				if ($i != "" && $i != "NA" ) {{
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
	
	# Clean up
	#rm -f {working_dir}/tmp/fst/{species_short}_fst_allpairs*.temp
	rm -f {working_dir}/tmp/fst/{os.path.basename(outputs['fst_file']).replace(".fst", "*.temp")}
    rm -f {working_dir}/tmp/fst/{os.path.basename(outputs['fst_file']).replace(".fst", "*.temp1")}
    rm -f {working_dir}/tmp/fst/{os.path.basename(outputs['fst_file']).replace(".fst", "*.tempPI")}
    rm -f {working_dir}/tmp/fst/{os.path.basename(outputs['fst_file']).replace(".fst", "*.tempPIavg1")}

	echo "END: $(date)"
	echo "$(jobinfo "$SLURM_JOBID")"
	"""
	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)







def fst_HIVERT_calc_from_AF_improved(allele_count_file: str, positions_file: str, working_dir: str, species_short: str, output_directory: str, landcover_type: str):
	"""
	Template: Hivert, using Poolfstats package, need allele count data. not frequencies.
	An inputfile with allele counts like this is needed:
		pop1 pop2 pop3 pop4
		5,1 1,1 4,0 0,4
		3,3 0,2 2,2 0,4
		1,5 0,2 2,2 1,3
	As well as a positions file like this (without header):
		chr_1 pos
	
	Fst is calculated like this:
		Fst = (Q1 - Q2) / (1 - Q2)
		Q1 = IIS (Identiti in state) probability for genes (snps?) sampled within pops/pools.
		Q2 = IIS (Identiti in state) probability for genes (snps?) sampled between pops/pools.
		
		The script from poolfstats implement two Fst measures:
		1: a decomposition of the total variance of allele or read count frequencies in an analysis-of-variance
			framework (Weir and Cockerham 1984) which is the default procedure of the functions (as specified
			by the argument method=“Anova”). The implemented estimators are derived in Weir (1996) (eq. 5.2)
			(see also Akey et al. 2002) for allele count data (i.e., countdata objects, see 2.1); and in Hivert et al.
			(2018) (eq. 9) for (Pool-Seq) read count data (i.e., pooldata objects, see 2.2).
		2: unbiased estimatorŝ Q1 and̂ Q2 of the IIS probabilities Q1 and Q2 (as specified by the method=“Identity”
			argument). For allele count data (i.e., countdata objects, see 2.1) this estimator actually correspond to
			the one used by Karlsson et al. (2007). For Pool-Seq read count data (i.e., pooldata objects, see 2.2),
			equations A39 and A43 in Hivert et al. (2018) Supplementary Materials describe the estimators for̂
			Q1 of thê Q2 respectively. By default, when using method=Identity, the overall̂ Q1 and pairwisê Q2 are
			computed as simple averages of all population-specifiĉ Q1 and pairwise population̂ Q2, respectively. For
			completion, when setting weightpid=TRUE, an alternative weighting scheme is performed, as described
			in eqs. A46 and A47 of Hivert et al. (2018) for PoolSeq data, and Rousset (2007) for allele count
			data.
		Note that multi-locus estimates (i.e., genome-wide estimates or sliding windows estimates) are derived as the
			sum of locus-specific numerators over the sum of locus-specific denominators of the different quantities (see,
			e.g., Hivert et al. 2018 and references therein).
		
	Template I/O::
	
		inputs = {}
		outputs = {}
	
	:param
	"""
	inputs = {'allele_count_file': allele_count_file,
		   		'positions_file': positions_file}
	outputs = {'fst_file': f'{output_directory}/{species_short}_{landcover_type}_FST_pairwise_poolfstat_blockJackn_WeirCocker.fst'}
	options = { 
		'cores': 1,
		'memory': '24g',
		'walltime': '11:00:00'
	}
	spec = f"""
	# Sources environment 										OBS EDIT:
	source /home/"$USER"/.bashrc
	conda activate ecogen_neutral_diversity_wf

	# info on job
	echo "START: $(date)"
	echo "JobID: $SLURM_JOBID"


	mkdir -p {working_dir}/tmp/fst/
	mkdir -p {os.path.dirname(outputs['fst_file'])}

	#################
	### Calc Fst  ###
	#################
	
	# implement R solution for Hivert Fst measure via poolfstat package
	# allele_count output from recalculate allele frequency template should be treemix format.
	# # make countdata object (via genotreemix2countdata)
	# see vignette: https://cran.r-project.org/web/packages/poolfstat/vignettes/vignette.pdf

	# make output directory:
	outdir_poolfstats={os.path.dirname(outputs['fst_file'])}/graphics
    mkdir -p $outdir_poolfstats
    outdir_poolfstats={os.path.dirname(outputs['fst_file'])}
    mkdir -p $outdir_poolfstats
    echo {inputs['allele_count_file']}
    # /faststorage/project/EcoGenetics/people/anneaa/derived_dat_scripts/neutral_diversity_pipeline/fst_pi_gwf_intermediate_steps/test/fst_pi/Collembola/conservation_agriculture/Entomobrya_nicoleti/intermediary_files/tmp/allele_freq/EntNic_Acount_allpops.count
	
    # get script
    if [[ -f ../../../../workflow_source/scripts/fst_poolfstat_from_allelecounts.r ]]; then
		script_to_run=../../../../workflow_source/scripts/fst_poolfstat_from_allelecounts.r
	elif [[ -f ../../workflow_source/scripts/fst_poolfstat_from_allelecounts.r ]]; then
		script_to_run=../../workflow_source/scripts/fst_poolfstat_from_allelecounts.r
	else
    	echo "script not found"
    fi
    
    # start R script
	Rscript $script_to_run {inputs['allele_count_file']} {inputs['positions_file']} $outdir_poolfstats {species_short} {landcover_type}
      #Rscript ../../../../workflow_source/scripts/fst_poolfstat_from_allelecounts.r {inputs['allele_count_file']} {inputs['positions_file']} $outdir_poolfstats {species_short} {landcover_type}
	# calculates fst
	#[1] "/faststorage/project/EcoGenetics/people/anneaa/derived_dat_scripts/neutral_diversity_pipeline/fst_pi_gwf_intermediate_steps/test/fst_pi/Collembola/conservation_agriculture/Entomobrya_nicoleti/intermediary_files/tmp/allele_freq/EntNic_Acount_allpops.count"
	#[1] "/faststorage/project/EcoGenetics/people/anneaa/derived_dat_scripts/neutral_diversity_pipeline/fst_pi_gwf_intermediate_steps/test/fst_pi/Collembola/conservation_agriculture/Entomobrya_nicoleti/intermediary_files/tmp/allele_freq/EntNic_positions.pos"
      
	echo {outputs['fst_file']}

	echo "END: $(date)"
	echo "$(jobinfo "$SLURM_JOBID")"
	"""
	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def fst_HUDSON_calc_from_AF_improved(allele_freq_file: str, working_dir: str, species_short: str, output_directory: str, landcover_type: str):
	"""
	Template: Calculates Hudson Fst from allele frequencies. 
	The mean Fst is calculated as Fst from mean divergence and pi measures across all sites.
	Fst_mean = (Avg(Dxy) - Avg(Pixy)) / Avg(Dxy)	 or as in the paper, the expected value:
	  		 = (E(Dxy) - E(Pixy)) / E(Dxy)

	Dxy(per site) = Px(1-Py) + Py(1-Px)
	Pixy (per site) = (Pix + Piy)/2
	Pix (per site) = 2px(1-px)
	Piy (per site) = 2py(1-py) 
	
	Implementation following: De Jong et al. 2024, bioarxive, doi: https://doi.org/10.1101/2024.09.24.614506 (https://www.biorxiv.org/content/10.1101/2024.09.24.614506v1.full#boxed-text-1)

	Template I/O::
	
		inputs = {}
		outputs = {}
	
	:param
	"""
	inputs = {'allele_freq_file': allele_freq_file}
	outputs = {'fst_file_hudson': f'{output_directory}/{species_short}_{landcover_type}_fst_allpairs_hudson.fst',
				'fst_file_mean_hudson': f'{output_directory}/{species_short}_{landcover_type}_fst_allpairs_mean_hudson.fst'}
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
	mkdir -p {os.path.dirname(outputs['fst_file_hudson'])}

	#################
	### Calc Fst  ###
	#################

	#NUMBER_MAINCOLS=7
	NUMBER_MAINCOLS=0 # New AF format
	echo input aF file: {inputs['allele_freq_file']}
	tot_cols=`awk 'NR == 1 {{print NF; exit}}' {inputs['allele_freq_file']}`
	echo total columns: $tot_cols
	head -n 1 {inputs['allele_freq_file']} 

	for columns in `seq $((NUMBER_MAINCOLS+1)) $tot_cols`
	do
		echo column $columns
		if [ $columns == $tot_cols ]
		then
			echo At last column - break
			break
		else
			for columns2 in `seq $((columns+1)) $tot_cols`
			do
				# the two first pops are now in sight
				echo versus column $columns2
				#popul_1=`awk -v idx=$columns '{{ print $idx; exit }}' {working_dir}/tmp/fst/AF_header.tmp1`
				#popul_2=`awk -v idx=$columns2 'NR=1{{ print $idx; exit }}' {working_dir}/tmp/fst/AF_header.tmp1`
				#echo $popul_1 - $popul_2

				# create output file name:
				#fst_output_file_mod={working_dir}/tmp/fst/{species_short}_fst_allpairs_hudson.fst
				fst_output_file_mod={working_dir}/tmp/fst/{os.path.basename(outputs['fst_file_hudson'])}
				fst_output_file_mod=${{fst_output_file_mod/.fst/"_pops_"$((columns-$NUMBER_MAINCOLS))"_"$((columns2-$NUMBER_MAINCOLS))".fst"}}
				#echo $fst_output_file_mod

				###################
				### Calculate fst
				###################
				echo ${{fst_output_file_mod/.fst/.temp}}
                rm -f ${{fst_output_file_mod/.fst/.temp}}
				awk -v firstpop_AF=$columns -v secondpop_AF=$columns2 '
					NR == 1{{print $firstpop_AF ":" $secondpop_AF }};
					NR > 1{{
					if ($firstpop_AF == "" || $secondpop_AF == "" || $firstpop_AF == "NA" || $secondpop_AF == "NA")  # if AF is na (aka empty field), put fst=NA
						{{ 	fst=NA; 
						print fst;
						first_pi=0; second_pi=0; pi_both=0; D_both=0; fst=NA; next	}}
					else
						# hudson calc pi per site
						{{ first_pi=2 * $firstpop_AF *(1 - $firstpop_AF);
						second_pi=2 * $secondpop_AF * (1 - $secondpop_AF);
						pi_both=(first_pi + second_pi) / 2;

						# hudson calc divergence per site
						D_both=$firstpop_AF * (1 - $secondpop_AF) + (1 - $firstpop_AF) * $secondpop_AF;
						
						# implement fst calc HUDSON
						if (D_both=="0")  # D_both will be 0, thus zero devision
							{{ 	fst=NA; 
								#print $firstpop_AF, $secondpop_AF, first_pi, second_pi, pi_both, D_both, fst;
								print fst;
								first_pi=0; second_pi=0; pi_both=0; D_both=0; fst=NA	}}
						else 
							{{ 	# sum D_both and Pi_both as well as implement a counter (these will be written at the bottom of file)
								sum_D_both += D_both;
								sum_Pi_both += pi_both;
								counter_nonzero ++;
								
								#print $firstpop_AF, $secondpop_AF, first_pi, second_pi, pi_both, D_both, fst;
								fst=(D_both - pi_both) / D_both;
								print fst;
								sum_fst += fst;
								first_pi=0; second_pi=0; pi_both=0; D_both=0; fst=NA  }}
						}};
					}}
					# After last line: calculate and Print mean pi, mean D and mean fst
					END {{ 
						###############
						## Calculate mean  and print
						##############
						D_mean = sum_D_both/counter_nonzero; 
						Pi_mean = sum_Pi_both/counter_nonzero;
						fst_mean = (D_mean - Pi_mean) / D_mean;
						fst_regular_mean = sum_fst/counter_nonzero;
						print D_mean;
						print Pi_mean; 
						print counter_nonzero;
						print fst_regular_mean;
						print fst_mean }}

					' {inputs['allele_freq_file']} >> ${{fst_output_file_mod/.fst/.temp}}			
			done
		fi
	done

	# paste files together
	paste -d'\t' {working_dir}/tmp/fst/{os.path.basename(outputs['fst_file_hudson']).replace(".fst", "*.temp")} > {outputs['fst_file_hudson']}.tmp
	
	# transfer mean calculations to new file    
    head -n 1 {outputs['fst_file_hudson']}.tmp > {outputs['fst_file_mean_hudson']}.tmp
	tail -n 5 {outputs['fst_file_hudson']}.tmp >> {outputs['fst_file_mean_hudson']}.tmp
	echo -e "name\nD_mean\npi_mean\nsitecount\nfst_regular_mean\nfst_mean" > {outputs['fst_file_mean_hudson']}.tmp1 
	paste -d'\t' {outputs['fst_file_mean_hudson']}.tmp1 {outputs['fst_file_mean_hudson']}.tmp > {outputs['fst_file_mean_hudson']}
	
	# clean up
	rm -f {outputs['fst_file_mean_hudson']}.tmp1
	rm -f {outputs['fst_file_mean_hudson']}.tmp

	# delete last five lines from per site file
	for num in `seq 1 5`; do 
		echo $num
		sed -i '$d' {outputs['fst_file_hudson']}.tmp
	done

	mv {outputs['fst_file_hudson']}.tmp {outputs['fst_file_hudson']}

	# Clean up
	rm -f {working_dir}/tmp/fst/{species_short}_fst_allpairs*.temp


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
	
	#echo {outputs['neutral_bed'].replace("neutral.","neutralTEMP1.")}
	cat {inputs['genes_bed']} {inputs['repeats_bed']} > {outputs['neutral_bed'].replace("neutral.","neutralTEMP1.")}
	##cat {inputs['repeats_bed']} >> {outputs['neutral_bed'].replace("neutral.","neutralTEMP.")}
	
	bedtools sort -i {outputs['neutral_bed'].replace("neutral.","neutralTEMP1.")} > {outputs['neutral_bed'].replace("neutral.","neutralTEMP.")}

	bedtools merge -i {outputs['neutral_bed'].replace("neutral.","neutralTEMP.")} > {outputs['neutral_bed'].replace("neutral.","neutralTEMP1.")}
	rm -f {outputs['neutral_bed'].replace("neutral.","neutralTEMP.")}
	
	bedtools sort -i {outputs['neutral_bed'].replace("neutral.","neutralTEMP1.")} \
		> {outputs['neutral_bed'].replace("neutral.","neutralTEMP.").replace(".bed", "_sort.bed")}
		
	# Now get the reverse complement, ie. only neutral and non repetitive parts.
	bedtools complement \
		-i {outputs['neutral_bed'].replace("neutral.","neutralTEMP.").replace(".bed", "_sort.bed")} \
		-g {working_dir}/tmp/genome_stuff/{os.path.basename(inputs['reference_genome_file']).replace(".fna", ".fna.fai")} \
		> {outputs['neutral_bed'].replace("neutral.","neutralTEMP.")}
			# when done test running, this should end up in same folder as genes_bed_file

		# only keep neutral file
	rm -f {outputs['neutral_bed'].replace("neutral.","neutralTEMP.").replace(".bed", "_sort.bed")}
	rm -f {outputs['neutral_bed'].replace("neutral.","neutralTEMP1.")}
	mv {outputs['neutral_bed'].replace("neutral.","neutralTEMP.")} {outputs['neutral_bed']}
	
	
	echo "END: $(date)"
	echo "$(jobinfo "$SLURM_JOBID")"
	"""
	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)






def make_neutral_bam_improved(bam_file: str, working_directory: str, neutral_bed: str):
	"""
	Template: Create netural bam for calculating pi using grenedalf.
	
	Template I/O::
	
		inputs = {}
		outputs = {}
	
	:param
	"""
	inputs = {'bam_file': bam_file,
		   	  'neutral_bed': neutral_bed}
	outputs = {'neutral_bam': f'{working_directory}/neutral_bam/{os.path.basename(bam_file).replace(".bam", "_neutral.bam")}',
				'neutral_bam_index': f'{working_directory}/neutral_bam/{os.path.basename(bam_file).replace(".bam.gz", "_neutral.bam.bai")}'}
	# change to just adjust vcf_file, if should be placed in same folder
	options = {
		'cores': 1,
		'memory': '5g',
		'walltime': '11:00:00'
	}
	spec = f"""
	# Sources environment 										OBS EDIT:
	source /home/"$USER"/.bashrc
	conda activate grenedalf_env

	echo "START: $(date)"
	echo "JobID: $SLURM_JOBID"

	# making directories
	mkdir -p {working_directory}
	mkdir -p {os.path.dirname(outputs['neutral_bam'])}

	# sort bed to make run quicker
	sort -k1,1 -k2,2n {inputs['neutral_bed']} > {working_directory}/{os.path.basename(inputs['neutral_bed']).replace(".bed", "_sorted.bed")}
	
	# make neutral bam
	samtools view -L {working_directory}/{os.path.basename(inputs['neutral_bed']).replace(".bed", "_sorted.bed")} \
		-b {inputs['bam_file']} > {outputs['neutral_bam'].replace(".bam", "_tmp.bam")}

	samtools index {outputs['neutral_bam'].replace(".bam", "_tmp.bam")}

	#rm -f {working_directory}/{os.path.basename(inputs['neutral_bed']).replace(".bed", "_sorted.bed")}
	mv {outputs['neutral_bam'].replace(".bam", "_tmp.bam")} {outputs['neutral_bam']}
	mv {outputs['neutral_bam'].replace(".bam", "_tmp.bam.bai")} {outputs['neutral_bam_index']}

	echo "END: $(date)"
	echo "$(jobinfo "$SLURM_JOBID")"
	"""
	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)





def pi_grenedalf_template(neutral_bams: dict, working_directory: str, FILT_MIN_COUNT: int, FILT_MIN_DEPTH: int, FILT_MAX_DEPTH: int, FILT_WINDOW: int, FILT_POOL_SIZE: int, ouput_dir: str, species_short: str, landcover_type: str):
	"""
	Template: Create netural bam for calculating pi using grenedalf.
	
	Template I/O::
	
		inputs = {}
		outputs = {}
	
	:param
	"""
	inputs = {'bam_files': neutral_bams}
	outputs = {'neutral_pi_greenedalf': f'{ouput_dir}/{species_short}_{landcover_type}_grenedalf_pi.csv',
			'neutral_pi_greenedalf_mean': f'{ouput_dir}/{species_short}_{landcover_type}_grenedalf_theta_pi_mean.tsv',
			'neutral_tajD_greenedalf_mean': f'{ouput_dir}/{species_short}_{landcover_type}_grenedalf_tajimas_d_mean.tsv',
			'neutral_watTet_greenedalf_mean': f'{ouput_dir}/{species_short}_{landcover_type}_grenedalf_theta_watterson_mean.tsv'}
	# change to just adjust vcf_file, if should be placed in same folder
	options = {
		'cores': 1,
		'memory': '1g',
		'walltime': '23:59:00'
	}
	spec = f"""
	# Sources environment 										OBS EDIT:
	source /home/"$USER"/.bashrc
	conda activate grenedalf_env

	echo "START: $(date)"
	echo "JobID: $SLURM_JOBID"

	# making directories
	mkdir -p {working_directory}
	mkdir -p {os.path.dirname(outputs['neutral_pi_greenedalf'])}


	grenedalf diversity --allow-file-overwriting --filter-sample-min-count {FILT_MIN_COUNT} --filter-sample-min-read-depth {FILT_MIN_DEPTH} --filter-sample-max-read-depth {FILT_MAX_DEPTH} --filter-total-only-biallelic-snps --window-type interval --window-interval-width {FILT_WINDOW} --pool-sizes {FILT_POOL_SIZE} --window-average-policy window-length --sam-path {get_dict_key_elements(inputs['bam_files'], "neutral_bams")} --out-dir {ouput_dir} --separator-char tab --file-prefix {species_short}_{landcover_type}_grenedalf_pi_tmp_
	## sbatch --account Ecogenetics --mem 16 --wrap="grenedalf diversity --filter-sample-min-count 2 --filter-sample-min-read-depth 200 --filter-sample-max-read-depth 660 --filter-total-only-biallelic-snps --window-type interval --window-interval-width 10000 --pool-sizes 100 --window-average-policy window-length --sam-path $bam_file --out-dir /home/anneaa/EcoGenetics/people/anneaa/tests/pi_bam_grenedalf/pi --separator-char tab --file-prefix EntNic_aaRJ_biallelic --allow-file-overwriting"
	## sbatch --account Ecogenetics --mem 16 --wrap="grenedalf diversity --filter-sample-min-count 2 --filter-sample-min-read-depth 300 --filter-sample-max-read-depth 600 --filter-total-only-biallelic-snps --window-type interval --window-interval-width 10000 --pool-sizes 100 --window-average-policy window-length --sam-path $bam_file --out-dir /home/anneaa/EcoGenetics/people/anneaa/tests/pi_bam_grenedalf/pi --separator-char tab --file-prefix EntNic_aaRJ_biallelic_300_600 --allow-file-overwriting"
	# callable within depth 
	# only SNPS - only BIALLELIC

	# citations
	grenedalf citation Czech2023-grenedalf Kofler2011-PoPoolation

	mv {outputs['neutral_pi_greenedalf'].replace("_pi.", "_pi_tmp_diversity.")} {outputs['neutral_pi_greenedalf']}
	# {species_short}_{landcover_type}_grenedalf_pi_tmp_
	
        # get script
    if [[ -f ../../../../workflow_source/scripts/calc_mean_pi.py ]]; then
		script_to_run=../../../../workflow_source/scripts/calc_mean_pi.py
	elif [[ -f ../../workflow_source/scripts/calc_mean_pi.py ]]; then
		script_to_run=../../workflow_source/scripts/calc_mean_pi.py
	else
    	echo "script not found"
    fi
    
    # calculate mean for each pop:
    if [[ -f {outputs['neutral_pi_greenedalf']} ]]; then
    
		$script_to_run {outputs['neutral_pi_greenedalf']} {species_short} {landcover_type} {ouput_dir}
            #../../../../workflow_source/scripts/calc_mean_pi.py {outputs['neutral_pi_greenedalf']} {species_short} {landcover_type} {ouput_dir}
        #../../../../calc_mean_pi.py {outputs['neutral_pi_greenedalf']} {species_short} {landcover_type} {ouput_dir}
            # calculates mean of all estimates -  pi, TajD and Wattersons theta.
    fi

	mv {outputs['neutral_pi_greenedalf_mean'].replace("mean.", "summary.")} {outputs['neutral_pi_greenedalf_mean']}
	mv {outputs['neutral_tajD_greenedalf_mean'].replace("mean.", "summary.")} {outputs['neutral_tajD_greenedalf_mean']}
	mv {outputs['neutral_watTet_greenedalf_mean'].replace("mean.", "summary.")} {outputs['neutral_watTet_greenedalf_mean']}
	

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
      #change to soft linking
	ln -s {inputs['vcf_file']}	{os.path.dirname(outputs['neutral_vcf'])}/{os.path.basename(inputs['vcf_file'])}
    #cp {inputs['vcf_file']}	{os.path.dirname(outputs['neutral_vcf'])}
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
	
	#echo May I remove: {os.path.dirname(outputs['neutral_vcf'])}/{os.path.basename(inputs['vcf_file'])}
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


def calculate_pi_template_improved(allele_freq_file: str, working_directory: str, count_file: str, species_short: str, output_directory: str, landcover_type: str, allele_count_file: str, indels_file:str):
	"""
	Template: Calculate pi from bed-style files with allele frequencies.
	They look like this: Scaff Start_0-type_position End_0-type_position Variant_type AlleleFrequency 
	positions_type: string to add to outputs filename. Is it all / common / another subset of positions?  eg. 'all' 'common'

	Outputs one file with sorted variant positions for all populations :
	pi per variant position: scaff pos pi
      
      OBS. IF AN (number of chromosomes) changes from 100 per populations. changes should be made.


	Template I/O::
	
		inputs = {}
		outputs = {}
	
	:param
	"""
	inputs = {'al_freq_file': allele_freq_file, 
				'count_file': count_file,
            	'allele_count_file': allele_count_file,
                'indels_file': indels_file }
	outputs = { 'pi': f'{output_directory}/{species_short}_{landcover_type}_pi.pi',
				'pi_mean': f'{output_directory}/{species_short}_{landcover_type}_pi_mean.pi'}
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
	mkdir -p {os.path.dirname(outputs['pi'])}	

   	## calculate number to divide by
		# read in count of sites
	echo get total count of sites from count file
    line_is=`grep "total" {inputs['count_file']} | grep "all" | grep "whole_genome" | grep "intergenic_excl_repeats"`
    echo Total count: $line_is
	total_count=`echo $line_is| rev | awk '{{print $1}}' | rev`
	echo Total sites: $total_count
      
    echo get count of sites within thresholds from count file
	line_is=`grep "within_threshold" {inputs['count_file']} | grep "all" | grep "whole_genome" | grep "intergenic_excl_repeats"`
	      ######## line_is=`grep "filter" {inputs['count_file']} | grep "all" | grep "whole_genome"
	echo within thresholds: $line_is
	within_threshold=`echo $line_is| rev | awk '{{print $1}}' | rev`
	echo sites within thresholds: $within_threshold
    
    	# calculate proportion of sites within threhshold
    #prop_within=within_threshold / total_count
	prop_within=`echo "$within_threshold / $total_count" | bc -l`
	
    
    	# read in count of indels:
    line_is=`grep "number of indels:" {inputs['indels_file']}|grep "SN"`
    nIndels=`echo $line_is| rev | awk '{{print $1}}' | rev`
        # calculate avg indel length using python function
	Avg_Ind_length={parse_indel_distribution(inputs['indels_file'])}
    
		## Calculate number of sites to devide by:
            # this is the logic:
			# ($total_count * $prop_within) - ($nIndels * $prop_within) *(10 + $Avg_Ind_length)
            # (sites within threshold) minus (amount of indels within tresholds) * (10 bases removed around indels + average length of indel, which is also being removed)
    sites_corrected_count=`echo "($total_count * $prop_within) - ($nIndels * $prop_within) *(10 + $Avg_Ind_length)" | bc -l`
    ind_removed=`echo "(($nIndels * $prop_within) *(10 + $Avg_Ind_length))" | bc -l`
    echo Number of sites removed around and including indels within thresholds: $ind_removed
    echo corrected sites count, after removal of sites around indels: $sites_corrected_count
    
          
	# initiate pi files
	echo initiate pi files
	rm -f {working_directory}/tmp/pi/{species_short}_pi_calc_allSites_*
	rm -f {working_directory}/tmp/pi/{species_short}_pi_calc_allSites_01aa.tmp

    # get variables for output information  
    short_landtype={classify_land_use(landcover_type)}_
    
    
	# built for loop, each iteration making a pop specific pi calc.
        
	count_fields=`awk 'NR == 1 {{print NF; exit}}' {inputs['al_freq_file']}`
	echo number of fields $count_fields
	echo run for loop for pi
	for number in `seq 1 $count_fields`
	do
		echo $number
		popul=`awk -v number=$number 'NR == 1 {{print $number}}' {inputs['al_freq_file']}`
		echo $popul
		#echo $popul > {working_directory}/tmp/pi/{species_short}_pi_calc_allSites_$popul.tmp

        # get pop AN
        AN_bcf=`awk -v number=$number 'NR==2 {{
				split($number, a, ",")
				sum = a[1] + a[2]
				print sum
			}}' {inputs['allele_count_file']}`
        echo "AN (Chr. number) is calculated as: $AN_bcf"
        
		# add pi calculated at every line
		awk -v short_landtype=$short_landtype -v number=$number -v AN_bcf=$AN_bcf '
			NR == 1 {{print short_landtype $number}}
			NR > 1 {{ AF=$number;
				#pi=1-(AF)^2-(1-AF)^2; #pi from popoolation script
				#pi=2*(AF)*(1-AF); #similar, but nei and Tajima
				#pi=2*(AF)*(1-AF)*(100/(100-1)); # with number of chromosomes added N/N-1
                pi=2*(AF)*(1-AF)*(AN_bcf/(AN_bcf-1)); # with number of chromosomes added N/N-1 via bcftools calc, if deviating from original 100 per sample.
				print pi
			}}' {inputs['al_freq_file']} >> {working_directory}/tmp/pi/{species_short}_pi_calc_allSites_$popul.tmp
		
		# calc mean pi
		awk -v sites_corrected_count=$sites_corrected_count '
				NR == 1{{print $1}} 
				NR > 1{{ sum += $1; n++ }} 
				END {{ if (sites_corrected_count > 0) print sum / sites_corrected_count }}
				' {working_directory}/tmp/pi/{species_short}_pi_calc_allSites_$popul.tmp > {working_directory}/tmp/pi/{species_short}_pi_calc_mean_$popul.tmp 
	done

	paste -d'\\t' {working_directory}/tmp/pi/{species_short}_pi_calc_allSites_*.tmp > {outputs['pi']}
	paste -d'\\t' {working_directory}/tmp/pi/{species_short}_pi_calc_mean*.tmp > {outputs['pi_mean']}

	# clean up 
	#rm -f {working_directory}/tmp/pi/{species_short}_pi_calc_allSites_*.tmp
	#rm -f {working_directory}/tmp/pi/{species_short}_pi_calc_mean*.tmp

	echo "END: $(date)"
	echo "$(jobinfo "$SLURM_JOBID")"
	"""
	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)



def add_pi_to_collection_file(mean_file: str, collection_output_file: str, working_directory: str, species_short: str, landcover_type: str, taxonomy: str, python_ignore_list: list):
	"""
	Template: Collects pi estimates and adds it to common file.
		one template for each diversity measure.


	Template I/O::
	
		inputs = {}
		outputs = {}
	
	:param
	"""
	inputs = {'file_mean': mean_file}
	outputs = {'dummy_file': f'{working_directory}/{taxonomy}_{species_short}_{landcover_type}_{os.path.splitext(os.path.basename(mean_file))[0]}_collection_done.txt'}
	options = {
		'cores': 1,
		'memory': '5g',
		'walltime': '11:00:00'
	}
	spec = f"""
	# Sources environment 										OBS EDIT:
	source /home/"$USER"/.bashrc
	conda activate ecogen_neutral_diversity_wf

	echo "START: $(date)"
	echo "JobID: $SLURM_JOBID"

	#get script path:
	#script_path=`find ../../../.. -maxdepth 4 -name 'collect_pi_estimates.py' -print -quit`
    #script_path=`find ../../../.. -maxdepth 4 -pattern 'collect_estimates.py' -print -quit`
	#echo $script_path
        # get script
    if [[ -f ../../../../workflow_source/scripts/update_collection_file1.py ]]; then
		script_to_run=../../../../workflow_source/scripts/update_collection_file1.py
	elif [[ -f ../../workflow_source/scripts/update_collection_file1.py ]]; then
		script_to_run=../../workflow_source/scripts/update_collection_file1.py
	else
    	echo "script not found"
    fi
    
	echo {collection_output_file}
	$script_to_run {inputs['file_mean']} {collection_output_file} {species_short}_{classify_land_use(landcover_type)} {collection_output_file.replace(".txt", ".txt.tmp")} {classify_land_use(landcover_type)} "{python_ignore_list}"
      #../../../../workflow_source/scripts/update_collection_file1.py {inputs['file_mean']} {collection_output_file} {species_short}_{classify_land_use(landcover_type)} {collection_output_file.replace(".txt", ".txt.tmp")} {classify_land_use(landcover_type)}
    #../../../../workflow_source/scripts/update_collection_file.py {inputs['file_mean']} {collection_output_file} {species_short}_{classify_land_use(landcover_type)} {collection_output_file.replace(".txt", ".txt.tmp")} {classify_land_use(landcover_type)}
     
     #/home/anneaa/EcoGenetics/general_workflows/population_genetics/fst_pi/workflow_source/scripts/update_collection_file.py 
     # 	/home/anneaa/EcoGenetics/general_workflow_outputs/population_genetics/pi/Collembola/conservation_agriculture/Entomobrya_nicoleti/EntNic_conservation_agriculture_grenedalf_theta_watterson_summary.tsv 
     #   /home/anneaa/EcoGenetics/people/anneaa/tests/tmp/tmp1/test_collection.tsv 
     #  CA_EntNic 
     #   /home/anneaa/EcoGenetics/people/anneaa/tests/tmp/tmp1/test_collection.tsv.tmp
  
        
    # script finished, now make copleting files appear
    #if [ -f {collection_output_file.replace(".txt", ".txt.tmp")} ]; then
    if [[ -s {collection_output_file.replace(".txt", ".txt.tmp")} ]]; then
      	mv {collection_output_file.replace(".txt", ".txt.tmp")} {collection_output_file}
    fi
	if [[ -s {collection_output_file.replace(".txt", ".txt.lock")} ]]; then
      	rm -f {collection_output_file.replace(".txt", ".txt.lock")}
    fi
    touch {outputs['dummy_file']}
    echo Delete this, if rerun is needed: {outputs['dummy_file']}
    
    
	echo "END: $(date)"
	echo "$(jobinfo "$SLURM_JOBID")"
	"""
	return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


