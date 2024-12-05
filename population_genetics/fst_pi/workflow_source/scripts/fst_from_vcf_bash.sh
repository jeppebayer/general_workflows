#!/usr/bin/env bash

## inputs:
vcf_file=$1
# vcf_file=/home/anneaa/EcoGenetics/general_workflow_outputs/museomics/vcf/Aglais_urticae/2022/sorted_Aglais_urticae_2022.vcf
# vcf_file=/home/anneaa/EcoGenetics/general_workflow_outputs/museomics/vcf/Gonepteryx_rhamni/GonRha.freebayes_n3_p100_minaltfrc0_minaltcnt2_singlecall.merged.norm.bcftoolsfilter_AF0_SnpGap5_typesnps_biallelic_DP200-dynamic_AO1.vcf.gz
# vcf_file=/home/anneaa/EcoGenetics/people/Sarah/data/EntNic_biallelic.vcf
# vcf_file=/home/anneaa/EcoGenetics/people/anneaa/tests/EntNic_biallelic_dotcorrected.vcf
working_directory=$2
# working_directory=/home/anneaa/EcoGenetics/people/anneaa/tests
fst_output_file=$3
# fst_output_file=fst_file_all_pairs
# conda activate ecogen_neutral_diversity_wf

# info on job
echo "START: $(date)"
echo "JobID: $SLURM_JOBID"

########################
### GET ALLELE FREQ  ###
########################

mkdir -p $working_directory
mkdir -p $working_directory/tmp/
#CHROM


if [[ $vcf_file =~ \.gz$ ]]; then
	pop_line=`zcat $vcf_file| grep --max-count 1 "^\#[A-Z]" `
else
    pop_line=`grep --max-count 1 "^\#[A-Z]" $vcf_file`
fi

echo $pop_line
pop_line=`echo $pop_line | awk -F' ' '{ for(i=10; i<=NF; i++) printf "%s%s", $i, (i<NF ? "\t" : "\n") }'`
echo $pop_line
#pop_line=`echo $pop_line |cut -d'-' -f2 | cut -d'/' -f1 | cut -d" " -f2`
#pop_line=`echo $pop_line | awk -F',' '{ for(i=1; i<=NF; i++) printf "%s%s", $i, (i<NF ? "\t" : "\n") }'`

# loop over each pop:
# while calculating allele frequencies for each population
for popname in $pop_line
do
	echo " "
	echo $popname
	# OBS search regular expression for col starting with /. and ending similarly, only replace dots within.
	# then recalculate using bcftools (maybe called recalcute?)
	bcftools view -Ou -s $popname $vcf_file | bcftools query -f '%CHROM\t%POS0\t%POS\t%TYPE\t%AC\t%AN\t%AF\n' | awk 'FS=OFS="\t" {AC=$5; AN=$6; if(AN==0) print $0, NA; else if (AN > 0) print $0, AC/AN}' > $working_directory/tmp/AF_$popname.tmp
done

################################################
### APPEND ALLELE FREQ TOGETHER  ###
################################################

# Append all population AF together
cat `ls $working_directory/tmp/AF_*.tmp|head -n1` > $working_directory/tmp/AF_allpops.tmp1
	# add one file to coming outfile
# make header file
echo -e "CHROM\tPOS0\tPOS\tTYPE\tAC\tAN\tAF" > $working_directory/tmp/AF_header.tmp1
for conseq_pops in `ls $working_directory/tmp/AF_*.tmp`
do 
	pop_is=`echo $conseq_pops | rev | cut -d"/" -f1|rev| cut -d"." -f 1`
	if [ $conseq_pops == `ls $working_directory/tmp/AF_*.tmp|head -n1` ]
	then
		echo "appending allele frequency (AC/AN) of pop: $conseq_pops"
		awk -v pop_is=$pop_is 'FS=OFS="\t" {print $0, pop_is}' $working_directory/tmp/AF_header.tmp1 > $working_directory/tmp/AF_header.tmp2
		continue
	else
		echo "appending allele frequency (AC/AN) of pop: $conseq_pops"
		awk '{print $8}' $conseq_pops > $conseq_pops.tmp2
			# make intermediary file with only AF

		paste -d'\t' $working_directory/tmp/AF_allpops.tmp1 $conseq_pops.tmp2 > $working_directory/tmp/AF_allpops.tmp2
		mv $working_directory/tmp/AF_allpops.tmp2 $working_directory/tmp/AF_allpops.tmp1
			# paste thw two files together (cbind)

		awk -v pop_is=$pop_is 'FS=OFS="\t" {print $0, pop_is}' $working_directory/tmp/AF_header.tmp2 > $working_directory/tmp/AF_header.tmp1
		mv $working_directory/tmp/AF_header.tmp1 $working_directory/tmp/AF_header.tmp2
			# add pop to header
	fi
done
mv $working_directory/tmp/AF_header.tmp2 $working_directory/tmp/AF_header.tmp1
cat $working_directory/tmp/AF_header.tmp1 $working_directory/tmp/AF_allpops.tmp1 > $working_directory/tmp/AF_allpops.frq





# calculate fst
NUMBER_MAINCOLS=7
tot_cols=`awk '{print NF}' $working_directory/tmp/AF_header.tmp1`
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
			#popul_1=`awk -v idx=$columns '{ print $idx; exit }' $working_directory/tmp/AF_header.tmp1`
			#popul_2=`awk -v idx=$columns2 'NR=1{ print $idx; exit }' $working_directory/tmp/AF_header.tmp1`
			echo $popul_1 - $popul_2

			# create output file name:
			fst_output_file_mod=$fst_output_file.fst
			fst_output_file_mod=${fst_output_file_mod/.fst/"_pops_"$((columns-$NUMBER_MAINCOLS))"_"$((columns2-$NUMBER_MAINCOLS))".fst"}
			###################
			### Calculate fst
			###################
			echo index: $columns-$columns2 > $working_directory/tmp/${fst_output_file_mod/.fst/.temp}
			echo pops_nr: $((columns-$NUMBER_MAINCOLS))-$((columns2-$NUMBER_MAINCOLS)) >> $working_directory/tmp/${fst_output_file_mod/.fst/.temp}
			echo pops: $popul_1-$popul_2 >> $working_directory/tmp/${fst_output_file_mod/.fst/.temp}
			
			awk -v firstpop=$columns -v secondpop=$columns2 '
				NR==1{print $firstpop "-" $secondpop };
				NR>1{
				if ($firstpop == "" || $secondpop == "" || $firstpop == "NA" || $secondpop == "NA")  # if AF is na (aka empty field), put fst=NA
					{ 	fst=NA; 
						print fst;
						first_pi=0; second_pi=0; pi_within=0; pi_total=0; fst=NA; next	}

				first_pi=1-($firstpop)^2-(1-$firstpop)^2;
				second_pi=1-($secondpop)^2-(1-$secondpop)^2;
				
				pi_within=(first_pi + second_pi)/2;

				pi_total=1-(($firstpop + $secondpop)/2)^2-(((1-$firstpop) + (1-$secondpop))/2)^2;
				
				if (pi_total=="0")  # pi tot will be 0, thus zero devision
					{ 	fst=NA; 
						#print $firstpop, $secondpop, first_pi, second_pi, pi_within, pi_total, fst;
						print fst;
						first_pi=0; second_pi=0; pi_within=0; pi_total=0; fst=NA	}
				else 
					{ 	fst=(pi_total-pi_within)/pi_total;
						#print $firstpop, $secondpop, first_pi, second_pi, pi_within, pi_total, fst;
						print fst;
						first_pi=0; second_pi=0; pi_within=0; pi_total=0; fst=NA  }
				}' $working_directory/tmp/AF_allpops.frq >> $working_directory/tmp/${fst_output_file_mod/.fst/.temp}

				mv $working_directory/tmp/${fst_output_file_mod/.fst/.temp} $working_directory/tmp/${fst_output_file_mod}			
			
		done
	fi
done

fst_output_file_new=allpairs_$fst_output_file.fst
paste -d'\t' $working_directory/tmp/${fst_output_file}*.fst	> $working_directory/$fst_output_file_new





#####################################
### Calculate mean
#####################################


# calculate mean:
awk 'BEGIN{ FS = OFS = "\t"};
	NR == 1 || NR==2 || NR==3 || NR==4 {print $0}
	NR > 2 {
		for (i=1; i<=NF; i++) {
			if ($i != "" || $i != "NA" ) {
				count[i]++
				sum[i] += $i };
		}
	}
	END {
		for (i=1; i<=NF; i++) {
			if (i == NF) {
				printf "%f", sum[i] / count[i]  # No tab for the last column
			} else {
				printf "%f\\t", sum[i] / count[i]  # Tab between columns
			}
		}
}' $working_directory/$fst_output_file_new > $working_directory/${fst_output_file_new/.fst/_mean.fst}