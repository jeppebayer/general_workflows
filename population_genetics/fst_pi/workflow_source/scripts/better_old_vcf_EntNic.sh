#!/usr/bin/env bash

# conda activate
## inputs:
vcf_file=$1
# vcf_file=/home/anneaa/EcoGenetics/general_workflow_outputs/museomics/vcf/Aglais_urticae/2022/sorted_Aglais_urticae_2022.vcf
# vcf_file=/home/anneaa/EcoGenetics/general_workflow_outputs/museomics/vcf/Gonepteryx_rhamni/GonRha.freebayes_n3_p100_minaltfrc0_minaltcnt2_singlecall.merged.norm.bcftoolsfilter_AF0_SnpGap5_typesnps_biallelic_DP200-dynamic_AO1.vcf.gz
# vcf_file=/home/anneaa/EcoGenetics/people/Sarah/data/EntNic_biallelic.vcf
working_directory=$2
# working_directory=/home/anneaa/EcoGenetics/people/anneaa/tests


########################################################################
### Replace dots in regex found fields with 0's  ###
########################################################################
new_file=`basename $vcf_file`
awk '
BEGIN { FS = OFS = "\t" }
{
    for (i = 1; i <= NF; i++) {
        # Match fields with repetitive "./" patterns followed by repetitive ":.:.:." patterns
        if ($i ~ /^(\.\/)*+(:\.)*/) {
            # Replace "./" with "0/"
            gsub(/\.\//, "0/", $i);
            # Replace "/." with "/0"
            gsub(/\/\./, "/0", $i);
        }
    }
    print
}' $vcf_file > $working_directory/${new_file/.vcf/_dotcorrected.vcf}

### Recalculate vcf fields (hopefully regain AC/AN fields)
#Probably does that when parsing columns later.
# bcftools view recalculates info fields


