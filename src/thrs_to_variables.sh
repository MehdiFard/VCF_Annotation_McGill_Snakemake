#!/bin/bash
QCTHR_IN1=$1
MVCFS_IN2=$2
MVCFSQC_OUT=$3

# A while loop to assign each threshold value to its name (assigns value of each field to field name)
nr=1
while IFS=, read -r -a ary; do
    if (( nr == 1 )); then
        col_names=("${ary[@]}")
    else
        for (( i = 0; i < ${#ary[@]}; i++ )); do
            printf -v "${col_names[i]}" "${ary[i]}"
        done
    fi
    (( nr++ ))
done < $QCTHR_IN1
echo "# Quality control ..."
echo "Thresholds:"
echo "Minor Allel Frequency: $MAF"
echo "Hardy-Weinberg test: $HWE"
echo "Missing genotype per site: $Missing_per_Mrks"
echo "Missing genotype per individual: $Missing_per_Inds"
vcftools --gzvcf $MVCFS_IN2 --missing-indv --out "intermediate/missDataPerInds"
awk '$5 > $Missing_per_Inds {print $1}' "intermediate/missDataPerInds.imiss" | sed '1d' > "intermediate/HighMissInds"
vcftools --gzvcf $MVCFS_IN2 --remove "intermediate/HighMissInds" --max-missing $Missing_per_Mrks --maf $MAF --hwe $HWE --recode --recode-INFO-all --out Merged_QCOK.vcf --stdout | bgzip -c > $MVCFSQC_OUT
