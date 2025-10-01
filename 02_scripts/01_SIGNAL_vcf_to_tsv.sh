#!/bin/bash


vcf=$1
output_tsv=$2


if [[ "$vcf" == *.gz ]]; then

        zcat $vcf | awk '!/^#/ { print $1 "\t" $2 "\t" $4 "\t" $5}' > $output_tsv

else
        awk '!/^#/ {print $1 "\t" $2 "\t" $4 "\t" $5}' $vcf > $output_tsv

fi
