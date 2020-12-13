#EYA4 eQTL colocalization (Bash)
# Step 1: extract all eQTL from all_association tissues
## names 
ls *.gz | cut -d_ -f 11,12 > names
sed -i 's/\.allpairs.txt.gz//g' names
## files_list 
ls *.gz
## Actual extraction
#!/bin/bash

i=$SLURM_ARRAY_TASK_ID
line=$(sed -n "${i}{p;}" < ./files_list)
name=$(sed -n "${i}{p;}" < ./names)

zgrep 'ENSG00000112319.17' ./${line} > ENSG00000112319.17_${name}.eqtl

# Step 2: Create a look_up table for GTEx ids and rsids, downloaded from GTEx website

awk '$2 == "chr6"' GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.lookup_table.txt > chr6_lookup.txt

# Step 3: Create a snp file with chrmosome 6 information only 

awk '$3 == 6' Biobank2-British-Bmd-As-C-Gwas-SumStats.txt > chr6_bmd_snps.txt

##The rest is done in R
