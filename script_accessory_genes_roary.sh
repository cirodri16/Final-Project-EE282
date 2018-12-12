#!/bin/bash

# Setting the working directory:
WORKDIR=/Users/cynthiarodriguez/Desktop/Final_project_informatics

# Take the header line and add to new file and extract only columns 1: gene name,4: number of isolates present in the cluster,5: number of sequences in the cluster,15: presence absence of genes in each sample.
head -n1 $WORKDIR/gene_presence_absence_copy.txt | cut -f1,4,5,15- > $WORKDIR/gene_presence_absence_sorted.txt

#Make single genes and all genes a varible so we can remove singletons and all core gene clusters from the presence absence file and only have the accessory genes.
singlegenes=1
allgenes=129

# Check that isolates and sequences with that gene are not 1 or all of them so we can be left with only the accessory genes and save the output to a new file.
awk -F"\t" -v var="$singlegenes" '(NR>1) && ($4 != var ) ' $WORKDIR/gene_presence_absence_copy.txt | \
awk -F"\t" -v var="$singlegenes" '(NR>1) && ($5 != var ) ' | \
awk -F"\t" -v var="$allgenes" '(NR>1) && ($4 != var ) ' | \
awk -F"\t" -v var="$allgenes" '(NR>1) && ($5 != var ) ' | \
awk -F"\t" '(NR>1) && ($4 == $5 ) ' | \
cut -f1,4,5,15- >> $WORKDIR/gene_presence_absence_sorted.txt

# Extract columns 1 and 4 from the previous output so we can visualize the accessory genes in RStudio.
cut -f1,4- $WORKDIR/gene_presence_absence_sorted.txt > $WORKDIR/gene_presence_absence_sorted.fix.txt