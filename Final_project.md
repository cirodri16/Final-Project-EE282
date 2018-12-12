# Final Project:
Cynthia I. Rodriguez
***
### Description:
I have a total of 129 bifidobacteria genomes downloaded from the NCBI and PATRIC databases.
#### Goals:
1. To annotate the genomes in an uniform way.
2. To get the pangenome of these strains: core genes and accessory genes.
3. To visualize the differences in accessory genes between the different strains.
***
#### - To address Goal 1 I used the program prokka to annotate the genomes:
https://github.com/tseemann/prokka 
##### Installation of prokka in the computer:
Bioconda: If you use Conda you can use the Bioconda channel:
```
conda install -c conda-forge -c bioconda prokka
```
##### For loop analysis and annotation use:
For .fna files:
```
for k in *.fna; do prokka $k --outdir "$k".prokka.output;echo $k;done
```
For .fasta files:
```
for k in *.fasta; do prokka $k --outdir "$k".prokka.output;echo $k;done
```
**Note1:** Since I had fasta files I used the second command and ran it in the folder where I had all the genomes.
**Note2:** Make sure your names inside the fasta file are shorter than 37 characters or else you will get an error message saying “Contig ID must be ,+37 characters long".

I copied all the .gff files from the prokka output into a folder and rename each file so each has a unique name with no spaces. Preferable rename them the accession number.

Now we are ready to use these annotated genomes in .gff format to run Roary.

***
#### - To address Goal 2 I used the program Roary to annotate the genomes:
https://sanger-pathogens.github.io/Roary/ 
##### Installation of roary in the computer:
```
Install conda. Then install bioconda and roary:
conda config --add channels r
conda config --add channels defaults
conda config --add channels conda-forge
conda config --add channels bioconda
conda install roary
```
##### How to use roary:
Get all the .gff files into one folder and rename each file to have the accession number on the name (no spaces)- last step from above.
Use the following commands while in the folder where the .gff files are:
```
roary --group_limit 1000000 -e --mafft -cd 99 -i 50 -p 32 *.gff
```
**Note1:** the -cd switch specifies the percentage of isolates a gene must be in to be core [99]; -i switch minimum percentage identity for blastp [50] I lowered the blastp percentage because the bifidobacteria genomes are very diverse and to ensure I got all possible orthologs; -p switch specifies the number of threads for the computer to use and run the analysis.

##### The two output files from Roary I cared about:
***summary_statistics.txt:*** Number of genes in the core and accessory. A text file with an overview of the genes and how frequently they occur in the input isolates.

***gene_presence_absence.csv:*** The gene presence and absence spreadsheet lists each gene and which samples it is present in.
From Roary page each column represents:
1. The gene name, which is the most frequently occurring gene name from the sequences in the cluster. If there is no gene name, then it is given a generic unique name group_XXX.
2. A non unique gene name, where sequences with the same gene name have ended up in different groups. It might be because of split genes, or miss annotation.
3. Functional annotation. The most frequently occurring functional annotation from the cluster is used.
4. Number of isolates represented in the cluster.
5. Number of sequences in the cluster.
6. Average number of sequences per isolate. This is normally 1. If this is greater than 1 then there is over clustering and the paralogs couldn’t be split.
7. Genome fragment, where there is evidence at the contig level that the genes are linked.
8. Order within fragment, combined with the genome fragment this gives an indication of the order of genes within the graph. In Excel, sort on Column G and H.
9. Accessory Fragment is where core genes are excluded and there is evidence at contig level that the genes are linked.
10. Accessory order with fragment, combined with the Accessory fragment this gives an indication of the order of genes within the accessory graph. In Excel, sort on columns I and J.
11. Comments on the quality of the cluster. Miss predictions are noted, as are single genes on single contigs, which can be evidence of low level contamination.
12. Minimum sequence length in nucleotides of the cluster.
13. Maximum sequence length in nucleotides of the cluster.
14. Average sequence length in nucleotides of the cluster.
15. Presence and absence of genes in each sample, with the corresponding source Gene ID.

The ***gene_presence_absence.csv*** file was converted to .txt to make downstream analysis easier.
***
#### - To address Goal 3 to visualize the differences in accessory genes between the different strains I ran the following script in the ***gene_presence_absence.txt*** file.
#### Script:
```
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
```
#### RStudio:
```
# Set working directory:
setwd('/Users/cynthiarodriguez/Desktop/Final_project_informatics')

#Read data output we obtained from the previous script"
tablein <- read.table('gene_presence_absence_sorted.fix.txt', 
                      row.names = 1,header = T, sep = '\t', quote = NULL, na.strings=c("","NA"))

# Load gplots and heatmap.plus to make a heatmap.
library(gplots)
library(heatmap.plus)

#Read the data as a dataframe.
tablenum <- as.data.frame(lapply(tablein, function(x) as.integer(x!="0")))
tablenum[is.na(tablenum)] <- 0

#Read the data as a matrix to make heatmap
rownames(tablenum) <- rownames(tablein)
input <-as.matrix(t(tablenum))

#Use the Euclidean distance for clustering of strains in the heatmap.
distfunc <- function(x) dist(x, method="euclidean")

#Make heatmap.
heatmap.2(input, main = "Accessory genes from 129 bifidobacteria strains", distfun=distfunc, 
          trace = "none", margins = c(10,13), cexRow=0.2, cexCol=0.3)
```
Save the image as a pdf.
