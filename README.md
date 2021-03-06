## dino_popgen
# **Final project for Environmental Bioinformatics, Fall 2021**

**Proposal Text**

In their study, Reich et al. (2021) use single nucleotide polymorphism (SNP) data to address whether strains of symbiotic dinoflagellates Symbiodinium 'fitti' are genetically differentiated across the two Caribbean acroporid coral species and their hybrid, with which S. 'fitti' is associated. Tissue samples were collected from 76 acroporids across their geographic distribution in the Caribbean, and resulting sequences were filtered, trimmed, aligned, and mapped to a draft assembly of the S. 'fitti' genome (available on dryad, https://doi.org/10.5061/dryad.xgxd254g8, ~326 MB). A subset of these samples were removed from downstream analysis due to co-infection with multiple S. ‘fitti’ strains. High-quality SNPs between the remaining sample genomes were identified and used to characterize population structure. Raw sequence data are publicly available on NCBI under SRA project PRJNA473816 (76 entries, ~340 GB). 

Code associated with this project, including all scripts used for data cleaning, analysis, and figure creation, has been uploaded to a github repository (https://github.com/hgreich/Sfitti). We propose to use the published raw data and code to repeat the following analysis steps: 1) sequence trimming and cleaning; 2) removal of samples with evidence of co-infection; 3) genome alignment and identification of high-quality SNPs; and 4) analysis of population structure using cleaned SNP data. If time allows, we may also undertake the final analysis of the paper, which identifies loci that showed strong evidence of selection. Our overall goal is to use manuscript resources, with modifications as necessary, to recreate Figure 2 (PCA analysis of S. ‘fitti’ genomic differentiation) and Figure 3 (dendogram displaying similarity between sample genomes) from Reich et al. (2021).

**Directory Contents**
 - README.md
 - notebook.txt
   * _text file for group members to log project activity/progress_
 - data/
   * project-PRJNA473816-metadata.csv
     - _project metadata from NCBI_
   * supp_table_1.csv
     - _Supplemental table 1 from Reich et al. (2021) paper, containing sample metadata_
   * seq_access_nums.txt
     - _List of SRA accession numbers for 47 sequences re-analyzed here_
   * sra/
     - _47 raw sequences (forward and reverse) from "shallow" reads, not included on github due to file size_
     - _download using load_seqs.txt script and seq_access_nums.txt file_
   * genome/
     - _draft genome for_ S. 'fitti' _from same study, not included on github due to file size_
     - _available for download at <https://doi.org/10.5061/dryad.xgxd254g8>_
 - envs/
   * trim_clean_call_snps.yml
     - _yaml file to load conda envs for seq trimming and aligning, and SNP ID and filtering_
   * analyze_snps.yml
     - _yaml file to load conda envs for SNP analysis_
 - figures
   * Fig2_JBW.pdf
     - _results of PCA_
   * Fig3_phylogeny_EF.pdf
     - _raw figure output from code_
   * Fig3_phylogeny_EF_edit.png
     - _prettied-up figure, to better match OG figure in paper_
   * Fig_5_manhattanplot_total.pdf
     - _results of selection analysis_
   * Fig2_original_Reich2021.png
     - _original figure 2 (PCA) from Reich et al. 2021_
   * Fig3_original_Reich2021.png
     - original figure 3 (phylogeny) from Reich et al 2021_
   * Fig5_original_Reich2021.png
     - _original figure 5 (selection anlaysis) from Reich et al. 2021_
 - jupyter-notebooks/
   * dino_popgen_final_writeup.ipynb
     - _the final comparison of reproduced results and original project_
   * trimlog_ext.ipynb
     - _code to extract statistics out of trimming logfile_
   * BayeScan_outliers.ipynb
     - _R code to pull outliers from BayeScan output_
   * PCAdapt_outliers.ipynb
     - _R code to find outliers in .vcf file using PCAdapt AND compare BayeScan to PCAdapt_
   * figure-2_JBW.Rmd
     - _R code for (re-)creation of figure 2_
   * figure-3_EF.ipynb
     - _python code used to format final treefile for visualization with iTOL_
   * figure-5_CLR.Rmd
     - _R code for (re-)creation of figure 5_
 - logs/
   * _log files from HPC sbatch runs_
 - scripts/
   * load_seqs.txt
     - _script to download raw seqs from NCBI, using data/seq_access_nums.txt_
   * zip.txt
     - _script to zip files so they take up less space_
   * trim.txt
     - _script to trim sequences using cutadapt_
   * align_seqs.txt
     - _script to align trimmed seqs to draft genome_
   * samtools_seqs.txt
     - _convert .sam files to .bam, remove PCR duplicates, and get alignment stats_
   * call_filter_snps/
     - mpileup.txt
       * _create mpileup file with all samples_
     - bcftools_run.txt
       * _call variants_
     - bcftools_view.txt
       * _identify SNPs_
     - bcftools_filter.txt
       * _filter out low and medium quality SNPs_
     - vcftools_filter.tx
       * _filter out sites with >20% missing data_
   * phylogeny/
     - vcftools_filter_phylo.txt
       * _rerun filtering step to subset only sites with NO missing data_
     - vcfkit_run.txt
       * _vcf to fasta alignment_
     - raxml_ng.txt
       * _run bootstrapping and create tree_
   * bayescan/
     - BayeScan2.1/
       * _run BayeScan from this folder_
     - bayescan.txt
       * _IDs SNP outliers based on population information_
     - vcf2bayescan.pl
       * _perl code to convert vcf file into a BayeScan file using a population map_
       * _taken from <https://github.com/santiagosnchez/vcf2bayescan/blob/master/vcf2bayescan.pl>_
     - _vcf2bayescan.pl output:_ host.txt, loc.txt, hostXloc.txt
   * supp_table_prep.R
     - _R script for subsetting and re-ordering supplemental table_
   * proc_for_vcf-seq-order.txt
     - _text file explaining bash and R code used to subset and re-order metadata table to match sample order in allhqSNPmm80.recode.vcf_
   * _code for PCA analysis and figure in jupyter-notebooks/figure-2_JBW.Rmd_
- output/
   * supp_table_subset.csv
     - _Supp table 1, subset with only the 47 samples re-analyzed here_
   * supp_table_ordered.csv
     - _Supp table 1, re-ordered to match vcf file_
   * proc_for_vcf-seq-order.txt
     - _file describing bash and R code used to convert supp_table_1.csv to supp_table_ordered.csv_
   * trim-data/
     - _trimmed sequences (after cutadapt run)_
   * align-data/
     - _aligned sequences (.bam files) (after bwa run)_
   * snp_calling/
     - allhqSNPmm80.recode.vcf
       * _final data file for high-quality SNPS, allowing for missing data_
     - _additoinal SNP data files (output from bcftools), not included on github due to file size_
   * phylogeny/
     - phylogeny_msa.fasta
       * _fasta alignment_
     - allhqSNPmm80_phylo.recode.vcf
       * _final data file for high-quality SNPS, allowing for NO missing data_
     - raxml_output/
       * _output from raxml run_
       * raxml_output.raxml.support 
         - file used in creation of Figure 3 (figure-3_EF.ipynb)
   * seq_outliers/
     - _outliers identified by BayeScan and PCAdapt and outlier statistics_

## Data Processing

### Data Download
Raw data from SRA Project PRJNA473816 was downloaded from NCBI (https://www.ncbi.nlm.nih.gov/bioproject/PRJNA473816/) using the fasterq-dump function from the sra-tools package (version 2.10.0). Only samples where colonization status of Symbiodinium was listed as "single" were selected for download, this was to prevent complications in data interpretation resulting from co-colonization in downstream analysis. These samples were identified by the "Colonization Status" field in the first supplemental table (S1) in the paper. The download script timed out on a large sequence (SRR7235983.fastq, 140 GB, total reads > 1.5 billion), this accession was excluded from our analysis due to size. All fastq files were zipped in order to save space on the HPC. The S. fitti reference genome was downloaded from dryad (https://doi.org/10.5061/dryad.xgxd254g8).  
```
#following commands run scripts/load_seqs.txt:
#input data/seq_access_nums.txt

for seq in $(cat ${1}) #input is a textfile with list of accession numbers
do
        fasterq-dump --split-3 --verbose -O data/sra/ $seq #retrieve SRA data
done
```
### Cleaning and Trimming
The cutadapt package (version 2.6) was used to trim and clean reads. This process was well-specified for the methods related to reference genome assembly in the paper, but not for individual samples; methods outlined for the reference genome were followed as closely as possible. Cutadapt 1.6 was specified in the paper but was not available through Bioconda, we decided to download the latest version.
```
#following commands run in scripts/trim.txt
for seq in $(cat ${1}) #input is a textfile with list of accession numbers
do
        cutadapt -m 50 -q 25 --trim-n -o data/trim-data/${seq}_1_out.fasta.gz -p data/trim-data/${seq}_2_out.fasta.gz data/sra/${seq}_1.fastq.gz data/sra/${seq}_2.fastq.gz #trim n's, phred score threshold of 25, minimum length of 50 bp
done
```
At this point it was discovered that one of the sequences (SRR7235991) did not download properly during the initial fasterq-dump; this sequence was redownloaded and trimmed separately. A python script was used to extract some relevant statistics from the trim log. The average percent of pairs written was 98.7%, and number of pairs written ranged from ~10 million to ~39 million.  
### Alignment
The Burrows-Wheeler Aligner (BWA) package (version 0.7.15) was used to index the reference genome and align the sequences. Alignment was run with the bwa mem function on 16 threads.
```
#following commands run in scripts/align_seqs.txt
bwa index /vortexfs1/omics/env-bio/collaboration/dino_popgen/data/genome/Sfitti_Apalm_v1.fa #index reference
#align each sample file, 16 threads
for seq in $(cat ${1})
do
        bwa mem -t 16 /vortexfs1/omics/env-bio/collaboration/dino_popgen/data/genome/Sfitti_Apalm_v1.fa /vortexfs1/omics/env-bio/collaboration/dino_popgen/data/trim-data/${seq}_1_out.fasta.gz /vortexfs1/omics/env-bio/collaboration/dino_popgen/data/trim-data/${seq}_2_out.fasta.gz > data/align-data/${seq}.sam 
done
```
Samtools and bcftools (both version 1.4.1) were installed and used to further process the sequences.  
Samtools view was used to convert sam files to bam format in order to save space on the HPC. Then, the .bam files were sorted and indexed so the PCR duplicates could be removed with samtools rmdup. Now that duplicates are removed, this .bam file was used in downstream analyses. Additionlly, the .bam files were indexed again and samtools flagstat used to calculate alignment statistics.
```
# following commands run in scripts/samtools_seqs.txt
# input is a textfile with list of accession numbers
for seq in $(cat ${1}) #for each sequence in provided text file
do
  	samtools view -bT /vortexfs1/omics/env-bio/collaboration/dino_popgen/da$
        #convert .sam to .bam files to save space on HPC
        #flag -T refers to reference files (genome), and -b asks for .bam output

        samtools sort /vortexfs1/omics/env-bio/collaboration/dino_popgen/output$
        samtools index -b /vortexfs1/omics/env-bio/collaboration/dino_popgen/ou$
        #sort and index .bam files, and give .bam output (-b); name sorted file$

        samtools rmdup /vortexfs1/omics/env-bio/collaboration/dino_popgen/outpu$
        #remove PCR duplicates and name output files as directed

        samtools index -b /vortexfs1/omics/env-bio/collaboration/dino_popgen/ou$
        #index no-PCR-duplicate .bam (-b) files

        samtools flagstat /vortexfs1/omics/env-bio/collaboration/dino_popgen/ou$
        #calculate alignment statistics and write to new text file
done
```
### SNP Calling and Filtering
The samtools mpileup function was used to create a single pileup bcf file from all bam files. The resulting bcf file was run through bcftools call in order to call variants and produce a vcf output file. All flags used for these functions were documented clearly in the methods section of the paper.  
```
#following commands run in snp_id.txt
#create a bcf pileup
samtools mpileup -ugAEf /vortexfs1/omics/env-bio/collaboration/dino_popgen/data/genome/Sfitti_Apalm_v1.fa -t AD,DP /vortexfs1/omics/env-bio/collaboration/dino_popgen/data/align-data/*.bam > /vortexfs1/omics/env-bio/collaboration/dino_popgen/data/snp_calling/symb_mpileup.bcf
#call variants
bcftools call -f GQ -vmO z --ploidy 1 -o symb_mpileup.vcf.gz /vortexfs1/omics/env-bio/collaboration/dino_popgen/data/snp_calling/symb_mpileup.bcf
```
At this point, we tried to run bcftools view to identify snps but the program threw an error when loading shared libraries. We chose to update bcftools to the current version (1.9) to proceed. The output of bcftools view was run with bcftools filter function in order to isolate only the high-quality variants (quality score >200). An additional filtering step was conducted with vcftools filter, which removed all sites with >20% missing data as well as indels. While flags related to bcftools were noted in the paper, the vcftools operation was not documented in the manuscript. The github repository included information on the vcftools operations, however the version noted in the github (0.1.13) was not available through bioconda so we opted to download the next available version (0.1.14). Interestingly, the github also included a bed file to be used with vcftools, but there was no documentation on how these files were created. After consultation, we decided to move forward without the bed file, as it seemed its primary purpose was to mask out coinfected S. fitti, which we had already subsetted out early in our analysis process.  
```
#following command run in bcftools_view.txt
#identify snps
bcftools view -O v --threads 8 -m2 -M2 -v snps -o /vortexfs1/omics/env-bio/collaboration/dino_popgen/data/snp_calling/snpview_symb_mpileup.vcf /vortexfs1/omics/env-bio/collaboration/dino_popgen/data/snp_calling/symb_mpileup.vcf
#following command run in bcftools_filter.txt
#filter out medium and low-quality snps
bcftools filter -i 'QUAL >= 200' --threads 8 -O v -o /vortexfs1/omics/env-bio/collaboration/dino_popgen/data/snp_calling/HQ_snps.vcf /vortexfs1/omics/env-bio/collaboration/dino_popgen/data/snp_calling/snpview_symb_mpileup.vcf
#following command run in vcftools_filter.txt
#filter out sites with >20% missing data and indels
vcftools --vcf /vortexfs1/omics/env-bio/collaboration/dino_popgen/data/snp_calling/HQ_snps_cp.vcf --max-missing 0.8 --recode --remove-indels --out allhqSNPmm80
```
The final vcf file contained 57,800 loci, closely matching the 58,000 that were identified in the paper. One reason that we may have identified less snps than the manuscript is that we chose to omit a single large file from our analysis. These final SNPs were used as inputs for the respective downstream analyses.

## Analyses

### PCA

Analysis was conducted on SNP data in allhqSNPmm80.recode.vcf in R, after reformatting the metadata file so that sample orders matched exactly between files. Metadata reformatting was achieved with a mix of bash and R code (described in scripts/proc_for_vcf-seq-order.txt and scripts/supp_table_prep.R, summarized here):

```
# bash
grep CHROM allhqSNPmm80.recode.vcf | tr [:space:] '\n' > [INTERMEDIATE-FILE]
        # pull out sample names from vcf file and put those in a list in new file

#manually remove top few lines, so that first line is sample name

cut -d '/' -f9 vcf_chrom_head.txt | cut -d '.' -f1 > vcf_seq_order.txt
        # cut out just SRR numbers and print those to new text file

# R
seq_order = read.csv("vcf_seq_order.txt",sep='\n',header=F) #read in text file
seq_order$index = seq(1:47) #add index column and numbers

Supp_tab_1 <- read.csv("supp_table_1.csv", na.strings = c("n.a.","")) #read in original supp. table file

Supp_tab_1_single <- subset.data.frame(Supp_tab_1, Supp_tab_1$Colonization_Status == "single")
  # remove any samples with colonization status of "multiple" or "incomplete"

Supp_tab_1_47 <- subset.data.frame(Supp_tab_1_single, Supp_tab_1_single$total_reads < 80000000)
  # remove the 2 samples with 10x the total reads

colnames(Supp_tab_1_47)[1] = "Sample_ID" #fix column name

write.csv(Supp_tab_1_47, file = "supp_table_subset.csv",row.names = F) #write subsetted table

poptab <- read.csv("supp_table_subset.csv") #load in the subsetted table

poptab$index = NA #create empty index column in subsetted metadata file
for (i in 1:47){ #for all rows in metadata
        for (k in 1:47){ #for all sequences in list
                if(poptab$NCBI_SRA_Accesion[i]==seq_order$V1[k]){ #if the sample IDs match...
                        poptab$index[i] = seq_order$index[k] #...assign the sequence index number
                }
	}
}

poptab_ordered = poptab[order(poptab$index),] #re-order csv by vcf-derived index number

write.csv(poptab_ordered, file = "supp_table_ordered.csv",row.names = F) #write subsetted table
```

After metadata table was reformatted, anlaysis was conducted in R. Figure code is reported in full in jupyter_notebooks/figure-2_JBW.Rmd, but below is a summary:

```
# load all necessary libraries -- versions not specified
library(reshape2)
library(pcadapt)
library(cowplot)
library(ggplot2)
library(ggpubr)
library(vcfR)

# import data
poptab <- read.csv("supp_table_ordered.csv")
path <- "allhqSNPmm80.recode.vcf"
vcf <- read.pcadapt(path, type = "vcf", type.out = "matrix", ploidy = 1)
	# NOTE: ploidy argument has been deprecated, but code still ran

# run pcadapt
acp <- pcadapt(input = vcf, K = 5)

# spit out % variance... JBW gets PC1 = 15.7% and PC2 = 12.2%
acp$singular.values

# wrangle into new dataframe
df1 <- data.frame(Species = poptab$Host_Species, Loc = poptab$Population, 
                  ID = poptab$VCF_ID, pca = pca$data)
all_pca <- as.data.frame(df1)

# assign colors to match paper
cols2 = c("#000000","#999999","#0072B2")

# make plot
ggplot(all_pca, aes(x=pca.PC_i, y=pca.PC_j)) + #JBW moved x and y into general ggplot aesthetic
  geom_point(aes(color = Species, 
                 fill = Species,
                 shape = Loc),
             size = 4, 
             alpha = 0.85) +
  scale_shape_manual(values = c(25,16,17,15,23), 
                     breaks = c("Bahamas", "Belize", "Curacao", "Florida", "US Virgin Islands"), 
                     labels = c("Bahamas", "Belize", "Curacao", "Florida", "US Virgin Islands"), 
                     name = "Location") +
  scale_color_manual(values = cols2, 
                     breaks = c("A. cervicornis", "Hybrid", "A. palmata"),
                     labels = c("A. cervicornis", "Hybrid", "A. palmata"),
                     name = "S. 'fitti' host") +
  theme_classic2() +
  ylim(0.4, -0.4) + #JBW had to reverse the y-axis to get the spatial layout to match
  theme(legend.position = "right",
        axis.title = element_text(face = "bold", size = 14),
        plot.title = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        panel.grid = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank()) + 
  stat_ellipse(geom = "polygon", type = "t", alpha = 0.25, 
               aes(fill= Species), 
               show.legend = FALSE) + 
 scale_fill_manual(values = cols2, 
                    breaks = c("A. cervicornis", "Hybrid", "A. palmata"),
                    labels = c("A. cervicornis", "Hybrid", "A. palmata"),
                    name = "S. 'fitti' host") +
  labs(title = "57,810 'genotyping' S. fitti SNPs", x="PC1 (15.7% variance)", y="PC2 (12.2% variance)")
```

### Phylogeny
SNPs were filtered with vcftools to subset only the locations with complete coverage across samples (allhqSNPmm80_phylo_recode.vcf, 5658 SNPs in total). In order to generate a fasta alignment for use in phylogenetic analysis, vcfkit phylo (version 0.2.9) was used to extract all SNPs from the vcf file and generate an alignment from the variant calls (phylogeny_msa.fasta). RAxML-NG (version 0.9.0) was run to identify the optimal maximum likelihood tree. Parameters used with RAxML (GTR+FO+G nucleotide model, 100 bootstraps)) were pulled directly from the Reich et al. paper github. Project metadata was used to match the NCBI identifiers to the actual sample names (documented in jupyter notebook figure-3_EF.ipynb), and then the best tree identified in RAxML was visualized using iTOL (https://itol.embl.de/) and exported in both pdf and newick text format.

```
#following code located in vcftools_filter_phylo.txt
#use vcftools to filter out only snps with no missing data from the high quality SNP file
vcftools --vcf /vortexfs1/omics/env-bio/collaboration/dino_popgen/data/snp_calling/HQ_snps_cp.vcf --max-missing 1.0 --recode --remove-indels --out allhqSNPmm80_phylo

#following code located in vcfkit_run.txt
#generate multiple sequence alignment using multi-sample vcf file
#Cory ran this since vcfkit wouldn't install in my environment
vk phylo fasta /vortexfs1/omics/env-bio/collaboration/dino_popgen/output/snp_calling/allhqSNPmm80_phylo.recode.vcf > /vortexfs1/omics/env-bio/collaboration/dino_popgen/output/snp_calling/test.fasta
#renamed test.fasta to phylogeny_msa.fasta

#following code located in raxml_ng.txt
#parameters for raxml were pulled directly from the paper github
raxml-ng --msa /vortexfs1/omics/env-bio/collaboration/dino_popgen/output/snp_calling/test.fasta --model GTR+FO+G  --opt-branches on --opt-model on --tree pars{50},rand{50} --all  --bs-trees 100 --force --threads 2 --prefix /vortexfs1/omics/env-bio/collaboration/dino_popgen/output/phylogeny/raxml_output
```

### SNPs under selection

Outliers were identified using two different programs: BayeScan 2.1 and PCAdapt 4.3.3. This was the version of BayeScan used in the paper, however PCAdapt 4.0.3 was used in the paper. Between possible errors in the linux download of 4.0.3 and the availability of 4.3.3 on bioconda, we opted to use 4.3.3. To utilize BayeScan 2.1, you need to download it directly (cmpg.unibe.ch/software/BayeScan/download.html) and use one of the executable file based on your operating system.  
#### PCAdapt
To use PCAdapt, plug the .vcf file containing only high quality snps (created in the last step by vcf tools) into the PCAdapt R code provided.
Outliers were defined as SNPs with qvalues less than 0.05. 
```
# taken from jupyter-notebooks/PCAdapt_outliers.ipynb 
# need entire code to run PCAdapt, this only shows the essential parts
pc <- pcadapt(input = vcf, K = 10)
qval <- qvalue(pc$pvalues)$qvalues
alpha <- 0.05
outliers_pcadapt <- which(qval < alpha)
length(outliers_pcadapt)
```
#### BayeScan
To create the input files need for BayeSxcan, I used a perl code from github (scripts/bayescan/vcf2bayescan.pl). _Note: this is *not* what is used in the paper to create this file._ Because we want the BayeScan file to take population into account when identifying selection outliers, both the .vcf file and a population map are required.
Outliers were defined as SNPs with false discovery rates (FDR) less than 0.05.
```
# usage of vcf2bayescan.pl
# -p is the population file (change this depending on what population info you want to include), -v is the .vcf file
perl vcf2bayescan.pl -p host.txt -v allhqSNPmm80.recode.vcf
``` 

Once the input files are created, run BayeScan one time for each population file. Here, we ran one where information on the host species was included, one where the location of the sample was included, and one where the interaction of the two was included. Remember to change the path to BayeScan to match where your files are downloaded.
```
# input host file as population info
scripts/bayescan/BayeScan2.1/binaries/BayeScan2.1_linux64bits scripts/bayescan/host.txt -snp -threads 16 -od output/seq_outliers/ -o host -out_pilot -out_freq
# input location file as population info
scripts/bayescan/BayeScan2.1/binaries/BayeScan2.1_linux64bits scripts/bayescan/loc.txt -snp -threads 16 -od output/seq_outliers/ -o loc -out_pilot -out_freq
# input interaction of host and location file as population info
scripts/bayescan/BayeScan2.1/binaries/BayeScan2.1_linux64bits scripts/bayescan/hostXloc.txt -snp -threads 16 -od output/seq_outliers/ -o hostXloc -out_pilot -out_freq
``` 

