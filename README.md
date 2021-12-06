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
   * Fig3_phylogeny_EF.pdf
   * Fig_5_manhattanplot_total.pdf
 - jupyter-notebooks/
   * data_download # DELETE
   * dino_popgen_final_writeup.ipynb
     - _a final comparison of reproduced results and original project_
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
   * load-seqs/
     - load_seqs.txt # MIGHT MOVE INTO MAIN FOLDER
     - seq_access_missing.txt #	MIGHT DELETE
     - seq_access_nums_short.txt # MIGHT DELETE
   * zip.txt
     - _script to zip files so they take up less space_
   * trim.txt
     - _script to trim sequences using cutadapt_
   * call_filter_snps/
     - mpileup.txt
       * _create a mpileup file with all samples_
     - bcftools_run.txt
       * _call variants_
     - bcftools_view.txt
       * _identify snps_
     - bcftools_filter.txt
       * _filter out low and medium quality SNPs_
     - vcftools_filter.tx
       * _filter out sites with >20% missing data_
   * samtools_seqs.txt
     - _convert .sam files to .bam, remove PCR duplicates, and get alignment stats_
   * snp_id.txt #MOVE TO CALL/FILTER AND DELETE MPILEUP.TXT AND BCFTOOLS_RUN.TXT
     - _code to create pileup file and call variants_
   * align_seqs.txt
     - _script to align trimmed seqs to draft genome_
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
     - _vcf2bayescan.pl output:_ host.txt, loc.txt, hostXloc.txt
   * proc_for_vcf-seq-order.txt
     - _text file explaining bash and R code used to subset and re-order metadata table to match sample order in allhqSNPmm80.recode.vcf_
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
     - _SNP data (output from bcftools)_
   * phylogeny/
     - _fasta alignment and raxml output_
     - [ALSO SHOULD THIS BE IN ANALYSIS OUTPUT?]
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
Samtools view was used to convert sam files to bam format in order to save space on the HPC. Then, the .bam files were sorted and indexed so the PCR duplicates could be removed with samtools rmdup. Now that duplicates are removed, this .bam file was used in downstream analyses. Additionlly, the .bam files were indexed again and samtools flagstat used to calculate alignment statistics. *Add part about python code to run through stats?*
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

### Phylogeny
*The short story:* SNPs were filtered with vcftools to subset only the locations with complete coverage across samples (allhqSNPmm80_phylo_recode.vcf, 5658 SNPs in total). In order to generate a fasta alignment for use in phylogenetic analysis, vcfkit phylo (version 0.2.9) was used to extract all SNPs from the vcf file and generate an alignment from the variant calls (phylogeny_msa.fasta). RAxML-NG (version 0.9.0) was run to identify the optimal maximum likelihood tree. Parameters used with RAxML (GTR+FO+G nucleotide model, 100 bootstraps)) were pulled directly from the Reich et al. paper github. Project metadata was used to match the NCBI identifiers to the actual sample names (documented in jupyter notebook figure-3_EF.ipynb), and then the best tree identified in RAxML was visualized using iTOL (https://itol.embl.de/) and exported in both pdf and newick text format.

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
*The long version:* Creating the phylogeny, while seemingly simple, ended up being an incredibly challenging process due to various issues and roadblocks. The step for transitioning from variant calls to multiple sequence alignment was not detailed in the manuscript or paper github, so we had to identify our own method for this. The package we used, vcfkit, was one of the only straightforward tools found online to achieve this, but was (for some reason) impossible to install on the HPC initially (days-long “Solving Environment” message, glibc error, etc.). Several other approaches were identified, but all proved unhelpful for various reasons:
* vcflib – Required individual-sample vcf files, output logfiles were all blank and unhelpful for troubleshooting. (https://github.com/vcflib/vcflib)
* vcf2msa – This was a wrapper file downloaded from github, ran in python but also yielded bugs and incomplete documentation made it challenging to troubleshoot. (https://github.com/tkchafin/vcf2msa.py)
* gatk4 – This approach was the most promising, as the FastaAlternateReferenceMaker function can take a vcf and create a fasta sequence for alignment. However, this approach required the final vcf to be split into individual sample vcfs and indexed. Significant time was spent getting gatk to work (initial installation trouble, then loaded the gatk module directly from HPC, then realized that we actually needed gatk4 to access FastaAlternateReferenceMaker which needed to be downloaded from bioconda). While we were able to generate final fasta files for each sample using this approach, the resulting files used the entire original reference, making them millions of base pairs long and not appropriate for use in alignment. 
* vcftools consensus – Similar to gatk, this output was not actually what we were looking for.  
Ultimately, vcfkit was installed successfully by an alternate HPC user (shoutout to Cory Berger!) and vcfkit was run to create the fasta file for alignment. 

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

