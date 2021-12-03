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
   * 
   * 
 - jupyter-notebooks/
   * data-download.ipynb
     - ...
   * dino_popgen_final_writeup.ipynb
     - ...
   * trimlog_ext.ipynb
     - ...
   * figure-2_JBW.Rmd
     - _R code for (re-)creation of figure 2_
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
     - bcftools_run.txt
     - bcftools_filter.txt
     - bcftools_view.txt
     - vcftools_filter.txt
   * samtools_seqs.txt
   * mpileup.txt
   * snp_id.txt
   * align_seqs.txt
     - _script to align trimmed seqs to draft genome_
   * phylogeny/
     - vcftools_filter_phylo.txt
     - vcftools_consensus.txt
     - vcf2msa/
       * vcf2msa.py
       * mpileup_all.txt
       * mpileup_all_*.log [SHOULD PROBS MOVE TO LOGS FOLDER]
       * vcf2msa_run.txt
   * bayescan/
     - ...
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
     - [EVIE SHOULD ENTER A BRIEF SUMMARY... ALSO SHOULD THIS BE IN ANALYSIS OUTPUT?]
