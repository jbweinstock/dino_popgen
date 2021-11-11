## dino_popgen
# **Final project for Environmental Bioinformatics, Fall 2021**

**Proposal Text**

In their study, Reich et al. (2021) use single nucleotide polymorphism (SNP) data to address whether strains of symbiotic dinoflagellates Symbiodinium 'fitti' are genetically differentiated across the two Caribbean acroporid coral species and their hybrid, with which S. 'fitti' is associated. Tissue samples were collected from 76 acroporids across their geographic distribution in the Caribbean, and resulting sequences were filtered, trimmed, aligned, and mapped to a draft assembly of the S. 'fitti' genome (available on dryad, https://doi.org/10.5061/dryad.xgxd254g8, ~326 MB). A subset of these samples were removed from downstream analysis due to co-infection with multiple S. ‘fitti’ strains. High-quality SNPs between the remaining sample genomes were identified and used to characterize population structure. Raw sequence data are publicly available on NCBI under SRA project PRJNA473816 (76 entries, ~340 GB). 

Code associated with this project, including all scripts used for data cleaning, analysis, and figure creation, has been uploaded to a github repository (https://github.com/hgreich/Sfitti). We propose to use the published raw data and code to repeat the following analysis steps: 1) sequence trimming and cleaning; 2) removal of samples with evidence of co-infection; 3) genome alignment and identification of high-quality SNPs; and 4) analysis of population structure using cleaned SNP data. If time allows, we may also undertake the final analysis of the paper, which identifies loci that showed strong evidence of selection. Our overall goal is to use manuscript resources, with modifications as necessary, to recreate Figure 2 (PCA analysis of S. ‘fitti’ genomic differentiation) and Figure 3 (dendogram displaying similarity between sample genomes) from Reich et al. (2021).

**Directory Contents**
 - README.md
 - notebook.txt
 - data/
   * SraRunTable_PRJNA473816.txt
     - project metadata (csv file)
   * sra/
     - 47 raw sequences from "shallow" reads
   * genome/
     - published draft genome for _S. 'fitti'_ from same study
   * trim-data/
     - trimmed sequences (after cutadapt run)
   * align-data/
     - ...
 - envs/
 - jupyter-notebooks/
 - logs/
 - output/
   * phylogeny/
   * popgen/
   * selection/
 - scripts/
   * load_seqs/
     - load_seqs.txt
     - one_num.txt
     - seq_access_missing.txt
     - seq_access_nums_short.txt
     - seq_access_nums.txt
   * zip.txt
   * trim.txt
   * align_seqs.txt
   * snp_id.txt
   * mpileup.txt
   * bcftools_run.txt

