# gRNA_environment

R scripts for the analysis of the environment around a genomic location based on a combination of datasets. 

cell-type: hTERT RPE-1 cells

Genome.version: hg19

#### Datasets included to characterize the environmment
- ChIP-sequencing data:
  - H3K4me1 (GEO: GSE163315)
  - H3K4me3 (GEO: GSE163315)
  - H3K27me3 (GEO: GSE163315)
  - H4K26me3 (GEO: GSE163315)
  - H3K9me3 (GEO: GSM3105086) 

- ChromHMM analysis of the 5 marks from the ChIP-sequencing data (GEO: GSE163315)
- Lamin B1 DAM-iD (4Dnucleome: 4DNESHGTQ73M)
- DNaseI hypersensitivity (ENCODE dataset: ENCSR000EON, file: ENCFF128BPC)
- RNA-sequencing (GEO: GSE163315) 


## Analysis part 1: Run the [gRNAinfo](https://github.com/eskoeleman/gRNA_environment/blob/30425e7ab38baeb6bb9d0dddc418f33d8e29f45d/single_cutter_gRNAinfo_final.R) script

Aim: For a table of gRNAs, this script will output tables containing the information of all marks included in the analysis within a user-defined window.

#### R packages required
rtracklayer, Gviz, ggplot2, reshape2, dplyr, ggsci, genomation, chromoMap, BSgenome.Hsapiens.UCSC.hg19, data.table, reshape2, plyr

#### Input
- List of genomic locations of interest (one column containing the number of the location (eg 123456), one column containing the chromosome (eg "chr1") loaded as a .txt file. For simplicity, it should be named gRNAloc.txt. 
- Link to the .txt files of the datasets used for analysis of the environment (ChIP-seq, DAMID, DNaseI, RNA-seq). See links to GEO and 4Dnucleome. 
- User-defined window of interest around the break-site (eg 4 kb) for which analysis will be performed.

#### Output
- Multiple tables, one for each characteristic, containing information for each gRNA from the input-file around the break site (long-format). 
  - "1.DNaseI_hypersensitivity.txt"
  - "2.DAMID.txt"
  - "3.chromHMM_score.txt"
  - "4.RNAseq.txt"
  - "5.H3K9me3_weightmean.txt"

The output tables from the gRNAinfo script can subsequently be used in the gRNA filter script to select locations which contain the characteristics of interest. 


## Analysis part 2: Run the [gRNAfilter](https://github.com/eskoeleman/gRNA_environment/blob/30425e7ab38baeb6bb9d0dddc418f33d8e29f45d/single_cutter_gRNAfilter_final.R) script
Aim: based on user-defined adaptations, filter gRNAs using the output file for gRNAinfo.

#### R packages required
rtracklayer, Gviz, ggplot2, reshape2, dplyr, ggsci, genomation, chromoMap, BSgenome.Hsapiens.UCSC.hg19, data.table, reshape2, dplyr, ggplot2, chromoMap

#### Input
- Load all files from the gRNAinfo script, they are used in the analysis.

#### Filtering steps
By altering the parameters within the code, gRNAs can be filtered on different characteristics step-by-step. After each step, the results are plotted for visualisation and comparison. 
- Step 1: Filter on ChromHMM state. In this part of the code, a specific ChromHMM mark can be chosen, in combination with a percentage of this marks which should be present around the break-site. 
- Step 2: Filter on DNaseI peak maximum. A threshold for maximal DNaseI peak height can be set. 
- Step 3: Filter on DAMID state by selecting which percentage of iLAD region is present around the break site.
- Step 4: Filter on RNAseq data, where gRNAs can be filtered on presence in a genic or non-genic region.
- Step 5: Filter on H3K9me3 data. Select the minimal weighted peak height in the window around the break site. Both the weighted average of the peaks and the maximum peak within the window can be visualised for each gRNA. 
Remark: For the analysis, it is not required to perform all filtering steps. If needed, any step of the gRNA filtering script can be skipped to obtain only the required filtering steps for gRNA selection. 

After filtering, this script also includes the option to visualise the gRNA locations of interest using a chromomap plot, where gRNA locations are plotted on the chromosomes.

#### Output
- Table with the locations of all filtered gRNAs that have been selected.
- Optional: tables with the data from all characteristics, containing only the data of the filtered gRNAs.
- Graphs for the gRNAs that are selected with each filtering step (should be saved separately in R during the analysis).
- Chromomap plot of the locations of the selected gRNAs.
