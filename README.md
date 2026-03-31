## Github Repository for
# Disturbance increases soil microbiome functional redundancy but decreases capacity for insurance via winnowed environmental responsiveness
## by Samuel Barnett and Ashley Shade
<i>This work is unpublished.</i>


### Data
Both the raw read data and metagenome assemblies for this study are available through NCBI under bioproject [PRJNA974462](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA974462/)

### To cite this work or code

Barnett, S.E., Shade, A., 2025. Disturbance increases soil microbiome functional redundancy but decreases capacity for insurance via winnowed environmental responsiveness. BioRxiv 2025.09.19.677300. doi:10.1101/2025.09.19.677300

### Abstract

Redundancy is the capacity of coexisting populations to perform similar functions, supporting stable outputs under disturbance. The insurance hypothesis proposes that diverse communities are more resilient because they can respond to varied environmental stressors. Both redundancy and insurance likely shape microbiome resilience, but their relative contributions and potential interactions remain unclear. Here, we examined functional redundancy and potential for insurance in soil microbial communities disturbed by an ongoing underground fire in Centralia, Pennsylvania, USA, using seven years of deep metagenome sequencing from heated and reference sites. Analyses of taxonomic and functional diversity showed that functional redundancy increased with disturbance intensity, driven by declines in both diversity and average genome size. Biogeochemistry-relevant metabolisms, such as carbohydrate and nitrogen pathways, showed greater per-genome investment with disturbance intensity, indicating enhanced redundancy. However, transcription factor evenness, which is linked to community environmental responsiveness, was lower in heated soils. Thus, prolonged disturbance that imposes a strong selection can increase redundancy and near-term stability while simultaneously diminishing insurance potential against future stress, revealing tradeoffs between these mechanisms in microbial community resilience.

### Contents

Code is split up into two directories: [Sequence_processing](https://github.com/ShadeLab/Centralia_metagenome_functional_redundancy_Barnett/tree/main/Sequence_processing) and [Analysis](https://github.com/ShadeLab/Centralia_metagenome_functional_redundancy_Barnett/tree/main/Analysis).

Data used within this code can be found in the [Data_sets](https://github.com/ShadeLab/Centralia_metagenome_functional_redundancy_Barnett/tree/main/Data_sets)

#### Sequence processing
Scripts were run using SLURM on the MSU HPCC using slurm batch files with suffix .sb and are numbered by their order in the processing workflow. Outputs such as logs, warnings, or errors if any, are designated by the suffix .out and named in accordence with the library, run number, and slurm batch file. 

#### Analysis
All analysis was run with R and code was run in Rmarkdown. In the analysis directory you'll find the raw Rmarkdown files (.Rmd), a github friendly markdown rendering (.md) and the associated figure files from the rendering in separate sub-directories. The analysis was broken down into multiple chunks in separate Rmarkdown files:
*  A1_Overall_stats.Rmd: Analysis of overall metagenome/community properties such as diversity, coverage, average genome size, taxonomy, and a map of the sites.
*  A2_Overall_annotations.Rmd: Broad scale analysis of annotated metagenome content including beta diversity measures, ordinations, and comparisons to 16S rRNA gene amplicon based analyses.
*  A3_MAGs_Functional_Redundancy.Rmd: Analysis of Metagenome Assembled Genomes (MAGs) including beta diversity, taxonomy, community functional redundancy. 
*  A4_Annotated_pathways.Rmd: Analysis of specific gene families and pathways from the assembled and annotated metagenomes, specifically nitrogen cycling pathways, carbohydrate metabolism pathways, and transcription factor diversity. 

#### Data sets
Data sets used for analysis can be found here. Data sets are:
xxxx

### Funding
This work was supported by the U.S. National Science Foundation CAREER award 1749544. This work was supported in part by Michigan State University through computational resources provided by the [Institute for Cyber-Enabled Research](https://icer.msu.edu/).

### More info
[ShadeLab](http://ashley17061.wixsite.com/shadelab/home)
