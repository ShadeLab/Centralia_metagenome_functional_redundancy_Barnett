## Github Repository for
# Disturbance increases soil microbiome functional redundancy but decreases capacity for insurance via winnowed environmental responsiveness
## by Samuel Barnett and Ashley Shade
<i>This work is unpublished.</i>


### Data
Both the raw read data and metagenome assemblies for this study are available through NCBI under bioproject [PRJNA974462](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA974462/)

### To cite this work or code

TBD

### Abstract

Redundancy is the capacity of coexisting populations to perform similar functions, supporting stable outputs under disturbance. The insurance hypothesis proposes that diverse communities are more resilient because they can respond to varied environmental stressors. Both redundancy and insurance likely shape resilience, but their relative contributions and potential interactions remain unclear. Here, we examined functional redundancy and potential for insurance in soil microbial communities disturbed by an ongoing underground fire in Centralia, Pennsylvania, USA, using seven years of deep metagenome sequencing from heated and reference sites. Analyses of taxonomic and functional diversity showed that functional redundancy increased with disturbance intensity, driven by declines in both diversity and average genome size. Biogeochemistry-relevant metabolisms, such as carbohydrate and nitrogen pathways, showed greater per-genome investment with disturbance intensity, indicating enhanced redundancy. However, transcription factor evenness – linked to community environmental responsiveness- was lower in heated soils. Thus, prolonged disturbance that imposes a strong selection can increase redundancy and near-term stability while simultaneously diminishing insurance potential against future stress, revealing tradeoffs between these mechanisms in microbial community resilience.

### Contents

Code is split up into two directories: [Sequence_processing](https://github.com/ShadeLab/Centralia_metagenome_functional_redundancy_Barnett/tree/main/Sequence_processing) and [Analysis](https://github.com/ShadeLab/Centralia_metagenome_functional_redundancy_Barnett/tree/main/Analysis).

#### Sequence processing
Code used for sequence processing including ... Scripts were run using SLURM on the MSU HPCC using slurm batch files with suffix .sb and are numbered by their order in the processing workflow. Outputs such as logs, warnings, or errors if any, are designated by the suffix .out and named in accordence with the library, run number, and slurm batch file. 

#### Analysis
Formal analysis can be found under ... All analysis was run with R and code was run in Rmarkdown. In the analysis directory you'll find the raw Rmarkdown files (.Rmd), a github friendly markdown rendering (.md) and the associated figure files from the rendering in separate sub-directories. The analysis was broken down into multiple chunks in separate Rmarkdown files:

### Funding
This work was supported by the U.S. National Science Foundation CAREER award 1749544. This work was supported in part by Michigan State University through computational resources provided by the [Institute for Cyber-Enabled Research](https://icer.msu.edu/).

### More info
[ShadeLab](http://ashley17061.wixsite.com/shadelab/home)
