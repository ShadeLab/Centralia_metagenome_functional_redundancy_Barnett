Analysis of Annotated Pathways
================
Sam Barnett
12 August, 2025

- [Introduction](#introduction)
  - [Librarys and global variables](#librarys-and-global-variables)
  - [Metadata](#metadata)
  - [KEGG annotations](#kegg-annotations)
  - [Find single copy genes](#find-single-copy-genes)
- [Nitrogen cycling](#nitrogen-cycling)
  - [Find the genes in each pathway](#find-the-genes-in-each-pathway)
  - [Nitrogen cycling investment over soil
    temperature](#nitrogen-cycling-investment-over-soil-temperature)
  - [Nitrogen cycling investment over soil nitrate and ammonium
    concentrations](#nitrogen-cycling-investment-over-soil-nitrate-and-ammonium-concentrations)
- [Carbon cycling](#carbon-cycling)
  - [Find the genes in each pathway](#find-the-genes-in-each-pathway-1)
  - [Carbohydrate metabolism investment over soil
    temperature](#carbohydrate-metabolism-investment-over-soil-temperature)
  - [CAZyme investment over soil
    temperature](#cazyme-investment-over-soil-temperature)
- [Environmental response genes](#environmental-response-genes)
  - [Find transcription factor genes](#find-transcription-factor-genes)
  - [Transcription factor investment over soil
    temperature](#transcription-factor-investment-over-soil-temperature)
  - [Transcription factor diversity](#transcription-factor-diversity)
- [Session info](#session-info)

# Introduction

In this notebook, lets look at some notable pathways, metabolisms, and
other genome features from the whole metagenome assembly. We will be
using this to back up and further examine some of the findings from the
MAGs in the previous analysis. For this analysis we will be using a
term, per-genome investment as a measure of gene/pathway abundance. This
measure is a normalized measure to account for differences in sequencing
depth across samples and calculated as the number of genes (KEGG
orthologes) present in in a sample divided by the number of single copy
genes present.

## Librarys and global variables

Here are some libraries used in this analysis and the global variables
that will be used throughout. Mostly variables for consistent plotting.

``` r
# Libraries for data
library(dplyr)
library(phyloseq)
library(ape)
library(readxl)
library(jsonlite)

# Libraries for analysis
library(vegan)
library(picante)
library(Nonpareil)
library(ecotraj)
library(adiv)
library(ggkegg)

# Libraries for plotting
library(ggplot2)
source("/Users/sambarnett/Documents/Misc_code/paul_tol_colors.R")

# Functon for extracting legends
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)} 

# Site lists
used_sites = c("Cen08", "Cen11", "Cen14", "Cen15", "Cen16", "Cen17", "Cen19", 
               "Cen21", "Cen22", "Cen23")

# Setting repeated plot aesthetics
## Sites
site.col = paultol_colors(length(used_sites))
names(site.col) = used_sites

site.shape = c(22, 24, 24, 24, 24, 22, 24, 24, 24, 22)
names(site.shape) = used_sites

## Fire Classifications
FC.col = c("FireAffected" = "red", "Reference" = "grey")
FC.shape = c("FireAffected" = 24, "Reference" = 22)

# Basic plotting theme so as not to continually repeat it
publication_theme = theme_bw() +
  theme(axis.text = element_text(size=6),
        axis.title = element_text(size=7),
        legend.text = element_text(size=6),
        legend.title = element_text(size=7, hjust=0.5),
        strip.text = element_text(size=7),
        plot.title = element_text(size=8, hjust=0.5))

present_theme = theme_bw() +
  theme(axis.text = element_text(size=10),
        axis.title = element_text(size=12),
        legend.text = element_text(size=10),
        legend.title = element_text(size=12, hjust=0.5),
        strip.text = element_text(size=10),
        plot.title = element_text(size=14, hjust=0.5))
```

## Metadata

Read in the metadata.

``` r
# Sample metadata
sample.meta = read_xlsx("/Users/sambarnett/Documents/Shade_lab/Centralia_project/Centralia_soil_metadata.xlsx", 
                        sheet = "Metagenomic_samples", na="NA") %>%
  filter(SampleID != "Cen08_07102019_R1") %>%
  arrange(SiteID, Year) %>%
  mutate(Seq_number = row_number()) %>%
  mutate(nonpareil_file = paste("/Users/sambarnett/Documents/Shade_lab/Centralia_project/Metagenomics/Data/nonpareil/", SampleID, "_S", Seq_number, ".npo", sep=""),
         SequenceID = paste(SampleID, Seq_number, sep="_S"))

# Contig mapped Read counts
mapped_reads.df = read.table("/Users/sambarnett/Documents/Shade_lab/Centralia_project/Metagenomics/Data/Mapped_read_totals.txt", 
                          header=TRUE, sep="\t", comment.char = "", quote = "")
```

## KEGG annotations

Read in the table of KEGG orthologue annotations.

``` r
KEGG_long.df = read.table("/Users/sambarnett/Documents/Shade_lab/Centralia_project/Metagenomics/Data/Annotations/KEGG_COG_annotations.txt",
                       header=TRUE, sep="\t", comment.char = "", quote = "") %>%
  filter(KEGG_ortho_kofamscan != "") %>%
  mutate(mapped_reads = Plus_reads + Minus_reads) %>%
  select(locus_tag, KEGG_ortho_kofamscan, SequenceID, mapped_reads) %>%
  arrange(locus_tag, KEGG_ortho_kofamscan) %>%
  left_join(mapped_reads.df, by = "SequenceID") %>%
  left_join(read.table("/Users/sambarnett/Documents/Shade_lab/Centralia_project/Metagenomics/Data/Annotations/prokka_CDS_annotations.txt",
                       header=TRUE, sep="\t", comment.char = "", quote = "") %>%
              select(locus_tag, SequenceID, length_bp)) %>%
  mutate(RPKM = mapped_reads/((length_bp/1000)*(total_mapped_reads/1000000))) %>%
  mutate(RPKM = ifelse(is.na(RPKM), 0, RPKM))
```

## Find single copy genes

To normalize gene presenece accounting for different sequencing depths
across samples, lets find the single copy genes that we will use as a
proxy for genome count. These single copy genes are those used by checkM
(Parks et al. 2015.)

First get a list of the single copy genes we will utilize identified by
their KEGG othologue ID (ko).

``` r
# Small subunit ribosomal proteins
KEGG_SSU_proteins.df = data.frame(KEGG_ortho_kofamscan = c("K02967", "K02982", "K02986", "K02988", "K02992", 
                                                           "K02994", "K02996", "K02946", "K02948", "K02950", 
                                                           "K02952", "K02956", "K02959", "K02961", "K02965", "K02968"),
                                  ribosomal_protein = c("S2", "S3", "S4", "S5", "S7", 
                                                        "S8", "S9", "S10", "S11", "S12", 
                                                        "S13", "S15", "S16", "S17", "S19", "S20"),
                                  ribosomal_subunit = "Small subunit")

# Large subunit ribosomal proteins
KEGG_LSU_proteins.df = data.frame(KEGG_ortho_kofamscan = c("K02886", "K02906", "K02926", "K02931", 
                                                           "K02864", "K02867", "K02871", "K02874", "K02876", "K02878", "K02879",
                                                           "K02881", "K02884", "K02887", "K02888", "K02890", "K02892", "K02895", "K02899"),
                                  ribosomal_protein = c("L2", "L3", "L4", "L5",
                                                        "L10", "L11", "L13", "L14", "L15", "L16", "L17",
                                                        "L18", "L19", "L20", "L21", "L22", "L23", "L24", "L27"),
                                  ribosomal_subunit = "Large subunit")

# Non-ribosomal proteins
KEGG_SC_proteins.df = data.frame(KEGG_ortho_kofamscan = c("K01873", "K03553", "K03076", "K01889", "K01869", "K01892",
                                                          "K02469", "K01876", "K01872"),
                                  ribosomal_protein = c("valS", "recA", "secY", "pheS", "leuS", "hisS", 
                                                        "gyrA", "aspS", "alaS"),
                                  ribosomal_subunit = "NonRibosomal")

# Join them together
KEGG_ribosomal_proteins.df = rbind(KEGG_SSU_proteins.df, KEGG_LSU_proteins.df, KEGG_SC_proteins.df)
```

``` r
#KEGG_SSU_proteins.df = data.frame(KEGG_ortho_kofamscan = c("K02967", "K02982", "K02986", "K02988", "K02992", 
#                                                           "K02994", "K02996", "K02946", "K02948", "K02950", 
#                                                           "K02952", "K02954", "K02956", "K02961", "K02965"),
#                                  ribosomal_protein = c("S2", "S3", "S4", "S5", "S7", 
#                                                        "S8", "S9", "S10", "S11", "S12", 
#                                                        "S13", "S14", "S15", "S17", "S19"),
#                                  ribosomal_subunit = "Small subunit")

#KEGG_LSU_proteins.df = data.frame(KEGG_ortho_kofamscan = c("K02863", "K02886", "K02906", "K02931", "K02933", 
#                                                           "K02864", "K02867", "K02871", "K02874", "K02876",
#                                                           "K02881", "K02890", "K02892", "K02895", "K02904",
#                                                           "K02907"),
#                                  ribosomal_protein = c("L1", "L2", "L3", "L5", "L6",
#                                                        "L10", "L11", "L13", "L14", "L15",
#                                                        "L18", "L22", "L23", "L24", "L29", "L30"),
#                                  ribosomal_subunit = "Large subunit")
```

Now pull out the the single copy genes from the metagenomes and
summarize within samples. We will be using the median count.

``` r
# Get the single copy genes
KEGG_SCG.df = filter(KEGG_long.df, KEGG_ortho_kofamscan %in% KEGG_ribosomal_proteins.df$KEGG_ortho_kofamscan) %>%
  group_by(SequenceID, locus_tag, KEGG_ortho_kofamscan) %>%
  summarize(n_hits = n()) %>%
  ungroup %>%
  arrange(-n_hits, locus_tag) %>%
  left_join(KEGG_ribosomal_proteins.df)

# Summarize counts within samples
KEGG_SCG.sum = KEGG_SCG.df %>%
  group_by(SequenceID, ribosomal_protein) %>%
  summarize(n_hits = n()) %>%
  ungroup %>%
  tidyr::spread(key="ribosomal_protein", value="n_hits") %>%
  tidyr::gather(key="ribosomal_protein", value="n_hits", -SequenceID) %>%
  mutate(n_hits = ifelse(is.na(n_hits), 0, n_hits)) %>%
  group_by(SequenceID) %>%
  summarize(median_RP_genes = median(n_hits),
            mean_RP_genes = mean(n_hits),
            sd_RP_genes = sd(n_hits),
            n_riboP = n()) %>%
  ungroup %>%
  mutate(SE_RP_genes = sd_RP_genes/sqrt(n_riboP))
```

Before moving on, lets take a look at the number of these single copy
genes we see across samples. Since metagenome coverage and diversity
changes across these samples I expect to see similar trends here.

``` r
# Linear mixed effects over temperature
SCG_temp.model = lme(mean_RP_genes ~ CoreTemp_C, random = ~1|SiteID, data=left_join(KEGG_SCG.sum, sample.meta))
SCG_temp.model.sum = data.frame(intercept = summary(SCG_temp.model)$tTable[1],
                                      slope = summary(SCG_temp.model)$tTable[2],
                                      p_value = summary(SCG_temp.model)$tTable[10])

SCG_temp.model.sum = SCG_temp.model.sum %>%
  mutate(sig = ifelse(p_value < 0.05, "p < 0.05", "p ≥ 0.05"))
  

SCG_temp.plot = ggplot(data=left_join(KEGG_SCG.sum, sample.meta), aes(x=CoreTemp_C, y=mean_RP_genes)) +
  geom_point(aes(shape=FireClassification, fill=SiteID), size=2) +
  geom_abline(data = filter(SCG_temp.model.sum, p_value < 0.05),
              aes(intercept = intercept, slope = slope),
              linetype = 2, linewidth=1, color="black") +
  scale_shape_manual(values=FC.shape) +
  scale_fill_manual(values=site.col) +
  labs(x="Soil temperature (˚C)", y="Median single copy genes") +
  publication_theme +
  guides(fill=guide_legend(override.aes=list(shape=site.shape), ncol = 2))
SCG_temp.plot
```

![](A4_Annotated_pathways_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

``` r
ggsave(SCG_temp.plot, file="/Users/sambarnett/Documents/Shade_lab/Centralia_project/Metagenomics/Manuscript/Figures/Supplemental/FigS11.tiff",
       device="tiff", width=5, height=3.5, units="in", bg = "white")
```

Call up the nonpareil data and compare the number of single copy genes
to both metagenome coverage and diversity.

``` r
# Read in the the nonpareil data. These are separate files.
nonpareil.full = Nonpareil.set(sample.meta$nonpareil_file)
```

![](A4_Annotated_pathways_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

``` r
nonpareil.sum = data.frame(summary(nonpareil.full)) %>%
  tibble::rownames_to_column(var="SequenceID") %>%
  left_join(KEGG_SCG.sum, by = "SequenceID") %>%
  left_join(sample.meta, by = "SequenceID")

# Linear relationship to coverage
SCG_coverage.model = lm(mean_RP_genes ~ C, data=nonpareil.sum)
SCG_coverage.model.sum = summary(SCG_coverage.model)
SCG_coverage.model.sum
```

    ## 
    ## Call:
    ## lm(formula = mean_RP_genes ~ C, data = nonpareil.sum)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -22.3519  -7.9366  -0.7976   7.1141  23.4875 
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)  -279.22      19.98  -13.97   <2e-16 ***
    ## C             377.79      23.08   16.37   <2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 10.49 on 67 degrees of freedom
    ## Multiple R-squared:    0.8,  Adjusted R-squared:  0.797 
    ## F-statistic:   268 on 1 and 67 DF,  p-value: < 2.2e-16

``` r
SCG_coverage.plot = ggplot(data=nonpareil.sum, aes(x=C, y=mean_RP_genes)) +
  geom_point(aes(shape=FireClassification, fill=SiteID), size=2) +
  geom_abline(intercept = SCG_coverage.model.sum$coefficients[1],
              slope = SCG_coverage.model.sum$coefficients[2],
              linetype = 2, linewidth=1, color="black") +
  scale_shape_manual(values=FC.shape) +
  scale_fill_manual(values=site.col) +
  labs(x="Metagenome coverage", y="Median single copy genes") +
  publication_theme +
  guides(fill=guide_legend(override.aes=list(shape=site.shape), ncol = 2))
SCG_coverage.plot
```

![](A4_Annotated_pathways_files/figure-gfm/unnamed-chunk-8-2.png)<!-- -->

``` r
# Linear relationship to diversity
SCG_diversity.model = lm(mean_RP_genes ~ diversity, data=nonpareil.sum)
SCG_diversity.model.sum = summary(SCG_diversity.model)
SCG_diversity.model.sum
```

    ## 
    ## Call:
    ## lm(formula = mean_RP_genes ~ diversity, data = nonpareil.sum)
    ## 
    ## Residuals:
    ##     Min      1Q  Median      3Q     Max 
    ## -36.307 -10.705   1.205  10.646  37.317 
    ## 
    ## Coefficients:
    ##             Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)   437.82      49.78   8.795 8.92e-13 ***
    ## diversity     -18.92       2.41  -7.853 4.43e-11 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 16.92 on 67 degrees of freedom
    ## Multiple R-squared:  0.4793, Adjusted R-squared:  0.4715 
    ## F-statistic: 61.67 on 1 and 67 DF,  p-value: 4.435e-11

``` r
SCG_diversity.plot = ggplot(data=nonpareil.sum, aes(x=diversity, y=mean_RP_genes)) +
  geom_point(aes(shape=FireClassification, fill=SiteID), size=2) +
  geom_abline(intercept = SCG_diversity.model.sum$coefficients[1],
              slope = SCG_diversity.model.sum$coefficients[2],
              linetype = 2, linewidth=1, color="black") +
  scale_shape_manual(values=FC.shape) +
  scale_fill_manual(values=site.col) +
  labs(x="Metagenome diversity", y="Median single copy genes") +
  publication_theme +
  guides(fill=guide_legend(override.aes=list(shape=site.shape), ncol = 2))
SCG_diversity.plot
```

![](A4_Annotated_pathways_files/figure-gfm/unnamed-chunk-8-3.png)<!-- -->

As you can see there is a strong correlation to coverage which makes
perfect sense. Less coverage means less single copy genes recovered.
Therefore single copy genes should be a good normalization for gene
counts.

# Nitrogen cycling

For the first pathway to look at, lets consider nitrogen cycling. For
this analysis lets look at genes for different nitrogen cycling pathways
as defined by KEGG modules.

## Find the genes in each pathway

First lets find the nitrogen cycling genes for each pathway and count
them.

``` r
# Define these pathways with the member KEGG orthologues
N_Cycling.def = rbind(data.frame(KEGG_ortho_kofamscan = module("M00175")@definitions[[1]]$definition_kos,
                                 N_path = "Nitrogen Fixation"),
                      data.frame(KEGG_ortho_kofamscan = module("M00531")@definitions[[1]]$definition_kos,
                                 N_path = "Assimilatory nitrate reduction"),
                      data.frame(KEGG_ortho_kofamscan = module("M00530")@definitions[[1]]$definition_kos,
                                 N_path = "Dissimilatory nitrate reduction"),
                      data.frame(KEGG_ortho_kofamscan = module("M00529")@definitions[[1]]$definition_kos,
                                 N_path = "Denitrification"),
                      data.frame(KEGG_ortho_kofamscan = module("M00528")@definitions[[1]]$definition_kos,
                                 N_path = "Nitrification"),
                      data.frame(KEGG_ortho_kofamscan = module("M00804")@definitions[[1]]$definition_kos,
                                 N_path = "comammox"),
                      data.frame(KEGG_ortho_kofamscan = module("M00973")@definitions[[1]]$definition_kos,
                                 N_path = "Anammox"))

# Now pull out genes for the pathways from our annotations and summarize
KEGG_N_Cycling.sum = KEGG_long.df %>%
  inner_join(N_Cycling.def) %>%
  group_by(SequenceID, N_path) %>%
  summarize(n_genes = n()) %>%
  ungroup %>%
  tidyr::spread(key="N_path", value="n_genes") %>%
  tidyr::gather(key="N_path", value="n_genes", -SequenceID) %>%
  mutate(n_genes = ifelse(is.na(n_genes), 0, n_genes)) %>%
  left_join(KEGG_SCG.sum) %>%
  mutate(pop_invest = n_genes/mean_RP_genes) %>%
  left_join(sample.meta)
```

## Nitrogen cycling investment over soil temperature

Now lets quantify the per-genome investment in these nitrogen cycling
pathways. First lets see how nitrogen cycling genes change over
temperature.

``` r
# Linear mixed effects over temperature
Ncyc_temp.model.sum = data.frame()
for(N_p in unique(N_Cycling.def$N_path)){
  sub.Ncyc_temp.model = lme(pop_invest ~ CoreTemp_C, random = ~1|SiteID, data=filter(KEGG_N_Cycling.sum, N_path == N_p))
  Ncyc_temp.model.sum = rbind(Ncyc_temp.model.sum,
                              data.frame(intercept = summary(sub.Ncyc_temp.model)$tTable[1],
                                         slope = summary(sub.Ncyc_temp.model)$tTable[2],
                                         p_value = summary(sub.Ncyc_temp.model)$tTable[10],
                                         N_path = N_p))
}


Ncyc_temp.model.sum = Ncyc_temp.model.sum %>%
  mutate(padj = p.adjust(p_value, method = "BH")) %>%
  mutate(sig = ifelse(padj < 0.05, "p < 0.05", "p ≥ 0.05"))
  
Ncyc_temp.plot = ggplot(data=KEGG_N_Cycling.sum, aes(x=CoreTemp_C, y=pop_invest)) +
  geom_point(aes(shape=FireClassification, fill=SiteID), size=2) +
  geom_abline(data = filter(Ncyc_temp.model.sum, p_value < 0.05),
              aes(intercept = intercept, slope = slope),
              linetype = 2, linewidth=1, color="black") +
  scale_shape_manual(values=FC.shape) +
  scale_fill_manual(values=site.col) +
  labs(x="Soil temperature (˚C)", y="Genes per genome") +
  publication_theme +
  facet_wrap(~N_path, scales = "free_y", nrow=2) +
  #theme(legend.position = "bottom",
  #     legend.direction = "vertical") +
  theme(legend.position = c(0.88,0.2)) +
  guides(fill=guide_legend(override.aes=list(shape=site.shape), ncol = 2),
         shape=guide_legend(ncol = 2))

Ncyc_temp.plot
```

![](A4_Annotated_pathways_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->

``` r
# Save plot for publication
ggsave(Ncyc_temp.plot, file="/Users/sambarnett/Documents/Shade_lab/Centralia_project/Metagenomics/Manuscript/Figures/Fig3.tiff",
       device="tiff", width=7, height=4.5, units="in", bg = "white")
```

## Nitrogen cycling investment over soil nitrate and ammonium concentrations

Next lets see how nitrogen cycling genes change over concentrations of
nitrate and ammonium in the soils.

``` r
# Linear mixed effects over nitrate levels
Ncyc_nitrate.model.sum = data.frame()
for(N_p in unique(N_Cycling.def$N_path)){
  sub.Ncyc_nitrate.model = lme(pop_invest ~ NO3N_ppm, random = ~1|SiteID, data=filter(KEGG_N_Cycling.sum, N_path == N_p))
  Ncyc_nitrate.model.sum = rbind(Ncyc_nitrate.model.sum,
                              data.frame(intercept = summary(sub.Ncyc_nitrate.model)$tTable[1],
                                         slope = summary(sub.Ncyc_nitrate.model)$tTable[2],
                                         p_value = summary(sub.Ncyc_nitrate.model)$tTable[10],
                                         N_path = N_p))
}
  
Ncyc_nitrate.model.sum = Ncyc_nitrate.model.sum %>%
  mutate(padj = p.adjust(p_value, method = "BH")) %>%
  mutate(sig = ifelse(padj < 0.05, "p < 0.05", "p ≥ 0.05"))
  
Ncyc_nitrate.plot = ggplot(data=KEGG_N_Cycling.sum, aes(x=NO3N_ppm, y=pop_invest)) +
  geom_point(aes(shape=FireClassification, fill=SiteID), size=2) +
  geom_abline(data = filter(Ncyc_nitrate.model.sum, padj < 0.05),
              aes(intercept = intercept, slope = slope),
              linetype = 2, linewidth=1, color="black") +
  scale_shape_manual(values=FC.shape) +
  scale_fill_manual(values=site.col) +
  labs(x="Nitrate nitrogen (ppm)", y="Genes per genome") +
  publication_theme +
  facet_wrap(~N_path, scales = "free_y", nrow=1) +
  theme(legend.position = "bottom",
        legend.direction = "vertical") +
  guides(fill=guide_legend(override.aes=list(shape=site.shape), nrow = 2))

# Linear mixed effects over ammonium levels
Ncyc_ammonium.model.sum = data.frame()
for(N_p in unique(N_Cycling.def$N_path)){
  sub.Ncyc_ammonium.model = lme(pop_invest ~ NH4N_ppm, random = ~1|SiteID, data=filter(KEGG_N_Cycling.sum, N_path == N_p))
  Ncyc_ammonium.model.sum = rbind(Ncyc_ammonium.model.sum,
                              data.frame(intercept = summary(sub.Ncyc_ammonium.model)$tTable[1],
                                         slope = summary(sub.Ncyc_ammonium.model)$tTable[2],
                                         p_value = summary(sub.Ncyc_ammonium.model)$tTable[10],
                                         N_path = N_p))
}
  
Ncyc_ammonium.model.sum = Ncyc_ammonium.model.sum %>%
  mutate(padj = p.adjust(p_value, method = "BH")) %>%
  mutate(sig = ifelse(padj < 0.05, "p < 0.05", "p ≥ 0.05"))
  
Ncyc_ammonium.plot = ggplot(data=KEGG_N_Cycling.sum, aes(x=NH4N_ppm, y=pop_invest)) +
  geom_point(aes(shape=FireClassification, fill=SiteID), size=2) +
  geom_abline(data = filter(Ncyc_ammonium.model.sum, padj < 0.05),
              aes(intercept = intercept, slope = slope),
              linetype = 2, linewidth=1, color="black") +
  scale_shape_manual(values=FC.shape) +
  scale_fill_manual(values=site.col) +
  labs(x="Ammonium nitrogen (ppm)", y="Genes per genome") +
  publication_theme +
  facet_wrap(~N_path, scales = "free_y", nrow=1) +
  theme(legend.position = "bottom",
        legend.direction = "vertical") +
  guides(fill=guide_legend(override.aes=list(shape=site.shape), nrow = 2))


# Plot together

cowplot::plot_grid(Ncyc_nitrate.plot + theme(legend.position = "none"), 
                   Ncyc_ammonium.plot + theme(legend.position = "none"), 
                   g_legend(Ncyc_nitrate.plot), rel_heights = c(1,1,0.3),
                   ncol=1, labels=c("A", "B"), label_size = 8)
```

![](A4_Annotated_pathways_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->

# Carbon cycling

Next lets take a look at carbon cycling. Specifically, lets look at
carbohydrate metabolism. For this analysis, lets look at level C
pathways under level B metabolism 09101 “Carbohydrate metabolism”.

## Find the genes in each pathway

First lets find the carbohydrate metabolism genes. To do this we need to
find all ko within each classification level.

``` r
# Get main metabolism dataset from KEGG.
ko00001.json = fromJSON("/Users/sambarnett/Documents/Shade_lab/Centralia_project/Metagenomics/Data/Annotations/ko00001.json", flatten = TRUE)

# Get the different levels of classification from the JSON file for each ko
ko00001_children = ko00001.json$children

ko00001.df = data.frame()
for (c1 in 1:nrow(ko00001_children)){
  LevelA_name = ko00001_children$name[[c1]]
  LevelA_children = ko00001_children$children[[c1]]
  for (c2 in 1:nrow(LevelA_children)){
    LevelB_name = LevelA_children$name[[c2]]
    LevelB_children = LevelA_children$children[[c2]]
    for (c3 in 1:nrow(LevelB_children)){
      LevelC_name = LevelB_children$name[[c3]]
      sub_BRITE.df = LevelB_children$children[[c3]]
      if (!is.null(sub_BRITE.df)){
        sub_BRITE.df = sub_BRITE.df %>%
          dplyr::rename(LevelD = name) %>%
          mutate(LevelC = LevelC_name,
                 LevelB = LevelB_name,
                 LevelA = LevelA_name) %>%
          select(LevelA, LevelB, LevelC, LevelD)
        ko00001.df = rbind(ko00001.df, sub_BRITE.df)
      }
    }
  }
}

ko00001.df = ko00001.df %>%
  mutate(KEGG_ortho_kofamscan = gsub(" .*", "", LevelD))
```

Pull out the carbohydrate metabolism genes from the dataset.

``` r
# Get the ko for carbohydrate metabolism and rename some of the levelC names so they fit better in figures.
Carbo_kegg.df = ko00001.df %>%
  filter(LevelB == "09101 Carbohydrate metabolism") %>%
  mutate(C_path = gsub(" \\[.*", "", LevelC)) %>%
  select(KEGG_ortho_kofamscan, C_path) %>%
  unique %>%
  mutate(C_path = gsub("metabolism", "met.", C_path)) %>%
  mutate(C_path = gsub("and", "&", C_path)) %>%
  mutate(C_path = gsub("interconversions", "interconv.", C_path))

# Get carbohydrate metabolism genes from the metagenomes and summarize across samples
KEGG_Carbo.sum = KEGG_long.df %>%
  inner_join(Carbo_kegg.df) %>%
  group_by(SequenceID, C_path) %>%
  summarize(n_genes = n()) %>%
  ungroup %>%
  tidyr::spread(key="C_path", value="n_genes") %>%
  tidyr::gather(key="C_path", value="n_genes", -SequenceID) %>%
  mutate(n_genes = ifelse(is.na(n_genes), 0, n_genes)) %>%
  left_join(KEGG_SCG.sum) %>%
  mutate(pop_invest = n_genes/mean_RP_genes) %>%
  #mutate(pop_invest = n_genes) %>%
  left_join(sample.meta)
```

## Carbohydrate metabolism investment over soil temperature

Now lets quantify the per-genome investment in these carbon cycling
pathways.

``` r
# Linear mixed effects over temperature
Ccyc_temp.model.sum = data.frame()
for(N_p in unique(Carbo_kegg.df$C_path)){
  sub.Ccyc_temp.model = lme(pop_invest ~ CoreTemp_C, random = ~1|SiteID, data=filter(KEGG_Carbo.sum, C_path == N_p))
  Ccyc_temp.model.sum = rbind(Ccyc_temp.model.sum,
                              data.frame(intercept = summary(sub.Ccyc_temp.model)$tTable[1],
                                         slope = summary(sub.Ccyc_temp.model)$tTable[2],
                                         p_value = summary(sub.Ccyc_temp.model)$tTable[10],
                                         C_path = N_p))
}
  
Ccyc_temp.model.sum = Ccyc_temp.model.sum %>%
  mutate(padj = p.adjust(p_value, method = "BH")) %>%
  mutate(sig = ifelse(padj < 0.05, "p < 0.05", "p ≥ 0.05"))
  

Ccyc_temp.plot = ggplot(data=KEGG_Carbo.sum, aes(x=CoreTemp_C, y=pop_invest)) +
  geom_point(aes(shape=FireClassification, fill=SiteID), size=2) +
  geom_abline(data = filter(Ccyc_temp.model.sum, p_value < 0.05),
              aes(intercept = intercept, slope = slope),
              linetype = 2, linewidth=1, color="black") +
  scale_shape_manual(values=FC.shape) +
  scale_fill_manual(values=site.col) +
  labs(x="Soil temperature (˚C)", y="Genes per genome") +
  publication_theme +
  facet_wrap(~C_path, scales = "free_y", ncol=3) +
  theme(legend.position = "bottom",
        legend.direction = "vertical") +
  guides(fill=guide_legend(override.aes=list(shape=site.shape), nrow = 2),
         shape=guide_legend(nrow = 2))

Ccyc_temp.plot
```

![](A4_Annotated_pathways_files/figure-gfm/unnamed-chunk-14-1.png)<!-- -->

``` r
# Save plot for publications
ggsave(Ccyc_temp.plot, file="/Users/sambarnett/Documents/Shade_lab/Centralia_project/Metagenomics/Manuscript/Figures/Fig4.tiff",
       device="tiff", width=7, height=7, units="in", bg = "white")
```

## CAZyme investment over soil temperature

Lets double check these findings by looking specifically at CAZymes
(carbohydrate active enzymes). Rather than using KEGG orthologues, we
will use genes annotated using the CAZy database. First call up the
CAZyme data

``` r
# Import CAZymes
CAZyme.df = read.table("/Users/sambarnett/Documents/Shade_lab/Centralia_project/Metagenomics/Data/Annotations/Specific_enzyme_annotations.txt", 
                       header=TRUE, sep="\t", comment.char = "", quote = "") %>%
  filter(Enzyme_type == "CAZyme")

# Filter to genes we are more confident are CAZymes.
CAZyme_filt.df = CAZyme.df %>%
  filter(alignment_probability > 0.8) %>%
  mutate(Target = gsub(".hmm", "", Target)) %>%
  group_by(SequenceID, locus_tag) %>%
  summarize(n_modules = n(),
            CAZy_modules = paste(Target, collapse = ";")) %>%
  ungroup

# Summarize per-genome investment in total CAZymes
CAZyme.sum = CAZyme_filt.df %>%
  group_by(SequenceID) %>%
  summarize(n_genes = n()) %>%
  ungroup %>%
  left_join(KEGG_SCG.sum) %>%
  mutate(pop_invest = n_genes/mean_RP_genes) %>%
  left_join(sample.meta, by = "SequenceID")

# Summarize per-genome investment in glycoside hydrolases
GH.sum = CAZyme_filt.df %>%
  filter(grepl("GH", CAZy_modules)) %>%
  group_by(SequenceID) %>%
  summarize(n_genes = n()) %>%
  ungroup %>%
  left_join(KEGG_SCG.sum) %>%
  mutate(pop_invest = n_genes/mean_RP_genes) %>%
  left_join(sample.meta, by = "SequenceID")
```

First lets take a look at all CAZymes over temperature.

``` r
# Wilcoxon test over fire classification
CAZyme_FC.wilcox = wilcox.test(x=filter(CAZyme.sum, FireClassification=="Reference")$pop_invest,
                               y=filter(CAZyme.sum, FireClassification=="FireAffected")$pop_invest,
                               conf.int=TRUE, conf.level=0.95)
CAZyme_FC.wilcox
```

    ## 
    ##  Wilcoxon rank sum exact test
    ## 
    ## data:  filter(CAZyme.sum, FireClassification == "Reference")$pop_invest and filter(CAZyme.sum, FireClassification == "FireAffected")$pop_invest
    ## W = 111, p-value = 5.109e-08
    ## alternative hypothesis: true location shift is not equal to 0
    ## 95 percent confidence interval:
    ##  -37.64438 -21.16768
    ## sample estimates:
    ## difference in location 
    ##              -28.79937

``` r
CAZyme_FC.plot = ggplot(data=CAZyme.sum, aes(x=FireClassification, y=pop_invest)) +
  geom_boxplot(outlier.shape=NA) +
  geom_jitter(aes(fill=SiteID, shape=FireClassification), size=2, width=0.25) +
  annotate("text", label="< 0.001", fontface="bold", x=1.5, y=max(CAZyme.sum$pop_invest), size=6*5/14) +
  #lims(y=c(NA, 22)) +
  scale_fill_manual(values=site.col) +
  scale_shape_manual(values=FC.shape) +
  labs(x="FireClassification", y="CAZyme genes per genome") +
  publication_theme +
  guides(fill=guide_legend(override.aes=list(shape=site.shape), ncol=2))
CAZyme_FC.plot
```

![](A4_Annotated_pathways_files/figure-gfm/unnamed-chunk-16-1.png)<!-- -->

``` r
# Linear mixed effects over temperature
CAZyme_temp.model = lme(pop_invest ~ CoreTemp_C, random = ~1|SiteID, data=CAZyme.sum)
summary(CAZyme_temp.model)
```

    ## Linear mixed-effects model fit by REML
    ##   Data: CAZyme.sum 
    ##        AIC      BIC    logLik
    ##   490.9005 499.7193 -241.4503
    ## 
    ## Random effects:
    ##  Formula: ~1 | SiteID
    ##         (Intercept) Residual
    ## StdDev:    15.10688 6.370932
    ## 
    ## Fixed effects:  pop_invest ~ CoreTemp_C 
    ##                Value Std.Error DF  t-value p-value
    ## (Intercept) 61.37431  6.232507 58 9.847452  0.0000
    ## CoreTemp_C   0.33548  0.169737 58 1.976478  0.0529
    ##  Correlation: 
    ##            (Intr)
    ## CoreTemp_C -0.63 
    ## 
    ## Standardized Within-Group Residuals:
    ##         Min          Q1         Med          Q3         Max 
    ## -2.17843967 -0.54008965 -0.05469466  0.49067615  2.52872552 
    ## 
    ## Number of Observations: 69
    ## Number of Groups: 10

``` r
CAZyme_temp.model.sum = data.frame(intercept = summary(CAZyme_temp.model)$tTable[1],
                                   slope = summary(CAZyme_temp.model)$tTable[2],
                                   p_value = summary(CAZyme_temp.model)$tTable[10])

CAZyme_temp.plot = ggplot(data=CAZyme.sum, aes(x=CoreTemp_C, y=pop_invest)) +
  geom_point(aes(shape=FireClassification, fill=SiteID), size=2) +
  geom_abline(data = filter(CAZyme_temp.model.sum, p_value < 0.05),
              aes(intercept = intercept, slope = slope),
              linetype = 2, linewidth=1, color="black") +
  scale_shape_manual(values=FC.shape) +
  scale_fill_manual(values=site.col) +
  labs(x="Soil temperature (˚C)", y="Genes per genome") +
  publication_theme +
  guides(fill=guide_legend(override.aes=list(shape=site.shape), ncol = 2))

CAZyme_temp.plot
```

![](A4_Annotated_pathways_files/figure-gfm/unnamed-chunk-16-2.png)<!-- -->

Next lets take a look at glycoside hydrolases over temperature.

``` r
# Wilcoxon test over fire classification
GH_FC.wilcox = wilcox.test(x=filter(GH.sum, FireClassification=="Reference")$pop_invest,
                               y=filter(GH.sum, FireClassification=="FireAffected")$pop_invest,
                               conf.int=TRUE, conf.level=0.95)
GH_FC.wilcox
```

    ## 
    ##  Wilcoxon rank sum exact test
    ## 
    ## data:  filter(GH.sum, FireClassification == "Reference")$pop_invest and filter(GH.sum, FireClassification == "FireAffected")$pop_invest
    ## W = 192, p-value = 4.102e-05
    ## alternative hypothesis: true location shift is not equal to 0
    ## 95 percent confidence interval:
    ##  -12.466814  -5.044215
    ## sample estimates:
    ## difference in location 
    ##              -9.439772

``` r
GH_FC.plot = ggplot(data=GH.sum, aes(x=FireClassification, y=pop_invest)) +
  geom_boxplot(outlier.shape=NA) +
  geom_jitter(aes(fill=SiteID, shape=FireClassification), size=2, width=0.25) +
  annotate("text", label="< 0.001", fontface="bold", x=1.5, y=max(GH.sum$pop_invest), size=6*5/14) +
  #lims(y=c(NA, 22)) +
  scale_fill_manual(values=site.col) +
  scale_shape_manual(values=FC.shape) +
  labs(x="FireClassification", y="GH genes per genome") +
  publication_theme +
  guides(fill=guide_legend(override.aes=list(shape=site.shape), ncol=2))
GH_FC.plot
```

![](A4_Annotated_pathways_files/figure-gfm/unnamed-chunk-17-1.png)<!-- -->

``` r
# Linear mixed effects over temperature
GH_temp.model = lme(pop_invest ~ CoreTemp_C, random = ~1|SiteID, data=GH.sum)
summary(GH_temp.model)
```

    ## Linear mixed-effects model fit by REML
    ##   Data: GH.sum 
    ##        AIC      BIC    logLik
    ##   402.4629 411.2817 -197.2314
    ## 
    ## Random effects:
    ##  Formula: ~1 | SiteID
    ##         (Intercept) Residual
    ## StdDev:    5.710726 3.447981
    ## 
    ## Fixed effects:  pop_invest ~ CoreTemp_C 
    ##                 Value Std.Error DF  t-value p-value
    ## (Intercept) 18.226307 2.7757982 58 6.566150  0.0000
    ## CoreTemp_C   0.093693 0.0892795 58 1.049439  0.2983
    ##  Correlation: 
    ##            (Intr)
    ## CoreTemp_C -0.745
    ## 
    ## Standardized Within-Group Residuals:
    ##        Min         Q1        Med         Q3        Max 
    ## -1.8541681 -0.5941598 -0.1686238  0.5116517  2.3908216 
    ## 
    ## Number of Observations: 69
    ## Number of Groups: 10

``` r
GH_temp.model.sum = data.frame(intercept = summary(GH_temp.model)$tTable[1],
                                   slope = summary(GH_temp.model)$tTable[2],
                                   p_value = summary(GH_temp.model)$tTable[10])

GH_temp.plot = ggplot(data=GH.sum, aes(x=CoreTemp_C, y=pop_invest)) +
  geom_point(aes(shape=FireClassification, fill=SiteID), size=2) +
  geom_abline(data = filter(GH_temp.model.sum, p_value < 0.05),
              aes(intercept = intercept, slope = slope),
              linetype = 2, linewidth=1, color="black") +
  scale_shape_manual(values=FC.shape) +
  scale_fill_manual(values=site.col) +
  labs(x="Soil temperature (˚C)", y="Genes per genome") +
  publication_theme +
  guides(fill=guide_legend(override.aes=list(shape=site.shape), ncol = 2))

GH_temp.plot
```

![](A4_Annotated_pathways_files/figure-gfm/unnamed-chunk-17-2.png)<!-- -->

Plot boxplot together for publication

``` r
CAZyme_GH_FC.plot = cowplot::plot_grid(cowplot::plot_grid(CAZyme_FC.plot + theme(legend.position = "none"), 
                                                          CAZyme_temp.plot + theme(legend.position = "none"),
                                                          nrow=1, rel_widths = c(0.5,1)),
                                       cowplot::plot_grid(GH_FC.plot + theme(legend.position = "none"), 
                                                          GH_temp.plot + theme(legend.position = "none"),
                                                          nrow=1, rel_widths = c(0.5,1)),
                                       g_legend(CAZyme_FC.plot + theme(legend.position = "bottom", legend.direction = "vertical") +
                                                  guides(fill=guide_legend(override.aes=list(shape=site.shape), nrow = 2))),
                                       nrow=3, rel_heights = c(1,1,0.3))
CAZyme_GH_FC.plot
```

![](A4_Annotated_pathways_files/figure-gfm/unnamed-chunk-18-1.png)<!-- -->

``` r
# Save plot for publications
ggsave(CAZyme_GH_FC.plot, file="/Users/sambarnett/Documents/Shade_lab/Centralia_project/Metagenomics/Manuscript/Figures/FigS12.tiff",
       device="tiff", width=5, height=5, units="in", bg = "white")
```

# Environmental response genes

Finally, lets take a look at genes related to environmental response.
Specifically, lets look at transcription factors. We will do this
similarly to above. The list of transcriptions factors used here are
taken from KEGG. We will be looking at transcription factors in
individual families.

## Find transcription factor genes

First we need to find the transcription factor genes and count them.

``` r
# Get the list of transcription factor ko
TFac.def = read_xlsx("/Users/sambarnett/Desktop/Transcription_factor_kos.xlsx") %>%
  rename(KEGG_ortho_kofamscan = ko, ER_path = family) %>%
  select(KEGG_ortho_kofamscan, ER_path) %>%
  mutate(ER_path = gsub(" transcriptional regulator", "", ER_path))

# Now find the transcription factor genes in our samples and summarize them within the families
KEGG_TFac.sum = KEGG_long.df %>%
  inner_join(TFac.def) %>%
  dplyr::group_by(SequenceID, ER_path) %>%
  dplyr::summarize(n_genes = n()) %>%
  ungroup %>%
  tidyr::spread(key="ER_path", value="n_genes") %>%
  tidyr::gather(key="ER_path", value="n_genes", -SequenceID) %>%
  mutate(n_genes = ifelse(is.na(n_genes), 0, n_genes)) %>%
  left_join(KEGG_SCG.sum) %>%
  mutate(pop_invest = n_genes/mean_RP_genes) %>%
  left_join(sample.meta)
```

## Transcription factor investment over soil temperature

Now lets see if the number of transcription factor genes-per-genome
varies across temperature for each of the different transcription factor
families

``` r
# Linear mixed effects over temperature
TFac_temp.model.sum = data.frame()
for(N_p in unique(TFac.def$ER_path)){
  if (nrow(filter(KEGG_TFac.sum, ER_path == N_p)) > 0){
    sub.TFac_temp.model = lme(pop_invest ~ CoreTemp_C, random = ~1|SiteID, data=filter(KEGG_TFac.sum, ER_path == N_p))
    TFac_temp.model.sum = rbind(TFac_temp.model.sum,
                                data.frame(intercept = summary(sub.TFac_temp.model)$tTable[1],
                                           slope = summary(sub.TFac_temp.model)$tTable[2],
                                           p_value = summary(sub.TFac_temp.model)$tTable[10],
                                           ER_path = N_p))
  } else{
    print(paste("No", N_p, "found"))
  }
}
```

    ## [1] "No SgrR family found"
    ## [1] "No HTH-type / antitoxin PezA found"
    ## [1] "No MucR family found"

``` r
TFac_temp.model.sum = TFac_temp.model.sum %>%
  mutate(padj = p.adjust(p_value, method = "BH")) %>%
  mutate(sig = ifelse(padj < 0.05, "p < 0.05", "p ≥ 0.05"))

sig_TFac_temp = filter(TFac_temp.model.sum, padj < 0.05)$ER_path
  

TFac_temp.plot = ggplot(data=filter(KEGG_TFac.sum, ER_path %in% sig_TFac_temp), aes(x=CoreTemp_C, y=pop_invest)) +
  geom_point(aes(shape=FireClassification, fill=SiteID), size=2) +
  geom_abline(data = filter(TFac_temp.model.sum, padj < 0.05),
              aes(intercept = intercept, slope = slope),
              linetype = 2, linewidth=1, color="black") +
  scale_shape_manual(values=FC.shape) +
  scale_fill_manual(values=site.col) +
  labs(x="Soil temperature (˚C)", y="Genes per Genome") +
  publication_theme +
  facet_wrap(~ER_path, scales = "free_y", nrow=4) +
  theme(legend.position = "bottom",
        legend.direction = "vertical") +
  guides(fill=guide_legend(override.aes=list(shape=site.shape), nrow = 2))

TFac_temp.plot
```

![](A4_Annotated_pathways_files/figure-gfm/unnamed-chunk-20-1.png)<!-- -->

``` r
ggsave(TFac_temp.plot, file="/Users/sambarnett/Documents/Shade_lab/Centralia_project/Metagenomics/Manuscript/Figures/Fig5.tiff",
       device="tiff", width=7, height=7, units="in", bg = "white")
```

## Transcription factor diversity

Now lets take a look at the diversity of the transcription factors found
across samples. Diversity may have something to do with the variety and
stability of the environments that these microbes live in.

``` r
# Get abundances of KEGG orthologues that are transcription factors. Abundance will be summed RPKM.
TFac.df = KEGG_long.df %>%
  inner_join(TFac.def) %>%
  group_by(SequenceID, ER_path, KEGG_ortho_kofamscan) %>%
  summarize(sum_RPKM = sum(RPKM)) %>%
  ungroup

# For each family of transcription factors, calculate alpha diversity (evenness especially) and see if there is a relationship to temperature.
TFac_diversity.df = data.frame()
TFac_evenness_temp.model.sum = data.frame()
for (family in unique(TFac.df$ER_path)){
  TFac.mat = TFac.df %>%
    filter(ER_path == family) %>%
    select(SequenceID, KEGG_ortho_kofamscan, sum_RPKM) %>%
    tidyr::spread(key=KEGG_ortho_kofamscan, value=sum_RPKM) %>%
    tibble::column_to_rownames(var="SequenceID") %>%
    as.matrix
  TFac.mat[is.na(TFac.mat)] = 0
  
  TFac.diversity = data.frame(richness = specnumber(TFac.mat),
                              shannon = vegan::diversity(TFac.mat)) %>%
    tibble::rownames_to_column(var="SequenceID") %>%
    mutate(ER_path = family,
           evenness = shannon/log(richness)) %>%
  left_join(sample.meta)
  
  if(max(TFac.diversity$richness >= 5)){
    sub.TFac_temp.model = lme(evenness ~ CoreTemp_C, random = ~1|SiteID, data=filter(TFac.diversity, !is.na(evenness) & !is.infinite(evenness)))
    TFac_evenness_temp.model.sum = rbind(TFac_evenness_temp.model.sum,
                                         data.frame(intercept = summary(sub.TFac_temp.model)$tTable[1],
                                                    slope = summary(sub.TFac_temp.model)$tTable[2],
                                                    p_value = summary(sub.TFac_temp.model)$tTable[10],
                                                    ER_path = family))
  }
  TFac_diversity.df = rbind(TFac_diversity.df, 
                            select(TFac.diversity, SequenceID, ER_path, richness, shannon, evenness))
}

# Add in metadata and adjust p-values for multiple comparisons.
TFac_diversity.df = TFac_diversity.df %>%
  left_join(sample.meta)
TFac_evenness_temp.model.sum = TFac_evenness_temp.model.sum %>%
  mutate(padj = p.adjust(p_value, method = "BH"))

# Get significant families
sig_TFac_evenness_temp = TFac_evenness_temp.model.sum %>%
  filter(padj < 0.05)

# Plot
TF_evenness.plot = ggplot(data=filter(TFac_diversity.df, ER_path %in% sig_TFac_evenness_temp$ER_path,
                                      !is.na(evenness), !is.infinite(evenness)), 
                          aes(x=CoreTemp_C, y=evenness)) +
  geom_point(aes(shape=FireClassification, fill=SiteID), size=2) +
  geom_abline(data = sig_TFac_evenness_temp,
              aes(intercept = intercept, slope = slope),
              linetype = 2, linewidth=1, color="black") +
  scale_shape_manual(values=FC.shape) +
  scale_fill_manual(values=site.col) +
  labs(x="Soil temperature (˚C)", y="Evenness") +
  publication_theme +
  facet_wrap(~ER_path, scales = "free_y", nrow=4) +
  theme(legend.position = "bottom",
        legend.direction = "vertical") +
  guides(fill=guide_legend(override.aes=list(shape=site.shape), nrow = 2))

TF_evenness.plot
```

![](A4_Annotated_pathways_files/figure-gfm/unnamed-chunk-21-1.png)<!-- -->

``` r
ggsave(TF_evenness.plot, file="/Users/sambarnett/Documents/Shade_lab/Centralia_project/Metagenomics/Manuscript/Figures/FigS13.tiff",
       device="tiff", width=5, height=5, units="in", bg = "white")
```

# Session info

``` r
sessionInfo()
```

    ## R version 4.4.1 (2024-06-14)
    ## Platform: aarch64-apple-darwin20
    ## Running under: macOS Ventura 13.0.1
    ## 
    ## Matrix products: default
    ## BLAS:   /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/lib/libRblas.0.dylib 
    ## LAPACK: /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.0
    ## 
    ## locale:
    ## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
    ## 
    ## time zone: America/Detroit
    ## tzcode source: internal
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ##  [1] ggkegg_1.3.4    tidygraph_1.3.1 igraph_2.0.3    XML_3.99-0.17  
    ##  [5] ggraph_2.2.1    ggplot2_3.5.2   adiv_2.2.1      ecotraj_1.1.0  
    ##  [9] Rcpp_1.0.13     Nonpareil_3.5.3 picante_1.8.2   nlme_3.1-166   
    ## [13] vegan_2.6-8     lattice_0.22-6  permute_0.9-7   jsonlite_1.8.8 
    ## [17] readxl_1.4.3    ape_5.8         phyloseq_1.48.0 dplyr_1.1.4    
    ## 
    ## loaded via a namespace (and not attached):
    ##   [1] RColorBrewer_1.1-3      rstudioapi_0.16.0       magrittr_2.0.3         
    ##   [4] magick_2.8.4            farver_2.1.2            rmarkdown_2.29         
    ##   [7] ragg_1.3.2              GlobalOptions_0.1.2     adegraphics_1.0-21     
    ##  [10] zlibbioc_1.50.0         vctrs_0.6.5             multtest_2.60.0        
    ##  [13] memoise_2.0.1           base64enc_0.1-3         htmltools_0.5.8.1      
    ##  [16] progress_1.2.3          curl_5.2.2              DEoptim_2.2-8          
    ##  [19] cellranger_1.1.0        Rhdf5lib_1.26.0         rhdf5_2.48.0           
    ##  [22] KernSmooth_2.23-24      htmlwidgets_1.6.4       plyr_1.8.9             
    ##  [25] cachem_1.1.0            uuid_1.2-1              lifecycle_1.0.4        
    ##  [28] iterators_1.0.14        pkgconfig_2.0.3         Matrix_1.7-0           
    ##  [31] R6_2.5.1                fastmap_1.2.0           GenomeInfoDbData_1.2.12
    ##  [34] digest_0.6.37           numDeriv_2016.8-1.1     colorspace_2.1-1       
    ##  [37] patchwork_1.2.0         AnnotationDbi_1.66.0    S4Vectors_0.42.1       
    ##  [40] phylobase_0.8.12        textshaping_0.4.0       RSQLite_2.3.7          
    ##  [43] org.Hs.eg.db_3.19.1     labeling_0.4.3          filelock_1.0.3         
    ##  [46] clusterGeneration_1.3.8 fansi_1.0.6             httr_1.4.7             
    ##  [49] polyclip_1.10-7         mgcv_1.9-1              compiler_4.4.1         
    ##  [52] bit64_4.0.5             withr_3.0.1             doParallel_1.0.17      
    ##  [55] optimParallel_1.0-2     DBI_1.2.3               viridis_0.6.5          
    ##  [58] highr_0.11              ggforce_0.4.2           maps_3.4.2             
    ##  [61] MASS_7.3-61             rjson_0.2.22            scatterplot3d_0.3-44   
    ##  [64] biomformat_1.32.0       tools_4.4.1             rncl_0.8.7             
    ##  [67] phytools_2.3-0          glue_1.7.0              quadprog_1.5-8         
    ##  [70] rhdf5filters_1.16.0     shadowtext_0.1.4        grid_4.4.1             
    ##  [73] cluster_2.1.6           reshape2_1.4.4          ade4_1.7-22            
    ##  [76] generics_0.1.3          lpSolve_5.6.20          gtable_0.3.5           
    ##  [79] tidyr_1.3.1             data.table_1.16.0       hms_1.1.3              
    ##  [82] sp_2.1-4                xml2_1.3.6              utf8_1.2.4             
    ##  [85] XVector_0.44.0          BiocGenerics_0.50.0     ggrepel_0.9.5          
    ##  [88] foreach_1.5.2           pillar_1.9.0            stringr_1.5.1          
    ##  [91] splines_4.4.1           Kendall_2.2.1           tweenr_2.0.3           
    ##  [94] BiocFileCache_2.12.0    bit_4.0.5               survival_3.7-0         
    ##  [97] deldir_2.0-4            tidyselect_1.2.1        Biostrings_2.72.1      
    ## [100] knitr_1.48              gridExtra_2.3           IRanges_2.38.1         
    ## [103] stats4_4.4.1            xfun_0.52               graphlayouts_1.1.1     
    ## [106] expm_1.0-0              Biobase_2.64.0          stringi_1.8.4          
    ## [109] UCSC.utils_1.0.0        yaml_2.3.10             boot_1.3-31            
    ## [112] evaluate_0.24.0         codetools_0.2-20        interp_1.1-6           
    ## [115] tibble_3.2.1            cli_3.6.3               systemfonts_1.1.0      
    ## [118] munsell_0.5.1           GenomeInfoDb_1.40.1     dbplyr_2.5.0           
    ## [121] coda_0.19-4.1           png_0.1-8               parallel_4.4.1         
    ## [124] blob_1.2.4              RNeXML_2.4.11           rgl_1.3.1              
    ## [127] prettyunits_1.2.0       latticeExtra_0.6-30     jpeg_0.1-10            
    ## [130] phangorn_2.11.1         viridisLite_0.4.2       scales_1.3.0           
    ## [133] purrr_1.0.2             crayon_1.5.3            combinat_0.0-8         
    ## [136] GetoptLong_1.0.5        rlang_1.1.4             cowplot_1.1.3          
    ## [139] KEGGREST_1.44.1         fastmatch_1.1-4         mnormt_2.1.1
