Analysis of MAGs and their Functional Redundancy
================
Sam Barnett
08 August, 2025

- [Introduction](#introduction)
  - [Librarys and global variables](#librarys-and-global-variables)
  - [Metadata](#metadata)
- [Amplicon analysis](#amplicon-analysis)
  - [Get OTU data](#get-otu-data)
  - [PERMANOVA](#permanova)
  - [Ordination](#ordination)
- [Initial MAG analysis](#initial-mag-analysis)
  - [Data import and quick look](#data-import-and-quick-look)
  - [How many reads map to MAGs](#how-many-reads-map-to-mags)
  - [Beta diversity](#beta-diversity)
- [Comparing the diversity of MAGs to that of
  OTUs](#comparing-the-diversity-of-mags-to-that-of-otus)
  - [Taxonomy](#taxonomy)
  - [Beta-diversity comparisons](#beta-diversity-comparisons)
- [Functional redundancy from KEGG
  orthologues](#functional-redundancy-from-kegg-orthologues)
  - [Get data](#get-data)
  - [Calculating functional
    redundancy](#calculating-functional-redundancy)
  - [Comparing functional redundancy](#comparing-functional-redundancy)
  - [Comparing functional diversity](#comparing-functional-diversity)
  - [Comparing taxonomic diversity](#comparing-taxonomic-diversity)
  - [Plot all together](#plot-all-together)
- [Functional redundancy from KEGG
  orthologues](#functional-redundancy-from-kegg-orthologues-1)
  - [Get data](#get-data-1)
  - [Calculating functional
    redundancy](#calculating-functional-redundancy-1)
  - [Comparing functional
    redundancy](#comparing-functional-redundancy-1)
  - [Comparing functional diversity](#comparing-functional-diversity-1)
  - [Comparing taxonomic diversity](#comparing-taxonomic-diversity-1)
  - [Plot all together](#plot-all-together-1)
- [Session info](#session-info)

# Introduction

In this notebook, lets look at the metagenome assembled genomes (MAGs)
and their annotations. In this analysis we will be looking at the
functional redundancy of the system using the MAGs as our measure of
species.

## Librarys and global variables

Here are some libraries used in this analysis and the global varaibles
that will be used throughout. Mostly variables for consistent plotting.

``` r
# Libraries for data
library(dplyr)
library(phyloseq)
library(ape)
library(readxl)

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

# Total read counts
readstats.df = read.table("/Users/sambarnett/Documents/Shade_lab/Centralia_project/Metagenomics/Data/Filtered_read_statistics.txt", 
                          header=TRUE, sep="\t", comment.char = "", quote = "") %>%
  mutate(total_reads = Forward_read_count + Reverse_read_count) %>%
  select(SampleID, total_reads) %>%
  rename(SequenceID = SampleID) 

# Contig mapped Read counts
mapped_reads.df = read.table("/Users/sambarnett/Documents/Shade_lab/Centralia_project/Metagenomics/Data/Mapped_read_totals.txt", 
                          header=TRUE, sep="\t", comment.char = "", quote = "")
```

# Amplicon analysis

First run basic beta diversity analyses on the amplicon data. This is
similar to what has been already published and found in the second
analyis notebook. We are essentially seeing if community composition
varies across our samples. We will use this compare to the MAG diveristy
later.

## Get OTU data

In order to compare to previous amplicon based data (OTUs) lets load up
and process that data like we did for that study.

For this data we need to remove the RNA sample set, the sample set from
2014 (not flash frozen and extracted with MoBio kits not the
Phenol-Chloroform method). We also want to remove any sites that don’t
have samples in less than 3 years so that we can get an analysis using
the timeseries.

``` r
# Import filtered phyloseq
DNA_RNA.physeq = readRDS(file="/Users/sambarnett/Documents/Shade_lab/Centralia_project/Multi_year_project/Data/RNA_DNA_physeq.RDS")

# Remove RNA samples
DNA.physeq = subset_samples(DNA_RNA.physeq, NucAcid_type == "DNA")
DNA.physeq = prune_taxa(taxa_sums(DNA.physeq) > 0, DNA.physeq)
DNA_RNA.physeq = NULL

# Samples just included in metagenome
meta.DNA.physeq = subset_samples(DNA.physeq, SampleID %in% sample.meta$SampleID)
meta.DNA.physeq = prune_taxa(taxa_sums(meta.DNA.physeq) > 0, meta.DNA.physeq)
DNA.physeq = NULL

## Rarefy dataset
set.seed(4242)
meta.rare.physeq = rarefy_even_depth(meta.DNA.physeq)
meta.DNA.physeq = NULL

rare_depth = mean(colSums(otu_table(meta.rare.physeq)))
print(paste("Rarifying to:", rare_depth))
```

    ## [1] "Rarifying to: 161171"

``` r
sample_data(meta.rare.physeq)$SiteID = factor(sample_data(meta.rare.physeq)$SiteID, levels=used_sites)
```

## PERMANOVA

Does fire class or sampling year explain any of the variation in
community OTU compositional differences?

``` r
# Set up the blocking design for the permanova. In this case since we are repeatedly sampling the same sites over multiple years I include SiteID as the block. This is similar to "strata" in the old version of adonis.
perm <- how(nperm = 999)
dat = data.frame(sample_data(meta.rare.physeq))
setBlocks(perm) <- with(dat, SiteID)

# Get the Bray-Curtis dissimilarity
OTU_BC.dist = vegdist(t(otu_table(meta.rare.physeq)), method="bray", binary=FALSE, diag=TRUE, upper=TRUE)

# Run adonis2 
set.seed(4242)
OTU_BC.adonis = adonis2(formula = OTU_BC.dist ~ FireClassification*as.factor(Year), 
                        permutations = perm, data = dat, by="terms")
OTU_BC.adonis
```

    ## Permutation test for adonis under reduced model
    ## Terms added sequentially (first to last)
    ## Blocks:  with(dat, SiteID) 
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## adonis2(formula = OTU_BC.dist ~ FireClassification * as.factor(Year), data = dat, permutations = perm, by = "terms")
    ##                                    Df SumOfSqs      R2       F Pr(>F)    
    ## FireClassification                  1   3.0389 0.16211 12.3303  0.001 ***
    ## as.factor(Year)                     6   1.3972 0.07453  0.9448  0.001 ***
    ## FireClassification:as.factor(Year)  6   0.7546 0.04026  0.5103  0.001 ***
    ## Residual                           55  13.5554 0.72310                   
    ## Total                              68  18.7461 1.00000                   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

## Ordination

Plot ordinations for the OTU based analysis. This will be later compared
to that from the MAGs

``` r
# Run PCoA
set.seed(4242)
OTU_BC.ord = pcoa(OTU_BC.dist)

# Get axes labels with percent variation
OTU_Xaxis = paste("PCo1 (", round(OTU_BC.ord$values[1,2]*100, digits=2), "%)", sep="")
OTU_Yaxis = paste("PCo2 (", round(OTU_BC.ord$values[2,2]*100, digits=2), "%)", sep="")

# Get PCoA points in dataframe
OTU_BC.ord.df = data.frame(OTU_BC.ord$vectors) %>%
  tibble::rownames_to_column(var="SampleID") %>%
  select(SampleID, Axis.1, Axis.2) %>%
  left_join(sample_data(meta.rare.physeq), by = "SampleID") %>%
  arrange(Year) %>%
  group_by(SiteID) %>%
  mutate(YearRank = row_number()) %>%
  ungroup %>%
  mutate(AirTemperature_C = as.numeric(AirTemperature_C)) %>%
  mutate(diff_temp = CoreTemp_C-AirTemperature_C)

## Now plot by different aesthetics to put together into one full figure.

# By site ID
OTU_BC_site.plot = ggplot(data=OTU_BC.ord.df, aes(x=Axis.1, y=Axis.2)) +
  geom_point(aes(fill=SiteID, shape=FireClassification), size=2) +
  scale_shape_manual(values=c("FireAffected" = 24, "Recovered" = 21, "Reference" = 22)) +
  scale_fill_manual(values=site.col) +
  labs(x=OTU_Xaxis, y=OTU_Yaxis) +
  present_theme +
  guides(shape = guide_legend(order = 1),
         fill = guide_legend(order = 2, override.aes=list(shape=site.shape)))

# By year
OTU_BC_time.plot = ggplot(data=OTU_BC.ord.df, aes(x=Axis.1, y=Axis.2)) +
  geom_point(aes(fill=Year, shape=FireClassification), size=2) +
  scale_shape_manual(values=c("FireAffected" = 24, "Recovered" = 21, "Reference" = 22)) +
  scale_fill_gradient(low="white", high="black") +
  labs(x=OTU_Xaxis, y=OTU_Yaxis) +
  present_theme +
  guides(shape = guide_legend(order = 1),
         fill = guide_legend(order = 2, override.aes=list(shape=22)))

# Now plot together
OTU_BC_FireClassYear.leg = g_legend(OTU_BC_time.plot + theme(legend.position = "bottom", legend.direction = "vertical") +
                                  guides(fill = guide_legend(ncol=2, override.aes=list(shape=22))))
OTU_BC_SiteID.leg = g_legend(OTU_BC_site.plot + guides(shape="none", fill = guide_legend(override.aes=list(shape=site.shape), nrow=4)))

cowplot::plot_grid(cowplot::plot_grid(OTU_BC_site.plot + theme(legend.position = "none"),
                                      OTU_BC_time.plot + theme(legend.position = "none"), nrow=1),
                   cowplot::plot_grid(OTU_BC_SiteID.leg, OTU_BC_FireClassYear.leg, nrow=1),
                   ncol=1, rel_heights = c(1,0.5))
```

![](A3_MAGs_Functional_Redundancy_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

# Initial MAG analysis

Now lets take a look at the MAGs generated for this study. Initially
lets take a look at how many there are, their quality, and their beta
diversity.

## Data import and quick look

First get the bin IDs and their qualities from checkM.

``` r
# Import table of dereplicated bins and their quality scores
centralia.metawrap_dRep_checkm.df = read.csv("/Users/sambarnett/Documents/Shade_lab/Centralia_project/Metagenomics/Data/metawrap_bins/Widb.csv",
                                     header=TRUE, comment.char = "", quote="") %>%
  select(genome, score, cluster, cluster_members, closest_cluster_member, furthest_cluster_member) %>%
  left_join(read.table("/Users/sambarnett/Documents/Shade_lab/Centralia_project/Metagenomics/Data/metawrap_bins/refined_bin_checkM.txt",
                       header=TRUE, comment.char = "", quote="")) %>%
  mutate(MIMAG_quality = ifelse(completeness >= 90 & contamination < 5, "High",
                          ifelse(completeness >= 50 & contamination < 10, "Medium",
                                 "Low"))) %>%
  rename(MAG_ID = genome) %>%
  mutate(MIMAG_quality = ifelse(is.na(MIMAG_quality), "Unknown", MIMAG_quality))

# Take a look at summaries
print(paste("There are", nrow(centralia.metawrap_dRep_checkm.df), "total dereplicated bins"))
```

    ## [1] "There are 1212 total dereplicated bins"

``` r
print(paste("There are", nrow(filter(centralia.metawrap_dRep_checkm.df, MIMAG_quality == "High")), "high quality bins"))
```

    ## [1] "There are 351 high quality bins"

``` r
print(paste("There are", nrow(filter(centralia.metawrap_dRep_checkm.df, MIMAG_quality == "Medium")), "medium quality bins"))
```

    ## [1] "There are 861 medium quality bins"

Now get their read mapping values across all samples. Since these MAGs
were dereplicated across samples we mapped reads to the dereplicated
bins. Also calculate abundance as RPKM.

``` r
# Import MAG read mapping file and calculate RPKM
centralia.MAG_coverages.df = read.table("/Users/sambarnett/Documents/Shade_lab/Centralia_project/Metagenomics/Data/metawrap_bins/MAG_contig_coverages.txt", 
                               header=TRUE, comment.char = "", quote="", sep="\t") %>%
  tidyr::gather(key="SampleID", value="Read_count", -MAG_num) %>%
  rename(SequenceID = SampleID) %>%
  group_by(SequenceID) %>%
  mutate(total_mapped_reads = sum(Read_count)) %>%
  ungroup %>%
  left_join(sample.meta) %>%
  left_join(read.table("/Users/sambarnett/Documents/Shade_lab/Centralia_project/Metagenomics/Data/metawrap_bins/MAG_contig_map.txt", 
                       header=TRUE, comment.char = "", quote="", sep="\t") %>%
              mutate(MAG_num = gsub("_.*", "", Mapped_ID)) %>%
              select(MAG_ID, MAG_num) %>%
              unique) %>%
  #left_join(centralia.metawrap_dRep_checkm.df) %>%
  inner_join(centralia.metawrap_dRep_checkm.df) %>%
  mutate(RPKM_m = Read_count/((size/1000)*(total_mapped_reads/1000000)))

# Make abundance matrix (similar to OTU table)
centralia.MAG_coverages.mat = centralia.MAG_coverages.df %>%
  select(SampleID, MAG_num, RPKM_m) %>%
  tidyr::spread(key="MAG_num", value="RPKM_m") %>%
  tibble::column_to_rownames(var="SampleID") %>%
  as.matrix

## Fill in NA's as 0 RPKM
centralia.MAG_coverages.mat[is.na(centralia.MAG_coverages.mat)] = 0

## Order the table same as the OTU table for directo comparisons later
centralia.MAG_coverages.mat = centralia.MAG_coverages.mat[sample_names(meta.rare.physeq),]
```

## How many reads map to MAGs

First lets see how many reads map to the MAGs from the samples. We will
compare this to the number of reads found in total from the samples

``` r
# Summarize the number of reads mapped
MAG_mapped.sum = centralia.MAG_coverages.df %>%
  group_by(SequenceID) %>%
  summarize(MAG_mapped_reads = sum(Read_count)) %>%
  ungroup %>%
  left_join(readstats.df) %>%
  mutate(percent_mapped = MAG_mapped_reads/total_reads*100) %>%
  left_join(sample.meta)

# Get min and max
print(paste("Minimim percent reads mapped to MAGs:", min(MAG_mapped.sum$percent_mapped)))
```

    ## [1] "Minimim percent reads mapped to MAGs: 11.2510594631662"

``` r
print(paste("Maximum percent reads mapped to MAGs:", max(MAG_mapped.sum$percent_mapped)))
```

    ## [1] "Maximum percent reads mapped to MAGs: 88.334965336685"

``` r
# How does MAG read mapping compare across fire classifications?
MAG_Mapping_FC.plot = ggplot(data=MAG_mapped.sum, aes(x=FireClassification, y=percent_mapped)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.25, height = 0, aes(fill=SiteID, shape=FireClassification)) +
  labs(x="Fire classification", y="Percent of reads mapped to MAGs (%)") +
  scale_shape_manual(values=c("FireAffected" = 24, "Reference" = 22)) +
  scale_fill_manual(values=site.col) +
  publication_theme +
  guides(fill=guide_legend(ncol=2, override.aes=list(shape=site.shape)))

# How does MAG read mapping compare over temperature?
MAG_Mapping_temp.plot = ggplot(data=MAG_mapped.sum, aes(x=CoreTemp_C, y=percent_mapped)) +
  geom_point(aes(fill=SiteID, shape=FireClassification)) +
  labs(x="Soil temperature (˚C)", y="Percent of reads mapped to MAGs (%)") +
  scale_shape_manual(values=c("FireAffected" = 24, "Reference" = 22)) +
  scale_fill_manual(values=site.col) +
  publication_theme +
  guides(fill=guide_legend(ncol=2, override.aes=list(shape=site.shape)))

# Plot together
MAG_Mapping.plot = cowplot::plot_grid(MAG_Mapping_FC.plot + theme(legend.position = "none"), 
                                      MAG_Mapping_temp.plot, rel_widths = c(0.5, 1), 
                                      nrow=1, labels=c("A", "B"), label_size = 8)
MAG_Mapping.plot
```

![](A3_MAGs_Functional_Redundancy_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

``` r
ggsave(MAG_Mapping.plot, file="/Users/sambarnett/Documents/Shade_lab/Centralia_project/Metagenomics/Manuscript/Figures/Supplemental/FigS6.tiff",
       device="tiff", width=7, height=3.5, units="in", bg = "white")
```

## Beta diversity

Now run the same beta diversity analysis you previously ran for the
OTUs. This includes PERMANOVA to see variation across fire class and
time and ordinations. We will compare this to the OTUs to see how
representative these MAGs are to OTU diversity.

First run the PERMANOVA

``` r
# Set up the blocking design for the permanova. In this case since we are repeatedly sampling the same sites over multiple years I include SiteID as the block. This is similar to "strata" in the old version of adonis.
perm <- how(nperm = 999)
dat = data.frame(sample_data(meta.rare.physeq))
setBlocks(perm) <- with(dat, SiteID)

# Get the Bray-Curtis dissimilarity
MAG_BC.dist = vegdist(centralia.MAG_coverages.mat, method="bray", binary=FALSE, diag=TRUE, upper=TRUE)

# Run adonis2 
set.seed(4242)
MAG_BC.adonis = adonis2(formula = MAG_BC.dist ~ FireClassification*as.factor(Year), 
                        permutations = perm, data = dat, by = "terms")
MAG_BC.adonis
```

    ## Permutation test for adonis under reduced model
    ## Terms added sequentially (first to last)
    ## Blocks:  with(dat, SiteID) 
    ## Permutation: free
    ## Number of permutations: 999
    ## 
    ## adonis2(formula = MAG_BC.dist ~ FireClassification * as.factor(Year), data = dat, permutations = perm, by = "terms")
    ##                                    Df SumOfSqs      R2       F Pr(>F)    
    ## FireClassification                  1   3.0127 0.14096 10.2032  0.001 ***
    ## as.factor(Year)                     6   1.3692 0.06406  0.7729  0.001 ***
    ## FireClassification:as.factor(Year)  6   0.7512 0.03515  0.4240  0.001 ***
    ## Residual                           55  16.2398 0.75983                   
    ## Total                              68  21.3729 1.00000                   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

Now make the ordinations

``` r
# Run PCoA
set.seed(4242)
MAG_BC.ord = pcoa(MAG_BC.dist)

# Get axes labels with percent variation
MAG_Xaxis = paste("PCo1 (", round(MAG_BC.ord$values[1,2]*100, digits=2), "%)", sep="")
MAG_Yaxis = paste("PCo2 (", round(MAG_BC.ord$values[2,2]*100, digits=2), "%)", sep="")

# Get PCoA points in dataframe
MAG_BC.ord.df = data.frame(MAG_BC.ord$vectors) %>%
  tibble::rownames_to_column(var="SampleID") %>%
  select(SampleID, Axis.1, Axis.2) %>%
  left_join(sample_data(meta.rare.physeq), by = "SampleID") %>%
  arrange(Year) %>%
  group_by(SiteID) %>%
  mutate(YearRank = row_number()) %>%
  ungroup %>%
  mutate(AirTemperature_C = as.numeric(AirTemperature_C)) %>%
  mutate(diff_temp = CoreTemp_C-AirTemperature_C)

## Now plot by different aesthetics to put together into one full figure.

# By site ID
MAG_BC_site.plot = ggplot(data=MAG_BC.ord.df, aes(x=Axis.1, y=Axis.2)) +
  geom_point(aes(fill=SiteID, shape=FireClassification), size=2) +
  scale_shape_manual(values=c("FireAffected" = 24, "Recovered" = 21, "Reference" = 22)) +
  scale_fill_manual(values=site.col) +
  labs(x=MAG_Xaxis, y=MAG_Yaxis) +
  present_theme +
  guides(shape = guide_legend(order = 1),
         fill = guide_legend(order = 2, override.aes=list(shape=site.shape)))

# By year
MAG_BC_time.plot = ggplot(data=MAG_BC.ord.df, aes(x=Axis.1, y=Axis.2)) +
  geom_point(aes(fill=Year, shape=FireClassification), size=2) +
  scale_shape_manual(values=c("FireAffected" = 24, "Recovered" = 21, "Reference" = 22)) +
  scale_fill_gradient(low="white", high="black") +
  labs(x=MAG_Xaxis, y=MAG_Yaxis) +
  present_theme +
  guides(shape = guide_legend(order = 1),
         fill = guide_legend(order = 2, override.aes=list(shape=22)))

# Now plot together
MAG_BC_FireClassYear.leg = g_legend(MAG_BC_time.plot + theme(legend.position = "bottom", legend.direction = "vertical") +
                                  guides(fill = guide_legend(ncol=2, override.aes=list(shape=22))))
MAG_BC_SiteID.leg = g_legend(MAG_BC_site.plot + guides(shape="none", fill = guide_legend(override.aes=list(shape=site.shape), ncol=5)))

cowplot::plot_grid(cowplot::plot_grid(MAG_BC_site.plot + theme(legend.position = "none"),
                                      MAG_BC_time.plot + theme(legend.position = "none"), nrow=1),
                   cowplot::plot_grid(MAG_BC_SiteID.leg, MAG_BC_FireClassYear.leg, nrow=1),
                   ncol=1, rel_heights = c(1,0.5))
```

![](A3_MAGs_Functional_Redundancy_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->

# Comparing the diversity of MAGs to that of OTUs

Now lets see how these MAGs compare in terms of diversity to what we see
from the OTUs. Ideally we will see the MAG diversity looking similar to
the OTUs thus we can use MAGs as a representative for community
diversity

## Taxonomy

First lets compare the taxonomic makeup of the OTUs and MAGs. We
previously did something like this with all reads and the contigs. MAG
taxonomy was determined using GTDB-TK2 while OTU taxonomy was determined
using SILVA138. In the previous publication, OTU taxonomy used an older
version of SILVA but there has been a renaming of a bunch of bacterial
taxa since then which is now included in SILVA138 and GTDB-TK so we
reran that taxonomic classifier. Still though, there might be some
taxonomy names that don’t quite match between the two databases.

``` r
# Load up the GTDB-TK data
MAG_taxonomy.df = rbind(read.table("/Users/sambarnett/Documents/Shade_lab/Centralia_project/Metagenomics/Data/metawrap_bins/gtdbtk.ar53.summary.tsv", 
                                   header=TRUE, comment.char = "", quote="", sep="\t") %>%
                          select(user_genome, classification),
                        read.table("/Users/sambarnett/Documents/Shade_lab/Centralia_project/Metagenomics/Data/metawrap_bins/gtdbtk.bac120.summary.tsv", 
                                   header=TRUE, comment.char = "", quote="", sep="\t") %>%
                          select(user_genome, classification)) %>%
  mutate(MAG_ID = paste(user_genome, ".fa", sep="")) %>%
  #inner_join(select(metawrap_dRep_checkm.df, MAG_ID, size, MIMAG_quality, N50), by = "MAG_ID") %>%
  mutate(classification = gsub(";.__", ";", gsub("d__", "", classification))) %>%
  tidyr::separate(classification, into=c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep=";")

# Get some quick summaries
print(paste("There are", nrow(filter(MAG_taxonomy.df, Domain == "Archaea")), "Archaeal MAGs"))
```

    ## [1] "There are 23 Archaeal MAGs"

``` r
print(paste("There are", nrow(filter(MAG_taxonomy.df, Domain == "Bacteria")), "Bacterial MAGs"))
```

    ## [1] "There are 1189 Bacterial MAGs"

``` r
print(paste("There are", length(unique(filter(MAG_taxonomy.df, Domain == "Archaea")$Phylum)), "Archaeal phyla represented"))
```

    ## [1] "There are 2 Archaeal phyla represented"

``` r
print(paste("There are", length(unique(filter(MAG_taxonomy.df, Domain == "Bacteria")$Phylum)), "Bacterial phyla represented"))
```

    ## [1] "There are 26 Bacterial phyla represented"

``` r
print(paste("There are", length(unique(filter(MAG_taxonomy.df, Domain == "Archaea")$Class)), "Archaeal classes represented"))
```

    ## [1] "There are 3 Archaeal classes represented"

``` r
print(paste("There are", length(unique(filter(MAG_taxonomy.df, Domain == "Bacteria")$Class)), "Bacterial classes represented"))
```

    ## [1] "There are 55 Bacterial classes represented"

``` r
# Add in the MAG abundances
centralia.MAG_coverages_Tax.df = MAG_taxonomy.df %>%
  left_join(centralia.MAG_coverages.df %>%
              select(SampleID, MAG_num, MAG_ID, RPKM_m))
```

First lets see the taxonomic breakdown of the MAGs (class level). For
this we will be converting RPKM into relative abundance.

``` r
# Get summaries of relative abundances for each taxonomic Class for the MAGs. We will cluster all classes that make up less than 10% of the population to make the chart readable.
MAG_RA.df = centralia.MAG_coverages_Tax.df %>%
  filter(RPKM_m > 0) %>%
  group_by(SampleID, Class) %>%
  summarize(PRPKM_m = sum(RPKM_m)) %>%
  ungroup %>%
  group_by(SampleID) %>%
  mutate(Total_RPKM_m = sum(PRPKM_m)) %>%
  ungroup %>%
  mutate(RA = PRPKM_m/Total_RPKM_m*100) %>%
  group_by(Class) %>%
  mutate(max_RA = max(RA)) %>%
  ungroup %>%
  mutate(Class = ifelse(max_RA < 10 | grepl("Unclassified", Class), "Less than 10% or Unclassified", Class)) %>%
  group_by(SampleID, Class) %>%
  summarize(RA = sum(RA)) %>%
  ungroup %>%
  left_join(sample.meta)

# Reorder the classes so that the "less than 10%" one is last. Then assign colors.
MAG_RA.df$Class = factor(MAG_RA.df$Class,
                          levels = c(sort(unique(filter(MAG_RA.df, Class != "Less than 10% or Unclassified")$Class)), "Less than 10% or Unclassified"))

MAG_class.col = c(paultol_colors(length(levels(MAG_RA.df$Class))-1), "#777777")
names(MAG_class.col) = levels(MAG_RA.df$Class)

# Plot!
MAG_class.plot = ggplot(data=MAG_RA.df, aes(x=as.factor(Year), y=RA)) +
  geom_bar(stat="identity", aes(fill=Class), color="black") +
  scale_fill_manual(values = MAG_class.col) +
  labs(x="Year", y="Relative abundance (%)", fill="Bacterial class", title="MAGs") +
  publication_theme +
  theme(axis.text.x = element_text(angle=90)) +
  facet_wrap(~FireClassification*SiteID, nrow=1)
MAG_class.plot
```

![](A3_MAGs_Functional_Redundancy_files/figure-gfm/unnamed-chunk-12-1.png)<!-- -->

Now do the same thing with the OTU taxonomies (at Phylum level). Note
that we are reading in a new taxonomy assignment using SILVA138 than
what was used in the phyloseq object.

``` r
# Read in new taxonomy and set labels
OTU_SILVA138_2.tax = read.table(file="/Users/sambarnett/Documents/Shade_lab/Centralia_project/Metagenomics/Data/all_samples_16S_OTU_SILVA138_2_tax.txt",
           sep="\t", header=FALSE, comment.char = "", quote = "") %>%
  select(V1, V4) %>%
  rename(OTU = V1) %>%
  mutate(Taxonomy = gsub("d:", "", gsub(",.:", ",", V4))) %>%
  select(-V4) %>%
  tidyr::separate(Taxonomy, into=c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep=",") %>%
  filter(OTU %in% rownames(data.frame(tax_table(meta.rare.physeq)))) %>%
  mutate(Domain = ifelse(Domain != "Bacteria", "Unclassified", Domain),
         Phylum = ifelse(Phylum == "Actinomycetota", "Actinobacteriota", Phylum))

# Manually rename a few Phyla to better match the GTDB-TK database. Chloroflexota class AD3 is now Phylum Dormibacterota.
OTU_phyla.df = OTU_SILVA138_2.tax %>%
  select(OTU, Domain, Phylum, Class) %>%
  mutate(Phylum = ifelse(Phylum == "Chloroflexota" & Class == "AD3", "Dormibacterota", Phylum)) %>%
  mutate(Taxa = ifelse(is.na(Phylum) | Phylum %in% c("uncultured", "uncultured bacterium", "metagenome", "Incertae Sedis"), "Unclassified", Phylum)) %>%
  mutate(Taxa = ifelse(Taxa == "Pseudomonadota", 
                       ifelse(is.na(Class) | Class %in% c("uncultured", "uncultured bacterium", "metagenome", "Incertae Sedis"), "Unclassified Pseudomonadota", Class), Taxa)) %>%
  mutate(Taxa = gsub("Candidatus ", "", Taxa))
```

Now compare the two taxonomic breakdowns

``` r
# Get summaries of relative abundances for each taxonomic phyla for the OTUs We will cluster all phyla that make up less than 10% of the population to make the chart readable.
OTU_RA.df = data.frame(otu_table(meta.rare.physeq)) %>%
  tibble::rownames_to_column(var="OTU") %>% 
  tidyr::gather(key="SampleID", value = "reads", -OTU) %>%
  filter(reads > 0) %>%
  left_join(OTU_phyla.df) %>%
  group_by(SampleID, Taxa) %>%
  summarize(Preads = sum(reads)) %>%
  ungroup %>%
  group_by(SampleID) %>%
  mutate(Total_reads = sum(Preads)) %>%
  ungroup %>%
  mutate(RA = Preads/Total_reads*100) %>%
  group_by(Taxa) %>%
  mutate(max_RA = max(RA)) %>%
  ungroup %>%
  mutate(Taxa = ifelse(max_RA < 10 | grepl("Unclassified", Taxa), "Less than 10% or Unclassified", Taxa)) %>%
  group_by(SampleID, Taxa) %>%
  summarize(RA = sum(RA)) %>%
  ungroup %>%
  mutate(Dataset = "OTU taxonomy")

# Get summaries of relative abundances for each taxonomic phyla for the MAGs We will cluster all phyla that make up less than 10% of the population to make the chart readable.
MAG_RA.df = centralia.MAG_coverages_Tax.df %>%
  filter(RPKM_m > 0) %>%
  filter(Domain == "Bacteria") %>%
  mutate(Taxa = ifelse(Phylum == "Proteobacteria", Class, Phylum)) %>%
  group_by(SampleID, Taxa) %>%
  summarize(PRPKM_m = sum(RPKM_m)) %>%
  ungroup %>%
  group_by(SampleID) %>%
  mutate(Total_RPKM_m = sum(PRPKM_m)) %>%
  ungroup %>%
  mutate(RA = PRPKM_m/Total_RPKM_m*100) %>%
  group_by(Taxa) %>%
  mutate(max_RA = max(RA)) %>%
  ungroup %>%
  mutate(Taxa = ifelse(max_RA < 10 | grepl("Unclassified", Taxa), "Less than 10% or Unclassified", Taxa)) %>%
  group_by(SampleID, Taxa) %>%
  summarize(RA = sum(RA)) %>%
  ungroup %>%
  mutate(Dataset = "MAG taxonomy")

# Combine the two datasets and rename Dormibacterota to Dormibacterota/Chloroflexota AD3
OTU_MAG_RA.df = rbind(OTU_RA.df, MAG_RA.df) %>%
  left_join(sample.meta) %>%
  mutate(Taxa = ifelse(Taxa == "Dormibacterota", "Dormibacterota/Chloroflexota AD3", Taxa))

# Set colors
Comb_taxa_list = sort(unique(OTU_MAG_RA.df$Taxa))
Comb_taxa_list = c(Comb_taxa_list[Comb_taxa_list != "Less than 10% or Unclassified"], "Less than 10% or Unclassified")
Comb_taxa.col = c(paultol_colors(length(Comb_taxa_list)-1), "#777777")
names(Comb_taxa.col) = Comb_taxa_list

OTU_MAG_RA.df$Taxa = factor(OTU_MAG_RA.df$Taxa, levels = Comb_taxa_list)


# Which abundant taxa are found in the MAGs but not in the OTUs?
print("Abundant taxa found in MAGs but not in OTUs")
```

    ## [1] "Abundant taxa found in MAGs but not in OTUs"

``` r
for (tax in Comb_taxa_list){
  if (tax %in% OTU_phyla.df$Phylum){
  } else{
    if(tax %in% OTU_phyla.df$Class){
    } else{
      print(tax)
    }
  }
}
```

    ## [1] "Desulfobacterota_B"
    ## [1] "Dormibacterota/Chloroflexota AD3"
    ## [1] "Less than 10% or Unclassified"

``` r
# Which abundant taxa are found in the OTUs but not in the MAGs?
print("Abundant taxa found in OTUs but not in MAGs")
```

    ## [1] "Abundant taxa found in OTUs but not in MAGs"

``` r
for (tax in Comb_taxa_list){
  if (tax %in% centralia.MAG_coverages_Tax.df$Phylum){
  } else{
    if(tax %in% centralia.MAG_coverages_Tax.df$Class){
    } else{
      print(tax)
    }
  }
}
```

    ## [1] "Bacillota"
    ## [1] "Dormibacterota/Chloroflexota AD3"
    ## [1] "Less than 10% or Unclassified"

Plot these taxonomic breakdowns together.

``` r
OTU_MAG_Taxa.plot = ggplot(data=OTU_MAG_RA.df, aes(x=as.factor(Year), y=RA)) +
  geom_bar(stat="identity", aes(fill=Taxa), color="black") +
  scale_fill_manual(values = Comb_taxa.col) +
  labs(x="Year", y="Relative abundance (%)", fill="Bacterial Phylum/Class") +
  publication_theme +
  theme(axis.text.x = element_text(angle=90),
        legend.position = "bottom",
        legend.direction = "vertical") +
  facet_grid(Dataset~FireClassification*SiteID) +
  guides(fill = guide_legend(nrow=3))
OTU_MAG_Taxa.plot
```

![](A3_MAGs_Functional_Redundancy_files/figure-gfm/unnamed-chunk-15-1.png)<!-- -->

``` r
ggsave(OTU_MAG_Taxa.plot, file="/Users/sambarnett/Documents/Shade_lab/Centralia_project/Metagenomics/Manuscript/Figures/Supplemental/FigS7.tiff",
       device="tiff", width=7, height=7, units="in", bg = "white")
```

## Beta-diversity comparisons

Now like we did for the overall annotations, lets see how well the
beta-diversity measures compare between the OTUs and MAGs. For this we
will use the Mantel and Procrustes analyses and also look at the
ordinations.

First run the Mantel test

``` r
# Full on mantel test between bacterial OTUs and MAGs
OTU_BC.dist.mat = as.matrix(OTU_BC.dist)
MAG_BC.dist.mat = as.matrix(MAG_BC.dist)

OTU_MAG.mantel = mantel(OTU_BC.dist.mat, MAG_BC.dist.mat, permutations = 9999)
OTU_MAG.mantel
```

    ## 
    ## Mantel statistic based on Pearson's product-moment correlation 
    ## 
    ## Call:
    ## mantel(xdis = OTU_BC.dist.mat, ydis = MAG_BC.dist.mat, permutations = 9999) 
    ## 
    ## Mantel statistic r: 0.8803 
    ##       Significance: 1e-04 
    ## 
    ## Upper quantiles of permutations (null model):
    ##    90%    95%  97.5%    99% 
    ## 0.0396 0.0519 0.0632 0.0779 
    ## Permutation: free
    ## Number of permutations: 9999

Now do the Procrustes analysis based on the PCoAs.

``` r
# Rerun the PCoA ordinations with the sampe ordered samples. Note this is important otherwise the comparison doesn't work right.
ordered_OTU_BC.ord = pcoa(as.dist(OTU_BC.dist.mat))
ordered_MAG_BC.ord = pcoa(as.dist(MAG_BC.dist.mat))

# Run the Procrustes analysis
protest(X = ordered_MAG_BC.ord$vectors, Y = ordered_OTU_BC.ord$vectors, permutations = 999, symmetric = TRUE)
```

    ## 
    ## Call:
    ## protest(X = ordered_MAG_BC.ord$vectors, Y = ordered_OTU_BC.ord$vectors,      permutations = 999, symmetric = TRUE) 
    ## 
    ## Procrustes Sum of Squares (m12 squared):        0.08699 
    ## Correlation in a symmetric Procrustes rotation: 0.9555 
    ## Significance:  0.001 
    ## 
    ## Permutation: free
    ## Number of permutations: 999

``` r
# Run the Procrustes analysis again and plot
MAG_OTU.procrustes = procrustes(X = ordered_MAG_BC.ord$vectors, Y = ordered_OTU_BC.ord$vectors, permutations = 999, symmetric = TRUE)
MAG_OTU.procrustes
```

    ## 
    ## Call:
    ## procrustes(X = ordered_MAG_BC.ord$vectors, Y = ordered_OTU_BC.ord$vectors,      symmetric = TRUE, permutations = 999) 
    ## 
    ## Procrustes sum of squares:
    ## 0.08699

``` r
plot(MAG_OTU.procrustes, kind=1)
```

![](A3_MAGs_Functional_Redundancy_files/figure-gfm/unnamed-chunk-17-1.png)<!-- -->

``` r
plot(MAG_OTU.procrustes, kind=2)
```

![](A3_MAGs_Functional_Redundancy_files/figure-gfm/unnamed-chunk-17-2.png)<!-- -->

``` r
# Save figure
tiff(filename = "/Users/sambarnett/Documents/Shade_lab/Centralia_project/Metagenomics/Manuscript/Figures/Supplemental/FigS8.tiff", 
     width = 5, height = 5, units = "in", res=300)
plot(MAG_OTU.procrustes, kind=1)
dev.off()
```

    ## quartz_off_screen 
    ##                 2

Now plot the two ordinations side by side like we did for the
annotations. We will recreate the ordination plots to include fire
class, siteID, and time all in one panel rather than in two panels above

``` r
# OTU ordination
OTU_BC.plot = ggplot(data=OTU_BC.ord.df, aes(x=Axis.1, y=Axis.2)) +
  geom_point(aes(fill=SiteID, shape=FireClassification, size=Year)) +
  scale_shape_manual(values=c("FireAffected" = 24, "Reference" = 22)) +
  scale_fill_manual(values=site.col) +
  scale_size_continuous(range=c(1,4)) +
  labs(x=OTU_Xaxis, y=OTU_Yaxis) +
  publication_theme +
  theme(legend.position = "bottom",
        legend.direction = "vertical") +
  guides(shape = guide_legend(order = 1, ncol=1),
         fill = guide_legend(order = 2, override.aes=list(shape=site.shape), ncol=2),
         size = guide_legend(order = 3, override.aes=list(shape=22), ncol=2))

# MAG ordination
MAG_BC.plot = ggplot(data=MAG_BC.ord.df, aes(x=Axis.1, y=Axis.2)) +
  geom_point(aes(fill=SiteID, shape=FireClassification, size=Year)) +
  scale_shape_manual(values=c("FireAffected" = 24, "Reference" = 22)) +
  scale_fill_manual(values=site.col) +
  scale_size_continuous(range=c(1,4)) +
  labs(x=MAG_Xaxis, y=MAG_Yaxis) +
  publication_theme +
  theme(legend.position = "bottom",
        legend.direction = "vertical") +
  guides(shape = guide_legend(order = 1, ncol=1),
         fill = guide_legend(order = 2, override.aes=list(shape=site.shape), ncol=2),
         size = guide_legend(order = 3, override.aes=list(shape=22), ncol=2))


# Plot together
MAG_OTU_BC.plot = cowplot::plot_grid(cowplot::plot_grid(OTU_BC.plot + ggtitle("Amplicon") + theme(legend.position = "none"),
                                                        MAG_BC.plot + ggtitle("Metagenome") + theme(legend.position = "none"), nrow=1),
                                       g_legend(MAG_BC.plot + theme(legend.position = "bottom", legend.direction="vertical") + 
                                                  guides(shape = guide_legend(order = 1, nrow=2),
                                                         fill = guide_legend(order = 2, override.aes=list(shape=site.shape), nrow=2),
                                                         size = guide_legend(order = 3, override.aes=list(shape=22), nrow=2))),
                                       ncol=1, rel_heights = c(1,0.3))
MAG_OTU_BC.plot
```

![](A3_MAGs_Functional_Redundancy_files/figure-gfm/unnamed-chunk-18-1.png)<!-- -->

``` r
ggsave(MAG_OTU_BC.plot, file="/Users/sambarnett/Documents/Shade_lab/Centralia_project/Metagenomics/Manuscript/Figures/Supplemental/FigS9.tiff",
       device="tiff", width=7, height=5, units="in", bg = "white")
```

# Functional redundancy from KEGG orthologues

Now that we see that MAGs seem to be representative of the community
structure, lets go ahead and quantify the functional redundancy in the
system using the MAGs found here as the species. To do this we will have
to call in the KEGG orthologue annotations we previously examined.

## Get data

First thing to do is read in those annotations and match them to the
MAGs. Note that these are big files so this might take a while.

``` r
# Get the MAG assignments for each of the contigs that are incorporated into MAGs
centralia.MAG_contig.map = read.table("/Users/sambarnett/Documents/Shade_lab/Centralia_project/Metagenomics/Data/metawrap_bins/MAG_contig_map.txt", 
                            header=TRUE, comment.char = "", quote="", sep="\t") %>%
  rename(ContigID_long = ContigID) %>%
  mutate(SequenceID = gsub("_bin.*", "", MAG_ID),
         node_number = gsub("NODE_", "", gsub("_length.*", "", ContigID_long)))

# Now read in genes recovered (ORFs actually) and, using the contig-MAG map join them to their MAG
MAG_annotations.df = read.table("/Users/sambarnett/Documents/Shade_lab/Centralia_project/Metagenomics/Data/Annotations/prokka_CDS_annotations.txt", 
                                header=TRUE, comment.char = "", quote="", sep="\t") %>%
  select(SequenceID, locus_tag, ContigID) %>%
  mutate(node_number = gsub(".*_", "", ContigID)) %>%
  inner_join(centralia.MAG_contig.map)

centralia.MAG_contig.map = NULL

# Now get the KEGG annotations and join them to the MAG genes.
MAG_KEGG.df = read.table("/Users/sambarnett/Documents/Shade_lab/Centralia_project/Metagenomics/Data/Annotations/KEGG_COG_annotations.txt",
                       header=TRUE, sep="\t", comment.char = "", quote = "") %>%
  filter(KEGG_ortho_kofamscan != "") %>%
  select(SequenceID, locus_tag, KEGG_ortho_kofamscan) %>%
  inner_join(MAG_annotations.df)

MAG_annotations.df = NULL
```

Now make a matrix of KEGG orthologues across all MAGs. Then use this to
get a measure of the dissimilarity across MAGs based on their KEGG ko
makeup. Basicially functional difference between MAGs.

``` r
# Make matrix of ko across all MAGs
MAG_KEGG.mat = MAG_KEGG.df %>%
  filter(MAG_ID %in% centralia.metawrap_dRep_checkm.df$MAG_ID) %>%
  mutate(MAG_ID_short = gsub("_.*", "", Mapped_ID)) %>%
  select(MAG_ID_short, KEGG_ortho_kofamscan) %>%
  unique %>%
  mutate(Presence = 1) %>%
  tidyr::spread(key="MAG_ID_short", value="Presence") %>%
  tibble::column_to_rownames(var="KEGG_ortho_kofamscan") %>%
  as.matrix
MAG_KEGG.mat[is.na(MAG_KEGG.mat)] = 0

# Get distance matrix across MAGs based on the presence and absence of ko
MAG_KEGG.dist = dist(t(MAG_KEGG.mat), method = "binary")
```

## Calculating functional redundancy

Now that we have MAG abundances across samples and differences in KEGG
orthologue makeup of MAGs lets calculat the functional redundancy as
defined by Ricotta et al., 2016.

``` r
# Filter MAG abundances to those MAGs with KEGG orthologue assignments
sub.centralia.MAG_coverages.mat = centralia.MAG_coverages.mat[,colnames(as.matrix(MAG_KEGG.dist))]

# Calculate functional diversity (Rao's quadradic entropy)
FD.df = QE(sub.centralia.MAG_coverages.mat, dis = MAG_KEGG.dist, formula="QE") %>%
  tibble::rownames_to_column(var="SampleID") %>%
  rename(FD = diversity)

# Calculate taxonomic diversity (Simpson's index)
TD.df = data.frame(TD = vegan::diversity(sub.centralia.MAG_coverages.mat, index="simpson")) %>%
  tibble::rownames_to_column(var="SampleID")

# Calculate functional redundancy
FR.df = inner_join(FD.df, TD.df) %>%
  #mutate(FR = TD-FD) %>%
  mutate(U = FD/TD) %>%
  mutate(FR = 1-U) %>%
  left_join(data.frame(sample_data(meta.rare.physeq)))
```

## Comparing functional redundancy

Now that you’ve calculated it, compare functional redundancy across
temp, time, and fire classification.

``` r
# Linear mixed effects over temperature
FR_temp.model = lme(FR ~ CoreTemp_C, random = ~1|SiteID, data=FR.df)
FR_temp.model.sum = summary(FR_temp.model)
FR_temp.model.sum
```

    ## Linear mixed-effects model fit by REML
    ##   Data: FR.df 
    ##         AIC       BIC   logLik
    ##   -337.9711 -329.1523 172.9855
    ## 
    ## Random effects:
    ##  Formula: ~1 | SiteID
    ##         (Intercept)   Residual
    ## StdDev:  0.02423539 0.01360915
    ## 
    ## Fixed effects:  FR ~ CoreTemp_C 
    ##                 Value   Std.Error DF   t-value p-value
    ## (Intercept) 0.2880428 0.011354561 58 25.368026  0.0000
    ## CoreTemp_C  0.0012216 0.000354938 58  3.441735  0.0011
    ##  Correlation: 
    ##            (Intr)
    ## CoreTemp_C -0.724
    ## 
    ## Standardized Within-Group Residuals:
    ##         Min          Q1         Med          Q3         Max 
    ## -2.11797622 -0.50304581  0.03312956  0.56975334  2.10739372 
    ## 
    ## Number of Observations: 69
    ## Number of Groups: 10

``` r
FR_temp.model.sig = ifelse(FR_temp.model.sum$tTable[10] < 0.05, 1, 2)

FR_temp.plot = ggplot(data=FR.df, aes(x=CoreTemp_C, y=FR)) +
  geom_point(aes(shape=FireClassification, fill=SiteID), size=2) +
  geom_abline(intercept = FR_temp.model.sum$tTable[1], 
              slope = FR_temp.model.sum$tTable[2],
              linetype = 2, linewidth=1, color="black") +
  scale_shape_manual(values=FC.shape) +
  scale_fill_manual(values=site.col) +
  labs(x="Soil temperature (˚C)", y="Functional redundancy") +
  publication_theme +
  theme(legend.position = "bottom",
        legend.direction = "vertical") +
  guides(fill=guide_legend(override.aes=list(shape=site.shape), nrow = 3))

# Wilcoxon test across fire classification
FR_FC.wilcox = wilcox.test(x=filter(FR.df, FireClassification=="Reference")$FR,
                                         y=filter(FR.df, FireClassification=="FireAffected")$FR,
                                         conf.int=TRUE, conf.level=0.95)
FR_FC.wilcox
```

    ## 
    ##  Wilcoxon rank sum exact test
    ## 
    ## data:  filter(FR.df, FireClassification == "Reference")$FR and filter(FR.df, FireClassification == "FireAffected")$FR
    ## W = 61, p-value = 1.141e-10
    ## alternative hypothesis: true location shift is not equal to 0
    ## 95 percent confidence interval:
    ##  -0.06705752 -0.04319431
    ## sample estimates:
    ## difference in location 
    ##            -0.05521131

``` r
FR_FC.wilcox.p = ifelse(FR_FC.wilcox$p.value < 0.001, "p < 0.001",
                  paste("p = ", round(FR_FC.wilcox$p.value, digits = 3), sep=""))

FR_FC.plot = ggplot(data=FR.df, aes(x=FireClassification, y=FR)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(shape=FireClassification, fill=SiteID), size=2,
              width = 0.25, height = 0) +
  annotate("text", x=1.5, y=max(FR.df$FR), size = 6*5/14,
           label=FR_FC.wilcox.p) +
  scale_shape_manual(values=FC.shape) +
  scale_fill_manual(values=site.col) +
  labs(x="Fire Classification", y="Functional redundancy") +
  publication_theme +
  guides(fill=guide_legend(override.aes=list(shape=site.shape)))

# Linear mixed effects over time
FR_time.model.df = data.frame()

for(FC in unique(FR.df$FireClassification)){
  sub_FR.model = lme(FR ~ Year, random = ~1|SiteID,
                                      data=filter(FR.df, FireClassification==FC))
  FR_time.model.df = rbind(FR_time.model.df,
                              data.frame(summary(sub_FR.model)$tTable) %>%
                                tibble::rownames_to_column(var="factor") %>%
                                mutate(FireClassification = FC))
}

FR_time.model.reg = FR_time.model.df %>%
  mutate(p_slope = ifelse(factor == "Year", p.value, 1),
         factor = ifelse(factor == "(Intercept)", "Intercept", factor)) %>%
  group_by(FireClassification) %>%
  mutate(p_slope = min(p_slope)) %>%
  ungroup %>%
  select(FireClassification, factor, Value, p_slope) %>%
  tidyr::spread(key="factor", value="Value") %>%
  mutate(sig = ifelse(p_slope < 0.05, "< 0.05", "≥ 0.05"))

FR_time.model.reg
```

    ## # A tibble: 2 × 5
    ##   FireClassification p_slope Intercept      Year sig   
    ##   <chr>                <dbl>     <dbl>     <dbl> <chr> 
    ## 1 FireAffected        0.0248     5.60  -0.00261  < 0.05
    ## 2 Reference           0.300     -0.990  0.000629 ≥ 0.05

``` r
FR_time.plot = ggplot(data=FR.df, aes(x=Year, y=FR)) +
  geom_line(aes(group=SiteID), color="black", linewidth=1) + 
  geom_line(aes(color=SiteID), linewidth=0.5) + 
  geom_point(aes(shape=FireClassification, fill=SiteID), size=2) +
  geom_abline(data=filter(FR_time.model.reg, p_slope < 0.05),
              aes(intercept = Intercept, slope = Year),
              linetype=2, linewidth=1, color="black") +
  scale_shape_manual(values=FC.shape) +
  scale_fill_manual(values=site.col) +
  scale_color_manual(values=site.col) +
  scale_linetype_manual(values=c("< 0.05" = 1, "≥ 0.05" = 2)) +
  labs(x="Year", y="Functional redundancy") +
  publication_theme +
  theme(legend.position = "bottom",
        legend.direction = "vertical") +
  guides(fill=guide_legend(override.aes=list(shape=site.shape), nrow = 2)) +
  facet_wrap(~FireClassification)

# Plot all comparisons together
## Get combined legend
FR.leg = cowplot::plot_grid(g_legend(FR_time.plot +
                                       guides(fill="none", color="none")),
                            g_legend(FR_time.plot +
                                       guides(linetype="none", shape="none",
                                              fill=guide_legend(override.aes=list(shape=site.shape), nrow = 3))),
                            ncol=1)

# Plot together
cowplot::plot_grid(FR_FC.plot + theme(legend.position = "none"), 
                   FR_temp.plot + theme(legend.position = "none"), 
                   FR_time.plot + theme(legend.position = "none"),
                   FR.leg, labels=c("A", "B", "C", ""), label_size = 8)
```

![](A3_MAGs_Functional_Redundancy_files/figure-gfm/unnamed-chunk-22-1.png)<!-- -->

## Comparing functional diversity

Now lets look at the first component of functional redundancy, compare
functional diversity across temp, time, and fire classification.

``` r
# Linear mixed effects over temperature
FD_temp.model = lme(FD ~ CoreTemp_C, random = ~1|SiteID, data=FR.df)
FD_temp.model.sum = summary(FD_temp.model)
FD_temp.model.sum
```

    ## Linear mixed-effects model fit by REML
    ##   Data: FR.df 
    ##         AIC       BIC   logLik
    ##   -319.9818 -311.1631 163.9909
    ## 
    ## Random effects:
    ##  Formula: ~1 | SiteID
    ##         (Intercept)   Residual
    ## StdDev:  0.02352864 0.01593522
    ## 
    ## Fixed effects:  FD ~ CoreTemp_C 
    ##                  Value   Std.Error DF  t-value p-value
    ## (Intercept)  0.7170471 0.012161540 58 58.96023       0
    ## CoreTemp_C  -0.0021786 0.000407179 58 -5.35057       0
    ##  Correlation: 
    ##            (Intr)
    ## CoreTemp_C -0.775
    ## 
    ## Standardized Within-Group Residuals:
    ##         Min          Q1         Med          Q3         Max 
    ## -2.58059492 -0.56482824 -0.04829687  0.73124518  2.31153077 
    ## 
    ## Number of Observations: 69
    ## Number of Groups: 10

``` r
FD_temp.model.sig = ifelse(FD_temp.model.sum$tTable[10] < 0.05, 1, 2)

FD_temp.plot = ggplot(data=FR.df, aes(x=CoreTemp_C, y=FD)) +
  geom_point(aes(shape=FireClassification, fill=SiteID), size=2) +
  geom_abline(intercept = FD_temp.model.sum$tTable[1], 
              slope = FD_temp.model.sum$tTable[2],
              linetype = 2, linewidth=1, color="black") +
  scale_shape_manual(values=FC.shape) +
  scale_fill_manual(values=site.col) +
  labs(x="Soil temperature (˚C)", y="Functional diversity (Rao's)") +
  publication_theme +
  theme(legend.position = "bottom",
        legend.direction = "vertical") +
  guides(fill=guide_legend(override.aes=list(shape=site.shape), nrow = 3))

# Wilcoxon test across fire classification
FD_FC.wilcox = wilcox.test(x=filter(FR.df, FireClassification=="Reference")$FD,
                                         y=filter(FR.df, FireClassification=="FireAffected")$FD,
                                         conf.int=TRUE, conf.level=0.95)
FD_FC.wilcox
```

    ## 
    ##  Wilcoxon rank sum exact test
    ## 
    ## data:  filter(FR.df, FireClassification == "Reference")$FD and filter(FR.df, FireClassification == "FireAffected")$FD
    ## W = 916, p-value = 1.766e-10
    ## alternative hypothesis: true location shift is not equal to 0
    ## 95 percent confidence interval:
    ##  0.04844364 0.07517236
    ## sample estimates:
    ## difference in location 
    ##              0.0593064

``` r
FD_FC.wilcox.p = ifelse(FD_FC.wilcox$p.value < 0.001, "p < 0.001",
                  paste("p = ", round(FD_FC.wilcox$p.value, digits = 3), sep=""))

FD_FC.plot = ggplot(data=FR.df, aes(x=FireClassification, y=FD)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(shape=FireClassification, fill=SiteID), size=2,
              width = 0.25, height = 0) +
  annotate("text", x=1.5, y=max(FR.df$FD), size = 6*5/14,
           label=FD_FC.wilcox.p) +
  scale_shape_manual(values=FC.shape) +
  scale_fill_manual(values=site.col) +
  labs(x="Fire Classification", y="Functional diversity (Rao's)") +
  publication_theme +
  guides(fill=guide_legend(override.aes=list(shape=site.shape)))

# Linear mixed effects over time
FD_time.model.df = data.frame()

for(FC in unique(FR.df$FireClassification)){
  sub_FD.model = lme(FD ~ Year, random = ~1|SiteID,
                                      data=filter(FR.df, FireClassification==FC))
  FD_time.model.df = rbind(FD_time.model.df,
                              data.frame(summary(sub_FD.model)$tTable) %>%
                                tibble::rownames_to_column(var="factor") %>%
                                mutate(FireClassification = FC))
}

FD_time.model.reg = FD_time.model.df %>%
  mutate(p_slope = ifelse(factor == "Year", p.value, 1),
         factor = ifelse(factor == "(Intercept)", "Intercept", factor)) %>%
  group_by(FireClassification) %>%
  mutate(p_slope = min(p_slope)) %>%
  ungroup %>%
  select(FireClassification, factor, Value, p_slope) %>%
  tidyr::spread(key="factor", value="Value") %>%
  mutate(sig = ifelse(p_slope < 0.05, "< 0.05", "≥ 0.05"))

FD_time.model.reg
```

    ## # A tibble: 2 × 5
    ##   FireClassification  p_slope Intercept     Year sig   
    ##   <chr>                 <dbl>     <dbl>    <dbl> <chr> 
    ## 1 FireAffected       0.000110   -10.2   0.00538  < 0.05
    ## 2 Reference          0.756        0.261 0.000222 ≥ 0.05

``` r
FD_time.plot = ggplot(data=FR.df, aes(x=Year, y=FD)) +
  geom_line(aes(group=SiteID), color="black", linewidth=1) + 
  geom_line(aes(color=SiteID), linewidth=0.5) + 
  geom_point(aes(shape=FireClassification, fill=SiteID), size=2) +
  geom_abline(data=filter(FD_time.model.reg, p_slope < 0.05),
              aes(intercept = Intercept, slope = Year),
              linetype=2, linewidth=1, color="black") +
  scale_shape_manual(values=FC.shape) +
  scale_fill_manual(values=site.col) +
  scale_color_manual(values=site.col) +
  labs(x="Year", y="Functional diversity (Rao's)") +
  publication_theme +
  theme(legend.position = "bottom",
        legend.direction = "vertical") +
  guides(fill=guide_legend(override.aes=list(shape=site.shape), nrow = 2)) +
  facet_wrap(~FireClassification)

# Plot all comparisons together
## Get combined legend
FD.leg = cowplot::plot_grid(g_legend(FD_time.plot +
                                       guides(fill="none", color="none")),
                            g_legend(FD_time.plot +
                                       guides(linetype="none", shape="none",
                                              fill=guide_legend(override.aes=list(shape=site.shape), nrow = 3))),
                            ncol=1)

cowplot::plot_grid(FD_FC.plot + theme(legend.position = "none"), 
                   FD_temp.plot + theme(legend.position = "none"), 
                   FD_time.plot + theme(legend.position = "none"),
                   FD.leg, labels=c("A", "B", "C", ""), label_size = 8)
```

![](A3_MAGs_Functional_Redundancy_files/figure-gfm/unnamed-chunk-23-1.png)<!-- -->

## Comparing taxonomic diversity

Now lets look at the second component of functional redundancy, compare
taxonomic diversity across temp, time, and fire classification.

``` r
# Linear mixed effects over temperature
TD_temp.model = lme(TD ~ CoreTemp_C, random = ~1|SiteID, data=FR.df)
TD_temp.model.sum = summary(TD_temp.model)
TD_temp.model.sum
```

    ## Linear mixed-effects model fit by REML
    ##   Data: FR.df 
    ##         AIC       BIC   logLik
    ##   -344.0034 -335.1846 176.0017
    ## 
    ## Random effects:
    ##  Formula: ~1 | SiteID
    ##         (Intercept)   Residual
    ## StdDev: 0.008148673 0.01483046
    ## 
    ## Fixed effects:  TD ~ CoreTemp_C 
    ##                  Value   Std.Error DF   t-value p-value
    ## (Intercept)  1.0054970 0.007699354 58 130.59499       0
    ## CoreTemp_C  -0.0013291 0.000303230 58  -4.38324       0
    ##  Correlation: 
    ##            (Intr)
    ## CoreTemp_C -0.913
    ## 
    ## Standardized Within-Group Residuals:
    ##        Min         Q1        Med         Q3        Max 
    ## -5.4804272 -0.3030917  0.1366075  0.4592193  1.5723345 
    ## 
    ## Number of Observations: 69
    ## Number of Groups: 10

``` r
TD_temp.model.sig = ifelse(TD_temp.model.sum$tTable[10] < 0.05, 1, 2)

TD_temp.plot = ggplot(data=FR.df, aes(x=CoreTemp_C, y=TD)) +
  geom_point(aes(shape=FireClassification, fill=SiteID), size=2) +
  geom_abline(intercept = TD_temp.model.sum$tTable[1], 
              slope = TD_temp.model.sum$tTable[2],
              linetype = 2, linewidth=1, color="black") +
  scale_shape_manual(values=FC.shape) +
  scale_fill_manual(values=site.col) +
  labs(x="Soil temperature (˚C)", y="Taxonomic diversity (Simpons)") +
  publication_theme +
  theme(legend.position = "bottom",
        legend.direction = "vertical") +
  guides(fill=guide_legend(override.aes=list(shape=site.shape), nrow = 3))

# Wilcoxon test across fire classification
TD_FC.wilcox = wilcox.test(x=filter(FR.df, FireClassification=="Reference")$TD,
                                         y=filter(FR.df, FireClassification=="FireAffected")$TD,
                                         conf.int=TRUE, conf.level=0.95)
TD_FC.wilcox
```

    ## 
    ##  Wilcoxon rank sum exact test
    ## 
    ## data:  filter(FR.df, FireClassification == "Reference")$TD and filter(FR.df, FireClassification == "FireAffected")$TD
    ## W = 725, p-value = 0.001546
    ## alternative hypothesis: true location shift is not equal to 0
    ## 95 percent confidence interval:
    ##  0.002023098 0.014000451
    ## sample estimates:
    ## difference in location 
    ##            0.007060666

``` r
TD_FC.wilcox.p = ifelse(TD_FC.wilcox$p.value < 0.001, "p < 0.001",
                  paste("p = ", round(TD_FC.wilcox$p.value, digits = 3), sep=""))

TD_FC.plot = ggplot(data=FR.df, aes(x=FireClassification, y=TD)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(shape=FireClassification, fill=SiteID), size=2,
              width = 0.25, height = 0) +
  annotate("text", x=1.5, y=max(FR.df$TD), size = 6*5/14,
           label=TD_FC.wilcox.p) +
  scale_shape_manual(values=FC.shape) +
  scale_fill_manual(values=site.col) +
  labs(x="Fire Classification", y="Taxonomic diversity (Simpons)") +
  publication_theme +
  guides(fill=guide_legend(override.aes=list(shape=site.shape)))

# Linear mixed effects over time
TD_time.model.df = data.frame()

for(FC in unique(FR.df$FireClassification)){
  sub_TD.model = lme(TD ~ Year, random = ~1|SiteID,
                                      data=filter(FR.df, FireClassification==FC))
  TD_time.model.df = rbind(TD_time.model.df,
                              data.frame(summary(sub_TD.model)$tTable) %>%
                                tibble::rownames_to_column(var="factor") %>%
                                mutate(FireClassification = FC))
}

TD_time.model.reg = TD_time.model.df %>%
  mutate(p_slope = ifelse(factor == "Year", p.value, 1),
         factor = ifelse(factor == "(Intercept)", "Intercept", factor)) %>%
  group_by(FireClassification) %>%
  mutate(p_slope = min(p_slope)) %>%
  ungroup %>%
  select(FireClassification, factor, Value, p_slope) %>%
  tidyr::spread(key="factor", value="Value") %>%
  mutate(sig = ifelse(p_slope < 0.05, "< 0.05", "≥ 0.05"))

TD_time.model.reg
```

    ## # A tibble: 2 × 5
    ##   FireClassification  p_slope Intercept    Year sig   
    ##   <chr>                 <dbl>     <dbl>   <dbl> <chr> 
    ## 1 FireAffected       0.000848     -7.55 0.00422 < 0.05
    ## 2 Reference          0.0853       -1.40 0.00118 ≥ 0.05

``` r
TD_time.plot = ggplot(data=FR.df, aes(x=Year, y=TD)) +
  geom_line(aes(group=SiteID), color="black", linewidth=1) + 
  geom_line(aes(color=SiteID), linewidth=0.5) + 
  geom_point(aes(shape=FireClassification, fill=SiteID), size=2) +
  geom_abline(data=filter(TD_time.model.reg, p_slope < 0.05),
              aes(intercept = Intercept, slope = Year),
              linetype=2, linewidth=1, color="black") +
  scale_shape_manual(values=FC.shape) +
  scale_fill_manual(values=site.col) +
  scale_color_manual(values=site.col) +
  labs(x="Year", y="Taxonomic diversity (Simpons)") +
  publication_theme +
  theme(legend.position = "bottom",
        legend.direction = "vertical") +
  guides(fill=guide_legend(override.aes=list(shape=site.shape), nrow = 2)) +
  facet_wrap(~FireClassification)

# Plot all comparisons together
## Get combined legend
TD.leg = cowplot::plot_grid(g_legend(TD_time.plot +
                                       guides(fill="none", color="none")),
                            g_legend(TD_time.plot +
                                       guides(linetype="none", shape="none",
                                              fill=guide_legend(override.aes=list(shape=site.shape), nrow = 3))),
                            ncol=1)

cowplot::plot_grid(TD_FC.plot + theme(legend.position = "none"), 
                   TD_temp.plot + theme(legend.position = "none"), 
                   TD_time.plot + theme(legend.position = "none"),
                   TD.leg, labels=c("A", "B", "C", ""), label_size = 8)
```

![](A3_MAGs_Functional_Redundancy_files/figure-gfm/unnamed-chunk-24-1.png)<!-- -->

## Plot all together

Now for the publication, lets plot all these figures together.

``` r
FR_TD_FD.plots = cowplot::plot_grid(cowplot::plot_grid(FR_FC.plot + labs(y="Functional\nRedundancy") + scale_y_continuous(breaks=c(0.25, 0.30, 0.35)) + 
                                                         publication_theme + theme(legend.position = "none"), 
                                      FR_temp.plot + publication_theme + theme(legend.position = "none", 
                                                                           axis.title.y = element_blank(),
                                                                           axis.text.y = element_blank()),
                                      nrow=1, align = "h", axis = "tb"),
                   cowplot::plot_grid(TD_FC.plot + labs(y="Taxonomic diversity\n(Simpson's index)") + publication_theme + theme(legend.position = "none"), 
                                      TD_temp.plot + publication_theme + theme(legend.position = "none", 
                                                                           axis.title.y = element_blank(),
                                                                           axis.text.y = element_blank()),
                                      nrow=1, align = "h", axis = "tb"),
                   cowplot::plot_grid(FD_FC.plot + labs(y="Functional diversity\n(Rao's entropy)") + publication_theme + theme(legend.position = "none"), 
                                      FD_temp.plot + publication_theme + theme(legend.position = "none", 
                                                                           axis.title.y = element_blank(),
                                                                           axis.text.y = element_blank()),
                                      nrow=1, align = "h", axis = "tb"),
                   g_legend(FR_temp.plot + publication_theme + 
                              theme(legend.direction = "vertical", legend.position = "bottom") +
                              guides(fill=guide_legend(override.aes=list(shape=site.shape), ncol = 3))),
                   labels=c("A", "B", "C", ""), label_size = 8, ncol=1, rel_heights = c(1,1,1,0.75))

FR_TD_FD.plots
```

![](A3_MAGs_Functional_Redundancy_files/figure-gfm/unnamed-chunk-25-1.png)<!-- -->

``` r
#ggsave(FR_TD_FD.plots, file="/Users/sambarnett/Documents/Shade_lab/Centralia_project/Metagenomics/Manuscript/Figures/Supplemental/Fig2.tiff",
#       device="tiff", width=3.5, height=7, units="in", bg = "white")
```

# Functional redundancy from KEGG orthologues

Because MAGs are not complete genomes we may be missing genes in some
MAGs that are present in others that may lead to the appearance of
greater functional difference than there actually is. To try to get
around this, lets double check these results using pathways that allow
for incompleteness. We used gapseq to annotate pathways.

## Get data

Load up the gapseq results. Note that parsing this file might take a
really long time. Its a bad format.

``` r
# Get list of MAGs used for gapseq
MAG_list = read.table("/Users/sambarnett/Documents/Shade_lab/Centralia_project/Metagenomics/Data/Annotations/MAG_list_gapseq.txt", sep = "\t", quote = "", comment.char = "")$V1

# Get gapseq output. This output is in a really annoying format. It is essentially just a concatinated output from separate gapseq runs. The MAGs were run in order as in the list above. Between each concatinated output will be a new set of rownames. To parse this file we have to go through the dataset and add MAG ID to each row, changing the MAG ID each time we come across a new rowname set. Basically each time the Predictions column has the entry also called "Predictions", we switch to the next MAG ID in the list. Sorry this is annoying!
MAG.path = read.table("/Users/sambarnett/Documents/Shade_lab/Centralia_project/Metagenomics/Data/Annotations/gapseq_full_output.txt",
          skip = 0, sep = "\t", quote = "", comment.char = "#")
colnames(MAG.path) = c("ID", "Name", "Prediction", "Completeness", "VagueReactions", "KeyReactions", "KeyReactionsFound", "ReactionsFound")

MAG.path$MAG_ID = "NA"

MAG.path = MAG.path %>%
  filter(Prediction %in% c("Prediction", "true"))

MAG_num = 0
for (i in 1:nrow(MAG.path)){
  if(MAG.path[i,]$Prediction == "Prediction"){
    MAG_num = MAG_num + 1
  }
  MAG.path[i,]$MAG_ID = MAG_list[MAG_num]
}

# Now filter out the repeated rowname rows and remove the ".fa" from the MAG IDs
MAG.path = MAG.path %>%
  filter(Prediction == "true") %>%
  inner_join(select(centralia.MAG_coverages.df, MAG_num, MAG_ID) %>%
              unique %>%
              mutate(MAG_ID = gsub(".fa", "", MAG_ID)))
```

Now make a matrix of pathways across all MAGs. Then use this to get a
measure of the dissimilarity across MAGs based on their gapseq pathway
makeup. Again, basically functional difference between MAGs.

``` r
MAG_path.mat = MAG.path %>%
  mutate(MAG_ID_short = MAG_num) %>%
  select(MAG_ID_short, ID) %>%
  unique %>%
  mutate(Presence = 1) %>%
  tidyr::spread(key="MAG_ID_short", value="Presence") %>%
  tibble::column_to_rownames(var="ID") %>%
  as.matrix

MAG_path.mat[is.na(MAG_path.mat)] = 0

MAG_path.dist = dist(t(MAG_path.mat), method = "binary")
MAG_path.mat[1:10, 1:10]
```

    ##                                          MAG1 MAG10 MAG100 MAG1000 MAG1001
    ## |12DICHLORETHDEG-PWY|                       0     0      0       0       0
    ## |1CMET2-PWY|                                0     0      0       0       0
    ## |2PHENDEG-PWY|                              0     0      0       0       0
    ## |3-HYDROXYPHENYLACETATE-DEGRADATION-PWY|    0     0      0       0       0
    ## |4-HYDROXYMANDELATE-DEGRADATION-PWY|        1     0      0       0       0
    ## |4TOLCARBDEG-PWY|                           0     0      0       0       0
    ## |ACETOACETATE-DEG-PWY|                      0     0      0       0       0
    ## |ADENOSYLHOMOCYSCAT-PWY|                    0     0      0       0       0
    ## |ALACAT2-PWY|                               0     0      0       0       0
    ## |ALADEG-PWY|                                0     0      0       0       0
    ##                                          MAG1002 MAG1003 MAG1004 MAG1005
    ## |12DICHLORETHDEG-PWY|                          0       0       0       0
    ## |1CMET2-PWY|                                   0       0       0       0
    ## |2PHENDEG-PWY|                                 0       0       0       0
    ## |3-HYDROXYPHENYLACETATE-DEGRADATION-PWY|       0       0       0       0
    ## |4-HYDROXYMANDELATE-DEGRADATION-PWY|           0       0       0       0
    ## |4TOLCARBDEG-PWY|                              0       0       0       0
    ## |ACETOACETATE-DEG-PWY|                         0       0       0       0
    ## |ADENOSYLHOMOCYSCAT-PWY|                       0       1       0       0
    ## |ALACAT2-PWY|                                  0       0       0       0
    ## |ALADEG-PWY|                                   0       0       0       0
    ##                                          MAG1006
    ## |12DICHLORETHDEG-PWY|                          1
    ## |1CMET2-PWY|                                   0
    ## |2PHENDEG-PWY|                                 0
    ## |3-HYDROXYPHENYLACETATE-DEGRADATION-PWY|       0
    ## |4-HYDROXYMANDELATE-DEGRADATION-PWY|           0
    ## |4TOLCARBDEG-PWY|                              0
    ## |ACETOACETATE-DEG-PWY|                         0
    ## |ADENOSYLHOMOCYSCAT-PWY|                       1
    ## |ALACAT2-PWY|                                  0
    ## |ALADEG-PWY|                                   0

## Calculating functional redundancy

Now that we have MAG abundances across samples and differences in
pathway makeup of MAGs lets calculate the functional redundancy as
defined by Ricotta et al., 2016.

``` r
# Filter MAG abundances to those MAGs with gapseq assignments
sub.centralia.MAG_coverages.mat = centralia.MAG_coverages.mat[,colnames(as.matrix(MAG_path.dist))]

# Calculate functional diversity (Rao's quadradic entropy)
path_FD.df = QE(sub.centralia.MAG_coverages.mat, dis = MAG_path.dist, formula="QE") %>%
  tibble::rownames_to_column(var="SampleID") %>%
  rename(FD = diversity)

# Calculate taxonomic diversity (Simpson's index)
path_TD.df = data.frame(TD = vegan::diversity(sub.centralia.MAG_coverages.mat, index="simpson")) %>%
  tibble::rownames_to_column(var="SampleID")

# Calculate functional redundancy
path_FR.df = inner_join(path_FD.df, path_TD.df) %>%
  mutate(U = FD/TD) %>%
  mutate(FR = 1-U) %>%
  left_join(data.frame(sample_data(meta.rare.physeq)))
```

## Comparing functional redundancy

Now that you’ve calculated it, compare functional redundancy across
temp, time, and fire classification.

``` r
# Linear mixed effects over temperature
path_FR_temp.model = lme(FR ~ CoreTemp_C, random = ~1|SiteID, data=path_FR.df)
path_FR_temp.model.sum = summary(path_FR_temp.model)
path_FR_temp.model.sum
```

    ## Linear mixed-effects model fit by REML
    ##   Data: path_FR.df 
    ##         AIC       BIC   logLik
    ##   -351.1707 -342.3519 179.5853
    ## 
    ## Random effects:
    ##  Formula: ~1 | SiteID
    ##         (Intercept)   Residual
    ## StdDev:  0.01865099 0.01262549
    ## 
    ## Fixed effects:  FR ~ CoreTemp_C 
    ##                  Value   Std.Error DF   t-value p-value
    ## (Intercept) 0.25234666 0.009637744 58 26.183167  0.0000
    ## CoreTemp_C  0.00095442 0.000322629 58  2.958249  0.0045
    ##  Correlation: 
    ##            (Intr)
    ## CoreTemp_C -0.775
    ## 
    ## Standardized Within-Group Residuals:
    ##         Min          Q1         Med          Q3         Max 
    ## -2.88677976 -0.47556104  0.04110005  0.50663625  1.83365400 
    ## 
    ## Number of Observations: 69
    ## Number of Groups: 10

``` r
path_FR_temp.model.sig = ifelse(path_FR_temp.model.sum$tTable[10] < 0.05, 1, 2)

path_FR_temp.plot = ggplot(data=path_FR.df, aes(x=CoreTemp_C, y=FR)) +
  geom_point(aes(shape=FireClassification, fill=SiteID), size=2) +
  geom_abline(intercept = path_FR_temp.model.sum$tTable[1], 
              slope = path_FR_temp.model.sum$tTable[2],
              linetype = 2, linewidth=1, color="black") +
  scale_shape_manual(values=FC.shape) +
  scale_fill_manual(values=site.col) +
  labs(x="Soil temperature (˚C)", y="Functional redundancy") +
  publication_theme +
  theme(legend.position = "bottom",
        legend.direction = "vertical") +
  guides(fill=guide_legend(override.aes=list(shape=site.shape), nrow = 3))

# Wilcoxon test across fire classification
path_FR_FC.wilcox = wilcox.test(x=filter(path_FR.df, FireClassification=="Reference")$FR,
                                         y=filter(path_FR.df, FireClassification=="FireAffected")$FR,
                                         conf.int=TRUE, conf.level=0.95)
path_FR_FC.wilcox
```

    ## 
    ##  Wilcoxon rank sum exact test
    ## 
    ## data:  filter(path_FR.df, FireClassification == "Reference")$FR and filter(path_FR.df, FireClassification == "FireAffected")$FR
    ## W = 197, p-value = 5.703e-05
    ## alternative hypothesis: true location shift is not equal to 0
    ## 95 percent confidence interval:
    ##  -0.03730002 -0.01288777
    ## sample estimates:
    ## difference in location 
    ##            -0.02532725

``` r
path_FR_FC.wilcox.p = ifelse(path_FR_FC.wilcox$p.value < 0.001, "p < 0.001",
                  paste("p = ", round(path_FR_FC.wilcox$p.value, digits = 3), sep=""))

path_FR_FC.plot = ggplot(data=path_FR.df, aes(x=FireClassification, y=FR)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(shape=FireClassification, fill=SiteID), size=2,
              width = 0.25, height = 0) +
  annotate("text", x=1.5, y=max(path_FR.df$FR), size = 6*5/14,
           label=path_FR_FC.wilcox.p) +
  scale_shape_manual(values=FC.shape) +
  scale_fill_manual(values=site.col) +
  labs(x="Fire Classification", y="Functional redundancy") +
  publication_theme +
  guides(fill=guide_legend(override.aes=list(shape=site.shape)))

# Linear mixed effects over time
path_FR_time.model.df = data.frame()

for(FC in unique(path_FR.df$FireClassification)){
  sub_FR.model = lme(FR ~ Year, random = ~1|SiteID,
                                      data=filter(path_FR.df, FireClassification==FC))
  path_FR_time.model.df = rbind(path_FR_time.model.df,
                              data.frame(summary(sub_FR.model)$tTable) %>%
                                tibble::rownames_to_column(var="factor") %>%
                                mutate(FireClassification = FC))
}

path_FR_time.model.reg = path_FR_time.model.df %>%
  mutate(p_slope = ifelse(factor == "Year", p.value, 1),
         factor = ifelse(factor == "(Intercept)", "Intercept", factor)) %>%
  group_by(FireClassification) %>%
  mutate(p_slope = min(p_slope)) %>%
  ungroup %>%
  select(FireClassification, factor, Value, p_slope) %>%
  tidyr::spread(key="factor", value="Value") %>%
  mutate(sig = ifelse(p_slope < 0.05, "< 0.05", "≥ 0.05"))

path_FR_time.model.reg
```

    ## # A tibble: 2 × 5
    ##   FireClassification p_slope Intercept      Year sig   
    ##   <chr>                <dbl>     <dbl>     <dbl> <chr> 
    ## 1 FireAffected        0.0862      4.05 -0.00187  ≥ 0.05
    ## 2 Reference           0.212      -1.38  0.000810 ≥ 0.05

``` r
path_FR_time.plot = ggplot(data=path_FR.df, aes(x=Year, y=FR)) +
  geom_line(aes(group=SiteID), color="black", linewidth=1) + 
  geom_line(aes(color=SiteID), linewidth=0.5) + 
  geom_point(aes(shape=FireClassification, fill=SiteID), size=2) +
  geom_abline(data=filter(path_FR_time.model.reg, p_slope < 0.05),
              aes(intercept = Intercept, slope = Year),
              linetype=2, linewidth=1, color="black") +
  scale_shape_manual(values=FC.shape) +
  scale_fill_manual(values=site.col) +
  scale_color_manual(values=site.col) +
  labs(x="Year", y="Functional redundancy") +
  publication_theme +
  theme(legend.position = "bottom",
        legend.direction = "vertical") +
  guides(fill=guide_legend(override.aes=list(shape=site.shape), nrow = 2)) +
  facet_wrap(~FireClassification)

# Plot all comparisons together
## Get combined legend
path_FR.leg = cowplot::plot_grid(g_legend(path_FR_time.plot +
                                       guides(fill="none", color="none")),
                            g_legend(path_FR_time.plot +
                                       guides(linetype="none", shape="none",
                                              fill=guide_legend(override.aes=list(shape=site.shape), nrow = 3))),
                            ncol=1)

cowplot::plot_grid(path_FR_FC.plot + theme(legend.position = "none"), 
                   path_FR_temp.plot + theme(legend.position = "none"), 
                   path_FR_time.plot + theme(legend.position = "none"),
                   path_FR.leg, labels=c("A", "B", "C", ""), label_size = 8)
```

![](A3_MAGs_Functional_Redundancy_files/figure-gfm/unnamed-chunk-29-1.png)<!-- -->

## Comparing functional diversity

Now lets look at the first component of functional redundancy, compare
functional diversity across temp, time, and fire classification.

``` r
# Linear mixed effects over temperature
path_FD_temp.model = lme(FD ~ CoreTemp_C, random = ~1|SiteID, data=path_FR.df)
path_FD_temp.model.sum = summary(path_FD_temp.model)
path_FD_temp.model.sum
```

    ## Linear mixed-effects model fit by REML
    ##   Data: path_FR.df 
    ##         AIC       BIC   logLik
    ##   -351.6201 -342.8014 179.8101
    ## 
    ## Random effects:
    ##  Formula: ~1 | SiteID
    ##         (Intercept)   Residual
    ## StdDev:  0.01556482 0.01289742
    ## 
    ## Fixed effects:  FD ~ CoreTemp_C 
    ##                  Value   Std.Error DF  t-value p-value
    ## (Intercept)  0.7518943 0.009030080 58 83.26553       0
    ## CoreTemp_C  -0.0019316 0.000319999 58 -6.03627       0
    ##  Correlation: 
    ##            (Intr)
    ## CoreTemp_C -0.821
    ## 
    ## Standardized Within-Group Residuals:
    ##         Min          Q1         Med          Q3         Max 
    ## -2.48423584 -0.58164040  0.04173726  0.61770437  1.93721034 
    ## 
    ## Number of Observations: 69
    ## Number of Groups: 10

``` r
path_FD_temp.model.sig = ifelse(path_FD_temp.model.sum$tTable[10] < 0.05, 1, 2)

path_FD_temp.plot = ggplot(data=path_FR.df, aes(x=CoreTemp_C, y=FD)) +
  geom_point(aes(shape=FireClassification, fill=SiteID), size=2) +
  geom_abline(intercept = path_FD_temp.model.sum$tTable[1], 
              slope = path_FD_temp.model.sum$tTable[2],
              linetype = 2, linewidth=1, color="black") +
  scale_shape_manual(values=FC.shape) +
  scale_fill_manual(values=site.col) +
  labs(x="Soil temperature (˚C)", y="Functional diversity (Rao's)") +
  publication_theme +
  theme(legend.position = "bottom",
        legend.direction = "vertical") +
  guides(fill=guide_legend(override.aes=list(shape=site.shape), nrow = 3))

# Wilcoxon test across fire classification
path_FD_FC.wilcox = wilcox.test(x=filter(path_FR.df, FireClassification=="Reference")$FD,
                                         y=filter(path_FR.df, FireClassification=="FireAffected")$FD,
                                         conf.int=TRUE, conf.level=0.95)
path_FD_FC.wilcox
```

    ## 
    ##  Wilcoxon rank sum exact test
    ## 
    ## data:  filter(path_FR.df, FireClassification == "Reference")$FD and filter(path_FR.df, FireClassification == "FireAffected")$FD
    ## W = 853, p-value = 2.436e-07
    ## alternative hypothesis: true location shift is not equal to 0
    ## 95 percent confidence interval:
    ##  0.02144829 0.04520443
    ## sample estimates:
    ## difference in location 
    ##             0.03391383

``` r
path_FD_FC.wilcox.p = ifelse(path_FD_FC.wilcox$p.value < 0.001, "p < 0.001",
                  paste("p = ", round(path_FD_FC.wilcox$p.value, digits = 3), sep=""))

path_FD_FC.plot = ggplot(data=path_FR.df, aes(x=FireClassification, y=FD)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(shape=FireClassification, fill=SiteID), size=2,
              width = 0.25, height = 0) +
  annotate("text", x=1.5, y=max(path_FR.df$FD), size = 6*5/14,
           label=path_FD_FC.wilcox.p) +
  scale_shape_manual(values=FC.shape) +
  scale_fill_manual(values=site.col) +
  labs(x="Fire Classification", y="Functional diversity (Rao's)") +
  publication_theme +
  guides(fill=guide_legend(override.aes=list(shape=site.shape)))

# Linear mixed effects over time
path_FD_time.model.df = data.frame()

for(FC in unique(path_FR.df$FireClassification)){
  sub_FD.model = lme(FD ~ Year, random = ~1|SiteID,
                                      data=filter(path_FR.df, FireClassification==FC))
  path_FD_time.model.df = rbind(path_FD_time.model.df,
                              data.frame(summary(sub_FD.model)$tTable) %>%
                                tibble::rownames_to_column(var="factor") %>%
                                mutate(FireClassification = FC))
}

path_FD_time.model.reg = path_FD_time.model.df %>%
  mutate(p_slope = ifelse(factor == "Year", p.value, 1),
         factor = ifelse(factor == "(Intercept)", "Intercept", factor)) %>%
  group_by(FireClassification) %>%
  mutate(p_slope = min(p_slope)) %>%
  ungroup %>%
  select(FireClassification, factor, Value, p_slope) %>%
  tidyr::spread(key="factor", value="Value") %>%
  mutate(sig = ifelse(p_slope < 0.05, "< 0.05", "≥ 0.05"))

path_FD_time.model.reg
```

    ## # A tibble: 2 × 5
    ##   FireClassification   p_slope Intercept      Year sig   
    ##   <chr>                  <dbl>     <dbl>     <dbl> <chr> 
    ## 1 FireAffected       0.0000364    -9.23  0.00492   < 0.05
    ## 2 Reference          0.887         0.572 0.0000783 ≥ 0.05

``` r
path_FD_time.plot = ggplot(data=path_FR.df, aes(x=Year, y=FD)) +
  geom_line(aes(group=SiteID), color="black", linewidth=1) + 
  geom_line(aes(color=SiteID), linewidth=0.5) + 
  geom_point(aes(shape=FireClassification, fill=SiteID), size=2) +
  geom_abline(data=filter(path_FD_time.model.reg, p_slope < 0.05),
              aes(intercept = Intercept, slope = Year),
              linetype=2, linewidth=1, color="black") +
  scale_shape_manual(values=FC.shape) +
  scale_fill_manual(values=site.col) +
  scale_color_manual(values=site.col) +
  labs(x="Year", y="Functional diversity (Rao's)") +
  publication_theme +
  theme(legend.position = "bottom",
        legend.direction = "vertical") +
  guides(fill=guide_legend(override.aes=list(shape=site.shape), nrow = 2)) +
  facet_wrap(~FireClassification)

# Plot all comparisons together
## Get combined legend
path_FD.leg = cowplot::plot_grid(g_legend(path_FD_time.plot +
                                       guides(fill="none", color="none")),
                            g_legend(path_FD_time.plot +
                                       guides(linetype="none", shape="none",
                                              fill=guide_legend(override.aes=list(shape=site.shape), nrow = 3))),
                            ncol=1)

cowplot::plot_grid(path_FD_FC.plot + theme(legend.position = "none"), 
                   path_FD_temp.plot + theme(legend.position = "none"), 
                   path_FD_time.plot + theme(legend.position = "none"),
                   path_FD.leg, labels=c("A", "B", "C", ""), label_size = 8)
```

![](A3_MAGs_Functional_Redundancy_files/figure-gfm/unnamed-chunk-30-1.png)<!-- -->

## Comparing taxonomic diversity

Now lets look at the second component of functional redundancy, compare
taxonomic diversity across temp, time, and fire classification.

``` r
# Linear mixed effects over temperature
path_TD_temp.model = lme(TD ~ CoreTemp_C, random = ~1|SiteID, data=path_FR.df)
path_TD_temp.model.sum = summary(path_TD_temp.model)
path_TD_temp.model.sum
```

    ## Linear mixed-effects model fit by REML
    ##   Data: path_FR.df 
    ##         AIC       BIC   logLik
    ##   -344.0034 -335.1846 176.0017
    ## 
    ## Random effects:
    ##  Formula: ~1 | SiteID
    ##         (Intercept)   Residual
    ## StdDev: 0.008148673 0.01483046
    ## 
    ## Fixed effects:  TD ~ CoreTemp_C 
    ##                  Value   Std.Error DF   t-value p-value
    ## (Intercept)  1.0054970 0.007699354 58 130.59499       0
    ## CoreTemp_C  -0.0013291 0.000303230 58  -4.38324       0
    ##  Correlation: 
    ##            (Intr)
    ## CoreTemp_C -0.913
    ## 
    ## Standardized Within-Group Residuals:
    ##        Min         Q1        Med         Q3        Max 
    ## -5.4804272 -0.3030917  0.1366075  0.4592193  1.5723345 
    ## 
    ## Number of Observations: 69
    ## Number of Groups: 10

``` r
path_TD_temp.model.sig = ifelse(path_TD_temp.model.sum$tTable[10] < 0.05, 1, 2)

path_TD_temp.plot = ggplot(data=path_FR.df, aes(x=CoreTemp_C, y=TD)) +
  geom_point(aes(shape=FireClassification, fill=SiteID), size=2) +
  geom_abline(intercept = path_TD_temp.model.sum$tTable[1], 
              slope = path_TD_temp.model.sum$tTable[2],
              linetype = 2, linewidth=1, color="black") +
  scale_shape_manual(values=FC.shape) +
  scale_fill_manual(values=site.col) +
  labs(x="Soil temperature (˚C)", y="Taxonomic diversity (Simpson's)") +
  publication_theme +
  theme(legend.position = "bottom",
        legend.direction = "vertical") +
  guides(fill=guide_legend(override.aes=list(shape=site.shape), nrow = 3))

# Wilcoxon test across fire classification
path_TD_FC.wilcox = wilcox.test(x=filter(path_FR.df, FireClassification=="Reference")$TD,
                                         y=filter(path_FR.df, FireClassification=="FireAffected")$TD,
                                         conf.int=TRUE, conf.level=0.95)
path_TD_FC.wilcox
```

    ## 
    ##  Wilcoxon rank sum exact test
    ## 
    ## data:  filter(path_FR.df, FireClassification == "Reference")$TD and filter(path_FR.df, FireClassification == "FireAffected")$TD
    ## W = 725, p-value = 0.001546
    ## alternative hypothesis: true location shift is not equal to 0
    ## 95 percent confidence interval:
    ##  0.002023098 0.014000451
    ## sample estimates:
    ## difference in location 
    ##            0.007060666

``` r
path_TD_FC.wilcox.p = ifelse(path_TD_FC.wilcox$p.value < 0.001, "p < 0.001",
                  paste("p = ", round(path_TD_FC.wilcox$p.value, digits = 3), sep=""))

path_TD_FC.plot = ggplot(data=path_FR.df, aes(x=FireClassification, y=TD)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(shape=FireClassification, fill=SiteID), size=2,
              width = 0.25, height = 0) +
  annotate("text", x=1.5, y=max(path_FR.df$TD), size = 6*5/14,
           label=path_TD_FC.wilcox.p) +
  scale_shape_manual(values=FC.shape) +
  scale_fill_manual(values=site.col) +
  labs(x="Fire Classification", y="Taxonomic diversity (Simpson's)") +
  publication_theme +
  guides(fill=guide_legend(override.aes=list(shape=site.shape)))

# Linear mixed effects over time
path_TD_time.model.df = data.frame()

for(FC in unique(path_FR.df$FireClassification)){
  sub_TD.model = lme(TD ~ Year, random = ~1|SiteID,
                                      data=filter(path_FR.df, FireClassification==FC))
  path_TD_time.model.df = rbind(path_TD_time.model.df,
                              data.frame(summary(sub_TD.model)$tTable) %>%
                                tibble::rownames_to_column(var="factor") %>%
                                mutate(FireClassification = FC))
}

path_TD_time.model.reg = path_TD_time.model.df %>%
  mutate(p_slope = ifelse(factor == "Year", p.value, 1),
         factor = ifelse(factor == "(Intercept)", "Intercept", factor)) %>%
  group_by(FireClassification) %>%
  mutate(p_slope = min(p_slope)) %>%
  ungroup %>%
  select(FireClassification, factor, Value, p_slope) %>%
  tidyr::spread(key="factor", value="Value") %>%
  mutate(sig = ifelse(p_slope < 0.05, "< 0.05", "≥ 0.05"))

path_TD_time.model.reg
```

    ## # A tibble: 2 × 5
    ##   FireClassification  p_slope Intercept    Year sig   
    ##   <chr>                 <dbl>     <dbl>   <dbl> <chr> 
    ## 1 FireAffected       0.000848     -7.55 0.00422 < 0.05
    ## 2 Reference          0.0853       -1.40 0.00118 ≥ 0.05

``` r
path_TD_time.plot = ggplot(data=path_FR.df, aes(x=Year, y=TD)) +
  geom_line(aes(group=SiteID), color="black", linewidth=1) + 
  geom_line(aes(color=SiteID), linewidth=0.5) + 
  geom_point(aes(shape=FireClassification, fill=SiteID), size=2) +
  geom_abline(data=filter(path_TD_time.model.reg, p_slope < 0.05),
              aes(intercept = Intercept, slope = Year),
              linetype=2, linewidth=1, color="black") +
  scale_shape_manual(values=FC.shape) +
  scale_fill_manual(values=site.col) +
  scale_color_manual(values=site.col) +
  labs(x="Year", y="Taxonomic diversity (Simpson's)") +
  publication_theme +
  theme(legend.position = "bottom",
        legend.direction = "vertical") +
  guides(fill=guide_legend(override.aes=list(shape=site.shape), nrow = 2)) +
  facet_wrap(~FireClassification)

# Plot all comparisons together
## Get combined legend
path_TD.leg = cowplot::plot_grid(g_legend(path_TD_time.plot +
                                       guides(fill="none", color="none")),
                            g_legend(path_TD_time.plot +
                                       guides(linetype="none", shape="none",
                                              fill=guide_legend(override.aes=list(shape=site.shape), nrow = 3))),
                            ncol=1)

cowplot::plot_grid(path_TD_FC.plot + theme(legend.position = "none"), 
                   path_TD_temp.plot + theme(legend.position = "none"), 
                   path_TD_time.plot + theme(legend.position = "none"),
                   path_TD.leg, labels=c("A", "B", "C", ""), label_size = 8)
```

![](A3_MAGs_Functional_Redundancy_files/figure-gfm/unnamed-chunk-31-1.png)<!-- -->

## Plot all together

Now for the publication, lets plot all these figures together.

``` r
path_FR_TD_FD.plots = cowplot::plot_grid(cowplot::plot_grid(path_FR_FC.plot + labs(y="Functional\nRedundancy") + scale_y_continuous(breaks=c(0.25, 0.30, 0.35)) + 
                                                         publication_theme + theme(legend.position = "none"), 
                                      path_FR_temp.plot + publication_theme + theme(legend.position = "none", 
                                                                           axis.title.y = element_blank(),
                                                                           axis.text.y = element_blank()),
                                      nrow=1, align = "h", axis = "tb"),
                   cowplot::plot_grid(path_TD_FC.plot + labs(y="Taxonomic diversity\n(Simpson's index)") + publication_theme + theme(legend.position = "none"), 
                                      path_TD_temp.plot + publication_theme + theme(legend.position = "none", 
                                                                           axis.title.y = element_blank(),
                                                                           axis.text.y = element_blank()),
                                      nrow=1, align = "h", axis = "tb"),
                   cowplot::plot_grid(path_FD_FC.plot + labs(y="Functional diversity\n(Rao's entropy)") + publication_theme + theme(legend.position = "none"), 
                                      path_FD_temp.plot + publication_theme + theme(legend.position = "none", 
                                                                           axis.title.y = element_blank(),
                                                                           axis.text.y = element_blank()),
                                      nrow=1, align = "h", axis = "tb"),
                   g_legend(path_FR_temp.plot + publication_theme + 
                              theme(legend.direction = "vertical", legend.position = "bottom") +
                              guides(fill=guide_legend(override.aes=list(shape=site.shape), ncol = 3))),
                   labels=c("A", "B", "C", ""), label_size = 8, ncol=1, rel_heights = c(1,1,1,0.75))

path_FR_TD_FD.plots
```

![](A3_MAGs_Functional_Redundancy_files/figure-gfm/unnamed-chunk-32-1.png)<!-- -->

``` r
ggsave(path_FR_TD_FD.plots, file="/Users/sambarnett/Documents/Shade_lab/Centralia_project/Metagenomics/Manuscript/Figures/Supplemental/FigS10.tiff",
       device="tiff", width=3.5, height=7, units="in", bg = "white")
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
    ## [13] vegan_2.6-8     lattice_0.22-6  permute_0.9-7   readxl_1.4.3   
    ## [17] ape_5.8         phyloseq_1.48.0 dplyr_1.1.4    
    ## 
    ## loaded via a namespace (and not attached):
    ##   [1] RColorBrewer_1.1-3      rstudioapi_0.16.0       jsonlite_1.8.8         
    ##   [4] magrittr_2.0.3          magick_2.8.4            farver_2.1.2           
    ##   [7] rmarkdown_2.29          ragg_1.3.2              GlobalOptions_0.1.2    
    ##  [10] adegraphics_1.0-21      zlibbioc_1.50.0         vctrs_0.6.5            
    ##  [13] multtest_2.60.0         memoise_2.0.1           base64enc_0.1-3        
    ##  [16] htmltools_0.5.8.1       progress_1.2.3          curl_5.2.2             
    ##  [19] DEoptim_2.2-8           cellranger_1.1.0        Rhdf5lib_1.26.0        
    ##  [22] rhdf5_2.48.0            KernSmooth_2.23-24      htmlwidgets_1.6.4      
    ##  [25] plyr_1.8.9              cachem_1.1.0            uuid_1.2-1             
    ##  [28] lifecycle_1.0.4         iterators_1.0.14        pkgconfig_2.0.3        
    ##  [31] Matrix_1.7-0            R6_2.5.1                fastmap_1.2.0          
    ##  [34] GenomeInfoDbData_1.2.12 digest_0.6.37           numDeriv_2016.8-1.1    
    ##  [37] colorspace_2.1-1        patchwork_1.2.0         AnnotationDbi_1.66.0   
    ##  [40] S4Vectors_0.42.1        phylobase_0.8.12        textshaping_0.4.0      
    ##  [43] RSQLite_2.3.7           org.Hs.eg.db_3.19.1     labeling_0.4.3         
    ##  [46] filelock_1.0.3          clusterGeneration_1.3.8 fansi_1.0.6            
    ##  [49] httr_1.4.7              polyclip_1.10-7         mgcv_1.9-1             
    ##  [52] compiler_4.4.1          bit64_4.0.5             withr_3.0.1            
    ##  [55] doParallel_1.0.17       optimParallel_1.0-2     DBI_1.2.3              
    ##  [58] viridis_0.6.5           highr_0.11              ggforce_0.4.2          
    ##  [61] maps_3.4.2              MASS_7.3-61             rjson_0.2.22           
    ##  [64] scatterplot3d_0.3-44    biomformat_1.32.0       tools_4.4.1            
    ##  [67] rncl_0.8.7              phytools_2.3-0          glue_1.7.0             
    ##  [70] quadprog_1.5-8          rhdf5filters_1.16.0     shadowtext_0.1.4       
    ##  [73] grid_4.4.1              cluster_2.1.6           reshape2_1.4.4         
    ##  [76] ade4_1.7-22             generics_0.1.3          lpSolve_5.6.20         
    ##  [79] gtable_0.3.5            tidyr_1.3.1             data.table_1.16.0      
    ##  [82] hms_1.1.3               sp_2.1-4                xml2_1.3.6             
    ##  [85] utf8_1.2.4              XVector_0.44.0          BiocGenerics_0.50.0    
    ##  [88] ggrepel_0.9.5           foreach_1.5.2           pillar_1.9.0           
    ##  [91] stringr_1.5.1           splines_4.4.1           Kendall_2.2.1          
    ##  [94] tweenr_2.0.3            BiocFileCache_2.12.0    bit_4.0.5              
    ##  [97] survival_3.7-0          deldir_2.0-4            tidyselect_1.2.1       
    ## [100] Biostrings_2.72.1       knitr_1.48              gridExtra_2.3          
    ## [103] IRanges_2.38.1          stats4_4.4.1            xfun_0.52              
    ## [106] graphlayouts_1.1.1      expm_1.0-0              Biobase_2.64.0         
    ## [109] stringi_1.8.4           UCSC.utils_1.0.0        yaml_2.3.10            
    ## [112] boot_1.3-31             evaluate_0.24.0         codetools_0.2-20       
    ## [115] interp_1.1-6            tibble_3.2.1            cli_3.6.3              
    ## [118] systemfonts_1.1.0       munsell_0.5.1           GenomeInfoDb_1.40.1    
    ## [121] dbplyr_2.5.0            coda_0.19-4.1           png_0.1-8              
    ## [124] parallel_4.4.1          blob_1.2.4              RNeXML_2.4.11          
    ## [127] rgl_1.3.1               prettyunits_1.2.0       latticeExtra_0.6-30    
    ## [130] jpeg_0.1-10             phangorn_2.11.1         viridisLite_0.4.2      
    ## [133] scales_1.3.0            purrr_1.0.2             crayon_1.5.3           
    ## [136] combinat_0.0-8          GetoptLong_1.0.5        rlang_1.1.4            
    ## [139] cowplot_1.1.3           KEGGREST_1.44.1         fastmatch_1.1-4        
    ## [142] mnormt_2.1.1
