# Genome Assembly

This section documents the end-to-end workflow used to produce chromosome-scale *Artocarpus* assemblies from PacBio HiFi (primary contig generation and polishing) and Hi-C (scaffolding and validation). Assemblies are generated in a consistent, reproducible way across all accessions using Puzzler, and then curated with a standardized post-processing suite: contaminant screening/removal (BlobToolKit/BlobTools), repeat annotation and harmonization across samples (EarlGrey + merged libraries), and assembly-level summaries for cross-accession comparison. The final outputs are additionally prepared for public deposition with NCBI-compliant sequence headers and upload packaging. 

The primary outputs of this section are assemblies.... and repeat figures. 

Assemblies & kmer histograms:

![genomes](/figures/20260409_Primary_Final_Contacts_Histos.png)

Repeat analyses:

![repeats](/figures/20260410_repeats.png))
___

## Puzzler Assembly

End-to-end assembly: [puzzler](https://github.com/merondun/puzzler) 

Here is the puzzler `samples.tsv`: 

| sample  | runtime   | container                                        | wd                                                           | hifi                                                         | hic_r1                                                       | hic_r2                                                       | num_chrs | reference                                                    | hom_cov | blob_database                                             | busco_lineage     | busco_database                                              |
| ------- | --------- | ------------------------------------------------ | ------------------------------------------------------------ | ------------------------------------------------------------ | ------------------------------------------------------------ | ------------------------------------------------------------ | -------- | ------------------------------------------------------------ | ------- | --------------------------------------------------------- | ----------------- | ----------------------------------------------------------- |
| HART001 | apptainer | /home/justin.merondun/apptainer/puzzler_v1.8.sif | /project/coffea_pangenome/Artocarpus/Assemblies/20250101_JustinAssemblies | /project/coffea_pangenome/Artocarpus/Concatenated_Reads/HART001.HiFi.fastq.gz | /project/coffea_pangenome/Artocarpus/Concatenated_Reads/HART001.HiC.R1.fastq.gz | /project/coffea_pangenome/Artocarpus/Concatenated_Reads/HART001.HiC.R2.fastq.gz | 28       | /project/coffea_pangenome/Artocarpus/Concatenated_Reads/ASM2540343.fa | 68      | /90daydata/coffea_pangenome/puzzler_trials/blob_downloads | embryophyta_odb10 | /project/coffea_pangenome/Software/Merondun/busco_downloads |
| HART027 | apptainer | /home/justin.merondun/apptainer/puzzler_v1.8.sif | /project/coffea_pangenome/Artocarpus/Assemblies/20250101_JustinAssemblies | /project/coffea_pangenome/Artocarpus/Concatenated_Reads/HART027.HiFi.fastq.gz | /project/coffea_pangenome/Artocarpus/Concatenated_Reads/HART027.HiC.R1.fastq.gz | /project/coffea_pangenome/Artocarpus/Concatenated_Reads/HART027.HiC.R2.fastq.gz | 28       | /project/coffea_pangenome/Artocarpus/Concatenated_Reads/ASM2540343.fa | 26      | /90daydata/coffea_pangenome/puzzler_trials/blob_downloads | embryophyta_odb10 | /project/coffea_pangenome/Software/Merondun/busco_downloads |
| HART058 | apptainer | /home/justin.merondun/apptainer/puzzler_v1.8.sif | /project/coffea_pangenome/Artocarpus/Assemblies/20250101_JustinAssemblies | /project/coffea_pangenome/Artocarpus/Concatenated_Reads/HART058.HiFi.fastq.gz | /project/coffea_pangenome/Artocarpus/Concatenated_Reads/HART058.HiC.R1.fastq.gz | /project/coffea_pangenome/Artocarpus/Concatenated_Reads/HART058.HiC.R2.fastq.gz | 28       | /project/coffea_pangenome/Artocarpus/Concatenated_Reads/ASM2540343.fa | 116     | /90daydata/coffea_pangenome/puzzler_trials/blob_downloads | embryophyta_odb10 | /project/coffea_pangenome/Software/Merondun/busco_downloads |
| HART060 | apptainer | /home/justin.merondun/apptainer/puzzler_v1.8.sif | /project/coffea_pangenome/Artocarpus/Assemblies/20250101_JustinAssemblies | /project/coffea_pangenome/Artocarpus/Concatenated_Reads/HART060.HiFi.fastq.gz | /project/coffea_pangenome/Artocarpus/Concatenated_Reads/HART060.HiC.R1.fastq.gz | /project/coffea_pangenome/Artocarpus/Concatenated_Reads/HART060.HiC.R2.fastq.gz | 28       | /project/coffea_pangenome/Artocarpus/Concatenated_Reads/ASM2540343.fa | 21      | /90daydata/coffea_pangenome/puzzler_trials/blob_downloads | embryophyta_odb10 | /project/coffea_pangenome/Software/Merondun/busco_downloads |
| HART061 | apptainer | /home/justin.merondun/apptainer/puzzler_v1.8.sif | /project/coffea_pangenome/Artocarpus/Assemblies/20250101_JustinAssemblies | /project/coffea_pangenome/Artocarpus/Concatenated_Reads/HART061.HiFi.fastq.gz | /project/coffea_pangenome/Artocarpus/Concatenated_Reads/HART061.HiC.R1.fastq.gz | /project/coffea_pangenome/Artocarpus/Concatenated_Reads/HART061.HiC.R2.fastq.gz | 28       | /project/coffea_pangenome/Artocarpus/Concatenated_Reads/ASM2540343.fa | 79      | /90daydata/coffea_pangenome/puzzler_trials/blob_downloads | embryophyta_odb10 | /project/coffea_pangenome/Software/Merondun/busco_downloads |
| HART062 | apptainer | /home/justin.merondun/apptainer/puzzler_v1.8.sif | /project/coffea_pangenome/Artocarpus/Assemblies/20250101_JustinAssemblies | /project/coffea_pangenome/Artocarpus/Concatenated_Reads/HART062.HiFi.fastq.gz | /project/coffea_pangenome/Artocarpus/Concatenated_Reads/HART062.HiC.R1.fastq.gz | /project/coffea_pangenome/Artocarpus/Concatenated_Reads/HART062.HiC.R2.fastq.gz | 28       | /project/coffea_pangenome/Artocarpus/Concatenated_Reads/ASM2540343.fa | 34      | /90daydata/coffea_pangenome/puzzler_trials/blob_downloads | embryophyta_odb10 | /project/coffea_pangenome/Software/Merondun/busco_downloads |
| HART063 | apptainer | /home/justin.merondun/apptainer/puzzler_v1.8.sif | /project/coffea_pangenome/Artocarpus/Assemblies/20250101_JustinAssemblies | /project/coffea_pangenome/Artocarpus/Concatenated_Reads/HART063.HiFi.fastq.gz | /project/coffea_pangenome/Artocarpus/Concatenated_Reads/HART063.HiC.R1.fastq.gz | /project/coffea_pangenome/Artocarpus/Concatenated_Reads/HART063.HiC.R2.fastq.gz | 28       | /project/coffea_pangenome/Artocarpus/Concatenated_Reads/ASM2540343.fa | 46      | /90daydata/coffea_pangenome/puzzler_trials/blob_downloads | embryophyta_odb10 | /project/coffea_pangenome/Software/Merondun/busco_downloads |
| HART067 | apptainer | /home/justin.merondun/apptainer/puzzler_v1.8.sif | /project/coffea_pangenome/Artocarpus/Assemblies/20250101_JustinAssemblies | /project/coffea_pangenome/Artocarpus/Concatenated_Reads/HART067.HiFi.fastq.gz | /project/coffea_pangenome/Artocarpus/Concatenated_Reads/HART067.HiC.R1.fastq.gz | /project/coffea_pangenome/Artocarpus/Concatenated_Reads/HART067.HiC.R2.fastq.gz | 28       | /project/coffea_pangenome/Artocarpus/Concatenated_Reads/ASM2540343.fa | 47      | /90daydata/coffea_pangenome/puzzler_trials/blob_downloads | embryophyta_odb10 | /project/coffea_pangenome/Software/Merondun/busco_downloads |
| HART068 | apptainer | /home/justin.merondun/apptainer/puzzler_v1.8.sif | /project/coffea_pangenome/Artocarpus/Assemblies/20250101_JustinAssemblies | /project/coffea_pangenome/Artocarpus/Concatenated_Reads/HART068.HiFi.fastq.gz | /project/coffea_pangenome/Artocarpus/Concatenated_Reads/HART068.HiC.R1.fastq.gz | /project/coffea_pangenome/Artocarpus/Concatenated_Reads/HART068.HiC.R2.fastq.gz | 28       | /project/coffea_pangenome/Artocarpus/Concatenated_Reads/ASM2540343.fa | 33      | /90daydata/coffea_pangenome/puzzler_trials/blob_downloads | embryophyta_odb10 | /project/coffea_pangenome/Software/Merondun/busco_downloads |
| N15_23  | apptainer | /home/justin.merondun/apptainer/puzzler_v1.8.sif | /project/coffea_pangenome/Artocarpus/Assemblies/20250101_JustinAssemblies | /project/coffea_pangenome/Artocarpus/Concatenated_Reads/N15_23.HiFi.fastq.gz | /project/coffea_pangenome/Artocarpus/Concatenated_Reads/N15_23.HiC.R1.fastq.gz | /project/coffea_pangenome/Artocarpus/Concatenated_Reads/N15_23.HiC.R2.fastq.gz | 14       | /project/coffea_pangenome/Artocarpus/Concatenated_Reads/ASM2540343.fa | 69      | /90daydata/coffea_pangenome/puzzler_trials/blob_downloads | embryophyta_odb10 | /project/coffea_pangenome/Software/Merondun/busco_downloads |
| N97_50  | apptainer | /home/justin.merondun/apptainer/puzzler_v1.8.sif | /project/coffea_pangenome/Artocarpus/Assemblies/20250101_JustinAssemblies | /project/coffea_pangenome/Artocarpus/Concatenated_Reads/N97_50.HiFi.fastq.gz | /project/coffea_pangenome/Artocarpus/Concatenated_Reads/N97_50.HiC.R1.fastq.gz | /project/coffea_pangenome/Artocarpus/Concatenated_Reads/N97_50.HiC.R2.fastq.gz | 28       | /project/coffea_pangenome/Artocarpus/Concatenated_Reads/ASM2540343.fa | 25      | /90daydata/coffea_pangenome/puzzler_trials/blob_downloads | embryophyta_odb10 | /project/coffea_pangenome/Software/Merondun/busco_downloads |

```bash
sbatch -J asm_HART001 puzzler -s HART001 -m samples.tsv --threads 48 --mem 384
```

## Remove contaminants

Blobtools identifies some scaffolds with bacteria, fungi, etc in some assemblies. I will remove all those scaffolds.

Remove with:

```bash
for file in $(ls *contaminants.txt); do 
	sample=$(basename ${file} .blob.contaminants.txt)
	#grep -v '#' ${file} | egrep -v 'no-hit|Streptophyta' | awk '{print $1}' > ${sample}.contam.scaffolds.list
	#grep -vwf ${sample}.contam.scaffolds.list ../${sample}.fa.fai | awk '{print $1}' > ${sample}.keep.scaffolds.list
	#samtools faidx ../${sample}.fa $(cat ${sample}.keep.scaffolds.list) > ../no_contaminants/${sample}.fa
	removed=$(cat ${sample}.contam.scaffolds.list | wc -l )
	seq_removed=$(grep -wf ${sample}.contam.scaffolds.list ../${sample}.fa.fai | awk '{print $2}' | datamash sum 1)
	echo "${sample} had $seq_removed bp removed across ${removed} scaffolds"
done
```

```bash
HART027 had  bp removed across 0 scaffolds
HART058 had 4457961 bp removed across 82 scaffolds
HART060 had  bp removed across 0 scaffolds
HART061 had 20535435 bp removed across 265 scaffolds
HART062 had  bp removed across 0 scaffolds
HART063 had  bp removed across 0 scaffolds
HART067 had 2043506 bp removed across 83 scaffolds
HART068 had 407849 bp removed across 14 scaffolds
N15_23 had 349922 bp removed across 11 scaffolds
N97_50 had  bp removed across 0 scaffolds
```

Some samples had almost 20mb of contaminant sequence! 

## Repeat Annotation

Annotate each genome with earlgrey v6.3.0

```bash
#!/bin/bash

#SBATCH --time=14-00:00:00   
#SBATCH --nodes=1  
#SBATCH --cpus-per-task=20
#SBATCH --mem=112Gb
#SBATCH --partition=ceres
#SBATCH --account=coffea_pangenome

t=20

#module load miniconda
#source activate earlgrey
#source activate eg6

SAMPLE=$1
WD=/project/coffea_pangenome/Artocarpus/Assemblies/20250101_JustinAssemblies/primary_asm

mkdir -p ${WD}/EarlGrey_SampleLibrary
cd ${WD}/EarlGrey_SampleLibrary

FASTA="${WD}/${SAMPLE}.fa" 

if [ ! -s ${WD}/EarlGrey_SampleLibrary/${SAMPLE}_EarlGrey/${SAMPLE}_summaryFiles/${SAMPLE}.softmasked.fasta ] && [ -s ${FASTA} ]; then
    echo -e "\e[43m~~~~ Starting repeat annotation for ${SAMPLE} ~~~~\e[0m"
    # Run earlgrey with eudicotyledons repeatmasker search time, generating soft-masked genome and run helitrons. 
    earlGrey -r eudicotyledons -d yes -e yes -t ${t} -g ${FASTA} -s ${SAMPLE} -o ${WD}/EarlGrey_SampleLibrary

else
    echo -e "\e[42m~~~~ Skipping repeat annotation for ${SAMPLE}, already exists ~~~~\e[0m"
fi 
```

And copy:

```bash
DIR=/90daydata/coffea_pangenome/scratch/repeats
REP=/project/coffea_pangenome/Artocarpus/Comparative_Paper/repeats/inhouse_genomes_only_compiled/

for SAMPLE in $(cat CompSamples.list) ; do 
cp ${DIR}/${SAMPLE}_EarlGrey/${SAMPLE}_summaryFiles/* ${REP}/
done 
```

Add accession to each output:

```bash
for SAMPLE in $(ls *.familyLevelCount.txt | sed 's/.familyLevelCount.txt//g'); do 
echo "${SAMPLE}"
awk -v s=${SAMPLE} '{OFS="\t"}{print $0, s}' ${SAMPLE}.familyLevelCount.txt > ${SAMPLE}.families.out
awk -v s=${SAMPLE} '{OFS="\t"}{print $0, s}' ${SAMPLE}.highLevelCount.txt > ${SAMPLE}.summary.out
awk -v s=${SAMPLE} '{OFS="\t"}{print $0, s}' ${SAMPLE}_divergence_summary_table.tsv > ${SAMPLE}.divergence.out
done 

mergem *families.out > Repeat_Families.txt
mergem *summary.out > Repeat_Summaries.txt
mergem *divergence.out > Divergence_Summaries.txt
```

### Repeat variation & plots

- Summarize high‑level TE composition per assembly (stacked proportions from EarlGrey / RepeatMasker).
- PGLS: `log(GSize) ~ log(Repeats)` (phylogeny‑corrected using species tree; λ estimated by ML).
- PGLS: `young_frac_LTR ~ log(GSize)` (young_frac = bp‑weighted fraction of LTR bp below chosen Kimura cutoff).

```R
setwd('/project/coffea_pangenome/Artocarpus/Comparative_Paper/repeats/')
library(tidyverse)
library(RColorBrewer)
library(ggpubr)
library(meRo) #devtools::install_github('merondun/meRo')
library(vegan)
library(broom)

# Add metadata information
md <- read_tsv('/project/coffea_pangenome/Artocarpus/Comparative_Paper/samples.txt')

##### High Level #####
high_level <- read_tsv('Repeat_Summaries.txt') 
names(high_level) <- c('Classification','Coverage','Count','Proportion','Gen','Distinct_Classifications','Accession')
non_repeat <- high_level %>% 
  group_by(Accession) %>% 
  summarize(Proportion = 1-(sum(Proportion)),
            Classification = "Non Repeat",
            Coverage=NA,Count=NA,Gen=NA,Distinct_Classifications=NA)
full_level <- rbind(high_level,non_repeat) %>% filter(!grepl('Anti|Ficus|Morus',Accession))
fl <- left_join(full_level,md)
ord <- md %>% dplyr::select(Accession,Group,`Accession Order`) %>% arrange(`Accession Order`)
grpord <- md %>% dplyr::select(Group,`Accession Order`) %>% arrange(`Accession Order`) %>% dplyr::select(Group) %>% distinct
fl$Accession <- factor(fl$Accession,levels=rev(ord$Accession))
fl$Group <- factor(fl$Group,levels=grpord$Group)
fl$Classification <- factor(fl$Classification,levels=c('Non Repeat','Unclassified','Other (Simple Repeat, Microsatellite, RNA)','DNA','Penelope','Rolling Circle','LTR','LINE','SINE'))
cols <- fl %>% dplyr::select(Classification) %>% distinct %>% mutate(Color = brewer.pal(9,'Paired'))
gcols <- md %>% dplyr::select(Group,Color,Shape) %>% distinct

# ltr labs
fl_labels <- fl %>%
  filter(Classification == "LTR") %>%
  mutate(
    label = paste0(round(Proportion * 100, 1), "%"),
    text_color = ifelse(Proportion > 0.08, "black", "black")
  )


# Plot landscape 
all <- fl %>% 
  mutate(Coverage = Coverage / 1e6) %>% 
  pivot_longer(c(Proportion)) %>%
  filter( !(name == 'Distinct_Classifications' & (Classification == 'Unclassified' | Classification == "Other (Simple Repeat, Microsatellite, RNA)")) ) %>% 
  ggplot(aes(y = Accession, x = value, fill = Classification)) +
  geom_bar(stat = 'identity', position = position_stack()) +
  # add LTR percent labels
  geom_text(
    data = fl_labels,
    aes(y = Accession, x = Proportion, label = label),
    position = position_stack(vjust = 0.5),
    color = fl_labels$text_color,
    size = 2.5
  ) +
  theme_bw() +
  facet_grid(Group ~ name, scales = 'free', space = 'free_y') +
  scale_fill_manual(values = cols$Color, breaks = cols$Classification) +
  theme(strip.text.y = element_text(angle = 0)) +
  ylab('') + xlab('Distinct Classifications') +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 3))

all
ggsave('~/symlinks/comp/figures/20260410_RepeatsHighLevelSummary.pdf',
       all,dpi=300,height=4,width=6.5)


#### 	PGLS ####
# Significance genome size ~ LTRs 
rep_df <- fl %>% filter(!grepl('Non Repeat',Classification)) %>% 
  dplyr::rename(GSize = `Genome Size (Assembly, Mb)`,
                phylo_order = `Accession Order`) %>% 
  group_by(Accession,Group,phylo_order,GSize) %>% 
  summarize(Repeats = sum(Coverage)/1e6)
rep_df

# Spearman correlation
cor_res <- cor.test(
  ~ Repeats + GSize,
  data = rep_df,
  method = "spearman"
) %>% tidy()
cor_res
# # A tibble: 1 × 5
# estimate statistic p.value method                          alternative
# <dbl>     <dbl>   <dbl> <chr>                           <chr>      
#   1    0.991      2.00       0 Spearman's rank correlation rho two.sided  

# pgls
library(caper)
pgls_in <- rep_df %>% ungroup %>% dplyr::select(Accession,Repeats,GSize) %>% mutate( Accession = gsub("_","", Accession) ) %>% as.data.frame
nwk <- read.tree('/project/coffea_pangenome/Artocarpus/Comparative_Paper/orthofinder/Liftover_N11/OrthoFinder/Results_Feb04/Species_Tree/SpeciesTree_rooted.txt')
nwk$node.label <- NULL
comp <- comparative.data(phy = nwk, data = pgls_in, names.col = Accession, vcv = TRUE, na.omit = FALSE )

pgls_model <- pgls( log(GSize) ~ log(Repeats) , data = comp, lambda = "ML" ) 
summary(pgls_model)

pgls_line <- data.frame(
  Repeats = pgls_in$Repeats,
  fitted = fitted(pgls_model)
)

profile <- pgls.profile(pgls_model) 
plot(profile)

# Call:
#   pgls(formula = log(GSize) ~ log(Repeats), data = comp, lambda = "ML")
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.18570 -0.06371  0.07070  0.14834  0.25227 
# 
# Branch length transformations:
#   
#   kappa  [Fix]  : 1.000
# lambda [ ML]  : 0.951
# lower bound : 0.000, p = 0.054709
# upper bound : 1.000, p = 0.59073
# 95.0% CI   : (NA, NA)
# delta  [Fix]  : 1.000
# 
# Coefficients:
#   Estimate Std. Error t value  Pr(>|t|)    
# (Intercept)  1.271753   0.336977   3.774  0.004389 ** 
#   log(Repeats) 0.854150   0.055247  15.461 8.672e-08 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.1539 on 9 degrees of freedom
# Multiple R-squared: 0.9637,	Adjusted R-squared: 0.9597 
# F-statistic:   239 on 1 and 9 DF,  p-value: 8.672e-08 

# add fitted values aligned to rows used by the model
pdat <- comp$data %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "Accession") %>%
  mutate(
    logG = log(GSize),
    logR = log(Repeats),
    fitted_logG = as.numeric(fitted(pgls_model))
  ) %>% 
  left_join(., 
            md %>% dplyr::select(Accession,Group,Color,Shape) %>% mutate(Accession = gsub('_','',Accession)))

# for a line, sort by x
pdat_line <- pdat %>% arrange(logR)

s <- summary(pgls_model)
beta <- s$coefficients["log(Repeats)", "Estimate"]
se   <- s$coefficients["log(Repeats)", "Std. Error"]
pval <- s$coefficients["log(Repeats)", "Pr(>|t|)"]
lam  <- pgls_model$param["lambda"]

label_pgls <- sprintf("PGLS: b = %.3f ± %.3f\np = %.2e\nl = %.2f", beta, se, pval, lam)


p_pgls <- ggplot(pdat, aes(x = logR, y = logG, fill = Group, shape = Group, label =Group)) +
  geom_point(size = 2) +
  geom_text_repel(size = 2, max.overlaps = Inf) +
  geom_line(data = pdat_line, aes(x = logR, y = fitted_logG), inherit.aes = FALSE,
            color = "blue", linewidth = 0.8) +
  theme_bw(base_size = 8) +
  scale_shape_manual(values=spcols$Shape,breaks=spcols$Group)+
  scale_fill_manual(values=spcols$Color,breaks=spcols$Group)+
  annotate("text", x = min(pdat$logR)+0.01, y = max(pdat$logG)-0.05, label = label_pgls,
           hjust = 0, vjust = 0, size = 1.2) +
  labs(
    x = "log(Repeat masked bp)",
    y = "log(Genome size (bp))"
  ) +
  theme(legend.position = "none")

p_pgls

ggsave('~/symlinks/comp/figures/20260410_RepeatsHighLevelSummary-Repeats-GenomeSize_PGLS.pdf',
       p_pgls,dpi=300,height=3,width=2.75)

###### Divergence Summaries ######
t <- read_tsv('inhouse_only/Divergence_Summaries.txt')
t <- t %>% dplyr::rename(Accession = HART001)
tm <- left_join(t,md)
tm$Accession <- factor(tm$Accession,levels=ord$Accession)
tm$Group <- factor(tm$Group,levels=rev(grpord$Group))

dom <- c("HART001","HART027")

# helper build div_traits + run pgls for a given cutoff since we run few models sensitivity
run_div_pgls <- function(young_cut_pct, tm, rep_df, nwk) {
  young_cut <- young_cut_pct / 100
  
  div_traits <- tm %>%
    filter(subclass %in% c("LTR","LINE","DNA")) %>%
    mutate(domest = if_else(as.character(Accession) %in% dom, "domesticated", "wild")) %>%
    group_by(Accession, domest, subclass) %>%
    summarise(
      te_mb = sum(total_bp, na.rm = TRUE) / 1e6,
      young_mb = sum(total_bp[mean_div <= young_cut], na.rm = TRUE) / 1e6,
      young_frac = young_mb / te_mb,
      wmean_div = weighted.mean(mean_div, w = total_bp, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    pivot_wider(
      names_from = subclass,
      values_from = c(te_mb, young_mb, young_frac, wmean_div),
      names_sep = "_"
    ) %>%
    mutate(Accession = gsub("_","", Accession)) %>%
    left_join(
      rep_df %>%
        dplyr::select(Accession, GSize, Repeats) %>%
        mutate(Accession = gsub("_","", Accession)),
      by = "Accession"
    ) %>%
    as.data.frame()
  
  comp_div <- comparative.data(
    phy = nwk,
    data = div_traits,
    names.col = "Accession",
    vcv = TRUE,
    na.omit = FALSE
  )
  
  # ensure domest is a factor
  comp_div$data$domest <- factor(comp_div$data$domest, levels = c("wild","domesticated"))
  
  m_dom <- pgls(wmean_div_LTR ~ domest, data = comp_div, lambda = "ML")
  m_div <- pgls(young_frac_LTR ~ log(GSize), data = comp_div, lambda = "ML")
  
  list(div_traits = div_traits, comp_div = comp_div, m_dom = m_dom, m_div = m_div)
}

res10 <- run_div_pgls(young_cut_pct = 10, tm = tm, rep_df = rep_df, nwk = nwk)
m_div <- res10$m_div

summary(res10$m_div)

# Call:
#   pgls(formula = young_frac_LTR ~ log(GSize), data = comp_div, 
#        lambda = "ML")
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.94304 -0.54048 -0.04824  0.37944  0.83333 
# 
# Branch length transformations:
#   
#   kappa  [Fix]  : 1.000
# lambda [ ML]  : 1.000
# lower bound : 0.000, p = 0.019737
# upper bound : 1.000, p = 1    
# 95.0% CI   : (0.393, NA)
# delta  [Fix]  : 1.000
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)  
# (Intercept) -4.11747    1.64560 -2.5021  0.03375 *
#   log(GSize)   0.66154    0.25404  2.6041  0.02855 *
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.6447 on 9 degrees of freedom
# Multiple R-squared: 0.4297,	Adjusted R-squared: 0.3663 
# F-statistic: 6.781 on 1 and 9 DF,  p-value: 0.02855 

pdat_young <- res10$comp_div$data %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "Accession") %>%
  mutate(
    logG = log(GSize),
    young = young_frac_LTR,
    fitted_young = as.numeric(fitted(m_div))
  )

pdat_young_line <- pdat_young %>% arrange(logG)

# label (β, p, λ)
s <- summary(m_div)
beta <- s$coefficients["log(GSize)", "Estimate"]
se   <- s$coefficients["log(GSize)", "Std. Error"]
pval <- s$coefficients["log(GSize)", "Pr(>|t|)"]
lam  <- m_div$param["lambda"]

label_pgls_young <- sprintf(
  "PGLS: b = %.3f ± %.3f\np = %.2e\nl = %.2f\nyoung cutoff = %d%%",
  beta, se, pval, lam, 10
)

# Plot 
p_young <- ggplot(
  pdat_young,
  aes(x = logG, y = young, fill = Group, shape = Group, label = Accession)
) +
  geom_point(size = 2) +
  geom_text_repel(size = 2, max.overlaps = Inf) +
  geom_line(
    data = pdat_young_line,
    aes(x = logG, y = fitted_young),
    inherit.aes = FALSE,
    color = "blue",
    linewidth = 0.8
  ) +
  theme_bw(base_size = 8) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 2)) +  
  scale_shape_manual(values = spcols$Shape, breaks = spcols$Group) +
  scale_fill_manual(values = spcols$Color, breaks = spcols$Group) +
  annotate(
    "text",
    x = min(pdat_young$logG, na.rm = TRUE) + 0.01,
    y = max(pdat_young$young, na.rm = TRUE) - 0.02,
    label = label_pgls_young,
    hjust = 0, vjust = 1, size = 1.2
  ) +
  labs(
    x = "log(Genome size (bp))",
    y = "LTR young fraction (bp-weighted)"
  ) +
  theme(legend.position = "none")

p_young

ggsave('~/symlinks/comp/figures/20260410_LTRYoungFrac_vs_GSize_PGLS.pdf',
       p_young, dpi=300, height=3, width=2.75)

# domesticates?
summary(res10$m_dom)

# Call:
#   pgls(formula = wmean_div_LTR ~ domest, data = comp_div, lambda = "ML")
# 
# Residuals:
#   Min        1Q    Median        3Q       Max 
# -0.109326 -0.050600  0.007363  0.046943  0.071145 
# 
# Branch length transformations:
#   
#   kappa  [Fix]  : 1.000
# lambda [ ML]  : 1.000
# lower bound : 0.000, p = 0.0017301
# upper bound : 1.000, p = 1    
# 95.0% CI   : (0.776, NA)
# delta  [Fix]  : 1.000
# 
# Coefficients:
#   Estimate Std. Error t value  Pr(>|t|)    
# (Intercept)         0.1559046  0.0086090 18.1095 2.176e-08 ***
#   domestdomesticated -0.0128897  0.0051162 -2.5194    0.0328 *  
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.06965 on 9 degrees of freedom
# Multiple R-squared: 0.4136,	Adjusted R-squared: 0.3484 
# F-statistic: 6.347 on 1 and 9 DF,  p-value: 0.0328

## Sensitivity loop over young_cut_pct values (5, 10, 15)
sens_vals <- c(5, 10, 15)
sens_out <- list()

for (sens in sens_vals) {
  
  res <- run_div_pgls(young_cut_pct = sens, tm = tm, rep_df = rep_df, nwk = nwk)
  
  # extract stats: genome size effect on young_frac_LTR
  sm_div <- summary(res$m_div)
  gs_beta <- sm_div$coefficients["log(GSize)", "Estimate"]
  gs_se   <- sm_div$coefficients["log(GSize)", "Std. Error"]
  gs_p    <- sm_div$coefficients["log(GSize)", "Pr(>|t|)"]
  gs_lam  <- res$m_div$param["lambda"]
  
  sens_out[[paste0(as.character(sens), "_2")]] <- tibble(
    young_cut_pct = sens,

    model = "young_frac_LTR ~ log(GSize)",
    beta = gs_beta,
    se = gs_se,
    p = gs_p,
    lambda = gs_lam
  )
}

sens_tbl <- bind_rows(sens_out) %>%
  arrange(model, young_cut_pct)

sens_tbl
# # A tibble: 3 × 6
# young_cut_pct model                        beta    se      p lambda
# <dbl> <chr>                       <dbl> <dbl>  <dbl>  <dbl>
#   1             5 young_frac_LTR ~ log(GSize) 0.374 0.132 0.0196      1
# 2            10 young_frac_LTR ~ log(GSize) 0.662 0.254 0.0285      1
# 3            15 young_frac_LTR ~ log(GSize) 0.388 0.224 0.117       1


```

## Gap & Telomere ID

Run telociraptor on each assembly to identify telomeres and summarize gap statistics:

```bash
#!/bin/bash

#SBATCH --time=2-00:00:00   
#SBATCH --nodes=1  
#SBATCH --cpus-per-task=4
#SBATCH --mem=16Gb
#SBATCH --partition=ceres
#SBATCH --account=coffea_pangenome

t=4

# Check if the correct number of arguments is provided
set -euo pipefail

module load miniconda
source activate chromsyn

SAMPLE="${1:?usage: $0 <SAMPLE>}"
WD=/project/coffea_pangenome/Artocarpus/Comparative_Paper/assemblies/telomere_ID
GENOME_DIR=/project/coffea_pangenome/Artocarpus/Comparative_Paper/assemblies/unmasked

# ID telomeres
TELO_DIR=/project/coffea_pangenome/Software/Merondun/telociraptor/code
if [ -f ${SAMPLE}.chr.telomeres.tdt]; then
        echo "Telociraptor output exists for ${SAMPLE} – skipping"
else
        python ${TELO_DIR}/telociraptor.py seqin=${GENOME_DIR}/${SAMPLE}.chr.fa basefile=${WD}/${SAMPLE} i=-1 tweak=F telonull=T
fi
```

Merge telomere:

```bash
out=all_samples.chr.telomeres.tsv
echo -e "sample\tSeqName\tSeqLen\tTel5\tTel3\tTel5Len\tTel3Len\tTrim5\tTrim3\tTelPerc" > "$out"

for f in *.chr.telomeres.tdt; do
  sample="${f%%.chr.telomeres.tdt}"
  awk -v s="$sample" 'BEGIN{FS=OFS="\t"} NR==1{next} {print s,$0}' "$f" >> "$out"
done
```

And gaps:

```bash
out=all_samples.chr.gaps.tsv

# write header once
echo -e "sample\tseqname\tstart\tend\tseqlen\tgaplen" > "$out"

for f in *.chr.gaps.tdt; do
  sample="${f%%.chr.gaps.tdt}"
  awk -v s="$sample" 'BEGIN{FS=OFS="\t"} NR==1{next} {print s,$1,$2,$3,$4,$5}' "$f" >> "$out"
done
```

in R:

```R
setwd('/project/coffea_pangenome/Artocarpus/Comparative_Paper/assemblies/telomere_ID')
library(readr)
library(dplyr)
library(tidyverse)

tel <- read_tsv("all_samples.chr.telomeres.tsv", show_col_types = FALSE)
tels <- tel %>% dplyr::select(sample, seqname = SeqName, chrlen = SeqLen, Tel5, Tel3)
gap <- read_tsv("all_samples.chr.gaps.tsv", show_col_types = FALSE)
gaps <- gap %>% group_by(sample,seqname) %>% summarize(seqlen = sum(seqlen),gaps =sum(gaplen))
md <- read_tsv('~/symlinks/comp/samples.txt') %>% dplyr::select(Accession,`Accession Order`, Group,Color)
ordg <- md %>% arrange(`Accession Order`) %>% distinct(Group)
orda <- md %>% arrange(`Accession Order`) %>% distinct(Accession)

df <- tels %>%
  left_join(gaps, by = c("sample", "SeqName")) %>%              # join keys (adjust if your chrom column has a different name)
  tidyr::replace_na(list(gaps = 0)) %>%                         # only fill the gaps column
  left_join(md, by = c("sample" = "Accession")) %>%
  mutate(
    t2t = Tel5 & Tel3,
    t2t_gapless = t2t & gaps == 0,
    Group = factor(Group, levels = ordg$Group),
    sample = factor(sample, levels = orda$Accession)
  )

dfs <- df %>%
  group_by(sample, Group) %>%
  summarise(
    t2t = sum(t2t, na.rm = TRUE),
    gapless = sum(gaps == 0, na.rm = TRUE),
    t2t_gapless = sum(t2t_gapless, na.rm = TRUE),
    .groups = "drop"
  )

tp <- dfs %>%
  mutate(y_lab = paste0(as.character(Group), " (", as.character(sample), ")"),
         y_lab = factor(y_lab, levels = rev(y_lab[order(Group, sample)]))) %>%
  pivot_longer(c(t2t, gapless, t2t_gapless), names_to = "name", values_to = "value") %>%
  ggplot(aes(y = y_lab, x = value, fill = Group)) +
  geom_col() +
  geom_text(aes(label = value), hjust = -0.5) +
  scale_fill_manual(values = md$Color, breaks = md$Group) +
  scale_x_continuous(
    breaks = scales::pretty_breaks(n = 3),
    labels = scales::label_number(accuracy = 1),
    expand = expansion(mult = c(0, 0.25))
  ) +
  coord_cartesian(clip = "off") +
  facet_grid(. ~ name, scales = "free_x") +
  theme_bw() +
  labs(y = "Species (accession)", x = NULL)+
  theme(legend.position = 'top')
tp

ggsave('~/symlinks/comp/figures/20260409_AssemblyStats.pdf',tp,height=5,width=6.5)  

```



## NCBI Formatting

Identify Min and Max gap sizes for upload:

```bash
for f in *.fa; do [ -f "$f" ] || continue; awk -v file="$f" 'BEGIN{RS=">"; FS="\n"} NR>1{seq=""; for(i=2;i<=NF;i++) seq=seq toupper($i); while(match(seq,/N+/)){ n=RLENGTH; if(min==0 || n<min) min=n; if(n>max) max=n; seq=substr(seq,RSTART+RLENGTH) }} END{printf "%s\tmin_N_stretch=%d\tmax_N_stretch=%d\n", file, min+0, max+0}' "$f"; done
```

Rename sequences into this format:

* Take 'Chr*' sequences and remove chr, add location=chr, and assign chr ID. Add accession/isolate lookup based on tab separated accession \t species txt file. 

```
>HART001_1 [organism=Artocarpus altilis] [isolate=HART001] [location=chromosome] [chromosome=1]
>HART001_2 [organism=Artocarpus altilis] [isolate=HART001] [location=chromosome] [chromosome=2]
>HART001_3 [organism=Artocarpus altilis] [isolate=HART001] [location=chromosome] [chromosome=3]
>HART001_4 [organism=Artocarpus altilis] [isolate=HART001] [location=chromosome] [chromosome=4]
>HART001_5 [organism=Artocarpus altilis] [isolate=HART001] [location=chromosome] [chromosome=5]
>HART001_6 [organism=Artocarpus altilis] [isolate=HART001] [location=chromosome] [chromosome=6]
>HART001_7 [organism=Artocarpus altilis] [isolate=HART001] [location=chromosome] [chromosome=7]
>HART001_8 [organism=Artocarpus altilis] [isolate=HART001] [location=chromosome] [chromosome=8]
>HART001_9 [organism=Artocarpus altilis] [isolate=HART001] [location=chromosome] [chromosome=9]
>HART001_10 [organism=Artocarpus altilis] [isolate=HART001] [location=chromosome] [chromosome=10]
>HART001_11 [organism=Artocarpus altilis] [isolate=HART001] [location=chromosome] [chromosome=11]
>HART001_12 [organism=Artocarpus altilis] [isolate=HART001] [location=chromosome] [chromosome=12]
>HART001_13 [organism=Artocarpus altilis] [isolate=HART001] [location=chromosome] [chromosome=13]
>HART001_14 [organism=Artocarpus altilis] [isolate=HART001] [location=chromosome] [chromosome=14]
>HART001_15 [organism=Artocarpus altilis] [isolate=HART001] [location=chromosome] [chromosome=15]
>HART001_16 [organism=Artocarpus altilis] [isolate=HART001] [location=chromosome] [chromosome=16]
>HART001_17 [organism=Artocarpus altilis] [isolate=HART001] [location=chromosome] [chromosome=17]
>HART001_18 [organism=Artocarpus altilis] [isolate=HART001] [location=chromosome] [chromosome=18]
>HART001_19 [organism=Artocarpus altilis] [isolate=HART001] [location=chromosome] [chromosome=19]
>HART001_20 [organism=Artocarpus altilis] [isolate=HART001] [location=chromosome] [chromosome=20]
>HART001_21 [organism=Artocarpus altilis] [isolate=HART001] [location=chromosome] [chromosome=21]
>HART001_22 [organism=Artocarpus altilis] [isolate=HART001] [location=chromosome] [chromosome=22]
>HART001_23 [organism=Artocarpus altilis] [isolate=HART001] [location=chromosome] [chromosome=23]
>HART001_24 [organism=Artocarpus altilis] [isolate=HART001] [location=chromosome] [chromosome=24]
>HART001_25 [organism=Artocarpus altilis] [isolate=HART001] [location=chromosome] [chromosome=25]
>HART001_26 [organism=Artocarpus altilis] [isolate=HART001] [location=chromosome] [chromosome=26]
>HART001_27 [organism=Artocarpus altilis] [isolate=HART001] [location=chromosome] [chromosome=27]
>HART001_28 [organism=Artocarpus altilis] [isolate=HART001] [location=chromosome] [chromosome=28]
>HART001_scaffold_29 [organism=Artocarpus altilis] [isolate=HART001]
>HART001_scaffold_30 [organism=Artocarpus altilis] [isolate=HART001]
>HART001_scaffold_31 [organism=Artocarpus altilis] [isolate=HART001]
>HART001_scaffold_32 [organism=Artocarpus altilis] [isolate=HART001]
```

Save this: 

```python
#!/usr/bin/env python3

import argparse
import os
import re
import sys


def read_lookup(path):
    lookup = {}
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 2:
                continue
            sample = parts[0]
            organism = " ".join(parts[1:])
            lookup[sample] = organism
    return lookup


def infer_sample_from_filename(fasta_path):
    base = os.path.basename(fasta_path)
    sample = re.sub(r'\.(fa|fasta|fna)(\.gz)?$', '', base, flags=re.IGNORECASE)
    return sample


def sanitize_seqid(text):
    # Keep only letters, digits, underscore, dot
    text = re.sub(r'[^A-Za-z0-9._]', '_', text)
    return text


def parse_chr_name(header_core):
    """
    Convert Chr01 -> 1
            Chr1  -> 1
            chr28 -> 28
    Return None if not a chromosome header.
    """
    m = re.fullmatch(r'Chr0*([1-9][0-9]*)', header_core, flags=re.IGNORECASE)
    if m:
        return str(int(m.group(1)))
    return None


def rename_header(orig_header, sample, organism):
    # take only first token from original header
    header_core = orig_header.strip().split()[0]

    chr_name = parse_chr_name(header_core)
    if chr_name is not None:
        seqid = sanitize_seqid(f"{sample}_{chr_name}")
        new_header = (
            f">{seqid} "
            f"[organism={organism}] "
            f"[isolate={sample}] "
            f"[location=chromosome] "
            f"[chromosome={chr_name}]"
        )
        return new_header

    # preserve scaffold/contig names as safe IDs
    safe_core = sanitize_seqid(header_core)
    seqid = sanitize_seqid(f"{sample}_{safe_core}")
    new_header = (
        f">{seqid} "
        f"[organism={organism}] "
        f"[isolate={sample}]"
    )
    return new_header


def process_fasta(infile, outfile, lookup, sample=None):
    if sample is None:
        sample = infer_sample_from_filename(infile)

    if sample not in lookup:
        sys.exit(f"ERROR: sample '{sample}' not found in lookup file")

    organism = lookup[sample]

    with open(infile) as fin, open(outfile, "w") as fout:
        for line in fin:
            if line.startswith(">"):
                orig_header = line[1:].rstrip("\n")
                fout.write(rename_header(orig_header, sample, organism) + "\n")
            else:
                fout.write(line)


def main():
    parser = argparse.ArgumentParser(
        description="Rename FASTA headers for NCBI genome submission."
    )
    parser.add_argument("-l", "--lookup", required=True,
                        help="Lookup table: sample_id<tab/space>organism name")
    parser.add_argument("-i", "--input", required=True,
                        help="Input FASTA")
    parser.add_argument("-o", "--output", required=True,
                        help="Output FASTA")
    parser.add_argument("-s", "--sample", default=None,
                        help="Sample ID override; default inferred from FASTA filename")
    args = parser.parse_args()

    lookup = read_lookup(args.lookup)
    process_fasta(args.input, args.output, lookup, sample=args.sample)


if __name__ == "__main__":
    main()
```

Run on all:

```bash
for i in $(ls *fa | sed 's/.fa//g'); do 
	python3 rename_ncbi_fasta.py   -l sp_lookup.tsv   -i ${i}.fa   -o ${i}.ncbi.fa
done
```

Compress

```
gzip -k *.ncbi.fa
```

Upload

```bash
#!/bin/bash

#SBATCH --time=48:00:00   
#SBATCH --cpus-per-task=1 
#SBATCH --partition=ceres
#SBATCH --account=coffea_pangenome
# source activate assembly
file_list=$1
output_dir=$2

ascp -i Keyfile2.ssh -QT -l500m -k1 --file-checksum md5 --overwrite diff --file-manifest text --file-manifest-path /project/coffea_pangenome/Artocarpus/RawData/Aspera_SRA_Logs -d $(cat ${file_list}) subasp@upload.ncbi.nlm.nih.gov:uploads/tropical_germplasm_genomics_outlook.com_Bb2CzrxE/${output_dir}

```

