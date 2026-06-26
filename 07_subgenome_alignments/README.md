# 07_subgenome_alignments/

Compare variation between Batocarpus and Artocarpus by delineating Artocarpus subgenomes, and then using comparative and phylogenetic approaches. 

This section will give:

![busco](/figures/panels/03_subgenome_delim/20260520_panel_AGORA_BUSCO.png)

___



## Delineate Artocarpus subgenomes

Based on the WGA, there are consistently 2 chrs from Artocarpus that align to Batocarpus. Identify those pairs:

```bash
R=/project/coffea_pangenome/Artocarpus/Comparative_Paper/assemblies/unmasked/N15_23.chr.fa
Q=/project/coffea_pangenome/Artocarpus/Comparative_Paper/assemblies/unmasked/HART063.chr.fa
mashmap -r ${R} -q ${Q} -t 4 -s 10000 --perc_identity 85 -o Arto_Bato.paf
samtools faidx ${R}
samtools faidx ${Q}
map_chromosomes --paf Arto_Bato.paf --fai ${Q}.fai --out map.txt --min_size 0.1
```

Assign subgenomes based on sequence similarity to Batocarpus:

```bash
awk '
BEGIN{FS=OFS="\t"}
{
  key=$2
  pct=$3; gsub(/%/,"",pct); pct+=0
  line[NR]=$0
  k[NR]=key
  p[NR]=pct
  if (!(key in best) || pct > best[key]) { best[key]=pct; besti[key]=NR }
}
END{
  for (i=1;i<=NR;i++){
    label = (i==besti[k[i]] ? "A" : "B")
    print line[i], label
  }
}
' map.txt > subgenome_assigned.txt
```

Outputs:

```bash
Chr01   Chr01   64.87%  21180220        +       A
Chr02   Chr01   41.30%  23630551        +       B
Chr03   Chr04   45.78%  23977082        +       B
Chr04   Chr04   63.44%  26232693        +       A
Chr05   Chr06   38.97%  22861709        +       B
Chr06   Chr06   56.79%  27820663        +       A
Chr07   Chr08   48.46%  29362978        +       B
Chr08   Chr08   67.99%  35213038        +       A
Chr09   Chr10   45.08%  27194823        +       B
Chr10   Chr10   63.20%  32340011        +       A
Chr11   Chr11   64.16%  33274442        +       A
Chr12   Chr11   46.75%  28575005        +       B
Chr13   Chr13   61.04%  28292573        +       A
Chr14   Chr13   43.43%  21252642        -       B
Chr15   Chr15   62.65%  24103480        +       A
Chr16   Chr15   36.11%  25476228        +       B
Chr17   Chr18   40.06%  27576369        -       B
Chr18   Chr18   62.49%  25124436        +       A
Chr19   Chr19   55.21%  30644708        +       A
Chr20   Chr19   38.86%  24241679        +       B
Chr21   Chr22   41.26%  28988368        +       B
Chr22   Chr22   65.39%  27677218        +       A
Chr23   Chr24   47.12%  29289060        +       B
Chr24   Chr24   69.63%  26712549        +       A
Chr25   Chr26   36.21%  25160240        +       B
Chr26   Chr26   52.99%  30441475        +       A
Chr27   Chr28   43.41%  27806804        +       B
Chr28   Chr28   57.10%  32679041        +       A
```

Apply proper stats excluding overlap from the PAF in R:

```R
library(data.table)
library(dplyr)
library(stringr)
library(meRo) #devtools::install_github('merondun/meRo')

paf <- fread('~/artocarpus_comparative_genomics/07_subgenome_alignments/subgenome_delineation/Arto_Bato.paf.gz', sep = "\t", header = FALSE, quote = "", fill = TRUE, data.table = FALSE)
min_block_bp <- 2500

colnames(paf)[1:12] <- c(
  "qname","qlen","qstart","qend","strand",
  "tname","tlen","tstart","tend",
  "nmatch","alen","mapq"
)

paf <- paf %>%
  mutate(
    qlen   = as.numeric(qlen),
    qstart = as.numeric(qstart),
    qend   = as.numeric(qend),
    tlen   = as.numeric(tlen),
    tstart = as.numeric(tstart),
    tend   = as.numeric(tend),
    alen   = as.numeric(alen)
  ) %>%
  filter(alen >= min_block_bp)

merge_intervals_len <- function(start, end) {
  if (length(start) == 0) return(0)
  o <- order(start, end)
  start <- start[o]
  end <- end[o]
  cs <- start[1]
  ce <- end[1]
  total <- 0
  if (length(start) > 1) {
    for (i in 2:length(start)) {
      if (start[i] <= ce) {
        ce <- max(ce, end[i])
      } else {
        total <- total + (ce - cs)
        cs <- start[i]
        ce <- end[i]
      }
    }
  }
  total + (ce - cs)
}

pair_stats <- paf %>%
  group_by(qname, tname) %>%
  summarise(
    qlen = first(qlen),
    tlen = first(tlen),
    q_cov_bp = merge_intervals_len(qstart, qend),
    t_cov_bp = merge_intervals_len(tstart, tend),
    .groups = "drop"
  ) %>%
  mutate(
    mean_len_bp = (qlen + tlen) / 2,
    q_cov_pct = 100 * q_cov_bp / qlen,
    t_cov_pct = 100 * t_cov_bp / tlen,
    reciprocal_synteny_pct = 100 * ((q_cov_bp + t_cov_bp) / 2) / mean_len_bp
  )

best_pairs <- bind_rows(
  pair_stats %>%
    group_by(qname) %>%
    slice_max(order_by = reciprocal_synteny_pct, n = 1, with_ties = FALSE) %>%
    ungroup(),
  pair_stats %>%
    group_by(tname) %>%
    slice_max(order_by = reciprocal_synteny_pct, n = 1, with_ties = FALSE) %>%
    ungroup()
) %>%
  distinct(qname, tname, .keep_all = TRUE) %>%
  arrange(qname, tname)

bp <- best_pairs %>%
         select(
           artocarpus_chr = qname,
           batocarpus_chr = tname,
           artocarpus_len_bp = qlen,
           batocarpus_len_bp = tlen,
           artocarpus_covered_bp = q_cov_bp,
           batocarpus_covered_bp = t_cov_bp,
           artocarpus_cov_pct = q_cov_pct,
           batocarpus_cov_pct = t_cov_pct,
           reciprocal_synteny_pct
         )

summary_tbl <- best_pairs %>%
  summarise(
    n_pairs = n(),
    mean_reciprocal_synteny_pct = mean(reciprocal_synteny_pct),
    median_reciprocal_synteny_pct = median(reciprocal_synteny_pct),
    min_reciprocal_synteny_pct = min(reciprocal_synteny_pct),
    max_reciprocal_synteny_pct = max(reciprocal_synteny_pct),
    mean_artocarpus_cov_pct = mean(q_cov_pct),
    mean_batocarpus_cov_pct = mean(t_cov_pct)
  )
summary_tbl
# # A tibble: 1 × 7
# n_pairs mean_reciprocal_synteny_pct median_reciprocal_synteny_pct min_reciprocal_synteny_pct max_reciprocal_synteny_pct mean_artocarpus_cov_pct mean_batocarpus_cov_pct
# <int>                       <dbl>                         <dbl>                      <dbl>                      <dbl>                   <dbl>                   <dbl>
#   1      28                        40.8                          37.4                       24.1                       55.4                    51.4                    33.3

# Plot 
best_pairs <- best_pairs %>%
  mutate(
    qname_num = as.integer(gsub("^Chr", "", qname)),
    tname_num = as.integer(gsub("^Chr", "", tname))
  ) %>%
  group_by(tname) %>%
  arrange(desc(reciprocal_synteny_pct), .by_group = TRUE) %>%
  mutate(
    subgenome = c("A", "B")[row_number()],
    subgenome = factor(subgenome, levels = c("A", "B"))
  ) %>%
  ungroup() %>%
  mutate(
    pair_label = paste0(qname, " \u2194 ", tname),
    pair_label = forcats::fct_reorder(pair_label, qname_num),
    recip_lab = sprintf("%.1f", reciprocal_synteny_pct)
  )
p <- ggplot(best_pairs, aes(x = pair_label, y = reciprocal_synteny_pct, fill = subgenome)) +
  geom_col(width = 0.8) +
  geom_text(aes(label = recip_lab), hjust = -0.1, size = 2.8) +
  coord_flip() +
  scale_fill_manual(values = c("A" = "#0072B2", "B" = "#E69F00")) +
  labs(
    x = NULL,
    y = "Reciprocal syntenic coverage (%)",
    fill = "Subgenome"
  ) +
  theme_bw(base_size = 10) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "top"
  ) +
  expand_limits(y = max(best_pairs$reciprocal_synteny_pct) + 5)
p

ggsave("~/artocarpus_comparative_genomics/figures/20260417_Arto_Bato_reciprocal_synteny_pairs.pdf", p, width = 5.5, height = 6.5)

best_pairs %>% sum_stats(reciprocal_synteny_pct)
# mean   min   max    sd    se median   iqr conf_low conf_high
# <dbl> <dbl> <dbl> <dbl> <dbl>  <dbl> <dbl>    <dbl>     <dbl>
#   1  40.8  24.1  55.4  9.32  1.76   37.4  15.8     37.2      44.5
best_pairs %>% group_by(subgenome) %>% sum_stats(reciprocal_synteny_pct)
# subgenome  mean   min   max    sd    se median   iqr conf_low conf_high
# <fct>     <dbl> <dbl> <dbl> <dbl> <dbl>  <dbl> <dbl>    <dbl>     <dbl>
#   1 A          49.1  38.6  55.4  4.37 1.17    50.0  2.32     46.6      51.6
# 2 B          32.6  24.1  36.2  3.72 0.993   34.0  3.98     30.4      34.7
```





## BUSCO-level Subgenome Synteny

This workflow splits each genome into haplotype A/B FASTAs (plus a Batocarpus reference), then runs chromsyn to place BUSCO genes onto chromosomes and summarize subgenome-scale synteny using BUSCO anchors. It generates plotting inputs (BUSCO tables, telomere tracks, and repeat/telomere-window scores), merges them into a chromsyn report (PDF/XLSX), and summarizes BUSCO counts and total syntenic block lengths/heatmaps in R.

```bash
WD=/project/coffea_pangenome/Artocarpus/Comparative_Paper/duplication_orthology/rigidus_busco_wga
GENOMES=/project/coffea_pangenome/Artocarpus/Comparative_Paper/assemblies/unmasked

for SAMPLE in HART063 HART027 HART058 HART060; do
	echo "Working on ${SAMPLE}"
	samtools faidx ${GENOMES}/${SAMPLE}.chr.fa $(cat hapA.list) > ${WD}/${SAMPLE}_A.fa
	samtools faidx ${GENOMES}/${SAMPLE}.chr.fa $(cat hapB.list) > ${WD}/${SAMPLE}_B.fa
done
samtools faidx ${GENOMES}/N15_23.chr.fa $(cat hapA.list) > ${WD}/Batocarpus.fa
```

Run the chromsyn workflow:

```bash
#!/bin/bash

#SBATCH --time=2-00:00:00   
#SBATCH --nodes=1  
#SBATCH --cpus-per-task=20
#SBATCH --mem=64Gb
#SBATCH --partition=ceres
#SBATCH --account=coffea_pangenome

t=20

# Check if the correct number of arguments is provided
set -euo pipefail

module load miniconda
source activate chromsyn

FASTA="${1:?usage: $0 <FASTA>}"
TARGET=$(basename ${FASTA} .fa)
FILE=$(realpath ${FASTA})
WD=/project/coffea_pangenome/Artocarpus/Comparative_Paper/duplication_orthology/rigidus_busco_wga

echo "Working on ${TARGET}, file ${FASTA}"
export PYTHONWARNINGS="ignore::SyntaxWarning"

# Prep busco db 
BUSCO_DB=/project/coffea_pangenome/Software/Merondun/busco_downloads
LINEAGE=embryophyta_odb12
if [ -d "${BUSCO_DB}/lineages/${LINEAGE}" ]; then
        echo "BUSCO db ${LINEAGE} already present at ${BUSCO_DB}/lineages/${LINEAGE} – skipping download"
else
        busco --download ${LINEAGE} --download_path ${BUSCO_DB}
fi

mkdir -p ${WD}/work ${WD}/plotting_inputs
cd ${WD}/work

# Generate inputs
TELO_DIR=/project/coffea_pangenome/Software/Merondun/telociraptor/code
if [ -f ${TARGET}.telomeres.tdt]; then
        echo "Telociraptor output exists for ${TARGET} – skipping"
else
        python ${TELO_DIR}/telociraptor.py seqin=${FILE} basefile=${FILE} i=-1 tweak=F telonull=T
fi

# busco 
if [ -f ${TARGET}.busco5.tsv]; then
        echo "BUSCO already ran on ${TARGET} - skipping"
else
        busco -f -o run_${TARGET} -i ${FILE} -l ${BUSCO_DB}/lineages/${LINEAGE} --cpu ${t} -m genome
        cp -v run_${TARGET}/run_${LINEAGE}/full_table.tsv ${TARGET}.busco5.tsv
        rm -rf run_${TARGET}*
fi

# repeat scores
if [ -f ${TARGET}.tidk.tsv]; then
        echo "TIDK already ran on ${TARGET} - skipping"
else
        tidk search --dir search --output ${TARGET} -s AACCCT ${FILE}
        cp -v search/${TARGET}_telomeric_repeat_windows.tsv ${TARGET}.tidk.tsv
fi

# Copy outputs
cp ${TARGET}.tidk.tsv ${TARGET}.busco5.tsv ${TARGET}.telomeres.tdt ${TARGET}.contigs.tdt ${WD}/plotting_inputs/
```

Merge the outputs and plot:

```bash
#!/bin/bash

#SBATCH --time=2-00:00:00   
#SBATCH --nodes=1  
#SBATCH --cpus-per-task=20
#SBATCH --mem=64Gb
#SBATCH --partition=ceres
#SBATCH --account=coffea_pangenome

WD=/project/coffea_pangenome/Artocarpus/Comparative_Paper/duplication_orthology/rigidus_busco_wga
cd ${WD}

> busco.fofn > gaps.fofn > sequences.fofn > tidk.fofn

for i in $(cat Samples.list); do 
    echo -e "${i} ${WD}/plotting_inputs/${i}.busco5.tsv" >> busco.fofn
    echo -e "${i} ${WD}/plotting_inputs/${i}.gaps.tdt" >> gaps.fofn
    echo -e "${i} ${WD}/plotting_inputs/${i}.telomeres.tdt" >> sequences.fofn
    echo -e "${i} ${WD}/plotting_inputs/${i}.tidk.tsv" >> tidk.fofn
done 

Rscript ~/symlinks/software/chromsyn/chromsyn.R labelsize=1.5 opacity=0.4
gs -sDEVICE=pdfwrite -dCompatibilityLevel=1.4 -dPDFSETTINGS=/printer -dNOPAUSE -dQUIET -dBATCH -sOutputFile=output.pdf chromsyn.pdf
```

BUSCO on all subgenomes:

```bash
WD=/project/coffea_pangenome/Artocarpus/Comparative_Paper/duplication_orthology/all_busco_wga
GENOMES=/project/coffea_pangenome/Artocarpus/Comparative_Paper/assemblies/unmasked
mkdir -p ${WD}

for SAMPLE in $(grep -v 'N15_23' CompSamples.list); do
	echo "Working on ${SAMPLE}"
	samtools faidx ${GENOMES}/${SAMPLE}.chr.fa $(cat A.haps) > ${WD}/${SAMPLE}_A.fa
	samtools faidx ${GENOMES}/${SAMPLE}.chr.fa $(cat B.haps) > ${WD}/${SAMPLE}_B.fa
done
samtools faidx ${GENOMES}/N15_23.chr.fa $(cat A.haps) > ${WD}/Batocarpus.fa
```

Estimate BUSCO:

```bash
#!/bin/bash

#SBATCH --time=2-00:00:00   
#SBATCH --nodes=1  
#SBATCH --cpus-per-task=10
#SBATCH --mem=64Gb
#SBATCH --partition=ceres
#SBATCH --account=coffea_pangenome

t=10

# Check if the correct number of arguments is provided
set -euo pipefail

#module load miniconda
#source activate chromsyn

FASTA="${1:?usage: $0 <FASTA>}"
TARGET=$(basename ${FASTA} .fa)
FILE=$(realpath ${FASTA})
WD=/project/coffea_pangenome/Artocarpus/Comparative_Paper/duplication_orthology/all_busco_wga

echo "Working on ${TARGET}, file ${FASTA}"
export PYTHONWARNINGS="ignore::SyntaxWarning"

# Prep busco db 
BUSCO_DB=/project/coffea_pangenome/Software/Merondun/busco_downloads
LINEAGE=embryophyta_odb12
if [ -d "${BUSCO_DB}/lineages/${LINEAGE}" ]; then
        echo "BUSCO db ${LINEAGE} already present at ${BUSCO_DB}/lineages/${LINEAGE} – skipping download"
else
        busco --download ${LINEAGE} --download_path ${BUSCO_DB}
fi

mkdir -p ${WD}/work ${WD}/plotting_inputs
cd ${WD}/work

# busco 
if [ -f ${TARGET}.busco5.tsv]; then
        echo "BUSCO already ran on ${TARGET} - skipping"
else
        busco -f -o run_${TARGET} -i ${FILE} -l ${BUSCO_DB}/lineages/${LINEAGE} --cpu ${t} -m genome
        cp -v run_${TARGET}/run_${LINEAGE}/full_table.tsv ${TARGET}.busco5.tsv
        rm -rf run_${TARGET}*
fi
```

Plot:

```R
setwd('/project/coffea_pangenome/Artocarpus/Comparative_Paper/duplication_orthology/rigidus_busco_wga')
library(tidyverse)
library(openxlsx)
library(ggpubr)
library(ggridges)

md <- read_tsv('~/artocarpus_comparative_genomics/samples.txt') %>% dplyr::select(Accession, ord = Accession_Order, Group, Color, Shape)
b <- read_tsv('/project/coffea_pangenome/Artocarpus/Comparative_Paper/duplication_orthology/all_busco_wga/work/20260415_BUSCO.tsv',col_names = F)
names(b) <- c('buscoid','status','chr','start','end','strand','s1','s2','s3','desc','Accession')
bs <- b %>% 
  drop_na(Accession) %>% 
  group_by(Accession) %>% 
  count(status) %>% 
  ungroup %>% 
  mutate(Accession = gsub('N97_50','N9750',Accession)) %>% 
  separate(Accession,into=c('Accession','hap')) %>% 
  replace_na(list(hap = 'Batocarpus')) 
# sanity on counts
bs %>% filter(status!='Complete') %>% ggplot(aes(y=Accession,col=hap,x=n,shape=status))+geom_point(size=3)+theme_bw() 
bc <- bs %>% filter(status == 'Complete') %>% 
  dplyr::select(Accession,hap,BUSCO=n) %>% 
  mutate(Accession = gsub('N9750','N97_50',Accession)) 

# add length
chrs <- read_tsv('/project/coffea_pangenome/Artocarpus/Comparative_Paper/assemblies/chr_lengths.tsv',col_names = T)
hapa <- read_tsv('/project/coffea_pangenome/Artocarpus/Comparative_Paper/assemblies/A.haps',col_names = F) %>% mutate(hap = 'A') %>% dplyr::rename(chr=X1)
hapb <- read_tsv('/project/coffea_pangenome/Artocarpus/Comparative_Paper/assemblies/B.haps',col_names = F) %>% mutate(hap = 'B') %>% dplyr::rename(chr=X1)
haps <- rbind(hapa,hapb)
chrinfo <- chrs %>% 
  filter(!grepl('Anti|Fic|Morus',ID)) %>% 
  left_join(haps) %>% 
  drop_na(hap) %>% 
  dplyr::rename(Accession=ID) %>% 
  left_join(md) %>% ungroup %>% 
  mutate(Accession = ifelse(Accession == 'N15_23','Batocarpus',Accession),
         hap = ifelse(Accession == 'Batocarpus','Batocarpus',hap),
         ID_md = paste0("*", Group, "*<br>(", Accession, ")"),
                #ID_md = factor(ID_md, levels = ID_md[order(ord)]),
                ID_md = factor(ID_md, levels = rev(unique(ID_md)))
  )
         
chrsum <- chrinfo %>% 
  group_by(Accession,ID_md,hap,Group,ord) %>% 
  summarize(lengths = mean(length)/1e6,
            lengthsd = sd(length)/1e6) %>% 
  ungroup %>% 
  left_join(bc)
lens <- chrsum %>% dplyr::select(ID_md,hap,Group,lengths,lengthsd)

# stats on diff
hap_totals <- chrinfo %>%
  filter(Accession != 'Batocarpus') %>% 
  group_by(Accession, hap) %>%
  summarise(total_bp = sum(length, na.rm = TRUE) / 1e6, .groups = "drop") %>%
  pivot_wider(names_from = hap, values_from = total_bp) %>%
  mutate(diff_mb = (A - B))

tt <- t.test(hap_totals$A, hap_totals$B, paired = TRUE)

mean_diff <- mean(hap_totals$diff_mb, na.rm = TRUE)
mean_diff
tt$p.value
tt

ann_df <- data.frame(
  label = sprintf(
    "Difference A~B\n %.1f (%.1f–%.1f)\nP=%s",
    unname(tt$estimate),
    tt$conf.int[1], tt$conf.int[2],
    format.pval(tt$p.value, digits = 2, eps = 1e-16)
  )
)

# just points
lp <- lens %>% 
  ggplot(aes(y = ID_md, x = lengths, xmin=lengths-lengthsd,xmax=lengths+lengthsd,fill = hap,shape=hap)) +
  geom_errorbar(width=0.5,position=position_dodge(width=0.6))+
  geom_point(size=3,width = 0.9,position=position_dodge(width=0.6)) +
  coord_cartesian(xlim=c(20,45))+
  scale_fill_manual(values = c('#f8766d','#00bf7d','black')) +
  scale_shape_manual(values=c(22,24,8))+
  theme_bw(base_size = 10) + xlab('')+ylab('Chromosome Length (Mb, Mean ± SD)')+
  theme(
    legend.position='top',
    axis.text.y = ggtext::element_markdown()
  ) +  
  geom_text(
    data = ann_df,
    aes(x = Inf, y = Inf, label = label),
    inherit.aes = FALSE,
    hjust = 1.05, vjust = 1.05,
    size = 3
  )
lp


hapalette <- c("#0072B2", "#E69F00", "black")
names(hapalette) <- c('A','B','Batocarpus')

cp <- chrinfo %>% 
  mutate(length = length / 1e6) %>% 
  ggplot(aes(x = length, y = hap, fill = hap)) +
  # half‑width density curves, nudged slightly up
  geom_density_ridges(
    rel_min_height = 0.01,
    scale = 0.2,
    alpha = 0.6,
    color = NA,
    position = position_nudge(y = 0.15)
  )+ 
  # thin boxplots, nudged slightly down
  geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.5,
               position = position_nudge(y = -0.2))+
  # jittered points, nudged up and jittered along y
  geom_jitter(aes(col=hap),height = 0.05, size = 1, stroke=0.01,alpha = 0.5,pch=21) +
  stat_summary(fun = mean, geom = "point",
               shape = 23, size = 2, fill = "white", colour = "black") +
  scale_fill_manual(values=hapalette)+
  scale_color_manual(values=hapalette)+
  theme_bw(base_size = 10) +
  theme(
    legend.position = "top",
    axis.text.y = ggtext::element_markdown(),
    panel.grid.major.y = element_blank(),      
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "gray60"),
    plot.title = element_text(face = "bold", hjust = 0.5)
  ) +
  geom_text(
    data = ann_df,
    aes(x = -Inf, y = Inf, label = label),
    inherit.aes = FALSE,
    hjust = -0.05, vjust = 1,
    size = 1
  )+
  scale_x_continuous(breaks = scales::pretty_breaks(n = 3)) +  
  xlab("Chromosome Length (Mb)")+ylab('')
cp

chrinfo %>%
  ggplot(aes(x = length, y = as.numeric(hap))) +
  # thin boxplots
  geom_boxplot(
    aes(y = as.numeric(hap) - 0.15, group = hap, fill = hap),
    width = 0.2, outlier.shape = NA, alpha = 0.5
  ) +
  # jittered points, nudged up and jittered along x
  geom_jitter(
    aes(y = as.numeric(hap) + 0.25, color = hap),
    width = 0.15, height = 0, size = 1, alpha = 0.7
  ) +
  # mean point (nudge down to sit with the box)
  stat_summary(fun = mean, geom = "point",
               shape = 23, size = 2, fill = "white", colour = "black") +
  # restore categorical y-axis labels
  #scale_y_continuous(breaks = seq_along(hap_levels), labels = hap_levels, expand = expansion(add = c(0.3, 0.3))) +
  scale_fill_manual(values = hapalette) +
  scale_color_manual(values = hapalette) +
  labs(
    x = "",
    y = "Chromosome Length (Mb)",
    title = "Chromosome lengths by haplotype"
  ) +
  coord_cartesian(clip = "off") + 
  theme_bw(base_size = 10) +
  theme(
    legend.position = "top",
    axis.text.y = ggtext::element_markdown(),
    panel.grid.major.y = element_blank(),      
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "gray60"),
    plot.title = element_text(face = "bold", hjust = 0.5)
  ) +
  geom_text(
    data = ann_df,
    aes(x = -Inf, y = Inf, label = label),
    inherit.aes = FALSE,
    hjust = -0.05, vjust = 1,
    size = 3
  )
lp
# Paired t-test
# 
# data:  hap_totals$A and hap_totals$B
# t = 11.305, df = 9, p-value = 1.277e-06
# alternative hypothesis: true mean difference is not equal to 0
# 95 percent confidence interval:
#   34.60258 51.91411
# sample estimates:
#   mean difference 
# 43.25834 


##### and buscos #####
# stats on diff
b_totals <- chrsum %>%
  filter(Accession != 'Batocarpus') %>% 
  dplyr::select(Accession,ID_md,Group,ord,hap,BUSCO) %>% 
  pivot_wider(names_from = hap, values_from = BUSCO) %>%
  mutate(diff_busco = (A - B))

ttb <- t.test(b_totals$A, b_totals$B, paired = TRUE)

mean_diffb <- mean(b_totals$diff_busco, na.rm = TRUE)
mean_diffb
ttb$p.value
ttb

ann_dfb <- data.frame(
  x = Inf,
  y = Inf,
  label = sprintf(
    "Difference A~B\n %.1f (%.1f–%.1f)\nP=%s",
    unname(ttb$estimate),
    ttb$conf.int[1], ttb$conf.int[2],
    format.pval(ttb$p.value, digits = 2, eps = 1e-16)
  )
)

bp <- chrsum %>% 
  ggplot(aes(x = BUSCO, y = hap, fill = hap)) +
  # half‑width density curves, nudged slightly up
  geom_density_ridges(
    rel_min_height = 0.01,
    scale = 0.2,
    alpha = 0.6,
    color = NA,
    position = position_nudge(y = 0.15)
  )+ 
  # thin boxplots, nudged slightly down
  geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.5,
               position = position_nudge(y = -0.2))+
  # jittered points, nudged up and jittered along y
  geom_jitter(aes(col=hap),height = 0.05, size = 1, alpha = 0.5) +
  stat_summary(fun = mean, geom = "point",
               shape = 23, size = 2, fill = "white", colour = "black") +
  scale_fill_manual(values=hapalette)+
  scale_color_manual(values=hapalette)+
  theme_bw(base_size = 10) +
  theme(
    legend.position = "top",
    axis.text.y = ggtext::element_markdown(),
    panel.grid.major.y = element_blank(),      
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "gray60"),
    plot.title = element_text(face = "bold", hjust = 0.5)
  ) +
  geom_text(
    data = ann_dfb,
    aes(x = -Inf, y = Inf, label = label),
    inherit.aes = FALSE,
    hjust = -0.05, vjust = 1,
    size = 1
  )+
  scale_x_continuous(breaks = scales::pretty_breaks(n = 3)) +  
  xlab("C&S BUSCOs")+ylab('')
bp

chrsum %>% 
  mutate(hap = ifelse(Group == 'Batocarpus spp.','Batocarpus',hap)) %>% 
  ggplot(aes(y = ID_md, x = BUSCO,fill = hap)) +
  geom_col(width=0.9,position=position_dodge(width=0.9))+
  #coord_cartesian(xlim=c(20,45))+
  scale_fill_manual(values = c('#f8766d','#00bf7d','black')) +
  theme_bw(base_size = 10) + xlab('')+ylab('Complete BUSCO')+
  theme(
    legend.position='top',
    axis.text.y = ggtext::element_markdown()
  ) +  
  geom_text(
    data = ann_dfb,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    hjust = 1.05, vjust = 1.05,
    size = 3
  )

both <- ggarrange(cp,bp,common.legend = TRUE)
both
ggsave('~/symlinks/comp/figures/20260416_SubgenomeChrSizeBUSCO_Boxes.pdf',both,height=2.5,width=4)


```

