# Inferring Whole Genome Duplication History

Compare variation between Batocarpus and Artocarpus by delineating Artocarpus subgenomes, and then using comparative and phylogenetic approaches. 

This section will give parts A & B:

![busco](/figures/20260318_panel_BUSCO_Synteny_SubTreeSimple.png)

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

Summarize total syntenic hit length and the number of buscos:

```R
setwd('/project/coffea_pangenome/Artocarpus/Comparative_Paper/duplication_orthology/rigidus_busco_wga')
library(tidyverse)
library(openxlsx)

md <- read_tsv('~/symlinks/comp/samples.txt')
cols <- md %>% select(Accession, Group, Color)
f <- read.xlsx('chromsyn.xlsx',sheet='BUSCO')
ord <- read_tsv('Samples.list',col_names = F)
f$Genome <- factor(f$Genome,levels=ord$X1)
f$Accession <- sub("_[AB]$", "", f$Genome)
counts <- f %>% count(Genome,Accession) 
counts
cp <- counts %>% 
  ggplot(aes(y = Genome, x = n, fill = Accession)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = n,
                x = n * 1.02),   #push labels slightly to the right
            size = 1.5, hjust = 0) +
  scale_fill_manual(values = cols$Color, breaks = cols$Accession) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.15))) + 
  theme_bw(base_size = 8) +
  theme(
    legend.position = "none",
    panel.grid.major.y = element_blank()
  ) +
  xlab("BUSCO Genes (C+S)") +
  ylab("")
cp

ggsave('~/symlinks/comp/figures/20260310_BUSCO_hits_subgenomes.pdf',cp,height=2.5,width=2)


#Heatmap
h <- read.xlsx('chromsyn.xlsx',sheet='Regions')
h %>% count(Genome,HitGenome)
h2 <- h %>% group_by(Genome,HitGenome) %>% summarize(len = sum(Length), .groups='drop')
h2$Genome <- factor(h2$Genome, levels = rev(ord$X1)) 
h2$HitGenome <- factor(h2$HitGenome, levels = rev(ord$X1))

h2_lower <- h2 %>%
  filter(as.integer(Genome) >= as.integer(HitGenome))

# Make a heatmap
hm <- ggplot(h2_lower, aes(x = Genome, y = HitGenome, fill = len)) +
  geom_tile(color = "white") +
  #scale_fill_viridis_c(option = "plasma", trans = "log10") +
  scale_fill_continuous(low='yellow',high='red', trans='log10') +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  ) +
  labs(
    x = "Hap A",
    y = "Hap B",
    fill = "Total Length",
  )
hm

ggsave('~/symlinks/comp/figures/20260310_HeatmapChromsynBuscoHits.pdf',hm,height=3,width=4.5)

```

