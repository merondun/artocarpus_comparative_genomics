# 02_genome_assembly/

This section documents the end-to-end workflow used to produce chromosome-scale *Artocarpus* assemblies from PacBio HiFi (primary contig generation and polishing) and Hi-C (scaffolding and validation). Assemblies are generated in a consistent, reproducible way across all accessions using Puzzler, and then curated with a standardized post-processing suite: contaminant screening/removal (BlobToolKit/BlobTools). Assembly-level summaries for cross-accession comparison. The final outputs are additionally prepared for public deposition with NCBI-compliant sequence headers and upload packaging. 

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
md <- read_tsv('~/artocarpus_comparative_genomics/samples.txt') %>% dplyr::select(Accession,Accession_Order, Group,Color)
ordg <- md %>% arrange(Accession_Order) %>% distinct(Group)
orda <- md %>% arrange(Accession_Order) %>% distinct(Accession)

df <- tels %>%
  left_join(gaps) %>%              # join keys (adjust if your chrom column has a different name)
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

## Merqury Completeness/Quality

Generate databases with HiFi and genome and estimate:

```bash
#!/bin/bash

#SBATCH --time=1-00:00:00   
#SBATCH --cpus-per-task=16
#SBATCH --mem=128Gb
#SBATCH --partition=ceres
#SBATCH --account=coffea_pangenome

module load miniconda
source activate puzzler192
module load merqury

SAMPLE=$1
WD=/project/coffea_pangenome/Artocarpus/Assemblies/20250101_JustinAssemblies
cd ${WD}

export K=21
export T=24

mkdir -p ${WD}/${SAMPLE}/09_merqury
cd ${WD}/${SAMPLE}/09_merqury

# Build read meryl db
HIFI=/project/coffea_pangenome/Artocarpus/Concatenated_Reads/${SAMPLE}.HiFi.fastq.gz
meryl k=${K} threads=${T} count output ${SAMPLE}.reads.meryl ${HIFI}

# BEFORE / AFTER
awk '/^S/{print ">"$2;print $3}' ${WD}/${SAMPLE}/${SAMPLE}.hic.p_ctg.gfa > ${WD}/${SAMPLE}/pri.init.fa
INIT="${WD}/${SAMPLE}/pri.init.fa"
GENOME="${WD}/primary_asm/${SAMPLE}.fa"
merqury.sh ${SAMPLE}.reads.meryl ${INIT} ${SAMPLE}.before
merqury.sh ${SAMPLE}.reads.meryl ${GENOME} ${SAMPLE}.after
```

Extract:

```bash
for i in $(cat CompSamples.list); do cp ${i}/09_merqury/*after*complete* ${i}/09_merqury/*after.qv merqury_spectra/; done
```

## Estimate heterozygosity

```bash
#!/bin/bash

#SBATCH --time=4-00:00:00    
#SBATCH --cpus-per-task=48
#SBATCH --mem=128Gb
#SBATCH --partition=ceres
#SBATCH --account=coffea_pangenome

module load apptainer
module load miniconda
#mamba create -n snp_array python bbtools minimap2 samtools bcftools glnexus seqkit r-tidyverse r-ape r-ggtree r-treeio bedtools iqtree admixture bioconductor-snprelate r-ranger r-randomforest r-tidymodels mosdepth 
source activate snp_array

SAMPLE=${1:?Missing SAMPLE argument}
WD=/project/coffea_pangenome/Artocarpus/Comparative_Paper
HIFI=/project/coffea_pangenome/Artocarpus/Concatenated_Reads/${SAMPLE}.HiFi.fastq.gz
GENOME=${WD}/assemblies/unmasked/${SAMPLE}.fa

t=48

# subsample based on dpeth to achieve 20x 
X=20
G=$(awk -v s=$SAMPLE '$1==s{print $2*1e6}' ${WD}/het_estimate/genome_sizes.tsv)  # Assembly Mb
BASES=$(printf "%.0f" "$(echo "$X * $G" | bc -l)")
mkdir -p \
  ${WD}/het_estimate/02_subset_reads \
  ${WD}/het_estimate/03_bams \
  ${WD}/het_estimate/04_vcfs \
  ${WD}/het_estimate/05_mosdepth \
  ${WD}/het_estimate/06_callable \
  ${WD}/het_estimate/07_stats

# # Subset reads so that we have roughly equal inputs
# echo "Subsetting Reads for: ${SAMPLE}"
bbduk.sh in=${HIFI} out=${WD}/het_estimate/02_subset_reads/${SAMPLE}.20gb.fastq.gz maxbasesout=${BASES}

# Align, with more relaxed parameters since jackfruit is ~20 MY diverged
echo "Aligning Reads for: ${SAMPLE}"
minimap2 -ax map-hifi -t ${t} \
    -R @RG\\tID:${SAMPLE}\\tPL:PACBIO\\tLB:${SAMPLE}\\tSM:${SAMPLE} ${GENOME} \
    ${WD}/het_estimate/02_subset_reads/${SAMPLE}.20gb.fastq.gz | \
    samtools view -F 4 -bS - | \
    samtools sort -@ ${t} -o ${WD}/het_estimate/03_bams/${SAMPLE}.sorted.bam
samtools index ${WD}/het_estimate/03_bams/${SAMPLE}.sorted.bam

# Call variants, just use diploid deep variant model since we will ignore dosage 
# apptainer pull deepvariant_latest.sif docker://google/deepvariant:latest
echo "Calling SNPs for: ${SAMPLE}"
DEEPVAR="/project/coffea_pangenome/Breadfruit_SNP_Array/containers/deepvariant_gh1060.sif"
apptainer exec \
    -B /project/coffea_pangenome:/project/coffea_pangenome \
    ${DEEPVAR} run_deepvariant \
    --make_examples_extra_args='small_model_call_multiallelics=false' \
    --model_type PACBIO \
    --ref ${GENOME} \
    --reads ${WD}/het_estimate/03_bams/${SAMPLE}.sorted.bam \
    --output_vcf ${WD}/het_estimate/04_vcfs/${SAMPLE}.pt.vcf.gz \
    --output_gvcf ${WD}/het_estimate/04_vcfs/${SAMPLE}.pt.gvcf.gz \
    --sample_name ${SAMPLE} \
    --num_shards ${t} \
    --postprocess_cpus ${t}
tabix -f -p vcf ${WD}/het_estimate/04_vcfs/${SAMPLE}.pt.vcf.gz

# identify callable regions MQ >=20 dp 0.5x-2x expected 
echo "Computing per-base depth (mosdepth) for: ${SAMPLE}"
BAM=${WD}/het_estimate/03_bams/${SAMPLE}.sorted.bam
MDIR=${WD}/het_estimate/05_mosdepth/${SAMPLE}
mkdir -p ${MDIR}

mosdepth -t ${t} -Q 20 ${MDIR}/${SAMPLE} ${BAM}

# estimate expected depth from mosdepth summary (total/mean)
MEAN_DP=$(awk '$1=="total"{print $4}' ${MDIR}/${SAMPLE}.mosdepth.summary.txt)
if [ -z "${MEAN_DP}" ]; then
  echo "ERROR: could not parse mean depth from mosdepth summary: ${MDIR}/${SAMPLE}.mosdepth.summary.txt" >&2
  exit 1
fi

MIN_DP=$(python3 - <<PY
m=float("${MEAN_DP}")
print(max(1, int(m*0.5)))
PY
)

MAX_DP=$(python3 - <<PY
m=float("${MEAN_DP}")
print(int(m*2.0))
PY
)

echo "Mean depth=${MEAN_DP}; callable depth range=[${MIN_DP},${MAX_DP}] (MAPQ>=20)" \
  | tee ${WD}/het_estimate/07_stats/${SAMPLE}.depth_thresholds.txt

# Build callable intervals from per-base depths:
CALLABLE_BED=${WD}/het_estimate/06_callable/${SAMPLE}.callable.MQ20.DP${MIN_DP}-${MAX_DP}.bed
zcat ${MDIR}/${SAMPLE}.per-base.bed.gz | \
  awk -v min=${MIN_DP} -v max=${MAX_DP} 'BEGIN{OFS="\t"} $4>=min && $4<=max {print $1,$2,$3}' | \
  bedtools merge -i - > ${CALLABLE_BED}

CALLABLE_BP=$(awk '{s+=$3-$2} END{print s+0}' ${CALLABLE_BED})
echo -e "${SAMPLE}\t${CALLABLE_BP}" > ${WD}/het_estimate/07_stats/${SAMPLE}.callable_bp.tsv

### heterozygosity time 
echo "Counting heterozygous SNPs in callable regions for: ${SAMPLE}"
VCF=${WD}/het_estimate/04_vcfs/${SAMPLE}.pt.vcf.gz
MIN_GQ=20
MIN_VCF_DP=10

HET_SNPS=$(bcftools view \
  -R ${CALLABLE_BED} \
  -f PASS \
  -v snps \
  -i "GT='het' && GQ>=${MIN_GQ} && DP>=${MIN_VCF_DP}" \
  ${VCF} | \
  bcftools view -H | wc -l)

# heterozygosity per bp (pi-like SNP density proxy)
HET_PER_BP=$(python3 - <<PY
het=int("${HET_SNPS}")
bp=int("${CALLABLE_BP}")
print("nan" if bp==0 else het/bp)
PY
)

printf "%s\t%s\t%s\t%s\n" \
  "${SAMPLE}" \
  "${CALLABLE_BP}" \
  "${HET_SNPS}" \
  "${HET_PER_BP}" \
  > ${WD}/het_estimate/07_stats/${SAMPLE}.het_density.tsv

echo "Done: ${SAMPLE} callable_bp=${CALLABLE_BP} het_snps=${HET_SNPS} het_per_bp=${HET_PER_BP}"
echo -e "${SAMPLE}\t${CALLABLE_BP}\t${HET_SNPS}\t${HET_PER_BP}" > ${WD}/het_estimate/07_stats/${SAMPLE}.output.tsv
```

Plot both the merqury estimates and the heterozygosity:

Plot:

```R
setwd('/project/coffea_pangenome/Artocarpus/Assemblies/20250101_JustinAssemblies/merqury_spectra')
library(tidyverse)
library(readr)
library(stringr)
library(patchwork)
library(ggtext)

md <- read_tsv('~/artocarpus_comparative_genomics/samples.txt') %>% dplyr::select(Accession,Group,Color,Shape,ord=Accession_Order)
qv_files <- list.files('.', pattern = "\\.after\\.qv$", full.names = TRUE)
comp_files <- list.files('.', pattern = "\\.after\\.completeness\\.stats$", full.names = TRUE)

qv <- map_dfr(qv_files, \(f) {
  x <- read_tsv(f, col_names = c("Accession","Errors","Bases","QV","ErrRate"), show_col_types = FALSE)
  x %>% select(Accession, QV)
})

comp <- map_dfr(comp_files, \(f) {
  txt <- read.table(f)
  tibble(Accession = txt$V1,
         completeness = txt$V5)
})

# import het
het <- read_tsv('~/symlinks/comp/het_estimate/20260415_Het_Estimates.tsv', col_names = F)
names(het) <- c('Accession','CallableBP','HeterozygousSNPs','HetRate')
df <- full_join(comp, qv, by = "Accession") %>%
  arrange(Accession) %>%
  left_join(het) %>% 
  left_join(md) %>% 
  mutate(
    ID_md = paste0("*", Group, "*<br>(", Accession, ")"),
    ID_md = factor(ID_md, levels = ID_md[order(ord)])
  ) %>% 
  dplyr::select(-CallableBP, -HeterozygousSNPs)

dp <- df %>% 
  pivot_longer(c(QV, HetRate, completeness)) %>% 
  mutate(
    label = case_when(
      name == "HetRate" ~ scales::label_number(accuracy = 0.0001)(value),
      TRUE              ~ scales::label_number(accuracy = 0.1)(value)
    )
  ) %>%
  ggplot(aes(x = ID_md, y = value, fill = Group, label = label)) +
  geom_col(width = 0.9) +
  geom_label(vjust = 0.5, alpha = 0.8, fill = "white",size=2.5) +
  scale_fill_manual(values = md$Color, breaks = md$Group) +
  facet_grid(name ~ ., scales = "free") +
  theme_bw(base_size = 10) + xlab('')+ylab('')+
  theme(
    legend.position='none',
    axis.text.x = ggtext::element_markdown(angle=45,vjust=1,hjust=1)
  )
dp

ggsave("~/symlinks/comp/figures/20260415_GenomeMerqury_Het.pdf", dp, width = 5, height = 5, dpi = 300)

# correlation?
pear <- cor.test(df$completeness, df$HetRate, method = "pearson")
pear

ct <- cor.test(df$completeness, df$HetRate, method = "pearson")
lab <- sprintf("Pearson r = %.3f\np = %.2g", unname(ct$estimate), ct$p.value)

dp2 <- ggplot(df, aes(x = completeness, y = HetRate)) +
  geom_smooth(aes(group = 1), method = "lm", se = TRUE,
              linewidth = 0.6, color = "grey30") +
  geom_point(aes(fill = Group, shape = Group), size = 2.5) +
  scale_fill_manual(values = md$Color, breaks = md$Group) +
  scale_shape_manual(values = md$Shape, breaks = md$Group) +
  annotate("text", x = Inf, y = Inf, label = lab, hjust = 1, vjust = 1.5) +
  coord_cartesian(clip = "off") +
  theme_bw(base_size = 10) +
  theme(plot.margin = margin(5.5, 25, 5.5, 5.5)) +
  labs(x = "Completeness (%)", y = "HetRate")
dp2
ggsave("~/symlinks/comp/figures/20260415_Correlation_MerqHet.pdf", dp2, width = 5.5, height = 3.5, dpi = 300)


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

