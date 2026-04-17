# Comparative Genomics: Artocarpus

This project builds a small comparative genomics resource for *Artocarpus* using 11 accessions. We generate and curate high-quality genome assemblies, perform whole-genome alignments, and produce consistent structural + functional annotations across accessions. Downstream, these data are integrated to infer a species tree, quantify gene family evolution (expansions/contractions), and support cross-accession comparisons of genome structure and gene content. 

# QC QA Basic Read Stats

This section captures QC/QA metrics used to validate inputs and track assembly readiness across all accessions. We summarize HiFi read length distributions, estimate genome-wide heterozygosity from GenomeScope2 outputs, and apply read deduplication when warranted. 

Read lengths:

```bash
#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=3
#SBATCH --partition=short

WD=/project/coffea_pangenome/Artocarpus/Assemblies/20241115_JustinAssemblies

SAMPLE=$1

seqkit stats -a ${WD}/${SAMPLE}/${SAMPLE}.HiFi.fastq.gz > ${SAMPLE}.stats
```

Afterwards:

```bash
mergem *stats | sed 's@.*\/@@g; s/.HiFi.fastq.gz//g; s/[,()]//g' > Read_Lengths.txt
```

And heterozygosity from genomescope2 (below):

```bash
for i in $(ls *summary.txt | sed 's/_summary.txt//g'); do grep 'Het' ${i}_summary.txt | sed 's/.*)//g' | sed 's/%.*//g' | awk -v i=${i} '{OFS="\t"}{print $0, i}' > ${i}.het; done

cat *het | sed 's/ //g' > Heterozygosity.txt
```

Deduplication if needed:

```bash
#!/bin/bash
#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=60Gb
#SBATCH --partition=short

FILE=$1
zcat ${FILE} | seqkit rmdup --threads 4 -n -o dedup/${FILE}
```

And Stats if needed:

```bash
#!/bin/bash
#SBATCH --time=24:00:00   
#SBATCH --nodes=1  
#SBATCH --cpus-per-task=1
#SBATCH --mem=30Gb
#SBATCH --partition=short

FILE="$1"

# Get stats on bases / reads, will output to slurm 
bbduk.sh in=${FILE}

# Extract 
BASES=$(grep 'Input:' slurm-$SLURM_JOB_ID.out | awk '{print $4}')
READS=$(grep 'Input:' slurm-$SLURM_JOB_ID.out | awk '{print $2}')
MD5SUM=$(md5sum ${FILE})
echo -e "$MD5SUM\t$BASES\t$READS" > ${FILE}.info
```



## K-mer Based Genome Size Estimation: Genomescope

```bash
#!/bin/bash

#SBATCH --time=08:00:00   # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=5   # 10 processor core(s) per node X 2 threads per core
#SBATCH --mem=75Gb
#SBATCH --partition=short    # standard node(s)

ID=$1
PLOIDY=$2

mkdir -p /project/coffea_pangenome/Artocarpus/GenomeSizeEstimates/20241106_JustinEstimates/Smudge/tmp
TMP=/project/coffea_pangenome/Artocarpus/GenomeSizeEstimates/20241106_JustinEstimates/Smudge/tmp
wd=/project/coffea_pangenome/Artocarpus/GenomeSizeEstimates/20241106_JustinEstimates/Smudge
FILES=/project/coffea_pangenome/Artocarpus/Assemblies/20241115_JustinAssemblies

mkdir -p ${wd}/work ${wd}/work/${ID} ${wd}/histograms 
cd ${wd}/work/${ID}

# Run FastK to create a k-mer database
FastK -v -t1 -k31 -P${TMP} -M70 -T3 -NFastK_Table_${ID} ${FILES}/${ID}/${ID}.HiFi.fastq.gz
Histex -G FastK_Table_${ID} > ${wd}/histograms/${ID}.hist 

rm tmp/${ID}.*fastq

# Genomescope
genomescope2 --input ${wd}/histograms/${ID}.hist --output ${wd}/histograms --ploidy ${PLOIDY} --kmer_length 31 --name_prefix ${ID} --max_kmercov 500
```

Extract info from genomescope: 

```bash
### Extract Stats
for i in $(ls *summary.txt | sed 's/_summary.txt//g'); do 

echo "Processing ${i}"

HETMIN=$(grep 'Het' ${i}_summary.txt | sed 's/.*)//g' | awk '{print $1}') 
HETMAX=$(grep 'Het' ${i}_summary.txt | sed 's/.*)//g' | awk '{print $2}') 
HAPMIN=$(grep 'Hap' ${i}_summary.txt | sed 's/.*)//g' | awk '{print $4}' | sed 's/,//g')
HAPMAX=$(grep 'Hap' ${i}_summary.txt | sed 's/.*)//g' | awk '{print $6}' | sed 's/,//g')
UNIMIN=$(grep 'Unique' ${i}_summary.txt | sed 's/.*)//g' | awk '{print $4}' | sed 's/,//g')
UNIMAX=$(grep 'Unique' ${i}_summary.txt | sed 's/.*)//g' | awk '{print $6}' | sed 's/,//g')

echo -e "Sample\tGS_HetMin\tGS_HetMax\tGS_SizeMin\tGS_SizeMax\tGS_UniqueMin\tGS_UniqueMax" > ${i}.info
echo -e "${i}\t${HETMIN}\t${HETMAX}\t${HAPMIN}\t${HAPMAX}\t${UNIMIN}\t${UNIMAX}" >> ${i}.info

done 
```

Plot:

```R
setwd('/project/coffea_pangenome/Artocarpus/GenomeSizeEstimates/20250101_JustinEstimates')
library(tidyverse)
library(RColorBrewer)
library(ggpubr)
library(ggrepel)
library(viridis)

md <- read_csv('/project/coffea_pangenome/Artocarpus/Assemblies/20250101_JustinAssemblies/samples.csv')
n15_22 <- md %>% filter(sample == 'N15_22') %>% mutate(sample = 'N15_22-NewLib')
md <- rbind(md,n15_22)

info_dat <- NULL
hist_dat <- NULL
files <- list.files('histograms/',pattern = '*info')
for (file in files) {
  id = gsub('.info','',file)
  cat('Processing: ',id,'\n')
  info <- read_tsv(paste0('histograms/',file))
  info_dat <- rbind(info_dat,info)
  hist <- read_tsv(paste0('histograms/',id,'.hist'),col_names=F)
  hist <- hist %>% mutate(Sample = id)
  hist_dat <- rbind(hist_dat,hist)
}

names(hist_dat) <- c('Bin','Coverage','Sample')
hm <- left_join(hist_dat,md %>% dplyr::rename(Sample = sample))
hm <- hm %>% mutate(label = paste0('Gb: ',HiFi,'\nLength: ',round(HiFiLength,0),'\nSize: ',Hist_Est_Size,'Mb'))
hm$Group <- factor(hm$Group,levels=c('A. altilis 2N','A. altilis hybrid 2N','A. altilis 3N','A. altilis hybrid 3N','A. altilis hybrid 4N','A. mariannensis','A. camansi','A. odoratissimus','A. rigidus','A. heterophyllus','A. dadah','A. lacucha','A. nitidus','Batocarpus','Treculia'))
order <- hm %>% arrange(Group) %>% select(Group) %>% unique 
hm$Group <- factor(hm$Group,levels=order$Group)

# Labels
labs <- hm %>% select(Sample,Group,label) %>% unique
hist_plot <- hm %>%  
  filter(Bin >= 5) %>% 
  ggplot(aes(x=Bin,y=Coverage,fill=Group))+
  geom_bar(stat='identity') +
  geom_vline(aes(xintercept=left_peak),lty=2,col='white')+
  geom_text(data=labs,aes(x=Inf,y=Inf,label=label),size=2,vjust=1.25,hjust=1.25,inherit.aes=F)+
  coord_cartesian(xlim=c(5,125))+ 
  facet_wrap(Group+Sample~.,scales='free')+
  #scale_fill_manual(values=brewer.pal(length(unique(hm$Group)),name='Paired'))+
  scale_fill_manual(values=viridis(15,option='turbo'))+
  theme_bw()+
  theme(axis.ticks.y=element_blank(),
        axis.text.y=element_blank())
hist_plot

ggsave('~/symlinks/figbf/20250310_Histograms_Breadfruit.pdf',hist_plot,height=9,width=9,dpi=300)

```

Color blindness checks

```R
library(colorBlindness)
library(viridis)
library(RColorBrewer)
displayAllColors(brewer.pal(10,''),color='white')
displayAllColors(viridis(10,option='turbo'),color='white')
```



# Genome Assembly

This section documents the end-to-end workflow used to produce chromosome-scale *Artocarpus* assemblies from PacBio HiFi (primary contig generation and polishing) and Hi-C (scaffolding and validation). Assemblies are generated in a consistent, reproducible way across all accessions using Puzzler, and then curated with a standardized post-processing suite: contaminant screening/removal (BlobToolKit/BlobTools), repeat annotation and harmonization across samples (EarlGrey + merged libraries), and assembly-level summaries for cross-accession comparison. The final outputs are additionally prepared for public deposition with NCBI-compliant sequence headers and upload packaging. 

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
library(ggrepel)
library(caper)

# Add metadata information
md <- read_tsv('~/artocarpus_comparative_genomics/samples.txt')

##### High Level #####
high_level <- read_tsv('~/artocarpus_comparative_genomics/03_annotation/repeat_annotation/Repeat_Summaries.txt') 
names(high_level) <- c('Classification','Coverage','Count','Proportion','Gen','Distinct_Classifications','Accession')
non_repeat <- high_level %>% 
  group_by(Accession) %>% 
  summarize(Proportion = 1-(sum(Proportion)),
            Classification = "Non Repeat",
            Coverage=NA,Count=NA,Gen=NA,Distinct_Classifications=NA)
full_level <- rbind(high_level,non_repeat) %>% filter(!grepl('Anti|Ficus|Morus',Accession))
fl <- left_join(full_level,md)
ord <- md %>% dplyr::select(Accession,Group,Accession_Order) %>% arrange(Accession_Order)
grpord <- md %>% dplyr::select(Group,Accession_Order) %>% arrange(Accession_Order) %>% dplyr::select(Group) %>% distinct
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

# fl %>% dplyr::select(Accession,Group,Classification,Proportion) %>% pivot_wider(names_from = Classification,values_from=Proportion)
# # A tibble: 11 × 11
# Accession Group               DNA    LINE   LTR `Other (Simple Repeat, Microsatellite, RNA)` `Rolling Circle`      SINE Unclassified Penelope `Non Repeat`
# <fct>     <fct>             <dbl>   <dbl> <dbl>                                        <dbl>            <dbl>     <dbl>        <dbl>    <dbl>        <dbl>
#   1 HART001   A. altilis       0.0385 0.00807 0.414                                       0.0199          0.00531   3.81e-4        0.200 NA              0.314
# 2 HART027   A. heterophyllus 0.0230 0.00449 0.466                                       0.0190          0.00320  NA              0.204  1.24e-4        0.280
# 3 HART058   A. rigidus       0.0342 0.00776 0.388                                       0.0198          0.00349   2.53e-4        0.232 NA              0.315
# 4 HART060   A. dadah         0.0386 0.00721 0.336                                       0.0322          0.00192  NA              0.257  1.71e-3        0.325
# 5 HART061   A. odoratissimus 0.0378 0.00347 0.394                                       0.0184          0.00681   2.42e-5        0.226 NA              0.313
# 6 HART062   A. odoratissimus 0.0381 0.00374 0.403                                       0.0175          0.00462  NA              0.221 NA              0.312
# 7 HART063   A. camansi       0.0467 0.00601 0.382                                       0.0194          0.00299  NA              0.208 NA              0.335
# 8 HART067   A. mariannensis  0.0418 0.00676 0.378                                       0.0222          0.00574   3.21e-4        0.209 NA              0.336
# 9 HART068   A. nitidus       0.0375 0.00495 0.364                                       0.0203          0.00438   4.30e-6        0.222  3.55e-4        0.346
# 10 N15_23    Batocarpus sp.   0.0593 0.00617 0.360                                       0.0379          0.00573   1.62e-4        0.214 NA              0.317
# 11 N97_50    A. lacucha       0.0437 0.00337 0.346                                       0.0186          0.00117   6.78e-4        0.230  6.51e-4        0.356
write.csv(fl %>% dplyr::select(Accession,Accession_Order,Group,Classification,Coverage,Count,Proportion,Distinct_Classifications),'~/artocarpus_comparative_genomics/03_annotation/repeat_annotation/Repeat_proportions_coverage_summarized.csv',quote = F,row.names = F)

#### 	PGLS ####
# Significance genome size ~ LTRs 
rep_df <- fl %>% filter(!grepl('Non Repeat',Classification)) %>% 
  dplyr::rename(GSize = Genome_Size_Assembly_Mb,
                phylo_order = Accession_Order) %>% 
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
pgls_in <- rep_df %>% ungroup %>% dplyr::select(Accession,Repeats,GSize) %>% mutate( Accession = gsub("_","", Accession) ) %>% as.data.frame
nwk <- read.tree('~/artocarpus_comparative_genomics/05_orthofinder/SpeciesTree_rooted_node_labels.txt')
nwk$node.label <- NULL
comp <- comparative.data(phy = nwk, data = pgls_in, names.col = Accession, vcv = TRUE, na.omit = FALSE )

pgls_model <- pgls( log(GSize) ~ log(Repeats) , data = comp, lambda = "ML" ) 
summary(pgls_model)
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


pgls_line <- data.frame(
  Repeats = pgls_in$Repeats,
  fitted = fitted(pgls_model)
)

profile <- pgls.profile(pgls_model) 
plot(profile)

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
  scale_shape_manual(values=md$Shape,breaks=md$Group)+
  scale_fill_manual(values=md$Color,breaks=md$Group)+
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
t <- read_tsv('~/artocarpus_comparative_genomics/03_annotation/repeat_annotation/Divergence_Summaries.txt')
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
  scale_shape_manual(values=md$Shape,breaks=md$Group)+
  scale_fill_manual(values=md$Color,breaks=md$Group)+
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





# Iso-seq Gene Annotation

We have Isoseq data for 2 samples (Artocarpus camansi (6 tissues) and Batocarpus sp. (2 tissues)).

These stpes will demultiplex Iso‑Seq reads, convert selected partitions to FASTA, and run eGAPx (with optional short reads just for sensitivity, it was found that Isoseq is sufficient) to generate gene/transcript models. Extract the longest isoform per gene and predict CDS/proteins with TransDecoder, then liftover reference annotations to other assemblies. Finally clean/standardize headers, produce mapping tables, and split CDS/proteomes by chromosome for downstream comparative analyses. 

Primary outputs after this massive codeblock: 

- eGAPx annotations: complete.genomic.gtf and run out/work directories.
- Predicted coding sequences/proteomes: .transdecoder.cds and .transdecoder.pep.
- Longest‑isoform exports: *.longest_transcript_per_gene.gtf and .fa.
- Liftoff transfers: per‑sample GFF3 liftover files.
- Cleaned inputs for downstream: proteomes/, cds/, gtf/ (clean headers) + genes.tsv mappings and chromosome‑split FASTA files.

## Demultiplexing & Input Prep

Split with [pbskera](https://skera.how/) into s-reads, and demultiplex with lima, using example [here](https://skera.how/examples.html). 

This uses the standard mas16 primer file:

```
cat mas16_primers.fasta 
>A
AGCTTACTTGTGAAGA
>B
ACTTGTAAGCTGTCTA
>C
ACTCTGTCAGGTCCGA
>D
ACCTCCTCCTCCAGAA
>E
AACCGGACACACTTAG
>F
AGAGTCCAATTCGCAG
>G
AATCAAGGCTTAACGG
>H
ATGTTGAATCCTAGCG
>I
AGTGCGTTGCGAATTG
>J
AATTGCGTAGTTGGCC
>K
ACACTTGGTCGCAATC
>L
AGTAAGCCTTCGTGTC
>M
ACCTAGATCAGAGCCT
>N
AGGTATGCCGGTTAAG
>O
AAGTCACCGGCACCTT
>P
ATGAAGTGGCTCGAGA
>Q
AGTAGCTGTGTGCA
```

And this isoseq barcode file for demultiplexing:

```
cat Barcodes.fa 
>IsoSeqX_bc01_5p
CTACACGACGCTCTTCCGATCTACTACACGCAATGAAGTCGCAGGGTTGGG
>IsoSeqX_bc02_5p
CTACACGACGCTCTTCCGATCTACTAGTAGCAATGAAGTCGCAGGGTTGGG
>IsoSeqX_bc03_5p
CTACACGACGCTCTTCCGATCTAGTGTACGCAATGAAGTCGCAGGGTTGGG
>IsoSeqX_bc04_5p
CTACACGACGCTCTTCCGATCTATCACTAGCAATGAAGTCGCAGGGTTGGG
>IsoSeqX_bc05_5p
CTACACGACGCTCTTCCGATCTCAGCTGTGCAATGAAGTCGCAGGGTTGGG
>IsoSeqX_bc06_5p
CTACACGACGCTCTTCCGATCTCAGTCACGCAATGAAGTCGCAGGGTTGGG
>IsoSeqX_bc07_5p
CTACACGACGCTCTTCCGATCTCATGTATGCAATGAAGTCGCAGGGTTGGG
>IsoSeqX_bc08_5p
CTACACGACGCTCTTCCGATCTCGTATGTGCAATGAAGTCGCAGGGTTGGG
>IsoSeqX_bc09_5p
CTACACGACGCTCTTCCGATCTGACATGTGCAATGAAGTCGCAGGGTTGGG
>IsoSeqX_bc10_5p
CTACACGACGCTCTTCCGATCTGAGTCTAGCAATGAAGTCGCAGGGTTGGG
>IsoSeqX_bc11_5p
CTACACGACGCTCTTCCGATCTGTAGATAGCAATGAAGTCGCAGGGTTGGG
>IsoSeqX_bc12_5p
CTACACGACGCTCTTCCGATCTGTATGACGCAATGAAGTCGCAGGGTTGGG
>IsoSeqX_3p
AAGCAGTGGTATCAACGCAGAGTAC
```

Run, submit positional for library:

```bash
#!/bin/bash

#SBATCH --time=48:00:00   
#SBATCH --nodes=1  
#SBATCH --cpus-per-task=10
#SBATCH --partition=ceres

module load miniconda
source activate isoseq_ann

#submit e.g. sbatch 00_Skera_Demul.sh m84125_250429_221012_s3.hifi_reads.bcM0002
if [ -z "$1" ]; then
    echo "Error: Library prefix file positional argument is required."
    exit 1
fi

LIBRARY=$1

echo -e "\e[43m~~~~ Demultiplexing isoseq file: ${LIBRARY} ~~~~\e[0m"
skera split ${LIBRARY}.bam mas16_primers.fasta ${LIBRARY}.skera.bam
lima ${LIBRARY}.skera.bam Barcodes.fa ${LIBRARY}.skera.bam --isoseq --overwrite-biosample-names
```

Copy those split and demultiplexed files for annotation:

```bash
#!/bin/bash

#SBATCH --time=2-00:00:00   
#SBATCH --nodes=1  
#SBATCH --ntasks-per-node=16
#SBATCH --partition=ceres
 
samtools fastq -@ 16 m84125_250429_221012_s3.hifi_reads.bcM0001.skera.IsoSeqX_bc01_5p--IsoSeqX_3p.bam | bgzip -c >  /project/coffea_pangenome/Artocarpus/Assemblies/20250101_JustinAssemblies/annotation/isoseq/raw_fastq_nopolyA/N15_23__YL.fastq.gz
samtools fastq -@ 16 m84125_250429_221012_s3.hifi_reads.bcM0001.skera.IsoSeqX_bc02_5p--IsoSeqX_3p.bam | bgzip -c >  /project/coffea_pangenome/Artocarpus/Assemblies/20250101_JustinAssemblies/annotation/isoseq/raw_fastq_nopolyA/N15_23__ML.fastq.gz
samtools fastq -@ 16 m84125_250429_221012_s3.hifi_reads.bcM0001.skera.IsoSeqX_bc03_5p--IsoSeqX_3p.bam | bgzip -c >  /project/coffea_pangenome/Artocarpus/Assemblies/20250101_JustinAssemblies/annotation/isoseq/raw_fastq_nopolyA/HART063__MF.fastq.gz
samtools fastq -@ 16 m84125_250429_221012_s3.hifi_reads.bcM0001.skera.IsoSeqX_bc04_5p--IsoSeqX_3p.bam | bgzip -c >  /project/coffea_pangenome/Artocarpus/Assemblies/20250101_JustinAssemblies/annotation/isoseq/raw_fastq_nopolyA/HART063__FF.fastq.gz
samtools fastq -@ 16 m84125_250429_221012_s3.hifi_reads.bcM0001.skera.IsoSeqX_bc05_5p--IsoSeqX_3p.bam | bgzip -c >  /project/coffea_pangenome/Artocarpus/Assemblies/20250101_JustinAssemblies/annotation/isoseq/raw_fastq_nopolyA/HART063__YL.fastq.gz
samtools fastq -@ 16 m84125_250429_221012_s3.hifi_reads.bcM0001.skera.IsoSeqX_bc06_5p--IsoSeqX_3p.bam | bgzip -c >  /project/coffea_pangenome/Artocarpus/Assemblies/20250101_JustinAssemblies/annotation/isoseq/raw_fastq_nopolyA/HART063__ML.fastq.gz
samtools fastq -@ 16 m84125_250429_221012_s3.hifi_reads.bcM0001.skera.IsoSeqX_bc07_5p--IsoSeqX_3p.bam | bgzip -c >  /project/coffea_pangenome/Artocarpus/Assemblies/20250101_JustinAssemblies/annotation/isoseq/raw_fastq_nopolyA/HART063__FR.fastq.gz
samtools fastq -@ 16 m84125_250429_221012_s3.hifi_reads.bcM0002.skera.IsoSeqX_bc08_5p--IsoSeqX_3p.bam | bgzip -c >  /project/coffea_pangenome/Artocarpus/Assemblies/20250101_JustinAssemblies/annotation/isoseq/raw_fastq_nopolyA/HART063__SD.fastq.gz
```

## eGAPx Annotation

See issues about long read data here: https://github.com/ncbi/egapx/issues/188 

To accomodate this, convert the long read data to fasta:

```bash
#!/bin/bash
#SBATCH --time=1-00:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --partition=ceres
#SBATCH --account=coffea_pangenome

# submit with 
# for f in *.fastq; do sbatch convert_fq_to_fa.sh "$f"; done

FILE="$1"
OUT="${FILE%.fastq}.fasta"

echo "Converting $FILE -> $OUT"
```

After:

```
lt *fasta
-rw-r-----. 1 justin.merondun proj-coffea_pangenome 6.1G Jan 21 18:42 HART063__FF.fasta
-rw-r-----. 1 justin.merondun proj-coffea_pangenome 6.1G Jan 21 18:42 HART063__SD.fasta
-rw-r-----. 1 justin.merondun proj-coffea_pangenome 7.7G Jan 21 18:42 HART063__ML.fasta
-rw-r-----. 1 justin.merondun proj-coffea_pangenome 7.4G Jan 21 18:42 HART063__FR.fasta
-rw-r-----. 1 justin.merondun proj-coffea_pangenome 7.9G Jan 21 18:42 HART063__YL.fasta
-rw-r-----. 1 justin.merondun proj-coffea_pangenome 8.4G Jan 21 18:42 N15_23__YL.fasta
-rw-r-----. 1 justin.merondun proj-coffea_pangenome  12G Jan 21 18:43 N15_23__ML.fasta
-rw-r-----. 1 justin.merondun proj-coffea_pangenome  15G Jan 21 18:43 HART063__MF.fasta
```

I will also try with some available bulk RNAseq from the previous genome:

```bash
total 187G
-rw-r-----. 1 justin.merondun proj-coffea_pangenome 5.3G Jan 20 17:54 external_SRR5997516_1.fastq
-rw-r-----. 1 justin.merondun proj-coffea_pangenome 5.3G Jan 20 17:55 external_SRR5997516_2.fastq
-rw-r-----. 1 justin.merondun proj-coffea_pangenome 5.6G Jan 20 17:55 external_SRR5997517_1.fastq
-rw-r-----. 1 justin.merondun proj-coffea_pangenome 5.6G Jan 20 17:55 external_SRR5997517_2.fastq
-rw-r-----. 1 justin.merondun proj-coffea_pangenome 5.3G Jan 20 17:55 external_SRR5997518_1.fastq
-rw-r-----. 1 justin.merondun proj-coffea_pangenome 5.3G Jan 20 17:55 external_SRR5997518_2.fastq
-rw-r-----. 1 justin.merondun proj-coffea_pangenome 2.2G Jan 20 17:55 external_SRR5997519_1.fastq
-rw-r-----. 1 justin.merondun proj-coffea_pangenome 2.2G Jan 20 17:55 external_SRR5997519_2.fastq
-rw-r-----. 1 justin.merondun proj-coffea_pangenome 2.9G Jan 20 17:55 external_SRR5997520_1.fastq
-rw-r-----. 1 justin.merondun proj-coffea_pangenome 2.9G Jan 20 17:55 external_SRR5997520_2.fastq
-rw-r-----. 1 justin.merondun proj-coffea_pangenome 6.3G Jan 20 17:55 external_SRR5997521_1.fastq
-rw-r-----. 1 justin.merondun proj-coffea_pangenome 6.3G Jan 20 17:55 external_SRR5997521_2.fastq
-rw-r-----. 1 justin.merondun proj-coffea_pangenome 7.4G Jan 20 17:56 external_SRR5997522_1.fastq
-rw-r-----. 1 justin.merondun proj-coffea_pangenome 7.4G Jan 20 17:56 external_SRR5997522_2.fastq
-rw-r-----. 1 justin.merondun proj-coffea_pangenome 5.8G Jan 20 17:56 external_SRR5997523_1.fastq
-rw-r-----. 1 justin.merondun proj-coffea_pangenome 5.8G Jan 20 17:56 external_SRR5997523_2.fastq
-rw-r-----. 1 justin.merondun proj-coffea_pangenome 6.5G Jan 20 17:56 external_SRR5997524_1.fastq
-rw-r-----. 1 justin.merondun proj-coffea_pangenome 6.5G Jan 20 17:56 external_SRR5997524_2.fastq
-rw-r-----. 1 justin.merondun proj-coffea_pangenome 7.8G Jan 20 17:56 external_SRR5997525_1.fastq
-rw-r-----. 1 justin.merondun proj-coffea_pangenome 7.8G Jan 20 17:56 external_SRR5997525_2.fastq
-rw-r-----. 1 justin.merondun proj-coffea_pangenome 7.6G Jan 20 17:56 external_SRR5997526_1.fastq
-rw-r-----. 1 justin.merondun proj-coffea_pangenome 7.6G Jan 20 17:57 external_SRR5997526_2.fastq
-rw-r-----. 1 justin.merondun proj-coffea_pangenome 3.6G Jan 20 17:57 external_SRR5997527_1.fastq
-rw-r-----. 1 justin.merondun proj-coffea_pangenome 3.6G Jan 20 17:57 external_SRR5997527_2.fastq
-rw-r-----. 1 justin.merondun proj-coffea_pangenome 7.6G Jan 20 17:57 external_SRR5997536_1.fastq
-rw-r-----. 1 justin.merondun proj-coffea_pangenome 7.6G Jan 20 17:57 external_SRR5997536_2.fastq
-rw-r-----. 1 justin.merondun proj-coffea_pangenome 8.2G Jan 20 17:57 external_SRR5997537_1.fastq
-rw-r-----. 1 justin.merondun proj-coffea_pangenome 8.2G Jan 20 17:57 external_SRR5997537_2.fastq
-rw-r-----. 1 justin.merondun proj-coffea_pangenome 6.8G Jan 20 17:57 external_SRR5997538_1.fastq
-rw-r-----. 1 justin.merondun proj-coffea_pangenome 6.8G Jan 20 17:57 external_SRR5997538_2.fastq
-rw-r-----. 1 justin.merondun proj-coffea_pangenome 5.0G Jan 20 17:58 external_SRR5997539_1.fastq
-rw-r-----. 1 justin.merondun proj-coffea_pangenome 5.0G Jan 20 17:58 external_SRR5997539_2.fastq
```

Have a long reads file with the paths to the long read fastas:

```bash
cat Hfastas.txt 
HFF     /90daydata/coffea_pangenome/scratch/egapx/isoseq_rename/HART063__FF.fasta
HFR     /90daydata/coffea_pangenome/scratch/egapx/isoseq_rename/HART063__FR.fasta
HMF     /90daydata/coffea_pangenome/scratch/egapx/isoseq_rename/HART063__MF.fasta
HML     /90daydata/coffea_pangenome/scratch/egapx/isoseq_rename/HART063__ML.fasta
HSD     /90daydata/coffea_pangenome/scratch/egapx/isoseq_rename/HART063__SD.fasta
HYL     /90daydata/coffea_pangenome/scratch/egapx/isoseq_rename/HART063__YL.fasta
```

and a short reads file:

```bash
cat Shortreads.txt 
S1      /90daydata/coffea_pangenome/scratch/egapx/rnaseq/external_SRR5997516_1.fastq
S1      /90daydata/coffea_pangenome/scratch/egapx/rnaseq/external_SRR5997516_2.fastq
S2      /90daydata/coffea_pangenome/scratch/egapx/rnaseq/external_SRR5997517_1.fastq
S2      /90daydata/coffea_pangenome/scratch/egapx/rnaseq/external_SRR5997517_2.fastq
S3      /90daydata/coffea_pangenome/scratch/egapx/rnaseq/external_SRR5997518_1.fastq
S3      /90daydata/coffea_pangenome/scratch/egapx/rnaseq/external_SRR5997518_2.fastq
S4      /90daydata/coffea_pangenome/scratch/egapx/rnaseq/external_SRR5997519_1.fastq
S4      /90daydata/coffea_pangenome/scratch/egapx/rnaseq/external_SRR5997519_2.fastq
S5      /90daydata/coffea_pangenome/scratch/egapx/rnaseq/external_SRR5997520_1.fastq
S5      /90daydata/coffea_pangenome/scratch/egapx/rnaseq/external_SRR5997520_2.fastq
S6      /90daydata/coffea_pangenome/scratch/egapx/rnaseq/external_SRR5997521_1.fastq
S6      /90daydata/coffea_pangenome/scratch/egapx/rnaseq/external_SRR5997521_2.fastq
S7      /90daydata/coffea_pangenome/scratch/egapx/rnaseq/external_SRR5997522_1.fastq
S7      /90daydata/coffea_pangenome/scratch/egapx/rnaseq/external_SRR5997522_2.fastq
S8      /90daydata/coffea_pangenome/scratch/egapx/rnaseq/external_SRR5997523_1.fastq
S8      /90daydata/coffea_pangenome/scratch/egapx/rnaseq/external_SRR5997523_2.fastq
S9      /90daydata/coffea_pangenome/scratch/egapx/rnaseq/external_SRR5997524_1.fastq
S9      /90daydata/coffea_pangenome/scratch/egapx/rnaseq/external_SRR5997524_2.fastq
S10     /90daydata/coffea_pangenome/scratch/egapx/rnaseq/external_SRR5997525_1.fastq
S10     /90daydata/coffea_pangenome/scratch/egapx/rnaseq/external_SRR5997525_2.fastq
S11     /90daydata/coffea_pangenome/scratch/egapx/rnaseq/external_SRR5997526_1.fastq
S11     /90daydata/coffea_pangenome/scratch/egapx/rnaseq/external_SRR5997526_2.fastq
S12     /90daydata/coffea_pangenome/scratch/egapx/rnaseq/external_SRR5997527_1.fastq
S12     /90daydata/coffea_pangenome/scratch/egapx/rnaseq/external_SRR5997527_2.fastq
S13     /90daydata/coffea_pangenome/scratch/egapx/rnaseq/external_SRR5997536_1.fastq
S13     /90daydata/coffea_pangenome/scratch/egapx/rnaseq/external_SRR5997536_2.fastq
S14     /90daydata/coffea_pangenome/scratch/egapx/rnaseq/external_SRR5997537_1.fastq
S14     /90daydata/coffea_pangenome/scratch/egapx/rnaseq/external_SRR5997537_2.fastq
S15     /90daydata/coffea_pangenome/scratch/egapx/rnaseq/external_SRR5997538_1.fastq
S15     /90daydata/coffea_pangenome/scratch/egapx/rnaseq/external_SRR5997538_2.fastq
S16     /90daydata/coffea_pangenome/scratch/egapx/rnaseq/external_SRR5997539_1.fastq
S16     /90daydata/coffea_pangenome/scratch/egapx/rnaseq/external_SRR5997539_2.fastq
```

Run egapx:

```bash
#!/bin/bash

#SBATCH --time=2-00:00:00   
#SBATCH --nodes=1  
#SBATCH --cpus-per-task=1
#SBATCH --mem=512Gb
#SBATCH --partition=ceres
#SBATCH --account=coffea_pangenome

module load miniconda
source activate egapx
module load nextflow
module load apptainer

# Variables for genome to annotate, and which isoseq reads to use 
ID=N15_23
ISOdata=ISOfastas
SRdata=None

# Static
RUN="${ID}_${ISOdata}_${SRdata}"
EGAPDIR=/project/coffea_pangenome/Software/Merondun/egapx
WD=/90daydata/coffea_pangenome/scratch/egapx
GENOMES=/project/coffea_pangenome/Artocarpus/Comparative_Paper/assemblies/unmasked
echo "WORKING ON ${RUN}"

# Create yaml for run 
YAML=${WD}/${RUN}.yaml
cat > "$YAML" <<EOF
genome: ${GENOMES}/${ID}.fa
taxid: 241863
long_reads: ${WD}/${ISOdata}.txt
cmsearch:
  enabled: false
trnascan:
  enabled: false
EOF

#long_reads: ${WD}/${ISOdata}.txt
#short_reads: ${WD}/Shortreads.txt
#241863 Bato  194251 Arto altilis 709039 camansi

# Run 
python3 ${EGAPDIR}/ui/egapx.py ${YAML} \
    -c ${EGAPDIR}/egapx_config/ \
    -e slurm \
    -w ${WD}/work/${RUN} -o ${WD}/out/${RUN}
```

Outputs

```bash
./N15_23_Nfastas_None/complete.genomic.gtf 21491
./HART063_Hfastas/complete.genomic.gtf  34139
./HART063_Hfastas_Sreads/complete.genomic.gtf   34914
```

For some analyses (Orthofinder, CAFE5), we only want a single transcript per gene. Use this to extract this, for the references:

```bash
#!/bin/bash

#SBATCH --time=2-00:00:00   
#SBATCH --nodes=1  
#SBATCH --cpus-per-task=6
#SBATCH --mem=32Gb
#SBATCH --partition=ceres
#SBATCH --account=coffea_pangenome

# module load miniconda
# source activate isoseq_ann
WD=/project/coffea_pangenome/Artocarpus/Comparative_Paper/annotation/egapx
mkdir -p ${WD}/only_longest_transcript_per_gene
cd ${WD}/only_longest_transcript_per_gene

for SAMPLE in HART063 N15_23; do 

	TARGET=/project/coffea_pangenome/Artocarpus/Comparative_Paper/assemblies/softmasked/${SAMPLE}.softmasked.fasta
	TDIR=${WD}/${SAMPLE}_IntraspecificISOseq_NoShortRead
	
	agat_sp_keep_longest_isoform.pl --gff ${TDIR}/complete.genomic.gtf -o ${SAMPLE}.longest_transcript_per_gene.gtf
	gffread ${SAMPLE}.longest_transcript_per_gene.gtf -o ${SAMPLE}.gff3 --keep-genes -O -g ${TARGET} -w ${SAMPLE}.fa
	TransDecoder.LongOrfs -t ${SAMPLE}.fa --output_dir .
	TransDecoder.Predict -t ${SAMPLE}.fa --output_dir . --single_best_only

done 
```

Sanity checks: count genes/transcripts/cds/proteins from the files:

```bash
awk 'BEGIN{genes=0; transcripts=0; cds=0} !/^#/{
  if($3=="gene") genes++;
  if($3=="transcript" || $3=="mRNA") transcripts++;
  if($3=="CDS") cds++;
  if(match($0,/protein_id "[^"]+"/)){
    pid[substr($0,RSTART+12,RLENGTH-13)]=1;
  }
} END{
  print "Genes:",genes,"Transcripts:",transcripts,"CDS:",cds,"Unique_Proteins:",length(pid)
}' complete.genomic.gtf
Genes: 34139 Transcripts: 55865 CDS: 304032 Unique_Proteins: 45455

```

## Liftover

Then run liftoff using that gtf onto the other genomes. Run this first once just on 1 sample to generate the database file. I will run using BOTH references for sensitivity, although we will only use liftover from HART063 for all Artocarpus, since Batocarpus is quite diverged. 

```bash
#!/bin/bash

#SBATCH --time=2-00:00:00   
#SBATCH --nodes=1  
#SBATCH --cpus-per-task=12
#SBATCH --mem=128Gb
#SBATCH --partition=ceres
#SBATCH --account=coffea_pangenome

# module load miniconda
# source activate isoseq_ann
t=12
SAMPLE=$1

TARGET=/project/coffea_pangenome/Artocarpus/Comparative_Paper/assemblies/softmasked/${SAMPLE}.softmasked.fasta
WD=/project/coffea_pangenome/Artocarpus/Comparative_Paper/annotation/egapx

for REF in HART063 N15_23; do 

    cd ${WD} 

    rm -rf ${WD}/liftoff_longestiso/${SAMPLE}
    mkdir -p ${WD}/liftoff_longestiso/${SAMPLE} ${WD}/liftoff_longestiso/ref_${REF}
    cd ${WD}/liftoff_longestiso/${SAMPLE}

    REFERENCE=/project/coffea_pangenome/Artocarpus/Comparative_Paper/assemblies/softmasked/${REF}.softmasked.fasta
    TDIR=${WD}/only_longest_transcript_per_gene/
    
    echo "Working on ${SAMPLE} with reference ${REF}"
    liftoff ${TARGET} ${REFERENCE} -g ${TDIR}/${REF}.longest_transcript_per_gene.gtf -o ${WD}/liftoff_longestiso/${SAMPLE}_ref${REF}.gff3 -p ${t}
    cd ${WD}
    rm -rf ${WD}/liftoff_longestiso/${SAMPLE}

done
```

Once the database is made, run for all samples:

```bash
#!/bin/bash

#SBATCH --time=2-00:00:00   
#SBATCH --nodes=1  
#SBATCH --cpus-per-task=12
#SBATCH --mem=128Gb
#SBATCH --partition=ceres
#SBATCH --account=coffea_pangenome

# module load miniconda
# source activate isoseq_ann
t=12
SAMPLE=$1

TARGET=/project/coffea_pangenome/Artocarpus/Comparative_Paper/assemblies/softmasked/${SAMPLE}.softmasked.fasta
WD=/project/coffea_pangenome/Artocarpus/Comparative_Paper/annotation/egapx

for REF in N15_23 HART063; do 

    cd ${WD} 

    rm -rf ${WD}/liftoff_longestiso/${SAMPLE}
    mkdir -p ${WD}/liftoff_longestiso/${SAMPLE} ${WD}/isoliftoff_longest_transcript_per_gene
    cd ${WD}/liftoff_longestiso/${SAMPLE}

    REFERENCE=/project/coffea_pangenome/Artocarpus/Comparative_Paper/assemblies/softmasked/${REF}.softmasked.fasta
    TDIR=${WD}/only_longest_transcript_per_gene/
    
    echo "Working on ${SAMPLE} with reference ${REF}"
    liftoff ${TARGET} ${REFERENCE} -db ${TDIR}/${REF}.longest_transcript_per_gene.gtf_db -o ${WD}/liftoff_longestiso/${SAMPLE}_ref${REF}.gff3 -p ${t} -exclude_partial
    cd ${WD}
    rm -rf ${WD}/liftoff_longestiso/${SAMPLE}
    
    # Extract protein fasta
    cd ${WD}/isoliftoff_longest_transcript_per_gene
    agat_sp_keep_longest_isoform.pl --gff ${WD}/liftoff_longestiso/${SAMPLE}_ref${REF}.gff3 -o ${SAMPLE}_ref${REF}.longest_transcript_per_gene.gtf
    gffread ${SAMPLE}_ref${REF}.longest_transcript_per_gene.gtf -g ${TARGET} -w ${SAMPLE}_ref${REF}.fa
    TransDecoder.LongOrfs -t ${SAMPLE}_ref${REF}.fa --output_dir .
    TransDecoder.Predict -t ${SAMPLE}_ref${REF}.fa --output_dir . --single_best_only
    rm -rf ${SAMPLE}_ref${REF}.fa.transdecoder_dir

done
```

## Formatting for Downstream

First, subset the files we want, and then create a tab delim file with gene info, including chromosome and function.

**For all downstream: I will use HART063 / N15_23 with their isoseq longest isoform annotations**

* **for outgroups (n=3), I will use the N15_23 liftover**
* **for Artocarpus, I will use the HART063 liftover**

```bash
WD=/project/coffea_pangenome/Artocarpus/Comparative_Paper/annotation/egapx/copies_isoliftoff_longest_transcript_per_gene
mkdir -p proteomes cds gtf
# runs.list has a list of all samples EXCEPT HART063 N15_23, because we use the isoseq peptide and not liftover peptide 
for SAMPLE in $(cat run.list); do 
    REF=HART063    
    cp ${SAMPLE}_ref${REF}.fa.transdecoder.pep proteomes/${SAMPLE}.fa;
    cp ${SAMPLE}_ref${REF}.fa.transdecoder.cds cds/${SAMPLE}.fa;
    cp ${SAMPLE}_ref${REF}.longest_transcript_per_gene.gtf gtf/${SAMPLE}.gtf;
done 

# outgroups 
for SAMPLE in Ficus_carica Morus_mongolica Antiaris_toxicaria; do 
	BASE=$(echo ${SAMPLE} | sed 's/_.*//g');
	echo "Copying ${BASE}"
	cp ${SAMPLE}_refN15_23.fa.transdecoder.pep proteomes/${BASE}.fa
	cp ${SAMPLE}_refN15_23.fa.transdecoder.cds cds/${BASE}.fa
	cp ${SAMPLE}_refN15_23.longest_transcript_per_gene.gtf gtf/${BASE}.gtf
done 

# and main targets
for SAMPLE in HART063 N15_23; do 
	cp ../only_longest_transcript_per_gene/${SAMPLE}.fa.transdecoder.pep proteomes/${SAMPLE}.fa
	cp ../only_longest_transcript_per_gene/${SAMPLE}.fa.transdecoder.cds cds/${SAMPLE}.fa
	cp ../only_longest_transcript_per_gene/${SAMPLE}.longest_transcript_per_gene.gtf gtf/${SAMPLE}.gtf	
done 

# and also for N=5 outgroup analysis, rename simply as artocarpus / batocarpus
cp proteomes/HART063.fa proteomes/Artocarpus.fa
cp gtf/HART063.gtf gtf/Artocarpus.gtf
cp cds/HART063.fa cds/Artocarpus.fa
cp proteomes/N15_23.fa proteomes/Batocarpus.fa
cp gtf/N15_23.gtf gtf/Batocarpus.gtf
cp cds/N15_23.fa cds/Batocarpus.fa

# remove underscores since taxon_geneid gets confusing with double underscore 
for dir in gtf cds proteomes; do
  for f in "$dir"/*; do
    mv "$f" "${f//_/}"
  done
done

```

For the N=3 outgroup samples, modify the chr labels for e.g. CM9823XXX becomes Chr01

```bash
cd gtf
for i in Ficus Morus Antiaris; do
    while read old new; do
    	echo "Fixing ${old} chr name to ${new} for ${i}"
        sed -i "s/$old/$new/g" "${i}.gtf"
    done < ../Outgroup_Chr_Map.tsv
done
```

**Format GTF** 

```bash
DIR=/project/coffea_pangenome/Artocarpus/Comparative_Paper/annotation/egapx/copies_isoliftoff_longest_transcript_per_gene/gtf
cd "$DIR"

for f in *.gtf; do
    sample="${f%.gtf}"
    echo "Working on ${sample}"

    awk -v sample="$sample" -F'\t' '
        $3 == "gene" {
            chr = $1
            start = $4
            end = $5
            strand = $7

            # Combine attributes (columns 9..end)
            attr = $9
            for (i = 10; i <= NF; i++) attr = attr " " $i

            # Defaults
            old_gene_id = "NA"
            desc = "NA"
            biotype = "NA"

            # Prefer gene_id, then copy_num_ID, then ID, then locus_tag
            if (match(attr, /gene_id=([^;]+)/, m)) old_gene_id = m[1]
            else if (match(attr, /copy_num_ID=([^;]+)/, m)) old_gene_id = m[1]
            else if (match(attr, /ID=([^;]+)/, m)) old_gene_id = m[1]
            else if (match(attr, /locus_tag=([^;]+)/, m)) old_gene_id = m[1]

            # Clean quotes if any
            gsub(/"/, "", old_gene_id)

            # Description and biotype (optional)
            if (match(attr, /description=([^;]+)/, m)) desc = m[1]
            if (match(attr, /gene_biotype=([^;]+)/, m)) biotype = m[1]

            # Derive gene_num: drop only leading egapxtmp_ prefix, keep copy suffix (_N)
            gene_num = old_gene_id
            sub(/^egapxtmp_/, "", gene_num)

            # First column: sample_geneNumber (ensure uniqueness with copy suffix)
            print sample "_" gene_num "\t" chr "\t" start "\t" end "\t" strand "\t" sample "\t" old_gene_id "\t" biotype "\t" desc
        }
    ' "$f" > "${sample}.genes.tsv"
done
```

Transcripts only:

```bash
for f in *.gtf; do

    sample="${f%.gtf}"
    echo "Working on ${sample}"

    awk -v sample="$sample" -F'\t' '
        $3 == "transcript" {
            chr   = $1; start = $4; end = $5; strand = $7

            # Merge attributes across fields to handle embedded spaces
            attr = $9; for (i = 10; i <= NF; i++) attr = attr " " $i

            # Defaults per record
            tid = old_gene_id = biotype = desc = "NA"

            # Parse key=value pairs exactly
            n = split(attr, a, ";")
            for (i = 1; i <= n; i++) {
                s = a[i]
                gsub(/^[[:space:]]+|[[:space:]]+$/, "", s)
                split(s, kv, "=")
                key = kv[1]; val = kv[2]
                if (key == "transcript_id")        { tid = val; old_gene_id = val }
                else if (key == "transcript_biotype") biotype = val
                else if (key == "product")  { desc = val; gsub(/%2C/, ",", desc) }
            }

            gene_num = tid; sub(/^egapxtmp_/, "", gene_num)

            print sample "_" gene_num "\t" chr "\t" start "\t" end "\t" strand "\t" sample "\t" old_gene_id "\t" biotype "\t" desc
        }
    ' "$f" > "${sample}.genes.tsv"
done
```

Should look like this:

```bash
head gtf/Artocarpus.genes.tsv 
Artocarpus_unassigned_transcript_1      Chr01   66885   69447   +       Artocarpus      unassigned_transcript_1 mRNA    U-box domain-containing protein 28 pseudogene
Artocarpus_004140-R1    Chr01   84232   88381   -       Artocarpus      egapxtmp_004140-R1      mRNA    root UVB sensitive-like protein (Protein of unknown function, DUF647)
Artocarpus_004211-R1    Chr01   88268   92161   -       Artocarpus      egapxtmp_004211-R1      mRNA    pentatricopeptide repeat-containing protein At2g27800, mitochondrial-like isoform 1-T1
Artocarpus_unassigned_transcript_4      Chr01   88276   91925   -       Artocarpus      unassigned_transcript_4 transcript      pentatricopeptide repeat-containing protein At2g27800, mitochondrial-like, transcript variant X7
Artocarpus_029154-R1    Chr01   132092  135486  +       Artocarpus      egapxtmp_029154-R1      mRNA    uridylate kinase PUMPKIN, chloroplastic-like
Artocarpus_014814-R1    Chr01   138641  140098  -       Artocarpus      egapxtmp_014814-R1      mRNA    protein MKS1-like
```

Then, clean up cds and proteins so they match: 

```bash
WD=/project/coffea_pangenome/Artocarpus/Comparative_Paper/annotation/egapx/copies_isoliftoff_longest_transcript_per_gene
cd ${WD}

for dir in cds proteomes; do 
	echo "Cleaning up ${dir}"
	mkdir -p ${WD}/${dir}/original
	cd ${WD}/${dir}
    for fa in *.fa; do
    	cp ${fa} ${WD}/${dir}/original/
        base=$(basename ${fa} .fa)
        echo "Working on ${base}"
        awk -v b=${base} '
            /^>/ {
                id=$1
                sub(/^>/, "", id)
                sub(/^rna-/, "", id)
                sub(/^egapxtmp_/, "", id)

                # remove everything after first "-" OR "."
                sub(/[.].*/, "", id)

                print ">" b "_" id
                next
            }
            { print }
        ' ${fa} > tmp && mv tmp ${fa}
    done

    # ensure there's no duplicates
    grep '^>' *.fa | sort | uniq -d
done 
```

Sanity:

```
head cds/Artocarpus.fa 
>Artocarpus_000001-R1
ATGATCATGTCTTCAAAAGGGTGTTTAGAGGAGATGGGAATATCTTCAACTAATATCAGT
GATGGTGGGAAAAATTGCTATAGAGGCCATTGGAGACCTGCGGAAGACGAGAAACTCCGA
CAACTCGTCGAACAATTCGGTCCTCAGAACTGGAATTTCATCGCCGAGCATCTACAAGGA

head proteomes/Artocarpus.fa 
>Artocarpus_000001-R1
MIMSSKGCLEEMGISSTNISDGGKNCYRGHWRPAEDEKLRQLVEQFGPQNWNFIAEHLQG
RSGKSCRLRWYNQLDPNINKKPFTEEEEERLLSAHRIYGNKWACIAKYFQGRTDNAVKNH
```

Sanity check to ensure the cds headers match the genes.tsv:

```
for fa in *.fa; do
  sample=${fa%.fa}
  tsv="../gtf/${sample}.genes.tsv"
  echo "Checking ${fa} vs ${tsv}"

  awk 'BEGIN{FS="\t"}
    # Read TSV col1 first
    NR==FNR { cnt[$1]++; tsv_total++; next }

    # Read FASTA headers
    /^>/ {
      h = substr($0, 2)
      sub(/[[:space:]]+$/, "", h)
      fasta_total++
      if (!(h in cnt)) { missing[h]=1; nmiss++ }
      else if (cnt[h] != 1) { bad[h]=cnt[h]; nbad++ }
      seen[h]=1
    }

    END {
      # Summary to stderr so the loop can still see failure via exit code
      print "FASTA headers:", fasta_total, " | TSV entries:", tsv_total > "/dev/stderr"

      for (h in missing) print "Missing in TSV:", h > "/dev/stderr"
      for (h in bad)     print "TSV count != 1:", h, "found", bad[h] > "/dev/stderr"
 

      if (nmiss + nbad > 0) exit 1
    }' "$tsv" "$fa" || echo "❌ Problem detected for ${sample}"
done

```

NOW, we need to use that map to divide the CDS into chromosome-specific files:

```bash
WD=/project/coffea_pangenome/Artocarpus/Comparative_Paper/annotation/egapx/copies_isoliftoff_longest_transcript_per_gene/cds
cd ${WD}

mkdir -p chrs
for fa in *.fa; do
    sample=${fa%.fa}
    echo "Splitting ${sample}"
    awk -v s="$sample" '
    NR==FNR { map[$1]=$2; next }
    /^>/ {
        id = substr($1, 2)
        chr = map[id]
        file = "chrs/" s "_" chr ".fa"
    }
    { if (chr != "") print > file }
    ' "${sample}.genes.tsv" "$fa"
done

```

Do the same for proteins:

```bash
WD=/project/coffea_pangenome/Artocarpus/Comparative_Paper/annotation/egapx/copies_isoliftoff_longest_transcript_per_gene/proteomes
cd ${WD}

mkdir -p chrs
for fa in *.fa; do
    sample=${fa%.fa}
    echo "Splitting ${sample}"
    awk -v s="$sample" '
    NR==FNR { map[$1]=$2; next }
    /^>/ {
        id = substr($1, 2)
        chr = map[id]
        file = "chrs/" s "_" chr ".fa"
    }
    { if (chr != "") print > file }
    ' "${sample}.genes.tsv" "$fa"
done
```

## Final Counts

Take the `.gtf` for each species and extract counts:

```R
setwd('/project/coffea_pangenome/Artocarpus/Comparative_Paper/annotation/egapx/copies_isoliftoff_longest_transcript_per_gene/gtf')
library(tidyverse)
library(data.table)
library(stringr)

read_gtf_like <- function(path) {
  dt <- data.table::fread(
    path, sep = "\t", skip = 3, header = FALSE, quote = "", fill = TRUE,
    data.table = TRUE, showProgress = FALSE
  )
  dt <- dt[!grepl("^#", V1)]
  if (ncol(dt) < 9) stop("File seems malformed (<9 columns): ", path)
  
  setnames(dt, c("seqid","source","type","start","end","score","strand","phase","attributes"))
  dt[, start := as.integer(start)]
  dt[, end   := as.integer(end)]
  dt
}

get_attr <- function(x, key) {
  str_match(x, paste0("(^|;)", key, "=([^;]+)"))[,3]
}

split_on_comma <- function(df, col) {
  df %>%
    mutate("{col}" := str_split(.data[[col]], ",")) %>%
    unnest(.data[[col]])
}

summarize_one_file <- function(path) {
  g <- read_gtf_like(path) %>%
    as_tibble() %>%
    mutate(
      ID     = get_attr(attributes, "ID"),
      Parent = get_attr(attributes, "Parent"),
      gene_biotype = get_attr(attributes, "gene_biotype"),
      transcript_biotype = get_attr(attributes, "transcript_biotype"),
      pseudo = get_attr(attributes, "pseudo")
    )
  
  genes <- g %>%
    filter(type == "gene", !is.na(ID)) %>%
    transmute(
      gene_id = ID,
      gene_biotype = gene_biotype,
      pseudo = pseudo,
      # normalize / coarsen biotypes for reporting
      gene_biotype2 = case_when(
        is.na(gene_biotype) | gene_biotype == "" ~ "unknown",
        gene_biotype %in% c("pseudogene", "transcribed_pseudogene") ~ "pseudogene",
        gene_biotype == "protein_coding" ~ "protein_coding",
        gene_biotype == "lncRNA" ~ "lncRNA",
        TRUE ~ gene_biotype
      )
    )
  
  tx <- g %>%
    filter(type %in% c("transcript","mRNA"), !is.na(ID)) %>%
    transmute(
      tx_id = ID,
      gene_id = Parent,
      tx_start = start,
      tx_end = end
    ) %>%
    split_on_comma("gene_id")
  
  exons <- g %>%
    filter(type == "exon", !is.na(Parent)) %>%
    transmute(
      tx_id = Parent,
      exon_len = (end - start + 1L)
    ) %>%
    split_on_comma("tx_id")
  
  cds <- g %>%
    filter(type == "CDS", !is.na(Parent)) %>%
    transmute(
      tx_id = Parent,
      cds_len = (end - start + 1L)
    ) %>%
    split_on_comma("tx_id")
  
  ex_by_tx <- exons %>%
    group_by(tx_id) %>%
    summarise(
      exon_count = n(),
      exon_bp = sum(exon_len),
      .groups = "drop"
    )
  
  cds_by_tx <- cds %>%
    group_by(tx_id) %>%
    summarise(
      cds_bp = sum(cds_len),
      .groups = "drop"
    )
  
  tx2 <- tx %>%
    left_join(ex_by_tx, by = "tx_id") %>%
    left_join(cds_by_tx, by = "tx_id") %>%
    mutate(
      exon_count = replace_na(exon_count, 0L),
      exon_bp    = replace_na(exon_bp, 0L),
      cds_bp     = replace_na(cds_bp, 0L),
      tx_len_bp = if_else(exon_bp > 0L, exon_bp, (tx_end - tx_start + 1L)),
      is_multiexon = exon_count >= 2L
    )
  
  # representative transcript per gene (longest)
  gene_rep <- tx2 %>%
    filter(!is.na(gene_id)) %>%
    group_by(gene_id) %>%
    slice_max(order_by = tx_len_bp, n = 1, with_ties = FALSE) %>%
    ungroup()
  
  gene_rep2 <- gene_rep %>%
    left_join(genes, by = "gene_id") %>%
    mutate(
      is_protein_coding = gene_biotype2 == "protein_coding",
      is_lncRNA         = gene_biotype2 == "lncRNA",
      is_pseudogene     = gene_biotype2 == "pseudogene"
    )
  
  tibble(
    Accession = str_replace(basename(path), "\\.gtf$", ""),
    file = basename(path),
    n_genes = n_distinct(genes$gene_id),
    n_transcripts = n_distinct(tx2$tx_id),
    
    # gene biotype counts (using the normalized gene_biotype2)
    n_protein_coding_genes = sum(genes$gene_biotype2 == "protein_coding", na.rm = TRUE),
    n_pseudogenes          = sum(genes$gene_biotype2 == "pseudogene", na.rm = TRUE),  # includes transcribed_pseudogene
    n_lncRNA_genes         = sum(genes$gene_biotype2 == "lncRNA", na.rm = TRUE),
    n_other_or_unknown_genes = sum(!genes$gene_biotype2 %in% c("protein_coding","pseudogene","lncRNA"), na.rm = TRUE),
    
    # structural stats (using representative transcript per gene)
    mean_tx_len_bp = mean(gene_rep2$tx_len_bp, na.rm = TRUE),
    mean_cds_len_bp = mean(gene_rep2$cds_bp[gene_rep2$is_protein_coding], na.rm = TRUE),
    mean_exons_per_gene = mean(gene_rep2$exon_count, na.rm = TRUE),
    pct_multiexon_genes = 100 * mean(gene_rep2$is_multiexon, na.rm = TRUE)
  )
}

gtfs <- list.files(pattern = "\\.gtf$", full.names = TRUE)

metrics <- purrr::map_dfr(gtfs, summarize_one_file) %>%
  mutate(
    Method = case_when(
      Accession %in% c("N1523","HART063") ~ "Iso-Seq annotation",
      TRUE ~ "Liftover"
    )
  )

metrics

#removes dups
mets <- metrics %>% filter(!grepl('Arto|Bato',Accession))
readr::write_tsv(metrics, "~/artocarpus_comparative_genomics/03_annotation/gene_annotation/AnnotationQC_fromGTF.tsv")

# What leads to the discrepancy between n_genes != pcgs + psuedo?
g <- data.table::fread("HART063.gtf", sep="\t", skip=3, header=FALSE, fill=TRUE, data.table=FALSE)
colnames(g)[1:9] <- c("seqid","source","type","start","end","score","strand","phase","attr")
genes <- g[g$type=="gene", ]
biotype <- stringr::str_match(genes$attr, "(^|;)gene_biotype=([^;]+)")[,3]
sort(table(biotype), decreasing=TRUE)[1:10]
# biotype
# protein_coding             pseudogene                 lncRNA transcribed_pseudogene                   <NA>                   <NA>                   <NA> 
#   27475                   5956                    703                      5                                                                      
# <NA>                   <NA>                   <NA> 

md <- read_tsv('~/artocarpus_comparative_genomics/samples.txt') %>% mutate(Accession = gsub('_','',Accession))
df <- md %>%
  left_join(mets, by = "Accession") %>%
  mutate(
    Species_short = gsub('Artocarpus','A.',Group),
    ylab = paste0(Species_short, " (", Accession, ")"),
    ylab = fct_reorder(ylab, Accession_Order,.desc = TRUE)
  )

plot_df <- df %>%
  transmute(
    ylab, Method,
    protein_coding = n_protein_coding_genes,
    pseudogene = n_pseudogenes,
    lncRNA = n_lncRNA_genes,
    other = n_other_or_unknown_genes
  ) %>%
  pivot_longer(cols = c(protein_coding, pseudogene, lncRNA, other),
               names_to = "biotype", values_to = "n") %>%
  mutate(
    biotype = factor(biotype, levels = c("protein_coding", "pseudogene", "lncRNA", "other"))
  )
lab_df <- plot_df %>%
  filter(biotype == "protein_coding") %>%
  group_by(ylab,Method) %>%
  summarise(n = sum(n), .groups = "drop") %>%
  mutate(ypos = n / 2)

p <- ggplot(plot_df, aes(x = ylab, y = n, fill = biotype)) +
  geom_col(width = 0.8) +
  geom_text(
    data = lab_df,
    aes(x = ylab, y = ypos, label = scales::comma(n)),
    inherit.aes = FALSE,
    color = "white",
    size = 2.8
  ) +
  coord_flip() +
  scale_fill_manual(values = c(
    protein_coding = "#377eb8",
    pseudogene = "#e41a1c",
    lncRNA = "#4daf4a",
    other = "grey70"
  )) +
  facet_grid(Method ~ ., scales = "free_y", space = "free_y") +
  theme_bw(base_size = 9) +
  theme(
    legend.position = "right",
    axis.title = element_blank()
  ) +
  labs(fill = NULL, y = "Genes")

p
ggsave('~/symlinks/comp/figures/20260413_Annotation_Counts.pdf',p,height=5,width=4.5)

```

And after running orthofinder later... return here to show orthogroup overlap by sample: 

```R
setwd('/project/coffea_pangenome/Artocarpus/Comparative_Paper/orthofinder/Subgenome_Divided/OrthoFinder/Results_Feb24/Orthogroups')
library(tidyverse)
library(stringr)
library(readr)

gc <- read_tsv("Orthogroups.GeneCount.tsv", show_col_types = FALSE)
md <- read_tsv('~/artocarpus_comparative_genomics/samples.txt') %>% mutate(Accession = gsub('_','',Accession)) %>% dplyr::select(Accession,Group,ord = Accession_Order)
md <- rbind(md, data.frame(Accession = c('Batocarpus','Morus'),Group = c('Batocarpus sp.','Morus mongolica'),ord = c(98,99)))

# First column is usually Orthogroup
og_col <- names(gc)[1]

long <- gc %>%
  pivot_longer(cols = -all_of(og_col), names_to = "Genome", values_to = "n_genes") %>%
  mutate(
    present = n_genes > 0,
    Accession = str_replace(Genome, "_[AB]$", ""),
    Subgenome = str_extract(Genome, "[AB]$") %>% replace_na("NA")
  )

missing_summary <- long %>%
  group_by(Genome, Accession, Subgenome) %>%
  summarise(
    n_orthogroups = n(),
    n_missing = sum(!present),
    pct_missing = 100 * mean(!present),
    .groups = "drop"
  )

# join to your metadata to color by Method / order nicely
missing_summary2 <- missing_summary %>%
  left_join(md %>% select(Accession, ord, Group), by = "Accession") %>% drop_na(Group) %>% 
  mutate(
    # make a nice label: Group (Accession_subgenome)
    ylab = if_else(Subgenome %in% c("A","B"),
                   paste0(Group, " (", Accession, "_", Subgenome, ")"),
                   paste0(Group, " (", Accession, ")")),
    ylab = fct_reorder(ylab, ord, .desc = TRUE)
  )

op <- ggplot(missing_summary2, aes(x = ylab, y = pct_missing, fill = Subgenome)) +
  geom_col(width = 0.8) +
  coord_flip() +
  theme_bw(base_size = 9) +
  scale_fill_manual(values=c('#f8766d','#00bf7d','black'))+
  labs(x = NULL, y = "% orthogroups missing", fill = "Subgenome")

ggsave('~/symlinks/comp/figures/20260413_Orthogroup_Missingness.pdf',op,height=5,width=4.5)

```





# Whole Genome Alignments

Per‑pair whole‑genome alignments are generated with nucmer, filtered by length/identity and merged into syntenic blocks (within a MAXGAP), then converted to BED/link files (.simple / .cols.simple). These files plus a layout/chromosome table are used by jcvi.graphics.karyotype to render karyotype/link plots that visualize intra‑ and inter‑chromosomal synteny across samples.

Run the alignment:

```bash
#!/bin/bash

#SBATCH --time=2-00:00:00   
#SBATCH --nodes=1  
#SBATCH --cpus-per-task=4
#SBATCH --mem=16Gb
#SBATCH --partition=ceres
#SBATCH --account=coffea_pangenome

# Check if the correct number of arguments is provided
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <first_pair> <second_pair>"
    exit 1
fi

#Submit as cat Pairs.list  | xargs -L 1 sbatch 01_WGA.sh
#module load miniconda
#source activate wga #for mummer
#source activate merothon #for last

# This comes from improvements on the parakeet synteny scripts from https://github.com/lurebgi/monkParakeet
# Use a file with $SAMPLE_1 \t $SAMPLE_2 within pairs.list, submit
target=$1  #reference
query=$2 #target
GENOME_DIR=/project/coffea_pangenome/Artocarpus/Comparative_Paper/wga
echo "Working on $target aligned to $query"

mkdir -p 1kb_merge1mb
cd 1kb_merge1mb
MAXGAP=1000000
MINLEN=2500
MINIDENT=95

# generating nucmer alignments
# nucmer -l 50 -b 400  -t 10 ${GENOME_DIR}/${target}.pri.chr.fa  ${GENOME_DIR}/${query}.pri.chr.fa -p ${target}-${query}
# delta-filter -1 -l ${MINLEN} ${target}-${query}.delta > ${target}-${query}.delta.filt
# show-coords -H -c -l -o -r -T ${target}-${query}.delta.filt > ${target}-${query}.delta.filt.coords

# Filter 
awk -v OFS="\t" -v minlen=$MINLEN -v minid=$MINIDENT '
{
  tchr=$12; tstart=$1; tend=$2
  qchr=$13; qstart=$3; qend=$4
  tlen=$5; qlen=$6
  if (tlen<minlen || qlen<minlen || $7<=minid) next

  torient = (tstart < tend ? "+" : "-")
  qorient = (qstart < qend ? "+" : "-")
  strand  = (torient == qorient ? "+" : "-")

  if (qorient == "-") { tmp=qstart; qstart=qend; qend=tmp }
  if (torient == "-") { tmp=tstart; tstart=tend; tend=tmp }

  print tchr, tstart, tend, qchr, qstart, qend, strand
}' ${target}-${query}.delta.filt.coords | sort -k1,1 -k2,2n -k4,4 -k5,5n > ${target}-${query}.delta.filtered.sorted.tsv

# merge any links that are within MAXGAP on the same chr, for both ref and query 
awk -v OFS="\t" -v maxgap=$MAXGAP '
function flush() {
    if (rt != "") {
        print rt, rs, re, qt, qs, qe, str
    }
}
NR==1 {
    rt=$1; rs=$2; re=$3;
    qt=$4; qs=$5; qe=$6;
    str=$7;
    next
}
{
    if ($1==rt && $4==qt && $7==str && $2 - re <= maxgap && $5 - qe <= maxgap) {
        # extend block
        re = ($3 > re ? $3 : re)
        qe = ($6 > qe ? $6 : qe)
    } else {
        flush()
        rt=$1; rs=$2; re=$3;
        qt=$4; qs=$5; qe=$6;
        str=$7;
    }
}
END { flush() }
' ${target}-${query}.delta.filtered.sorted.tsv > ${target}-${query}.merged.tsv

# Reference and query BEDs
> ${target}-${query}_ref.bed
> ${target}-${query}_qry.bed
> ${target}-${query}.simple

awk -v OFS="\t" -v prefix="${target}-${query}" '
{
    # reference anchors
    print $1,$2,$2+1,$1"__"$2,0,"+" >> (prefix "_ref.bed")
    print $1,$3,$3+1,$1"__"$3,0,"+" >> (prefix "_ref.bed")
    # query anchors
    print $4,$5,$5+1,$4"__"$5,0,"+" >> (prefix "_qry.bed")
    print $4,$6,$6+1,$4"__"$6,0,"+" >> (prefix "_qry.bed")
}
' ${target}-${query}.merged.tsv

# Simple link file
awk -v OFS="\t" -v prefix="${target}-${query}" '
{
    print $1"__"$2,$1"__"$3,$4"__"$5,$4"__"$6,($3-$2),$7
}
' ${target}-${query}.merged.tsv > ${target}-${query}.simple

# Color inter-chromosomal hits
awk '{
  split($1,a,"__"); split($3,b,"__");
  if(a[1] == b[1]) print $0;
  else print "g*"$0;
}' ${target}-${query}.simple > ${target}-${query}.cols.simple

raw=$(cat ${target}-${query}.delta.filt.coords | wc -l)
filtered=$(cat ${target}-${query}.delta.filtered.sorted.tsv | wc -l)
merged=$(cat ${target}-${query}.simple | wc -l)
echo "${raw} raw links; ${filtered} filtered links; ${merged} merged blocks"
```

Submit with:

```
sed '1d' Pairs.list  | xargs -L 1 sbatch 02_Mummer.sh
```

Combine the outputs: 

Create a layout file:

```
cat chr_layout.txt 
#       y,      xstart, xend,   rotation,       color,  label,  va,     bed
        0.92,   0.2,    0.95,   0,      ,       HART001,        top,    HART001.bed
        0.83,   0.2,    0.95,   0,      ,       HART067,        top,    HART067.bed
        0.75,   0.2,    0.95,   0,      ,       HART063,        top,    HART063.bed
        0.67,   0.2,    0.95,   0,      ,       HART061,        top,    HART061.bed
        0.58,   0.2,    0.95,   0,      ,       HART062,        top,    HART062.bed
        0.5,    0.2,    0.95,   0,      ,       HART058,        top,    HART058.bed
        0.42,   0.2,    0.95,   0,      ,       HART027,        top,    HART027.bed
        0.33,   0.2,    0.95,   0,      ,       HART060,        top,    HART060.bed
        0.25,   0.2,    0.95,   0,      ,       N97_50, top,    N97_50.bed
        0.17,   0.2,    0.95,   0,      ,       HART068,        top,    HART068.bed
        0.08,   0.2,    0.95,   0,      ,       N15_23, bottom, N15_23.bed
#       edges
e,      0,      1,      HART001-HART067.cols.simple
e,      1,      2,      HART067-HART063.cols.simple
e,      2,      3,      HART063-HART061.cols.simple
e,      3,      4,      HART061-HART062.cols.simple
e,      4,      5,      HART062-HART058.cols.simple
e,      5,      6,      HART058-HART027.cols.simple
e,      6,      7,      HART027-HART060.cols.simple
e,      7,      8,      HART060-N97_50.cols.simple
e,      8,      9,      N97_50-HART068.cols.simple
e,      9,      10,     HART068-N15_23.cols.simple
```

And the chrs.txt file:

```
cat chrs.txt 
Chr01,Chr02,Chr03,Chr04,Chr05,Chr06,Chr07,Chr08,Chr09,Chr10,Chr11,Chr12,Chr13,Chr14,Chr15,Chr16,Chr17,Chr18,Chr19,Chr20,Chr21,Chr22,Chr23,Chr24,Chr25,Chr26,Chr27,Chr28
Chr01,Chr02,Chr03,Chr04,Chr05,Chr06,Chr07,Chr08,Chr09,Chr10,Chr11,Chr12,Chr13,Chr14,Chr15,Chr16,Chr17,Chr18,Chr19,Chr20,Chr21,Chr22,Chr23,Chr24,Chr25,Chr26,Chr27,Chr28
Chr01,Chr02,Chr03,Chr04,Chr05,Chr06,Chr07,Chr08,Chr09,Chr10,Chr11,Chr12,Chr13,Chr14,Chr15,Chr16,Chr17,Chr18,Chr19,Chr20,Chr21,Chr22,Chr23,Chr24,Chr25,Chr26,Chr27,Chr28
Chr01,Chr02,Chr03,Chr04,Chr05,Chr06,Chr07,Chr08,Chr09,Chr10,Chr11,Chr12,Chr13,Chr14,Chr15,Chr16,Chr17,Chr18,Chr19,Chr20,Chr21,Chr22,Chr23,Chr24,Chr25,Chr26,Chr27,Chr28
Chr01,Chr02,Chr03,Chr04,Chr05,Chr06,Chr07,Chr08,Chr09,Chr10,Chr11,Chr12,Chr13,Chr14,Chr15,Chr16,Chr17,Chr18,Chr19,Chr20,Chr21,Chr22,Chr23,Chr24,Chr25,Chr26,Chr27,Chr28
Chr01,Chr02,Chr03,Chr04,Chr05,Chr06,Chr07,Chr08,Chr09,Chr10,Chr11,Chr12,Chr13,Chr14,Chr15,Chr16,Chr17,Chr18,Chr19,Chr20,Chr21,Chr22,Chr23,Chr24,Chr25,Chr26,Chr27,Chr28
Chr01,Chr02,Chr03,Chr04,Chr05,Chr06,Chr07,Chr08,Chr09,Chr10,Chr11,Chr12,Chr13,Chr14,Chr15,Chr16,Chr17,Chr18,Chr19,Chr20,Chr21,Chr22,Chr23,Chr24,Chr25,Chr26,Chr27,Chr28
Chr01,Chr02,Chr03,Chr04,Chr05,Chr06,Chr07,Chr08,Chr09,Chr10,Chr11,Chr12,Chr13,Chr14,Chr15,Chr16,Chr17,Chr18,Chr19,Chr20,Chr21,Chr22,Chr23,Chr24,Chr25,Chr26,Chr27,Chr28
Chr01,Chr02,Chr03,Chr04,Chr05,Chr06,Chr07,Chr08,Chr09,Chr10,Chr11,Chr12,Chr13,Chr14,Chr15,Chr16,Chr17,Chr18,Chr19,Chr20,Chr21,Chr22,Chr23,Chr24,Chr25,Chr26,Chr27,Chr28
Chr01,Chr02,Chr03,Chr04,Chr05,Chr06,Chr07,Chr08,Chr09,Chr10,Chr11,Chr12,Chr13,Chr14,Chr15,Chr16,Chr17,Chr18,Chr19,Chr20,Chr21,Chr22,Chr23,Chr24,Chr25,Chr26,Chr27,Chr28
Chr01,Chr04,Chr06,Chr08,Chr10,Chr11,Chr13,Chr15,Chr18,Chr19,Chr22,Chr24,Chr26,Chr28
```

And create karyoplot:

```bash
for i in $(cat Samples.list ); do cat ${i}.bed*-* > ${i}.bed; done

# And add a color (e.g. g=green) for inter-chromosomal translocations
for i in $(ls *simple | sed 's/.simple//g'); do 
awk '{
  split($1,a,"__"); split($3,b,"__");
  if(a[1] == b[1]) print $0;
  else print "g*"$0;
}' ${i}.simple > ${i}.cols.simple
done 

python -m jcvi.graphics.karyotype chrs.txt chr_layout.txt
```

Count the green (interchr) and gray (no *g) aligments:

```bash
for f in *.cols.simple; do   awk -v OFS="\t" -v file="$f" '
  function pos(tok,   t) {
    # tok looks like: "Chr01__45905" or "g*Chr27__12441166"
    gsub(/^g\*/, "", tok)
    split(tok, t, /__/)
    return t[2] + 0
  }
  BEGIN{
    gq=gr=aq=ar=0
  }
  {
    # columns: qstart qend rstart rend ...
    qs = pos($1); qe = pos($2)
    rs = pos($3); re = pos($4)

    qlen = (qe - qs); if (qlen < 0) qlen = -qlen
    rlen = (re - rs); if (rlen < 0) rlen = -rlen

    if ($1 ~ /^g\*/) { gq += qlen; gr += rlen }
    else             { aq += qlen; ar += rlen }
  }
  END{
    print file, gq, gr, aq, ar
  }' "$f"; done
```



## Dotplots against HART063

Using HART063 (A. camansi) as the reference, also produce dotplots:

```bash
#!/bin/bash

#SBATCH --time=1-00:00:00    
#SBATCH --cpus-per-task=16
#SBATCH --partition=ceres
#SBATCH --account=coffea_pangenome

SAMPLE=$1

module load apptainer
WD=/project/coffea_pangenome/Artocarpus/Comparative_Paper/wga/2604_HART063_Ref
GENOMES=/project/coffea_pangenome/Artocarpus/Comparative_Paper/assemblies/unmasked
mkdir -p pafs

echo "WGA for ${SAMPLE}"
mashmap -r ${GENOMES}/HART063.chr.fa -q ${GENOMES}/${SAMPLE}.chr.fa -t 10 -s 5000 --perc_identity 90 -o pafs/${SAMPLE}.paf 2> pafs/${SAMPLE}.mashmap.log
Rscript ~/apptainer/paf2dotplot.R pafs/${SAMPLE}.paf -r 1e6 -m 1e4 -p 4 -c 1 --sort-by-refid --identity-lower-color 100

```

Submit: 

```bash
cat Samples.list  | xargs -I {} sbatch -J wga_{} 01_WGA_HART063Ref.sh {}
zip wga_dotplots_hart063.zip pafs/*pdf
```



# Orthofinder

Run orthofinder on our n=10 Artocarpus accessions and our Batocarpus sample, using the protein data from above. 

```bash
#!/bin/bash

#SBATCH --time=2-00:00:00   
#SBATCH --nodes=1  
#SBATCH --ntasks-per-node=20
#SBATCH --mem=128Gb
#SBATCH --partition=ceres
#SBATCH --account=coffea_pangenome

#module load miniconda
#source activate orthofinder

WD=/project/coffea_pangenome/Artocarpus/Comparative_Paper/orthofinder

t=20
RUN=$1
cd ${WD} 

orthofinder -f ${RUN} -a ${t}

# Afterwards, clean up for long term storage since this generates thousdands of files 
#find proteomes/${REF} -type d \( -name "MultipleSequenceAlignments" -o -name "Gene_Trees" -o -name "WorkingDirectory" -o -name "Orthogroup_Sequences" -o -name "Orthologues" -o -name "Resolved_Gene_Trees" \) -exec rm -rf {} +
```



# Ancestral Reconstruction

This builds an N=5 comparative dataset (Artocarpus, Batocarpus, Morus, Ficus, Antiaris) by repeat-masking outgroups, lifting over gene annotations, and running OrthoFinder to infer orthogroups and a rooted species tree. Those outputs are then converted into AGORA inputs to reconstruct ancestral gene order/karyotypes (ancestral chromosomes/CARs) and generate karyotype-style plots plus CAR→extant chromosome mapping summaries.

| Genus              | Genes | Size (Mb) | Chrs |
| ------------------ | ----- | --------- | ---- |
| Artocarpus         | 27052 | 770       | 28   |
| Batocarpus         | 16863 | 532       | 14   |
| Morus mongolica    | 15920 | 342       | 14   |
| Ficus carica       | 13612 | 318       | 12   |
| Antiaris toxicaria | 13349 | 682       | 13   |

## Download and annotate outgroup repeats 

Repeat the earlGrey repeat masking step on the downloaded outgroup genomes:

```bash
#!/bin/bash

#SBATCH --time=14-00:00:00   
#SBATCH --nodes=1  
#SBATCH --cpus-per-task=20
#SBATCH --mem=256Gb
#SBATCH --partition=ceres
#SBATCH --account=coffea_pangenome

t=40

#module load miniconda
#source activate earlgrey
#source activate eg6

SAMPLE=$1
FASTA="/project/coffea_pangenome/Artocarpus/Comparative_Paper/outgroups/${SAMPLE}.fa"

earlGrey -r eudicotyledons -d yes -e yes -t ${t} -g ${FASTA} -s ${SAMPLE} -o earlgrey
```

## Annotate genes with liftover

Using the same script as for the Artocarpus genomes from above. 

```bash
#!/bin/bash

#SBATCH --time=2-00:00:00   
#SBATCH --nodes=1  
#SBATCH --cpus-per-task=12
#SBATCH --mem=128Gb
#SBATCH --partition=ceres
#SBATCH --account=coffea_pangenome

t=12

SAMPLE=$1

TARGET=/project/coffea_pangenome/Artocarpus/Comparative_Paper/assemblies/softmasked/${SAMPLE}.softmasked.fasta
WD=/project/coffea_pangenome/Artocarpus/Comparative_Paper/annotation
REPEATS=/project/coffea_pangenome/Artocarpus/Comparative_Paper/repeats

for REF in N15_23 HART063; do 

    cd ${WD} 

    rm -rf ${WD}/liftoff/${SAMPLE}
    mkdir -p ${WD}/liftoff/${SAMPLE}
    cd ${WD}/liftoff/${SAMPLE}

    REFERENCE=/project/coffea_pangenome/Artocarpus/Comparative_Paper/assemblies/softmasked/${REF}.softmasked.fasta
	TDIR=/project/coffea_pangenome/Artocarpus/Comparative_Paper/annotation/egapx/${REF}_IntraspecificISOseq_NoShortRead
	
    echo "Working on ${SAMPLE} with reference ${REF}"
    liftoff ${TARGET} ${REFERENCE} -db ${TDIR}/complete.genomic.gff_db -o ${WD}/liftoff/${SAMPLE}_ref${REF}.gff3 -p ${t} -exclude_partial
    cd ${WD}
    rm -rf ${WD}/liftoff/${SAMPLE}
	
    # Extract protein fasta
    cd ${WD}/only_longest_transcript_per_gene
    agat_sp_keep_longest_isoform.pl --gff ${WD}/liftoff/${SAMPLE}_ref${REF}.gff3 -o ${SAMPLE}_ref${REF}.longest_transcript_per_gene.gtf
    gffread ${SAMPLE}_ref${REF}.longest_transcript_per_gene.gtf -g ${TARGET} -w ${SAMPLE}_ref${REF}.fa
    TransDecoder.LongOrfs -t ${SAMPLE}_ref${REF}.fa --output_dir .
    TransDecoder.Predict -t ${SAMPLE}_ref${REF}.fa --output_dir . --single_best_only
    rm -rf ${SAMPLE}_ref${REF}.fa.transdecoder_dir

done
```

## Orthofinder on N=5

The rest of these steps need to be performed within  `Comparative_Paper/orthofinder/AGORA_N5`

* **Use the proteins and `.genes.tsv` lookup files from the annotation section which have cleaned IDs for longest isoform only.**

There are also many transcripts which are part of an 'unassigned_transcript' from Gnomon, remove those:

```bash
for fa in *.fa; do
    echo "Cleaning $fa ..."
	mkdir -p original
	cp ${fa} original/
    awk '
        /^>/ {
            # check if header contains "unass"
            if ($0 ~ /unass/) {
                skip = 1
            } else {
                skip = 0
                print
            }
            next
        }
        {
            if (!skip) print
        }
    ' "$fa" > tmp && mv tmp "$fa"
done

for i in $(ls *fa | sed 's/.fa//g'); do 
	raw=$(grep '>' original/${i}.fa | wc -l)
	after=$(grep '>' ${i}.fa | wc -l)
	echo -e "${i}\t${raw}\t${after}"
done 
```

Run orthofinder `sbatch 01_Orthofinder.sh AGORA_N5`: 

```bash
#!/bin/bash

#SBATCH --time=2-00:00:00   
#SBATCH --cpus-per-task=20
#SBATCH --mem=128Gb
#SBATCH --partition=ceres
#SBATCH --account=coffea_pangenome

#module load miniconda
#source activate orthofinder

WD=/project/coffea_pangenome/Artocarpus/Comparative_Paper/orthofinder

t=40
RUN=$1
cd ${WD} 

orthofinder -f ${RUN} -a ${t}
```

## Prep AGORA Files

**Gene files** 

For AGORA, we need to prep a few files. Here is the first one with gene information:

Prep the extant genes files, using the `genes.tsv` files from the annotation step above. 

```bash
for tsv in *.genes.tsv; do
    species=$(basename "$tsv" .genes.tsv)

    echo "Processing $species ..."

    awk -v sp="$species" '
        BEGIN { OFS="\t" }
        {
            gene  = $1
            chr   = $2
            start = $3
            end   = $4
            strand = ($5 == "+" ? 1 : -1)

            print chr, start, end, strand, gene
        }
    ' "$tsv" | bzip2 > genes.${species}.list.bz2
done

mkdir -p /project/coffea_pangenome/Artocarpus/Comparative_Paper/orthofinder/AGORA_N5/Agora_data/{sub,genes_extant}
cp *bz2 ../Agora_data/genes_extant
# just for checking manually
bzcat genes.Antiaris.list.bz2 | head 
Chr01   42648   43256   1       Antiaris_007359
Chr01   48539   54236   -1      Antiaris_006245
Chr01   70301   73934   -1      Antiaris_006319
```

**Prep Ancestral Gene Files**

This will copy all the necessary inputs:

* species tree
* genes at the ancestral nodes using the HOG script
* gene coordinates

``` bash
cd /project/coffea_pangenome/Artocarpus/Comparative_Paper/orthofinder/AGORA_N5/OrthoFinder/Results_Feb04/

mkdir -p /project/coffea_pangenome/Artocarpus/Comparative_Paper/orthofinder/AGORA_N5/hogs/project/coffea_pangenome/Artocarpus/Comparative_Paper/orthofinder/AGORA_N5/Agora_data/{sub,genes_extant}

# Copy species tree
cp Species_Tree/SpeciesTree_rooted_node_labels.txt ../../Agora_data/

# Convert
python /project/coffea_pangenome/Artocarpus/Comparative_Paper/orthofinder/convert_hogs_sp.py -of_hogs Phylogenetic_Hierarchical_Orthogroups/ -outdir /project/coffea_pangenome/Artocarpus/Comparative_Paper/orthofinder/AGORA_N5/hogs
cd /project/coffea_pangenome/Artocarpus/Comparative_Paper/orthofinder/AGORA_N5/hogs

# Copy over
for i in $(ls *list); do cat ${i} | bzip2 > ../Agora_data/sub/${i}.bz2; done
cd /project/coffea_pangenome/Artocarpus/Comparative_Paper/orthofinder/AGORA_N5

# ensure there's no duplicate hogs
for n in 1 2 3; do
  echo "--- Checking N" $n " ---"
  bzcat Agora_data/sub/orthologyGroups.N${n}.list.bz2 \
  | awk '{
    delete seen
    dup=0
    for (i=1;i<=NF;i++) {
      g=$i; sub(/^[[:space:]]+/,"",g); sub(/[[:space:]]+$/,"",g)
      if (seen[g]++) dup=1
    }
    if (dup) print "Line", NR, "has duplicates:", $0
  }'
done
```

Ensure they match with gene files:

```bash
[Agora_data]$ tree
.
├── genes_extant
│   ├── genes.Antiaris.list.bz2
│   ├── genes.Artocarpus.list.bz2
│   ├── genes.Batocarpus.list.bz2
│   ├── genes.Ficus.list.bz2
│   └── genes.Morus.list.bz2
├── SpeciesTree_rooted_node_labels.txt
└── sub
    ├── orthologyGroups.N0.list.list.bz2
    ├── orthologyGroups.N1.list.list.bz2
    ├── orthologyGroups.N2.list.list.bz2
    └── orthologyGroups.N3.list.list.bz2

2 directories, 10 files

bzcat Agora_data/sub/*bz2 | head 
Artocarpus_unassigned_transcript_4244 Ficus_005562-R1  Ficus_015044-R1 
Antiaris_unassigned_transcript_2420 Artocarpus_unassigned_transcript_5716  Artocarpus_unassigned_transcript_5722  Artocarpus_unassigned_transcript_8548 Batocarpus_unassigned_transcript_1669 Ficus_013835-R1 Morus_012309-R1
Ficus_016763-R1  Ficus_unassigned_transcript_4215 

bzcat Agora_data/genes_extant/* | head
Chr01   42648   43256   1       Antiaris_007359-R1
Chr01   48539   54236   -1      Antiaris_006245-R1
Chr01   70301   73934   -1      Antiaris_006319-R1
Chr01   97852   100097  1       Antiaris_006364-R1
Chr01   113618  122025  -1      Antiaris_006746-R1
```

## Run AGORA

```bash
#!/bin/bash

#SBATCH --time=10-00:00:00   
#SBATCH --nodes=1  
#SBATCH --cpus-per-task=5
#SBATCH --mem=16Gb
#SBATCH --partition=ceres
#SBATCH --account=coffea_pangenome

AGORA="/project/coffea_pangenome/Software/Merondun/Agora/src/agora-plants.py"
${AGORA} Agora_data/SpeciesTree_rooted_node_labels.txt Agora_data/sub/orthologyGroups.%s.list.bz2 Agora_data/genes_extant/genes.%s.list.bz2 -workingDir=results -nbThreads=4
```

**Plot AGORA**

Plot using [AGORA](https://github.com/DyogenIBENS/Agora/issues/18):

*1. plot internal ancestral nodes N0-N3 to draw the reference colors (1 bloc / 1 colour):*

```bash
COMPARE="/project/coffea_pangenome/Software/Merondun/Agora/src/misc.compareGenomes.py"
DIR="/project/coffea_pangenome/Artocarpus/Comparative_Paper/orthofinder/AGORA_N5/results"
for NODE in N0 N1 N2 N3; do 

	python ${COMPARE} ${DIR}/ancGenomes/plants-workflow/ancGenome.${NODE}.list.bz2 \
		${DIR}/ancGenomes/plants-workflow/ancGenome.N0.list.bz2 \
		${DIR}/ancGenes/all/ancGenes.N0.list.bz2 \
		-minChrSize=20 -mode=drawKaryotype +sortBySize > ${NODE}.ps
	ps2pdf ${NODE}.ps ${NODE}.pdf
	rm ${NODE}.ps
done 
```

*2. plot species according to N0 colors of blocs:*

```bash
COMPARE="/project/coffea_pangenome/Software/Merondun/Agora/src/misc.compareGenomes.py"
DIR="/project/coffea_pangenome/Artocarpus/Comparative_Paper/orthofinder/AGORA_N5/results"

for SPECIES in Artocarpus Batocarpus Ficus Antiaris Morus; do 

	python ${COMPARE} ${DIR}/genes/genes.${SPECIES}.list.bz2 ${DIR}/ancGenomes/plants-workflow/ancGenome.N0.list.bz2 ${DIR}/ancGenes/all/ancGenes.N0.list.bz2 -mode=drawKaryotype -minChrSize=20 > ${SPECIES}.ps
	ps2pdf ${SPECIES}.ps ${SPECIES}.pdf
	rm ${SPECIES}.ps
done 
zip karyotypes_chrsize20_20260218.zip *pdf*
```

**Plot CARS**

This is useful to identify which chromosomes are syntenic across different assemblies. 

```R
setwd('/project/coffea_pangenome/Artocarpus/Comparative_Paper/orthofinder/AGORA_N5/')
library(tidyverse)
library(viridis)
library(ggrepel)
library(RColorBrewer)
library(ggpubr)
library(data.table)
library(gggenes)
library(ggtree)
library(tidytree)
library(treeio)

# need to bzip2 -d -k *bz2 on results/ancGenomes/plants-workflow/ and results/genes

# concatenated Tree 
tr <- read.tree('Agora_data/SpeciesTree_rooted_node_labels.txt')
t <- ggtree(tr, layout = "rectangular") + 
  geom_tiplab()+
  geom_nodelab(aes(label=label))
t
ggsave('~/symlinks/comp/figures/20260203_AGORA_N5_Species_Tree.pdf',t,height=4,width=3)

anc_file <- c('results/ancGenomes/plants-workflow/ancGenome.N0.list')

# Load extant coordinates for all species
extant_files <- list.files('results/genes/', pattern = "list$", full.names = TRUE)

extant_df <- map_df(extant_files, \(f){
  sp <- str_match(basename(f), "^genes\\.(.+)\\.list$")[,2]
  read_tsv(f, col_names = c("chr","start","end","strand","transcript"),
           col_types = "cdiic") %>%
    mutate(species = sp,
           transcript = sub('^[^.]*\\.','',transcript))
})
extant_df

# Parse AGORA ancestral gene lists
agora <- fread(anc_file)
af <- read_tsv(anc_file, col_names = FALSE, col_types = cols(.default = "c")) %>%
  rename(CAR = X1, anc_i0 = X2, anc_i1 = X3, orient = X4, raw = X5) %>%
  mutate(
    anc_i0 = as.integer(anc_i0),
    anc_i1 = as.integer(anc_i1),
    orient = as.integer(orient)
  ) %>%
  # Split raw column into anc_gene and extant tokens
  mutate(raw_split = str_split(raw, "\\s+")) %>%
  mutate(
    anc_gene = map_chr(raw_split, 1),
    extant_tokens = purrr::map(raw_split, ~ .x[-1])
  ) %>%
  select(-raw, -raw_split) %>%
  unnest(extant_tokens) %>%
  rename(token = extant_tokens) %>%
  separate(token, into = c("species", "transcript_id"), sep = "\\.", extra = "merge", fill = "right") %>%
  mutate(
    node = str_match(basename(anc_file), "ancGenome\\.(N\\d+)\\.list$")[,2]
  ) %>%
  select(node, CAR, anc_i0, anc_i1, orient, anc_gene, species, transcript_id) %>%
  distinct()

# merge with gene info to get extant chrs
car <- left_join(af,extant_df %>% dplyr::rename(transcript_id = transcript)) %>% 
  mutate(chr = gsub('Chr','',chr))
ord <- car %>% 
  mutate(CARord = as.numeric(gsub('CAR_','',CAR)),
         CHRord = as.numeric(gsub('b','',gsub('a','',chr))))
carord <- ord %>% select(CAR,CARord) %>% distinct %>% arrange(CARord)  
chrord <- ord %>% select(chr,CHRord) %>% distinct %>% arrange(CHRord)
car <- car %>% 
  mutate(chr = factor(chr,levels=chrord$chr),
         CAR = factor(CAR,levels=carord$CAR),
         species = factor(species,levels=c('Artocarpus','Batocarpus','Morus','Antiaris','Ficus')))

# top 20 cars
keep <- car %>% count(CAR) %>% top_n(14) %>% pull(CAR)

default_colors <- c(
  "red3", "chartreuse4", "turquoise2", "gold", "coral2", "olivedrab2", "orange",
  "darkseagreen4", "firebrick4", "royalblue4", "lightsalmon", "darkviolet", "magenta2",
  "mediumaquamarine", "salmon", "paleturquoise2", "darkseagreen1", "khaki1", "thistle2",
  "peachpuff", "lightblue", "skyblue1", "lightgoldenrod3", "wheat1", "darkolivegreen1",
  "lavender", "darkorange", "royalblue4"
)

car_cols <- data.frame(CAR=keep,col=default_colors[1:14])

car_plot <- car %>%
  filter(CAR %in% keep) %>%
  count(species, CAR, chr) %>%
  filter(n > 50) %>%
  ggplot(aes(x = n, y = CAR)) +
  geom_hline(aes(yintercept = CAR), color = "grey70", linewidth = 0.3) +
  geom_point(aes(color = CAR), size = 2) +
  ggrepel::geom_text_repel(aes(label = chr), size = 3, max.overlaps = Inf) +
  scale_color_manual(values=car_cols$col,breaks=car_cols$CAR) +
  facet_grid(. ~ species, scales = "free") +
  theme_bw(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    strip.background = element_rect(fill = "grey95"),
    legend.position = "none"
  ) +
  labs(
    x = "Number of genes mapping to CAR",
    y = "Top 20 ancestral CARs"
  )

car_plot
ggsave('~/symlinks/comp/figures/20260218_CAR_Ancestral_Sequences.pdf',
       height=4,width=7,dpi=300,car_plot)
```



# Inferring Whole Genome Duplication History

Compare variation between Batocarpus and Artocarpus by delineating Artocarpus subgenomes, and then using comparative and phylogenetic approaches. 

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

## Subgenome-divided Orthofinder

This section continues the subgenome-specific evolution analyses, dealing with Artocarpus A/B (HART067), Batocarpus, and Morus.

For most analyses, I will only run with Morus as the outgroup since that's deep enough. However, for proper rooting, I will first run an orthofinder run to get a species tree that is properly rooted with Ficus:

```R
library(ape)
library(ggtree)
library(ggpubr)
library(RColorBrewer)
library(stringr)
library(tidyverse)

# Read in and prune ficus 
t <- read.tree('/project/coffea_pangenome/Artocarpus/Comparative_Paper/orthofinder/Ficus_Subgenome_Divided/OrthoFinder/Results_Mar04/Species_Tree/SpeciesTree_rooted_node_labels.txt')
plot(t)
t2 <- drop.tip(t,'Ficus')
plot(t2)
is.rooted(t2)
is.binary(t2)
is.ultrametric(t2)

# Drop the duplicate HART061
tr_22 <- root(as.phylo(drop.tip(t2,c('HART061_A','HART061_B'))),outgroup = 'Morus',resolve.root = TRUE)
tr_22$tip.label <- gsub('_','',tr_22$tip.label)
is.rooted(tr_22)
is.binary(tr_22)
is.ultrametric(tr_22) # shouldn't be! 

# Subset both A and B subgenome-only trees:
Atips <- grep("B$", tr_22$tip.label, value = TRUE)
Btips <- grep("A$", tr_22$tip.label, value = TRUE)
tr_A <- drop.tip(tr_22, Atips)
tr_B <- drop.tip(tr_22, Btips)

# Reordering so they are all similar: Desired B order
desired_B_order <- tr_A$tip.label %>% str_replace("(A)$", "B") 

# Sanity check: ensure all desired labels exist in B
missing_in_B <- setdiff(desired_B_order, tr_B$tip.label)
missing_in_B

# Constrain rotation (rotate internal nodes) to match A's tip order for B
tr_B_aligned <- rotateConstr(tr_B, desired_B_order)

### also fix n=22 subgenome B 
desired_22_order <- c(str_subset(tr_A$tip.label, "Bato|Morus", negate = TRUE),desired_B_order)
tr_22_aligned <- rotateConstr(tr_22, desired_22_order)

# label nodes for HyPhy
label_internal_nodes <- function(tr) {
  tr$node.label <- paste0("Node", seq.int(from = length(tr$tip.label) + 1,
                                          length.out = tr$Nnode))
  tr
}

tr_A_named  <- label_internal_nodes(tr_A)
tr_B_named  <- label_internal_nodes(tr_B_aligned)
tr_22_named <- label_internal_nodes(tr_22_aligned)

# Verify the plotted order now matches desired_B_order (for present labels)
xt <- 0.1
a1 <- ggtree(tr_A_named)+geom_tiplab()+xlim(0,xt)+geom_nodelab(aes(label=label))
b1 <- ggtree(tr_B)+geom_tiplab()+xlim(0,xt)
b2 <- ggtree(tr_B_named)+geom_tiplab()+xlim(0,xt)+geom_nodelab(aes(label=label))
c1 <- ggtree(tr_22)+geom_tiplab()+xlim(0,xt)
c2 <- ggtree(tr_22_named)+geom_tiplab()+xlim(0,xt)+geom_nodelab(aes(label=label))
ggarrange(b1,b2,a1,c1,c2,labels=c('initB','B','A','initN22','N22'))

# for aligning:
ggarrange(a1,b2,c2,labels=c('A','B','22'),nrow=1)

tr_A_named$node.label <- c('A','B','C','F','G','D','E','H','I','J'); ggtree(tr_A_named)+geom_tiplab()+xlim(0,xt)+geom_nodelab(aes(label=label))
tr_B_named$node.label <- c('A','B','C','F','G','D','E','I','J','H'); ggarrange(ggtree(tr_A_named)+geom_tiplab()+xlim(0,xt)+geom_nodelab(aes(label=label)),ggtree(tr_B_named)+geom_tiplab()+xlim(0,xt)+geom_nodelab(aes(label=label)),nrow=2)
tr_22_named$node.label <- c('A','X','B','C1','F1','G1','D1','E1','H1','I1','J1','C2','F2','G2','D2','E2','I2','J2','H2'); ggarrange(ggtree(tr_22_named)+geom_tiplab()+xlim(0,xt)+geom_nodelab(aes(label=label)),ggtree(tr_A_named)+geom_tiplab()+xlim(0,xt)+geom_nodelab(aes(label=label)),nrow=2)

#Afterwards, confirm!
a1 <- ggtree(tr_A_named)+geom_tiplab()+xlim(0,xt)+geom_nodelab(aes(label=label),geom='label')
b2 <- ggtree(tr_B_named)+geom_tiplab()+xlim(0,xt)+geom_nodelab(aes(label=label),geom='label')
c2 <- ggtree(tr_22_named)+geom_tiplab()+xlim(0,xt)+geom_nodelab(aes(label=label),geom='label')

ggarrange(a1,b2,c2,labels=c('A','B','22'),nrow=1)

# convert to time tree
ggtree(tr_22_named) + geom_tiplab() + geom_nodelab(aes(label=node))

write.tree(tr_A_named,'/project/coffea_pangenome/Artocarpus/Comparative_Paper/subgenome_divided_dnds/tree_A.nwk')
write.tree(tr_B_named,'/project/coffea_pangenome/Artocarpus/Comparative_Paper/subgenome_divided_dnds/tree_B.nwk')
write.tree(tr_22_named,'/project/coffea_pangenome/Artocarpus/Comparative_Paper/subgenome_divided_dnds/tree_22.nwk')

# For cafe5, also convert to ultrametric time tree, use same prior as beast
cal <- data.frame(
  node        = 21,
  age.min     = 74.85,  # Ma
  age.max     = 92.65,  # Ma
  soft.bounds = FALSE
)

dated <- chronos(
  phy        = tr_22_named,
  model      = "correlated",
  lambda     = 1,
  calibration = cal,
  quiet      = FALSE
)

is.rooted(dated)
is.ultrametric(dated)
summary(branching.times(dated))

#### Plot subgenome tree #####
cols <- brewer.pal(5,'Set2')[c(2,1,3,4,5)]
mdraw <-  read_tsv('~/artocarpus_comparative_genomics/samples.txt')
md1 <- mdraw %>% 
  mutate(Accession = case_when(
    Accession == "N97_50" ~ "N9750",
    Accession == "N15_23" ~ "Batocarpus",
    TRUE ~ Accession
  )) %>% 
  dplyr::select(Accession,Group) %>% 
  rbind(.,
        data.frame(
          Accession = 'Morus',
          Group = 'Morus mongolica'
        ))
hapa <- md1 %>% filter(grepl('A. ',Group)) %>% 
  mutate(Accession = paste0(Accession,'A'),
         Haplotype='A')
hapb <- md1 %>% filter(grepl('A. ',Group)) %>% 
  mutate(Accession = paste0(Accession,'B'),
         Haplotype='B')
ogs <- md1 %>% filter(!grepl('A. ',Group)) %>% mutate(Haplotype=Accession)
md <- rbind(hapa,hapb,ogs)
tp <- ggtree(dated, layout = "rectangular")  %<+% md
tp$data <- tp$data %>% mutate(label = ifelse(isTip == TRUE, gsub(' spp.','',gsub('A. ','',Group)), label))
sp_tree <- tp +
  geom_tiplab(hjust = -0.1,size=2)+
  geom_tippoint(aes(fill = Haplotype, shape = Haplotype), size=2)+
  scale_fill_manual(values=cols)+
  scale_shape_manual(values=c(21,21,4,8))+
  xlim(0,max(tp$data$x)*1.3)+
  theme(legend.text = element_text(size = 5),legend.title = element_text(size = 6),
        legend.key.size = unit(0.03, "cm"),    legend.position = 'top')+
  theme_tree2()
sp_tree

ggsave('~/symlinks/comp/figures/20260318_species-time-tree-orthofinder-Subgenome.pdf',sp_tree,height=3,width=4.25)
write.tree(dated,'/project/coffea_pangenome/Artocarpus/Comparative_Paper/cafe5/Subgenome_Divided/resolved_ultrametric_binary_tree.nwk')

```

![image-20260318133526578](C:\Users\justin.merondun\AppData\Roaming\Typora\typora-user-images\image-20260318133526578.png)

Dividing into A and B trees, and aligning labels with node rotation:

![image-20260318134009062](C:\Users\justin.merondun\AppData\Roaming\Typora\typora-user-images\image-20260318134009062.png)

After node alignment:

![image-20260318134937047](C:\Users\justin.merondun\AppData\Roaming\Typora\typora-user-images\image-20260318134937047.png)



## Subgenome-divided dNdS

This workflow splits CDS FASTAs into subgenome (A/B) partitions using chromosome haplotype lists, then identifies orthologous CDS sets via reciprocal best‑hit BLAST using *Morus* as the anchor reference. Per-gene multi-sample CDS alignments are built, filtered, and pruned to matching taxa, and HyPhy (MG94) is run on each gene to estimate branch-specific dN/dS across the tree.

Primary outputs

- Subgenome CDS FASTAs: `*A.fa`, `*B.fa` (plus outgroup FASTAs in `cds_files/`).
- RBH ortholog pairs: `blast/RBH_Morus_<SAMPLE>.txt` and per-gene sequence folders: `genes/<Morus_gene>/*.fa`.
- Per-gene multi-FASTA alignments: `raw/<gene>.fa`.
- HyPhy per-gene results: `hyphy_out/<gene>.tsv` and `hyphy_out/<gene>.tree.nwk` (+ `*.skip.txt` for filtered genes).
- Compiled table: `Node_dNdS_20260406.tsv`.

```bash
# Config
CHR_DIR="chrs"            # dir containing per-chromosome FASTAs
OUTDIR="."                # where to write SAMPLE_A.fa, SAMPLE_B.fa
HAPS=("A.haps" "B.haps")  # haplotype lists to process
SAMPLE_LIST="CompSamples.list"

# Sanity checks
for hap in "${HAPS[@]}"; do
  [[ -s "$hap" ]] || { echo "ERROR: Hap file '$hap' missing or empty." >&2; exit 1; }
done
[[ -d "$CHR_DIR" ]] || { echo "ERROR: Chromosome directory '$CHR_DIR' not found." >&2; exit 1; }
[[ -s "$SAMPLE_LIST" ]] || { echo "ERROR: Sample list '$SAMPLE_LIST' missing or empty." >&2; exit 1; }

# Process each sample
while IFS=$'\r' read -r SAMPLE || [[ -n "$SAMPLE" ]]; do
  # Skip blanks and comments
  [[ -z "$SAMPLE" || "$SAMPLE" =~ ^# ]] && continue

  for hap in "${HAPS[@]}"; do
    # Derive suffix "A" or "B" from filename (before first dot)
    suffix="${hap%%.*}"
    out="${OUTDIR}/${SAMPLE}${suffix}.fa"
    # Truncate/initialize output
    : > "$out"

    echo "Building ${out} from ${hap}…"

    # Read chromosomes in order from hap file
    while IFS=$'\r' read -r CHR || [[ -n "$CHR" ]]; do
      # Skip blanks/comments
      [[ -z "$CHR" || "$CHR" =~ ^# ]] && continue

      chr_file="${CHR_DIR}/${SAMPLE}_${CHR}.fa"
      if [[ -s "$chr_file" ]]; then

          # Rewrite headers so `>SAMPLE_...` becomes `>SAMPLESUFFIX_...` (e.g., HART001A_...)
          awk -v s="$SAMPLE" -v suf="$suffix" '
            BEGIN { OFS="" }
            /^>/ {
              # Only modify headers that begin with the exact sample ID
              if ($0 ~ "^>" s) {
                sub("^>" s, ">" s suf);
              }
              print; next
            }
            { print }
          ' "$chr_file" >> "$out"

      else
        echo "WARNING: Missing file: ${chr_file} (skipping)" >&2
      fi
    done < "$hap"

    echo "Done: ${out}"
  done
  
done < "$SAMPLE_LIST"
```

Copy the files:

```bash
cp *A.fa *B.fa Batocarpus.fa Morus.fa ~/symlinks/comp/subgenome_divided_dnds/cds_files
```

This will output:

```
ls cds_files/*fa
cds_files/Batocarpus.fa  cds_files/HART027A.fa  cds_files/HART058B.fa  cds_files/HART061A.fa  cds_files/HART062B.fa  cds_files/HART067A.fa  cds_files/HART068B.fa  cds_files/N9750B.fa
cds_files/HART001A.fa    cds_files/HART027B.fa  cds_files/HART060A.fa  cds_files/HART061B.fa  cds_files/HART063A.fa  cds_files/HART067B.fa  cds_files/Morus.fa
cds_files/HART001B.fa    cds_files/HART058A.fa  cds_files/HART060B.fa  cds_files/HART062A.fa  cds_files/HART063B.fa  cds_files/HART068A.fa  cds_files/N9750A.fa

head cds_files/HART061A.fa
>HART061A_000001-R1
ATGATCATGTCTTCAAAAGGGTGTTTAGAGGAGATGGGAATATCTTCAACTAATATCAGT
GATGGTGGGAAAAATTGCTATAGAGGCCATTGGAGACCTGCGGAAGACGAGAAACTCCGA
```

1. This takes a reference sample (Morus), performs reciprocal best‑hit BLAST searches against all other samples, extracts matching CDS sequences and outputs them into a directory named after the Morus gene in `/genes/` 

```bash
#!/bin/bash

#SBATCH --time=2-00:00:00   
#SBATCH --nodes=1  
#SBATCH --cpus-per-task=16
#SBATCH --mem=32Gb
#SBATCH --partition=ceres
#SBATCH --account=coffea_pangenome

#module load miniconda
#source activate isoseq_ann

JOBS=16
WD=/project/coffea_pangenome/Artocarpus/Comparative_Paper/subgenome_divided_dnds
cd ${WD}

# Submit with sample 
B="${1:?usage: $0 <SAMPLE>}"
first=0
[ "$B" = "Batocarpus" ] && first=1

# dNdS Across tree 
A="Morus"
Af="cds_files/${A}.fa"

# output dirs
mkdir -p blast db genes

echo "Blasting ${A} against ${B}"
Bf="cds_files/${B}.fa"

# Make blast dbs
if [ ! -f "db/${B}.nhr" ]; then
makeblastdb -in $Bf -dbtype nucl -out db/${B}
fi

# blast each against the other
blastn -query $Af -db db/${B} -out blast/${A}_vs_${B}.tsv -outfmt "6 qseqid sseqid pident length evalue bitscore" -max_target_seqs 1 -evalue 1e-5
blastn -query $Bf -db db/${A} -out blast/${B}_vs_${A}.tsv -outfmt "6 qseqid sseqid pident length evalue bitscore" -max_target_seqs 1 -evalue 1e-5

# sort for best hits 
sort -k1,1 -k6,6nr blast/${A}_vs_${B}.tsv \
  | awk -F'\t' '!seen[$1]++ {print $1"\t"$2}' > blast/best_${A}_to_${B}.txt
sort -k1,1 -k6,6nr blast/${B}_vs_${A}.tsv \
  | awk -F'\t' '!seen[$1]++ {print $1"\t"$2}' > blast/best_${B}_to_${A}.txt
awk 'NR==FNR {a[$1]=$2; next} {if (a[$2]==$1) print $2"\t"$1}' \
  blast/best_${A}_to_${B}.txt \
  blast/best_${B}_to_${A}.txt \
  > blast/RBH_${A}_${B}.txt

# Export for parallel subshells
export WD Af Bf A B first

# SELF per ida, only if first==1; run once per unique ida in parallel ---
if [ "${first:-0}" -eq 1 ]; then
  cut -f1 "blast/RBH_${A}_${B}.txt" | sort -u | \
  parallel --jobs ${JOBS} '
    ida={}
    mkdir -p ${WD}/genes/${ida}
    self="${WD}/genes/${ida}/${A}.fa"
    # overwrite or create; trimming terminal stop codon
    samtools faidx "'"$Af"'" "$ida" | sed -E "s/(TAA|TAG|TGA)$//" > "$self"
  '
fi

# Extract each idb from Bf in parallel, one file per pair ---
parallel --jobs ${JOBS} --colsep '\t' '
  ida={1}; idb={2}
  mkdir -p ${WD}/genes/${ida}
  out="${WD}/genes/${ida}/${B}.fa"
  samtools faidx "'"$Bf"'" "$idb" | sed -E "s/(TAA|TAG|TGA)$//" > "$out"
' :::: "blast/RBH_${A}_${B}.txt"
```

Sanity, check number of genes per sample:

```bash
find genes -maxdepth 2 -type f -printf "%f\n" \
    | grep -o -f <(sed 's/$/\\.fa/' AllSamples.list) \
    | sort | uniq -c
    
  13487 Batocarpus.fa
  11431 HART001A.fa
  10163 HART001B.fa
  11106 HART027A.fa
   9849 HART027B.fa
  11408 HART058A.fa
   9644 HART058B.fa
  10814 HART060A.fa
   9305 HART060B.fa
  11090 HART062A.fa
   9781 HART062B.fa
  11508 HART063A.fa
  10541 HART063B.fa
  11289 HART067A.fa
  10117 HART067B.fa
  10650 HART068A.fa
   8839 HART068B.fa
  13487 Morus.fa
  10745 N9750A.fa
   9242 N9750B.fa
```

Concatenate the alignents:

```bash
#!/bin/bash

#SBATCH --time=2-00:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=64Gb
#SBATCH --partition=ceres
#SBATCH --account=coffea_pangenome

set -euo pipefail

WD=/project/coffea_pangenome/Artocarpus/Comparative_Paper/subgenome_divided_dnds
SRC="${WD}/genes"
OUT="${WD}/raw"
JOBS=48

mkdir -p "$OUT"

export SRC OUT

find "$SRC" -mindepth 1 -maxdepth 1 -type d | \
parallel -j ${JOBS} '
  d={}
  gene=$(basename "$d")
  out="${OUT}/${gene}.fa"
  tmp=$(mktemp "${OUT}/.${gene}.XXXXXX")

  # require Morus.fa and Batocarpus.fa to exist and be non-empty
  [ -s "${d}/Morus.fa" ] || exit 0
  [ -s "${d}/Batocarpus.fa" ] || exit 0

  # concatenate all non-empty fasta files in sorted filename order
  find "$d" -maxdepth 1 -type f -name "*.fa" -size +0c -printf "%f\n" \
    | LC_ALL=C sort \
    | while read -r f; do
        cat "${d}/${f}" >> "$tmp"
      done

  # only keep if something was written
  if [ -s "$tmp" ]; then
    mv "$tmp" "$out"
    echo "Built ${gene}.fa"
  else
    rm -f "$tmp"
  fi
'
```

Subset those tips from the tree, and run hyphy for dnds estimation:

```bash
find raw -type f -print0 | xargs -0 realpath > genes_all.list
split -n l/8 -d genes_all.list genes_all.list_
ls *list_* | xargs -I {} echo sbatch -J dnds_{} 03_NodeSpecific_dNdS.sh {}
ls *list_* | grep -v '00' | xargs -I {} echo sbatch -J dnds_{} 03_NodeSpecific_dNdS.sh {}
```

Run this:

```bash
#!/bin/bash

#SBATCH --time=3-00:00:00    
#SBATCH --cpus-per-task=48
#SBATCH --mem=64Gb
#SBATCH --partition=ceres
#SBATCH --account=coffea_pangenome

# module load miniconda
# source activate isoseq_ann

JOBS=44
LIST="${1:?usage: $0 <genes_list_file>}"
WD=/project/coffea_pangenome/Artocarpus/Comparative_Paper/subgenome_divided_dnds
TREE="${WD}/tree.nwk"
GENEDIR="${WD}/raw"

mkdir -p ${WD}/hyphy ${WD}/hyphy_out

cd ${WD}

echo "Working on ${LIST}"

# Ensure required tools exist
for tool in macse pal2nal.pl hyphy jq parallel; do
  command -v "$tool" >/dev/null 2>&1 || { echo "ERROR: $tool not found in PATH"; exit 1; }
done

SCR=${SLURM_TMPDIR:-/tmp}
export WD TREE RUN GENEDIR SCR

process_gene() {
  local fa="$1"           
  local gene
  gene="$(basename "${fa%.fa}")" 

  # Work in a per-gene dir
  (
    outdir="${SCR}/hyphy/${gene}"
    mkdir -p "$outdir"
    cd "$outdir"
    echo "Working on ${gene}"

    # 1) Align with MACSE
    macse -prog alignSequences \
      -seq "${fa}" \
      -out_NT aln_NT.fasta \
      -out_AA aln_AA.fasta > macse.log 2>&1

    # 2) Clean alignment
    macse -prog exportAlignment \
      -align aln_NT.fasta \
      -codonForExternalFS --- \
      -codonForFinalStop --- \
      -codonForInternalFS --- \
      -codonForInternalStop --- \
      -charForRemainingFS - \
      -out_NT aln_NT.clean.fasta 2>&1

    # 3) Normalization, ensure underscores stripped
    sed 's/_.*//g' aln_NT.clean.fasta > aln_NT.clean.fa

    # 4) Trim gappy codon positions
    clipkit aln_NT.clean.fa -o aln_NT.clipkit.fa --codon -m kpic > clipkit.log 2>&1

    # 5) Extract present taxa from fasta headers
    grep '^>' aln_NT.clipkit.fa | sed 's/^>//' > tips.txt

    ntips=$(wc -l < tips.txt)
    if [ "$ntips" -lt 8 ]; then
      echo "Skipping ${gene}: fewer than 8 taxa after filtering" > skip.txt
      cp skip.txt "${WD}/hyphy_out/${gene}.skip.txt"
      exit 0
    fi

    # 6) Prune master tree to present taxa
    cp "${TREE}" tree.full.nwk
    nw_prune -v tree.full.nwk $(tr '\n' ' ' < tips.txt) > tree.pruned.nwk 2> prune.log

    # make sure pruning succeeded and tree is non-empty
    if [ ! -s tree.pruned.nwk ]; then
      echo "Skipping ${gene}: pruned tree is empty" > skip.txt
      cp skip.txt "${WD}/hyphy_out/${gene}.skip.txt"
      exit 0
    fi

    # 7) Run HyPhy on pruned tree
    hyphy "${WD}/FitMG94.bf" \
      --alignment aln_NT.clipkit.fa \
      --tree tree.pruned.nwk \
      --output fitmg94.json \
      --type local \
      --lrt Yes > hyphy.log 2>&1

    # 8) Extract branch results
    jq -r '
      .input["number of sites"] as $L
      | ["Branch","dN","dS","N_exp","S_exp","LB","MLE","UB","LRT","FDR"],
        (.["branch attributes"]["0"]
          | to_entries
          | map([
              .key,
              (.value["dN"] // "NA"),
              (.value["dS"] // "NA"),
              ((.value["dN"] // 0) * $L),
              ((.value["dS"] // 0) * $L),
              (.value["Confidence Intervals"]["LB"] // "NA"),
              (.value["Confidence Intervals"]["MLE"] // "NA"),
              (.value["Confidence Intervals"]["UB"] // "NA"),
              (try .value["LRT"]["LRT"] catch "NA"),
              (try .value["LRT"]["FDR"] catch "NA")
            ])
        )[] | @tsv
    ' fitmg94.json > fitmg94.tsv

    # 9) Build mapping of terminal branches
    grep '^>' aln_NT.clipkit.fa \
      | sed 's/^>//' \
      | awk 'BEGIN{OFS="\t"} {print $1, $1}' > mapping.tsv

    grep '^Node' fitmg94.tsv | awk 'BEGIN{OFS="\t"} {print $1, "NA"}' >> mapping.tsv || true
    sort -k1,1 mapping.tsv > mapping.sorted.tsv

    # 10) Merge metadata into final results
    awk -F'\t' -v OFS='\t' -v gene=${gene} '
      NR==FNR {map[$1]=$2; next}
      FNR==1 {print $0, "Gene", "Mapping"; next}
      {print $0, gene, (map[$1] ? map[$1] : "NA")}
    ' mapping.sorted.tsv fitmg94.tsv > merged.tsv

    cp merged.tsv "${WD}/hyphy_out/${gene}.tsv"
    cp tree.pruned.nwk "${WD}/hyphy_out/${gene}.tree.nwk"
  )
}

export -f process_gene

# Run in parallel
parallel --will-cite -j ${JOBS} process_gene :::: "${LIST}"
```

Check after:

```
find hyphy_out/ -type f -name '*tsv' | wc -l
find hyphy_out/ -type f -name '*skip*' | wc -l
```

Compile:

```bash
# Compile results From the directory containing hyphy_out/
awk -v OFS='\t' '
  FNR==1 {
    # Extract gene from filename (strip directory and .tsv)
    gene = FILENAME
    sub(/.*\//, "", gene)
    sub(/\.tsv$/, "", gene)
  }

  # Print header once, from the first file, and add Gene column
  FNR==1 && NR==1 { print $0, "Gene"; next }

  # Skip headers in subsequent files
  FNR==1 { next }

  # Print data line + gene ID
  NF { print $0, gene }
' hyphy_out/*.tsv > Node_dNdS_20260406.tsv
```

#### Summarize & Plot Tree 

Treads the per-branch HyPhy dN/dS results, assigns branches to subgenome A/B (shared vs unique gene sets), filters unstable estimates, and then compares selection signals between A and B using Fisher tests and tree-based visualizations.

Outputs

- Filtered results table: `Node_dNdS_20260406-RInput.tsv`
- Purifying-selection contrast plot: `20260406_PurifiedBranches_Log2Odds.pdf`
- Tree plots (purifying selection): `20260406_Subgenome_dNdS_BranchColorsPURIFYING.pdf`
- Candidate gene lists: `TopGenes_Puri_20260406.tsv`

```R
# Subgenome-specific dNdS along the tree 
setwd('/project/coffea_pangenome/Artocarpus/Comparative_Paper/subgenome_divided_dnds')
library(tidyverse)
library(scales)
library(stringr)
library(data.table)
library(ggrepel)
library(ape)
library(ggtree)
library(ggpubr)
library(viridis)
library(meRo)

# read in dnds
dnds <- read_tsv('Node_dNdS_20260406.tsv')

# Read in trees
tr <- read.tree('tree.nwk')
tr_A  <- read.tree('tree_A.nwk')
tr_B  <- read.tree('tree_B.nwk')
tr_22 <- read.tree('tree_22.nwk')

gene_subset_tips <- dnds %>%
  mutate(
    tip_side = case_when(
      str_detect(Branch, "^(HART\\d+A|N9750A)$") ~ "A",
      str_detect(Branch, "^(HART\\d+B|N9750B)$") ~ "B",
      TRUE ~ NA_character_
    )
  ) %>%
  group_by(Gene) %>%
  summarise(
    has_A = any(tip_side == "A", na.rm = TRUE),
    has_B = any(tip_side == "B", na.rm = TRUE),
    n_A_tips = n_distinct(Branch[tip_side == "A"],na.rm=TRUE),
    n_B_tips = n_distinct(Branch[tip_side == "B"],na.rm=TRUE),
    Subset = case_when(
      has_A & has_B ~ "22",
      has_A & !has_B ~ "A",
      !has_A & has_B ~ "B",
      TRUE ~ NA_character_
    ),
    .groups = "drop"
  )
gene_subset_tips

# merge
wstats <- left_join(dnds,gene_subset_tips) %>% 
  mutate(
    ID = # single-letter 1/2 -> keep letter; non-single-letter ending A/B -> strip tail
      case_when(
        str_detect(Branch, "^[A-Z][12]$") ~ str_replace(Branch, "[12]$", ""),
        str_detect(Branch, "^[A-Za-z]+\\d+[AB]$") ~ str_replace(Branch, "[AB]$", ""),
        TRUE ~ Branch
      ),
    # Subg: A if ends in A or 1; B if ends in B or 2; else NA
    Subg = case_when(
      str_detect(Branch, "A$|1$") ~ "A",
      str_detect(Branch, "B$|2$") ~ "B",
      TRUE ~ NA_character_
    ),
    Localization = ifelse(
      Subset == "22",
      paste0('Shared_',Subg),
      paste0('Unique_',Subg)
    ),
    Localization = factor(Localization,
                          levels = c("Shared_A", "Shared_B", "Unique_A","Unique_B"))
  )

# Now apply some sanity thresholds like # syn mutations 
omega_cap <- 5    # ceiling for dnds
min_syn <- 1; max_syn <- 5000     # require at least 1 syn mutation for stability 
max_non <- 5000
max_ub <- 100 # maximum omega upper bound, above indicates model instability - these will be REMOVED
wstats_clean <- wstats %>%
  # Remove extrmee 
  filter(S_exp >= min_syn & S_exp < max_syn & N_exp < max_non & UB < max_ub) %>% 
  mutate(omega = pmin(MLE, omega_cap),
         pos = FDR < 0.10 & MLE > 1,
         puri = FDR < 0.10 & MLE < 1) %>% 
  filter(!grepl('^A$|^B$|Morus|Batocarpus',Branch)) 
hist(wstats_clean$omega)
write_tsv(wstats_clean,'Node_dNdS_20260406-RInput.tsv')
wstats_clean <- read_tsv('Node_dNdS_20260406-RInput.tsv')

# First, aggregate: compare # purified genes to total for each branch 
# Build per-ID counts
counts <- wstats_clean %>%
  mutate(
    pos  = as.integer(FDR < 0.10 & MLE > 1),
    puri = as.integer(FDR < 0.10 & MLE < 1)
  ) %>%
  group_by(ID, Localization) %>%
  summarise(pos = sum(pos), puri = sum(puri), total = n(), .groups = "drop")
counts

res_df <- tibble()
pairs_list <- list(
  c("Shared_A","Shared_B"),
  c("Shared_A","Unique_A"),
  c("Shared_B","Unique_B"),
  c("Unique_A","Unique_B")
)


for (ID_i in unique(counts$ID)) {
  for (pair in pairs_list) {
    cat("Working on:", ID_i, "contrast:", paste(pair, collapse = " vs "), "\n")
    
    # Extract rows for this ID and the two localizations; fill missing with zeros
    df_id <- counts %>%
      filter(ID == ID_i, Localization %in% pair) %>%
      complete(ID = ID_i, Localization = pair,
               fill = list(pos = 0L, puri = 0L, total = 0L)) %>%
      arrange(Localization)
    L1 <- pair[1]; L2 <- pair[2]
    
    ## PURIFYING
    puri1 <- df_id$puri[df_id$Localization == L1]
    tot1  <- df_id$total[df_id$Localization == L1]
    puri2 <- df_id$puri[df_id$Localization == L2]
    tot2  <- df_id$total[df_id$Localization == L2]
    
    not1  <- tot1 - puri1
    not2  <- tot2 - puri2
    
    # Build 2x2 table
    puri_mat <- matrix(c(puri1, not1, puri2, not2), nrow = 2, byrow = TRUE,
                       dimnames = list(Localization = c(L1, L2), Outcome = c("puri", "not")))
    ft_puri <- fisher.test(puri_mat)
    puri_p <- unname(ft_puri$p.value)
    puri_or <- unname(ft_puri$estimate) # odds ratio
    
    res_df <- bind_rows(res_df, tibble(
      ID = ID_i, L1 = L1, L2 = L2, outcome = "puri",
      a = puri1, b = not1, c = puri2, d = not2,
      odds_ratio = puri_or, p = puri_p
    ))
    
    ## POSITIVE
    pos1 <- df_id$pos[df_id$Localization == L1]
    tot1 <- df_id$total[df_id$Localization == L1]  # re-use variables
    pos2 <- df_id$pos[df_id$Localization == L2]
    tot2 <- df_id$total[df_id$Localization == L2]
    
    not1 <- tot1 - pos1
    not2 <- tot2 - pos2
    
    pos_mat <- matrix(c(pos1, not1, pos2, not2), nrow = 2, byrow = TRUE,
                      dimnames = list(Localization = c(L1, L2), Outcome = c("pos", "not")))
    ft_pos <- fisher.test(pos_mat)
    pos_p <- unname(ft_pos$p.value)
    pos_or <- unname(ft_pos$estimate) # odds ratio
    
    res_df <- bind_rows(res_df, tibble(
      ID = ID_i, L1 = L1, L2 = L2, outcome = "pos",
      a = pos1, b = not1, c = pos2, d = not2,
      odds_ratio = pos_or, p = pos_p
    ))
  }
}

# Add Holm-adjusted p-values and significance per contrast & outcome
res_df <- res_df %>%
  group_by(outcome, L1, L2) %>%
  mutate(p_adj = p.adjust(p, method = "holm"),
         sig   = !is.na(p_adj) & p_adj < 0.05) %>%
  ungroup()

idord <- c('HART001','HART063','HART067','HART058','HART062','HART027','HART068','N9750','HART060','J','I','H','G','F','E','D','C')
idord <- c('A. altilis','A. mariannensis','A. camansi','A. rigidus','A. odoratissimus','A. heterophyllus','A. nitidus','A. lacucha','A. dadah','J','I','H','G','F','E','D','C')

# import md
md <-  read_tsv('~/artocarpus_comparative_genomics/samples.txt') %>% dplyr::select(ID=Accession,Group) %>% mutate(ID = gsub('_','',ID))

# Add Fisher exact 95% CI for purifying rows (no continuity correction)
res_puri <- res_df %>%
  filter(outcome == "puri") %>%
  rowwise() %>%
  mutate(
    ft = list(fisher.test(matrix(c(a, b, c, d), nrow = 2, byrow = TRUE))),  # exact
    ci_lo = ft$conf.int[1],
    ci_hi = ft$conf.int[2]
  ) %>%
  ungroup() %>%
  dplyr::select(-ft) %>% 
  left_join(.,md) %>% 
  mutate(Group = factor(ifelse(is.na(Group),ID,Group),levels=rev((idord))),
         contrast = paste(L1, "vs", L2))

puri_plot <- res_puri %>% 
  ggplot(aes(x = odds_ratio, y = Group)) +
  geom_errorbarh(aes(xmin = ci_lo, xmax = ci_hi), height = 0.25, color = "grey55") +
  geom_point(aes(color = sig), size = 2.4) +
  geom_text(aes(label = ifelse(sig, "*", "")), nudge_x = 0.05, size = 5) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  scale_color_manual(values = c(`TRUE` = "#D95F02", `FALSE` = "#1B9E77"), guide = "none") +
  #scale_x_log10() +
  facet_grid(.~ contrast) +
  labs(x = "Odds ratio",y = NULL) +
  #coord_cartesian(xlim=c(-1,1))+
  theme_bw(base_size = 12) +
  theme(strip.text = element_text(face = "bold"))
puri_plot

## log2 scale
plt_df <- res_puri %>%
  mutate(
    x     = log2(odds_ratio),
    xmin  = log2(ci_lo),
    xmax  = log2(ci_hi)
  ) 

# Axis breaks and labels: -2,-1,0,1,2 => L2×4, L2×2, equal, L1×2, L1×4
brks <- c(-2, -1, 0, 1, 2)
labs <- c("×4","×2","equal","×2","×4")

puri_plot <- ggplot(plt_df, aes(x = x, y = Group)) +
  geom_errorbarh(aes(xmin = xmin, xmax = xmax), height = 0.25, color = "grey55") +
  geom_point(aes(color = sig), size = 2.8) +
  geom_text(aes(label = ifelse(sig, "*", "")), nudge_x = 0.15, size = 5) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_color_manual(values = c(`TRUE` = "#D95F02", `FALSE` = "#1B9E77"), guide = "none") +
  scale_x_continuous(breaks = brks, labels = labs) +
  labs(x = "log2(OR)",y = NULL) +
  facet_grid(.~ contrast) +
  theme_bw(base_size = 8)
puri_plot
ggsave('~/symlinks/comp/figures/20260406_PurifiedBranches_Log2Odds.pdf',puri_plot,height=3,width=5)

##### Second, compare exact homeologs #####
homeolog_diff <- wstats_clean %>% 
  # filtered for matched IDs (both homeologs) from the same genes 
  group_by(Gene,ID) %>% 
  mutate(count = n()) %>% 
  filter(count > 1) %>% 
  dplyr::select(Gene,ID,Subset,Subg,omega,pos,puri) %>% 
  # calc log2 diff 
  pivot_wider(names_from = Subg, values_from = c(omega,pos,puri)) %>% 
  mutate(log2_omega_ratio = log2((omega_A + 1e-6)/(omega_B + 1e-6)))

# filter for very small omega
# drop near-zeros
min_omega <- 1e-3 
homeolog_diff2 <- homeolog_diff %>%
  filter(omega_A >= min_omega, omega_B >= min_omega) %>%
  mutate(log2_omega_ratio = log2(omega_B/omega_A))

# summary stats 
idord <- c('altilis','mariannensis','camansi','rigidus','odoratissimus','heterophyllus','nitidus','lacucha','dadah','J','I','H','G','F','E','D','C')
homeo_puri <- homeolog_diff2 %>% 
  group_by(ID) %>% 
  sum_stats(log2_omega_ratio) %>%   
  left_join(.,md %>% mutate(Group = gsub('A. ','',Group))) %>%
  mutate(Group = factor(ifelse(is.na(Group),ID,Group),levels=rev((idord))))

homeo_puri_plot <- homeo_puri %>% 
  ggplot(aes(x = mean, y = Group)) +
  geom_errorbarh(aes(xmin = conf_low, xmax = conf_high), height = 0.25, color = "grey55") +
  geom_point() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  coord_cartesian(xlim=c(-0.65,0.65))+
  xlab('log2(B/A) w')+ylab('')+
  theme_bw(base_size = 12) +
  theme(strip.text = element_text(face = "bold"))
homeo_puri_plot
ggsave('~/symlinks/comp/figures/20260414_PurifiedBranches.pdf',homeo_puri_plot,height=3,width=2.5)

# log odds ratio & stats
min_omega <- 1e-3

homeolog_bin <- homeolog_diff %>%
  filter(omega_A >= min_omega, omega_B >= min_omega) %>%
  mutate(
    B_more_constrained = omega_B < omega_A,
    A_more_constrained = omega_A < omega_B
  )

log2or_by_ID <- homeolog_bin %>%
  group_by(ID) %>%
  summarise(
    b = sum(B_more_constrained, na.rm = TRUE),
    a = sum(A_more_constrained, na.rm = TRUE),
    ties = sum(omega_A == omega_B, na.rm = TRUE),
    n_discord = a + b,
    .groups = "drop"
  ) %>%
  mutate(
    OR = (b + 0.5) / (a + 0.5),
    log2OR = log2(OR),
    se = sqrt(1/(b+0.5) + 1/(a+0.5)),
    ci_lo = log2OR - 1.96 * se,
    ci_hi = log2OR + 1.96 * se
  ) %>%
  rowwise() %>%
  mutate(
    p = if (n_discord > 0) binom.test(b, n_discord, p = 0.5)$p.value else NA_real_
  ) %>%
  ungroup() %>%
  mutate(
    p_adj = p.adjust(p, method = "BH"),
    sig = !is.na(p_adj) & p_adj < 0.05
  ) %>%
  left_join(md %>% mutate(Group = gsub('A. ','',Group)), by = "ID") %>%
  mutate(Group = factor(ifelse(is.na(Group), ID, Group), levels = rev(idord)))

brks <- log2(c(0.5, 2/3, 1, 1.5, 2))
labs <- c("0.5× odds", "0.67× odds", "equal", "1.5× odds", "2× odds")
homeo_puri_plot_or <- ggplot(log2or_by_ID, aes(x = log2OR, y = Group)) +
  geom_errorbarh(aes(xmin = ci_lo, xmax = ci_hi), height = 0.25, color = "grey55") +
  geom_point(aes(color = sig), size = 2.8) +
  geom_text(aes(label = ifelse(sig, "*", "")), nudge_x = 0.10, size = 4) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_color_manual(values = c(`TRUE` = 'darkviolet', `FALSE` = 'grey30'), guide = "none") +
  labs(x = "log2(OR) (B more constrained vs A)", y = NULL) +
  theme_bw(base_size = 8)+
  scale_x_continuous(breaks = brks, labels = labs,limits=c(-1.1,1.1))
homeo_puri_plot_or
ggsave('~/symlinks/comp/figures/20260414_logOR_Homeolog_PurifiedBranches.pdf',homeo_puri_plot,height=3,width=2.5)


# extract top 
puri_cands <- homeolog_diff2 %>% 
  filter(ID == 'C' | ID == 'J') %>% 
  ungroup %>% 
  mutate(cand = ifelse(ID == 'C' & log2_omega_ratio > 0.5,'candidate',
                       ifelse(ID == 'J' & log2_omega_ratio < 0.5, 'candidate','other'))) %>% 
  filter(cand != 'other') %>% filter(puri_A == TRUE | puri_B == TRUE) %>% mutate(cand = 'purifying')
pos_cands <- homeolog_diff2 %>% 
  filter(ID == 'C' | ID == 'J') %>% 
  ungroup %>% 
  mutate(cand = ifelse(ID == 'C' & log2_omega_ratio > 0.5,'candidate',
                       ifelse(ID == 'J' & log2_omega_ratio < 0.5, 'candidate','other'))) %>% 
  filter(cand != 'other') %>% filter(pos_A == TRUE | pos_B == TRUE) %>% mutate(cand = 'positive')
homeo_cands <- rbind(puri_cands,pos_cands)
write.table(homeo_cands,file='TopHomeologs_20260414.tsv',quote=F,sep='\t',row.names=F)

##### Plot per-branch variation #####
# per-branch summaries
summarize_branch_metrics <- function(df) {
  df %>%
    mutate(
      ci_w = pmax(0, UB - LB),
      pos = FDR < 0.10 & MLE > 1,
      puri = FDR < 0.10 & MLE < 1,
    ) %>%
    group_by(Subset, Subg, ID, Branch) %>%
    summarise(
      med_omega = median(omega, na.rm = TRUE),
      pos_frac  = mean(pos,  na.rm = TRUE),
      total_pos = sum(pos, na.rm = TRUE),
      puri_frac = mean(puri, na.rm = TRUE),
      total_puri = sum(puri, na.rm = TRUE),
      ci_w_med  = median(ci_w, na.rm = TRUE),
      n_genes   = n(),
      .groups   = "drop"
    ) 
}

# apply 
branch_summ <- summarize_branch_metrics(wstats_clean)
branch_summ 

##### Plot dNdS on the tree: lwd and color POSITIVE #####
# Limits for the plot 
omega_limits <- range(branch_summ$med_omega, na.rm = TRUE)
pos_limits <- range(branch_summ$pos_frac, na.rm = TRUE)
probs <- c(0.10, 0.30, 0.70, 0.90)
breaks <- round(quantile(branch_summ$med_omega, probs = probs, na.rm = TRUE, names = FALSE),2)
colors <- viridis(length(breaks) + 1, option = "plasma")

# now throw that on the trees
plot_branch_summary_tree <- function(tr, branch_summ, subset_id, title = NULL, expand_x = 1.3) {
  # Build tree data and join branch summaries
  td <- fortify(tr) %>% 
    dplyr::left_join(
      branch_summ %>% 
        dplyr::filter(Subset == subset_id) %>% 
        mutate(Branch = ifelse(Subset != '22',gsub('[1|2]$','',Branch),Branch)),
      by = c("label" = "Branch"))
  
  # Base tree
  p <- ggtree(tr) %<+% td +
    geom_tree(aes(color = med_omega, linewidth = pos_frac)) +
    geom_tiplab(size = 3.0) +
    theme_tree() +
    labs(title = paste0("Subset ", subset_id)) +
    scale_color_stepsn(
      name    = "Median dN/dS",
      colors  = colors,
      breaks  = breaks,
      limits  = omega_limits,
      na.value = "grey80"
    ) +
    scale_linewidth_continuous(
      name  = "Frac FDR<0.1 & dnds>1",
      range = c(0.5, 2.5),limits=pos_limits)
  
  # Expand x-axis to create room for tip labels
  xmax <- max(p$data$x, na.rm = TRUE)
  p <- p + ggplot2::xlim(0, xmax * expand_x)
  
  return(p)
}


pA  <- plot_branch_summary_tree(tr_A,  branch_summ, subset_id = "A")
pB  <- plot_branch_summary_tree(tr_B,  branch_summ, subset_id = "B")
p22 <- plot_branch_summary_tree(tr_22, branch_summ, subset_id = "22")
allplot <- ggarrange(pA,pB,p22,nrow=1,common.legend = TRUE)
allplot
ggsave('~/symlinks/comp/figures/20260319_Subgenome_dNdS_BranchColors.pdf',
       allplot,height=4,width=8)


##### Plot dNdS on the tree: lwd and color PURIFYING #####

# read in md
md <- read_tsv('~/artocarpus_comparative_genomics/samples.txt') %>% dplyr::select(ID=Accession,Species,Group) %>% mutate(ID = gsub('_','',ID),Group = gsub('A. ','',Group))

# Limits for the plot, in case we want to do discrete bins  
pur_limits <- range(branch_summ$puri_frac, na.rm = TRUE)
probs <- c(0.05, 0.30, 0.70, 0.95)
breaks <- round(quantile(branch_summ$puri_frac, probs = probs, na.rm = TRUE, names = FALSE),2)
colors <- viridis(length(breaks) + 1, option = "plasma")

# now throw that on the trees
plot_branch_summary_tree <- function(tr, branch_summ, subset_id, title = NULL, expand_x = 1.3) {
  # Build tree data and join branch summaries
  td <- fortify(tr) %>% 
    dplyr::left_join(
      branch_summ %>% 
        dplyr::filter(Subset == subset_id) %>% 
        mutate(Branch = ifelse(Subset != '22',gsub('[1|2]$','',Branch),Branch)),
      by = c("label" = "Branch"))
  
  # Base tree
  trm <- left_join(td,md)
  p <- ggtree(tr) %<+% trm +
    geom_tree(aes(color = puri_frac, linewidth = puri_frac)) +
    geom_tiplab(aes(label=Group),size = 3.0) +
    theme_tree() +
    labs(title = paste0("Subset ", subset_id)) +
    scale_color_viridis(
      option = 'plasma',
      name = "Fraction Purifying"
    ) +
    scale_linewidth_continuous(
      name  = "Frac FDR<0.1 & dnds<1",
      range = c(0.5, 2.5),limits=pur_limits)
  
  # Expand x-axis to create room for tip labels
  xmax <- max(p$data$x, na.rm = TRUE)
  p <- p + ggplot2::xlim(0, xmax * expand_x)
  
  return(p)
}


pA  <- plot_branch_summary_tree(tr_A,  branch_summ, subset_id = "A")
pB  <- plot_branch_summary_tree(tr_B,  branch_summ, subset_id = "B")
p22 <- plot_branch_summary_tree(tr_22, branch_summ, subset_id = "22")
allplot <- ggarrange(pA,pB,p22,nrow=1,common.legend = TRUE)
allplot
ggsave('~/symlinks/comp/figures/20260406_Subgenome_dNdS_BranchColorsPURIFYING.pdf',
       allplot,height=2.75,width=6)

# confirm
library(meRo)
wstats_clean %>% group_by(Subset,Branch) %>% sum_stats(puri) %>% filter(grepl('C',Branch))

# Subset Branch  mean   min   max    sd      se median   iqr conf_low conf_high
# <chr>  <chr>  <dbl> <int> <int> <dbl>   <dbl>  <dbl> <dbl>    <dbl>     <dbl>
#   1 22     C1     0.321     0     1 0.467 0.00543      0     1    0.310     0.331
# 2 22     C2     0.663     0     1 0.473 0.00504      1     1    0.653     0.673
# 3 A      C1     0.385     0     1 0.487 0.0122       0     1    0.361     0.409
# 4 B      C2     0.699     0     1 0.459 0.0154       1     1    0.668     0.729

pA
pB
ggtree(tr_A)+geom_tiplab()+geom_nodelab(aes(label=label),geom='label')
td <- fortify(tr_A)
td %>% filter(label %in% branch_summ$Branch) %>% count(isTip)

##### A vs B: branch‑wise comparison plots #####
#### Positive selection A/B
subg_sel <- counts %>%
  left_join(md) %>% 
  mutate(pos_frac = pos/total, puri_frac = puri/total) %>% 
  separate(Localization, into = c('Subset','Subg'),remove=TRUE) %>% 
  pivot_wider(names_from=Subg, id_cols = c(ID,Group,Subset), values_from = puri_frac) %>% 
  filter(Subset == 'Shared') %>% 
  mutate(Group = ifelse(is.na(Group),ID,Group)) %>% 
  ggplot(aes(x = A, y = B, col= B - A) ) +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted") +
  geom_point(aalpha = 0.95) +
  geom_text_repel(aes(label = Group),col='black', size = 3, min.segment.length = 3, seed =111) +
  scale_size(name = "Number of\nGenes") +
  scale_color_viridis(
    option = 'mako',
    name = "Fraction Selected"
  ) +
  labs(x = "Fraction A", y = "Fraction B") +
  theme_bw(base_size = 8)
subg_sel
ggsave('~/symlinks/comp/figures/20260406_SubgenomeFraction_Purifying.pdf',subg_sel,height=5,width=7)

##### OF those branches (C) with high A/B variation: what are the top candidates? ######
gene_pair <- wstats_clean %>%
  filter(Localization %in% c("Shared_A", "Shared_B"), ID == "C") %>%
  group_by(Gene, Localization) %>%
  summarise(
    omega_med = median(omega, na.rm = TRUE),
    puri_frac = mean(puri, na.rm = TRUE),
    pos_frac = mean(pos, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = Localization,
    values_from = c(omega_med, puri_frac, pos_frac, n),
    names_sep = "."
  ) %>%
  mutate(
    log2_omega_ratio = log2((omega_med.Shared_A + 1e-6)/(omega_med.Shared_B + 1e-6)),
    delta_puri = puri_frac.Shared_A - puri_frac.Shared_B,
    category = case_when(
      puri_frac.Shared_A >= 0.5 & puri_frac.Shared_B <= 0.1 &
        pos_frac.Shared_B <= 0.05 & log2_omega_ratio < -1 ~ "A_biased_purifying",
      puri_frac.Shared_B >= 0.5 & puri_frac.Shared_A <= 0.1 &
        pos_frac.Shared_A <= 0.05 & log2_omega_ratio > 1 ~ "B_biased_purifying",
      TRUE ~ "other"
    )
  )
top_puri_genes <- gene_pair %>% dplyr::select(Gene,category) %>% filter(category != 'other')
top_puri_uniq <- wstats_clean %>% filter(ID == 'C' & puri == TRUE & grepl('Unique',Localization))  %>% dplyr::select(Gene,category=Localization)
top_puri <- rbind(top_puri_genes,top_puri_uniq)
write.table(top_puri,file='TopGenes_Puri_20260406.tsv',quote=F,sep='\t',row.names=F)

# also save full file with the gene pair dnds 
savedf <- gene_pair %>% dplyr::select(Gene, omega_med.Shared_A, omega_med.Shared_B, log2_omega_ratio, delta_puri, category)
write.table(savedf,file='AllGenes_Puri_20260406.tsv',quote=F,sep='\t',row.names=F)

```

#### Plot Functions

First, use interproscan to identify the function of the Morus genes (which are used as anchors):

```bash
#!/bin/bash

#SBATCH --time=2-00:00:00   
#SBATCH --nodes=1  
#SBATCH --cpus-per-task=8
#SBATCH --mem=52Gb
#SBATCH --partition=ceres
#SBATCH --account=coffea_pangenome

#module load miniconda
#mamba create -n interpro -c conda-forge -c bioconda interproscan nextflow>=25.04.6 -y
#source activate interpro

RUN="${1:?usage: $0 <FASTA>}"
DATA_DIR=/project/coffea_pangenome/Artocarpus/Comparative_Paper/annotation/interproscan
BASE=$(basename ${RUN})
echo "Running on ${BASE}"

interproscan.sh \
  -i ${RUN} \
  -f tsv \
  --goterms --pathways \
  --cpu 8
```



This maps candidate subgenome-biased purifying-selection genes to GO terms (from InterProScan on *Morus*) and runs Fisher-test GO enrichment (vs a genome-wide universe) to highlight overrepresented functions.

Outputs

- GO enrichment plot: `20260406_PurifyingSubgenomeA_Enrichment.pdf`
- Enrichment tables in-session: `mf_sA`, `mf_sB`, `mf_uA`, `mf_uB` (top terms subset in `plot_df`)
- Candidate gene <-> GO table for highlighted terms: `TopGenes_Puri_SubgenomeA_20260406.tsv`

```R
setwd("/project/coffea_pangenome/Artocarpus/Comparative_Paper/subgenome_divided_dnds/")
# Subgenome-specific selection: functional annotation of genes 
library(tidyverse)
library(scales)
library(stringr)
library(GO.db)
library(AnnotationDbi)
library(data.table)	
library(ggrepel)

# Read in genes 
genes <- read_tsv("TopGenes_Puri_20260406.tsv")
wstats_clean <- read_tsv("Node_dNdS_20260406-RInput.tsv")
homeo <- read_tsv("TopHomeologs_20260414.tsv")

id <- readr::read_tsv(
  "/project/coffea_pangenome/Artocarpus/Comparative_Paper/orthofinder/Subgenome_Divided/interproscan/Morus.tsv",
  col_names = FALSE, show_col_types = FALSE
) %>%
  transmute(Gene = X1, source_db = X4, interpro_id = X12, go_raw = X14) %>%
  distinct()

gene2go_long <- id %>%
  filter(!is.na(go_raw), go_raw != "-", str_detect(go_raw, "GO:")) %>%
  # split on |
  separate_rows(go_raw, sep = "\\|") %>%
  # extract the GO:######## part only 
  mutate(GO = str_extract(go_raw, "GO:\\d{7}")) %>%
  filter(!is.na(GO)) %>%
  distinct(Gene, GO)

# All genes 
universe_genes <- sort(unique(wstats_clean$Gene))
universe_genes <- intersect(universe_genes, unique(gene2go_long$Gene))
length(universe_genes)

##### PURIFYING SELECTION #####
# Candidates 
sA <- unique(genes$Gene[genes$category == "A_biased_purifying"])
sB <- unique(genes$Gene[genes$category == "B_biased_purifying"])
uA <- unique(genes$Gene[genes$category == "Unique_A"])
uB <- unique(genes$Gene[genes$category == "Unique_B"])
hC <- homeo %>% filter(ID == 'C' & cand == 'purifying') %>% pull(Gene)
hJ <- homeo %>% filter(ID == 'J' & cand == 'purifying') %>% pull(Gene)
pos <- homeo %>% filter(cand == 'positive') %>% pull(Gene)

sA <- intersect(sA, universe_genes)
sB <- intersect(sB, universe_genes)
uA <- intersect(uA, universe_genes)
uB <- intersect(uB, universe_genes)
hC <- intersect(hC, universe_genes)
hJ <- intersect(hJ, universe_genes)
pos <- intersect(pos, universe_genes)

# Map ontology
go_ontology <- AnnotationDbi::select(
  GO.db,
  keys = unique(gene2go_long$GO),
  columns = c("ONTOLOGY", "TERM"),
  keytype = "GOID"
) %>%
  as_tibble() %>%
  dplyr::rename(GO = GOID) %>%
  distinct(GO, ONTOLOGY, TERM)

gene2go_annot <- gene2go_long %>%
  inner_join(go_ontology, by = "GO") %>%
  filter(!is.na(ONTOLOGY))

# Enrichment
enrich_go_fisher <- function(cand_genes, universe_genes, gene2go_annot, ontology = "MF") {
  # gene > GO restricted to universe + ontology
  g2 <- gene2go_annot %>%
    filter(ONTOLOGY == ontology, Gene %in% universe_genes) %>%
    distinct(Gene, GO, TERM)
  
  # term counts in universe
  bg <- g2 %>%
    group_by(GO, TERM) %>%
    summarise(bg_with_term = n_distinct(Gene), .groups = "drop")
  
  # term counts in candidates
  fg <- g2 %>%
    filter(Gene %in% cand_genes) %>%
    group_by(GO, TERM) %>%
    summarise(cand_with_term = n_distinct(Gene), .groups = "drop")
  
  res <- bg %>%
    left_join(fg, by = c("GO", "TERM")) %>%
    mutate(
      cand_with_term = replace_na(cand_with_term, 0L),
      cand_total = length(unique(cand_genes)),
      bg_total   = length(unique(universe_genes)),
      # 2x2 Fisher: in candidates vs not; has term vs not
      a = cand_with_term,
      b = cand_total - cand_with_term,
      c = bg_with_term - cand_with_term,
      d = (bg_total - bg_with_term) - (cand_total - cand_with_term)
    ) %>%
    rowwise() %>%
    mutate(
      p = fisher.test(matrix(c(a,b,c,d), nrow=2, byrow=TRUE), alternative="greater")$p.value,
      odds_ratio = fisher.test(matrix(c(a,b,c,d), nrow=2, byrow=TRUE), alternative="greater")$estimate[[1]]
    ) %>%
    ungroup() %>%
    mutate(
      p_adj = p.adjust(p, method = "BH"),
      fold_enrich = (a / cand_total) / (bg_with_term / bg_total)
    ) %>%
    arrange(p_adj, p)
  
  res
}

mf_sA <- enrich_go_fisher(sA, universe_genes, gene2go_annot, ontology = "MF")
mf_sB <- enrich_go_fisher(sB, universe_genes, gene2go_annot, ontology = "MF")
mf_uA <- enrich_go_fisher(uA, universe_genes, gene2go_annot, ontology = "MF")
mf_uB <- enrich_go_fisher(uB, universe_genes, gene2go_annot, ontology = "MF")
mf_hC <- enrich_go_fisher(hC, universe_genes, gene2go_annot, ontology = "MF")
mf_hJ <- enrich_go_fisher(hJ, universe_genes, gene2go_annot, ontology = "MF")
mf_pos <- enrich_go_fisher(pos, universe_genes, gene2go_annot, ontology = "MF")

plot_df <- bind_rows(
  mf_sA %>% mutate(Set = "Shared A"),
  mf_sB %>% mutate(Set = "Shared B"),
  mf_uA %>% mutate(Set = "Unique A"),
  mf_uB %>% mutate(Set = "Unique B"),
  mf_hC %>% mutate(Set = "Homeologs C"),
  mf_hJ %>% mutate(Set = "Homeologs J"),
  mf_pos %>% mutate(Set = "Positive Selection")
) %>%
  mutate(
    sig10 = !is.na(p_adj) & p_adj < 0.10,
    label = str_trunc(TERM, 60),
    fold_enrich = as.numeric(fold_enrich)
  ) %>%
  group_by(Set) %>%
  arrange(p_adj, p, desc(fold_enrich)) %>%
  slice_head(n = 15) %>%
  ungroup()

plot_df %>% filter(Set == 'Homeologs C') %>% dplyr::select(TERM,p_adj,fold_enrich,sig10)
plot_df %>% filter(Set == 'Homeologs J') %>% dplyr::select(TERM,p_adj,fold_enrich,sig10)
plot_df %>% filter(Set == 'Positive Selection') %>% dplyr::select(TERM,p_adj,fold_enrich,sig10)

plot_df %>% filter(Set == 'Unique A') %>% dplyr::select(TERM,p_adj,fold_enrich,sig10)
plot_df %>% filter(Set == 'Unique A' & grepl('manno|carb',TERM)) %>% dplyr::select(GO,bg_with_term,cand_with_term,TERM,p_adj,fold_enrich,sig10)
# GO         bg_with_term cand_with_term TERM                         p_adj fold_enrich sig10
# <chr>             <int>          <int> <chr>                        <dbl>       <dbl> <lgl>
#   1 GO:0030246           61             13 carbohydrate binding       0.00957        4.03 TRUE 
# 2 GO:0004559            5              4 alpha-mannosidase activity 0.0171        15.1  TRUE 

# plot
enrich_plot <- plot_df %>% 
  filter(Set == 'Unique A' | Set == 'Homeologs C') %>% 
  ggplot(aes(x = fold_enrich, y = reorder(label, fold_enrich))) +
  geom_point(aes(shape = sig10, size = cand_with_term, color = fold_enrich)) +
  facet_wrap(~ Set, scales = "free_y") +
  scale_x_log10(
    labels = label_number(accuracy = 0.1)
  ) +
  scale_shape_manual(
    values = c(`FALSE` = 1, `TRUE` = 16),
    name = "FDR < 10%"
  ) +
  scale_color_gradientn(
    colours = c("#2b83ba", "#abdda4", "#fdae61", "#d7191c"),
    trans = "log10",
    name = "Fold\nenrichment"
  ) +
  labs(
    x = "Fold enrichment (log10 scale)",
    y = NULL
  ) +
  scale_y_discrete(labels = function(x) str_wrap(x, width = 25)) +
  theme_bw(base_size = 8) +
  theme(
    legend.position = "right",
    strip.background = element_rect(fill = "grey95")
  )
enrich_plot

ggsave('~/symlinks/comp/figures/20260414_GO_Enrichment.pdf',enrich_plot,height=2.5,width=5.5)

# Extract those genes
cands <- gene2go_annot %>% filter(grepl('0030246|0004559',GO)) %>% filter(Gene %in% uA) %>% data.frame
write_tsv(cands,file='TopGenes_Puri_SubgenomeA_20260406.tsv')

```



## Subgenome-specific Expression

This section quantifies RNA-seq/Iso-Seq read expression against the combined HART063 A+B CDS transcriptome with Salmon, maps transcripts to *Morus* orthologs via RBH lookups, and summarizes per-gene subgenome A vs B expression bias across dN/dS-based gene categories.

Outputs

- RBH lookup table: `Morus_Lookup.tsv`
- Expression-bias figure: `20260406_Expression_Bias_Subgenome.pdf`

```bash
#!/bin/bash

#SBATCH --time=3-00:00:00    
#SBATCH --cpus-per-task=8
#SBATCH --mem=24Gb
#SBATCH --partition=ceres
#SBATCH --account=coffea_pangenome

module load miniconda
source activate isoseq_ann

READS=${1:?Missing READS argument}
SAMPLE=${2:?Missing SAMPLE argument}
WD=/project/coffea_pangenome/Artocarpus/Comparative_Paper/expression_subgenome
CDS=/project/coffea_pangenome/Artocarpus/Comparative_Paper/subgenome_divided_dnds/cds_files
BASE=$(basename ${READS} .fastq.gz)

cd ${WD}

# Merge HART063 subgenomes into a single transcriptome...
# cat ../subgenome_divided_dnds/cds_files/HART063A.fa ../subgenome_divided_dnds/cds_files/HART063B.fa > HART063.fa 

if [ -d "${SAMPLE}_index" ] && [ -f "${SAMPLE}_index/info.json" ]; then
  echo "Index exists for ${SAMPLE}, skipping"
else
  salmon index -t "${SAMPLE}.fa" -i "${SAMPLE}_index" -k 31
fi

salmon quant -i ${SAMPLE}_index -l U -r ${READS} --validateMappings -o ${BASE}_quant
```

Extract the Morus ~ Artocarpus/Batocarpus gene RBH pairs from the subgenome-divided dnds dir: 

```bash
cat blast/RBH_Morus_HART063A.txt blast/RBH_Morus_HART063B.txt blast/RBH_Morus_Batocarpus.txt > ../expression_subgenome/Morus_Lookup.tsv
```

Plot subgenome expression bias across dnds categories:

```R
setwd('/project/coffea_pangenome/Artocarpus/Comparative_Paper/expression_subgenome') 
library(tidyverse)
library(tximport)
library(DESeq2)
library(pheatmap)
library(readr)
library(ggrepel)

# Read in top puri genes from scan 
genes <- read_tsv("/project/coffea_pangenome/Artocarpus/Comparative_Paper/subgenome_divided_dnds/TopGenes_Puri_SubgenomeA_20260406.tsv")

# import dnds for shared copies of all paired genes 
savedf <- read_tsv('/project/coffea_pangenome/Artocarpus/Comparative_Paper/subgenome_divided_dnds/AllGenes_Puri_20260406.tsv')

# read in morus ~ arto/bato gene ID lookups
ids <- read_tsv('Morus_Lookup.tsv',col_names = F)
names(ids) <- c('Morus','RBH')

# read in isoseq data, first for artocarpus
dirs <- list.files(path = ".", pattern = "_quant$", full.names = FALSE)

mat <- str_match(dirs, "^(.*?)__(.*?)_quant$")

samples <- tibble(dir = dirs) %>%
  mutate(
    prefix = mat[,2],
    tissue = mat[,3],
    species = case_when(
      prefix == "HART063" ~ "Artocarpus",
      prefix == "N15_23" ~ "Batocarpus",
      TRUE ~ prefix
    ),
    sample = paste0(species, "_", tissue),
    quant = file.path(dir, "quant.sf")
  ) %>%
  dplyr::select(sample, species, tissue, prefix, quant)

samples

# import separately
samples_hart <- samples %>% filter(prefix == "HART063")
files_hart <- samples_hart$quant
names(files_hart) <- samples_hart$sample
txi_hart <- tximport(files_hart, type = "salmon", txOut = TRUE, ignoreTxVersion = TRUE)

# convert to tpm martrix
hart_tpm <- as_tibble(txi_hart$abundance, rownames = "RBH") %>%
  pivot_longer(-RBH, names_to = "sample", values_to = "TPM") %>%
  left_join(samples_hart, by = "sample")

# map to hart
hart_tpm2 <- hart_tpm %>%
  left_join(ids, by = "RBH") %>%
  mutate(
    subgenome = case_when(
      str_detect(RBH, "^HART063A_") ~ "A",
      str_detect(RBH, "^HART063B_") ~ "B",
      TRUE ~ NA_character_
    )
  )

# sumarize Arto A/B to Morus
hart_AB <- hart_tpm2 %>%
  filter(!is.na(Morus), !is.na(subgenome)) %>%
  group_by(Morus, tissue, subgenome) %>%
  summarize(TPM = sum(TPM, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = subgenome, values_from = TPM, values_fill = 0) %>%
  mutate(
    total_AB = A + B,
    log2_BA = log2((B + 1) / (A + 1))
  )
write_tsv(hart_AB,   "hart_AB.tsv")

# Plot log2(B/A) expression across purifying categories
exp_omega <- left_join(hart_AB %>% dplyr::rename(Gene = Morus),savedf %>% mutate(omega_diff = omega_med.Shared_A - omega_med.Shared_B)) %>% drop_na(category)

# reduce to gene level 
exp_gene <- exp_omega %>%
  group_by(Gene, category) %>%
  summarise(mean_log2_BA = mean(log2_BA, na.rm = TRUE), .groups = "drop")
pairwise.wilcox.test(exp_gene$mean_log2_BA, exp_gene$category, p.adjust.method = "holm")

# plot it 
plot_sum <- exp_gene %>%
  mutate(category = factor(category,
                           levels = c("A_biased_purifying", "other", "B_biased_purifying"))) %>%
  group_by(category) %>%
  summarise(
    mean = mean(mean_log2_BA, na.rm = TRUE),
    sd = sd(mean_log2_BA, na.rm = TRUE),
    n = n(),
    se = sd / sqrt(n),
    conf_low = mean - 1.96 * se,
    conf_high = mean + 1.96 * se,
    .groups = "drop"
  )

exp_plot <- ggplot(plot_sum, aes(x = category, y = mean, color = category)) +
  geom_hline(yintercept = 0, linetype = 2, color = "grey55") +
  geom_errorbar(aes(ymin = conf_low, ymax = conf_high), width = 0.10, linewidth = 0.7) +
  geom_point(size = 2.8) +
  scale_color_manual(values = c(
    "A_biased_purifying" = "#1b9e77",
    "other" = "grey55",
    "B_biased_purifying" = "#d95f02"
  )) +
  labs(x = NULL, y = "Mean log2(B/A) expression per gene") +
  theme_bw(base_size = 9) +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(face = "bold")
  )

ggsave('~/symlinks/comp/figures/20260406_Expression_Bias_Subgenome.pdf',exp_plot,height=2.5,width=1.5)

```



## kSRates: Ks distributions

This prepares cleaned CDS FASTAs and the ARtocarpus GFF3, then runs kSRates (in manual mode) to estimate paralog and ortholog Ks distributions across Artocarpus–Batocarpus–Morus and perform rate adjustment before generating plots. 

Outputs

- Ortholog peak / Ks list databases: `ortholog_peak_db.tsv`, `ortholog_ks_list_db.tsv`
- Figures: ortholog Ks distribution plots, adjusted mixed paralog–ortholog Ks plot(s), inferred WGD component plot(s), and the Ks-scaled phylogram/tree plot

I'm running it in [manual mode](https://ksrates.readthedocs.io/en/latest/usage.html#run-example-case-as-a-nextflow-pipeline-recommended) , takes around 1 day.

From the genome annotation, take the transdecoder CDS and the genome-coordinate .gff3 file from Artocarpus:

```bash
WD=/90daydata/coffea_pangenome/scratch/gene_work/ksrates/raw_manual/sequences
sed '/^>/ s/[-.].*//' HART063.fa.transdecoder.cds > ${WD}/Artocarpus.fa
sed '/^>/ s/[-.].*//'  N15_23.fa.transdecoder.cds > ${WD}/Batocarpus.fa
sed '/^>/ s/[-.].*//' Morus_mongolica_refN15_23.fa.transdecoder.cds > ${WD}/Morus.fa
cp HART063.gff3 ${WD}/Artocarpus.gff3
```

Sanity check to confirm that the fasta headers exist in the gff3:

```bash
awk '
  ## Pass 1: FASTA (unique) ##
  FNR==NR {
    if ($0 ~ /^>/) {
      id = $0
      sub(/^>/, "", id)
      if (!(id in fasta)) { fasta[id] = 1; total_fa_unique++ }
    }
    next
  }

  ## Pass 2: GFF3 transcripts, unique matching ##
  $0 ~ /^#/ { next }
  $3 == "transcript" {
    tid=""
    if (match($9, /ID=([^;]+)/, m)) tid=m[1]
    base=tid; sub(/-.*/, "", base)
    if (!(base in fasta)) {
      if (match($9, /Parent=([^;]+)/, p)) {
        base=p[1]; sub(/-.*/, "", base)
      }
    }
    if (base in fasta) matched_ids[base]=1
  }

  END {
    for (k in matched_ids) matched_unique++
    not_matched = total_fa_unique - matched_unique

    printf("Unique FASTA IDs:              %d\n", total_fa_unique)
    printf("Unique FASTA IDs matched:      %d\n", matched_unique)
    printf("Unique FASTA IDs not matched:  %d\n", not_matched)
  }
' Artocarpus.fa Artocarpus.gff3
```

Great, moving on, config file `config_files/config_arto.txt`:

```bash
cat config_arto.txt 
[SPECIES]
focal_species = arto
# informal name of the focal species from the input tree

newick_tree = ((arto, bato), morus);
# input phylogenetic tree in newick format; use the informal names

latin_names =       arto    : Artocarpus altilis, 
                    bato     : Batocarpus spp,
                    morus : Morus spp
# informal names associated to their latin name through a colon and separated by comma

fasta_filenames =   arto    : sequences/Artocarpus.fa, 
                    bato     : sequences/Batocarpus.fa, 
                    morus : sequences/Morus.fa
gff_filename = sequences/Artocarpus.gff3
# informal names associated to their filename/path through a colon and separated by comma

peak_database_path = ortholog_peak_db.tsv
ks_list_database_path = ortholog_ks_list_db.tsv
# filenames/paths of the ortholog data databases


[ANALYSIS SETTING]
paranome = yes
collinearity = yes
reciprocal_retention = no
# analysis type for paralog data; allowed values: 'yes' or 'no'

gff_feature = transcript
# keyword to parse the sequence type from the gff file (column 3); can be 'gene', 'mrna'...

gff_attribute = ID
# keyword to parse gene id from the gff file (column 9); can be 'id', 'name'...

max_number_outgroups = 4
# maximum number of outspecies/trios selected to correct each divergent species pair (default: 4)

consensus_mode_for_multiple_outgroups = mean among outgroups
# allowed values: 'mean among outgroups' or 'best outgroup' (default: 'mean among outgroups')


[PARAMETERS]
x_axis_max_limit_paralogs_plot = 5
# highest value of the x axis in the mixed distribution plot (default: 5)

bin_width_paralogs = 0.1
# bin width in paralog ks histograms (default: 0.1, ten bins per unit)

y_axis_max_limit_paralogs_plot = None
# highest value of the y axis in the mixed distribution plot  (default: none)

num_bootstrap_iterations = 200
# number of bootstrap iterations for ortholog peak estimate

divergence_colors =  Red, MediumBlue, Goldenrod, Crimson, ForestGreen, Gray, SaddleBrown, Black
# color of the divergence lines drawn in correspondence of the ortholog peaks
# use color names/codes separated by comma and use at least as many colors as the number of divergence nodes

x_axis_max_limit_orthologs_plots = 5
# highest value of the x axis in the ortholog distribution plots (default: 5)

bin_width_orthologs = 0.1
# bin width in ortholog ks histograms (default: 0.1, ten bins per unit)

max_ks_paralogs = 5
# maximum paralog ks value accepted from ks data table (default: 5)

max_ks_orthologs = 10
# maximum ortholog ks value accepted from ks data table (default: 10)
```

Run ksRates:

```bash
#!/bin/bash

#SBATCH --time=2-00:00:00   
#SBATCH --cpus-per-task=16
#SBATCH --mem=32Gb
#SBATCH --partition=ceres
#SBATCH --account=coffea_pangenome

module load apptainer

# Host working directory (where config_files/ and rate_adjustment/ live)
WD="/90daydata/coffea_pangenome/scratch/gene_work/ksrates/raw_manual"
MPLCACHE="/tmp/mpl-$$"

# Bind WD to /work, start in /work so relative paths resolve
ks="apptainer exec \
  -B ${WD}:/work \
  --pwd /work \
  --no-home --cleanenv \
  --env MPLCONFIGDIR=${MPLCACHE} \
  /90daydata/coffea_pangenome/scratch/containers/vibpsb-ksrates-latest.img"

mkdir -p "${MPLCACHE}"

# Estimate the mode and associated standard deviation for each ortholog KS distribution:
${ks} ksrates orthologs-analysis --expert config_files/config_expert.txt config_files/config_arto.txt

# Plot the ortholog KS distributions for each focal species–other species pair (and each of their trios):
${ks} ksrates plot-orthologs config_files/config_arto.txt --expert config_files/config_expert.txt

# perform the rate-adjustment.
${ks} ksrates orthologs-adjustment config_files/config_arto.txt --expert config_files/config_expert.txt

# Plot the adjusted mixed paralog–ortholog KS distribution plot
${ks} ksrates plot-paralogs config_files/config_arto.txt --expert config_files/config_expert.txt

# Plot the phylogram based on the input phylogenetic tree with branch lengths equal to the KS distances estimated from the ortholog KS distirbutions
${ks} ksrates plot-tree config_files/config_arto.txt --expert config_files/config_expert.txt

# Plot the adjusted mixed paralog–ortholog KS distribution with inferred WGD components:
${ks} ksrates paralogs-analyses config_files/config_arto.txt --expert config_files/config_expert.txt

cp ./rate_adjustment/arto/mixed_arto_adjusted_anchors.pdf . # syntenic anchor pairs kS rate-adjusted distribution 
cp ./rate_adjustment/arto/mixed_arto_adjusted_paranome.pdf . # rate-adjusted paralog kS distribution, incorporating correction for lineage rate differences
cp ./rate_adjustment/arto/mixed_arto_lmm_paranome.pdf . # LMM version of the adjusted paranome 
cp ./rate_adjustment/arto/orthologs_arto_bato.pdf . # rate-adjusted ortholog distribution 
cp ./rate_adjustment/arto/tree_arto_distances.pdf . # kS scaled branch length tree 
```



## Quartet Gene Discordance & BEAST

This builds 4-taxon ortholog quartets (Morus–Batocarpus–ArtocarpusA–ArtocarpusB) via reciprocal best-hit BLAST, aligns each quartet, extracts shared 4-fold degenerate sites, and infers gene trees (IQ-TREE) and a quartet species tree (ASTRAL), and summarizes gene-tree/species-tree discordance by topology counting, and concatenates 4-fold sites for BEAST divergence dating.

Outputs

- Input quartets: `quadruplets/Morus_Bato_ArtoA_ArtoB.unique.tsv`
- Per-quartet FASTAs: `fastas/*.fa` and 4-fold-only FASTAs: `fastas_4fold/*.fa`
- Clean alignments: `alignments/*.aln.fa`
- Gene trees: `trees/*.treefile` and concatenated set: `astral/all_gene_trees.tre` + `astral/species_tree.tre`
- Topology counts: `topology_counts.tsv` and discordance figure: `20260304_WGD_Topology_Discordance.pdf`
- BEAST inputs/outputs (optional): `FourFoldDegenerate.nex`, `*.xml`, `*.ann`, and `20260225_BEAST_Divergence_Dating_All.pdf`

Taking the chr specific CDS files from the annotation section, divide into ArtoA / ArtoB/ Bato / Morus cds fastas. 

Compare the Batocarpus / Artocarpus A / Artocarpus B subgenome maps:

```bash
OD=/project/coffea_pangenome/Artocarpus/Comparative_Paper/duplication_orthology/copies_quartet
awk '{print $1}' Chr_Map.txt | xargs -I {} sh -c 'cat chrs/Artocarpus*{}*.fa' > ${OD}/ArtoA.fa
awk '{print $2}' Chr_Map.txt | xargs -I {} sh -c 'cat chrs/Artocarpus*{}*.fa' > ${OD}/ArtoB.fa
cat chrs/Batocarpus*fa > ${OD}/Bato.fa
cat chrs/Morus*fa > ${OD}/Morus.fa 
```

And then to ensure each CDS is unique, ensure the haps for Artocarpus are indicated in CDS headers: 

```bash
sed -i 's/Artocarpus/ArtoA/g' ArtoA.fa
sed -i 's/Artocarpus/ArtoB/g' ArtoB.fa
sed -i 's/Batocarpus/Bato/g' Bato.fa
```

```bash
#!/bin/bash

#SBATCH --job-name=RBH4_Morus_anchor
#SBATCH --time=8-00:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --partition=ceres
#SBATCH --account=coffea_pangenome

# module load miniconda
# source activate beast 
set -euo pipefail

# CONFIG
WD=/project/coffea_pangenome/Artocarpus/Comparative_Paper/duplication_orthology/copies_quartet
cd "${WD}"

# Input CDS FASTAs
MORUS_FA=${WD}/Morus.fa
BATO_FA=${WD}/Bato.fa
ARTOA_FA=${WD}/ArtoA.fa
ARTOB_FA=${WD}/ArtoB.fa

# Output dirs
mkdir -p db blast maps rbh quadruplets fastas work

# Threads for BLAST
THREADS=8

# 1. MAKE BLAST DATABASES
echo "[`date`] Making BLAST databases..."

# Anchor (Morus) DB for reciprocal BLAST
if [ ! -f db/Morus.nhr ]; then
    makeblastdb -in "${MORUS_FA}" -dbtype nucl -out db/Morus
fi

# Target DBs: Bato, ArtoA, ArtoB
if [ ! -f db/Bato.nhr ]; then
    makeblastdb -in "${BATO_FA}" -dbtype nucl -out db/Bato
fi

if [ ! -f db/ArtoA.nhr ]; then
    makeblastdb -in "${ARTOA_FA}" -dbtype nucl -out db/ArtoA
fi

if [ ! -f db/ArtoB.nhr ]; then
    makeblastdb -in "${ARTOB_FA}" -dbtype nucl -out db/ArtoB
fi

# 2. BLAST: MORUS -> OTHERS
echo "[`date`] BLAST: Morus -> Bato / ArtoA / ArtoB..."

blastn -query "${MORUS_FA}" -db db/Bato \
    -out blast/Morus_vs_Bato.tsv \
    -outfmt "6 qseqid sseqid pident length evalue bitscore" \
    -max_target_seqs 5 -evalue 1e-5 -num_threads ${THREADS}

blastn -query "${MORUS_FA}" -db db/ArtoA \
    -out blast/Morus_vs_ArtoA.tsv \
    -outfmt "6 qseqid sseqid pident length evalue bitscore" \
    -max_target_seqs 5 -evalue 1e-5 -num_threads ${THREADS}

blastn -query "${MORUS_FA}" -db db/ArtoB \
    -out blast/Morus_vs_ArtoB.tsv \
    -outfmt "6 qseqid sseqid pident length evalue bitscore" \
    -max_target_seqs 5 -evalue 1e-5 -num_threads ${THREADS}

# 3. BLAST: OTHERS -> MORUS
echo "[`date`] BLAST: Bato / ArtoA / ArtoB -> Morus..."

blastn -query "${BATO_FA}" -db db/Morus \
    -out blast/Bato_vs_Morus.tsv \
    -outfmt "6 qseqid sseqid pident length evalue bitscore" \
    -max_target_seqs 5 -evalue 1e-5 -num_threads ${THREADS}

blastn -query "${ARTOA_FA}" -db db/Morus \
    -out blast/ArtoA_vs_Morus.tsv \
    -outfmt "6 qseqid sseqid pident length evalue bitscore" \
    -max_target_seqs 5 -evalue 1e-5 -num_threads ${THREADS}

blastn -query "${ARTOB_FA}" -db db/Morus \
    -out blast/ArtoB_vs_Morus.tsv \
    -outfmt "6 qseqid sseqid pident length evalue bitscore" \
    -max_target_seqs 5 -evalue 1e-5 -num_threads ${THREADS}

# 4. GET BEST HITS (PER QUERY)
echo "[`date`] Selecting best hits per query..."

# Sort by query, then by descending bitscore, keep first per query
best_hit() {
    in=$1
    out=$2
    sort -k1,1 -k6,6nr "${in}" \
      | awk -F'\t' '!seen[$1]++ {print $1"\t"$2}' > "${out}"
}

best_hit blast/Morus_vs_Bato.tsv   maps/Morus2Bato.best
best_hit blast/Morus_vs_ArtoA.tsv  maps/Morus2ArtoA.best
best_hit blast/Morus_vs_ArtoB.tsv  maps/Morus2ArtoB.best

best_hit blast/Bato_vs_Morus.tsv   maps/Bato2Morus.best
best_hit blast/ArtoA_vs_Morus.tsv  maps/ArtoA2Morus.best
best_hit blast/ArtoB_vs_Morus.tsv  maps/ArtoB2Morus.best

# 5. RECIPROCAL BEST HITS (PAIRWISE)
echo "[`date`] Computing pairwise RBH with Morus anchor..."
# Morus <-> Bato
# Morus2Bato: MorusID  BatoID
# Bato2Morus: BatoID   MorusID
# We want MorusID  BatoID where both agree
join -t $'\t' -1 1 -2 2 \
    <(sort -k1,1 maps/Morus2Bato.best) \
    <(sort -k2,2 maps/Bato2Morus.best) \
    | awk -F'\t' '{print $1"\t"$2}' \
    > rbh/Morus_Bato.rbh

# Morus <-> ArtoA
join -t $'\t' -1 1 -2 2 \
    <(sort -k1,1 maps/Morus2ArtoA.best) \
    <(sort -k2,2 maps/ArtoA2Morus.best) \
    | awk -F'\t' '{print $1"\t"$2}' \
    > rbh/Morus_ArtoA.rbh

# Morus <-> ArtoB
join -t $'\t' -1 1 -2 2 \
    <(sort -k1,1 maps/Morus2ArtoB.best) \
    <(sort -k2,2 maps/ArtoB2Morus.best) \
    | awk -F'\t' '{print $1"\t"$2}' \
    > rbh/Morus_ArtoB.rbh

# Each RBH file now has:
# MorusID   BatoID
# MorusID   ArtoAID
# MorusID   ArtoBID

# 6. INTERSECT TO GET 1:1:1:1 SETS
echo "[`date`] Intersecting RBHs to get Morus–Bato–ArtoA–ArtoB 1:1:1:1 sets..."

# Join on Morus ID across the three RBH maps
# Step 1: Morus–Bato with Morus–ArtoA
join -t $'\t' -1 1 -2 1 \
    <(sort -k1,1 rbh/Morus_Bato.rbh) \
    <(sort -k1,1 rbh/Morus_ArtoA.rbh) \
    > quadruplets/tmp_Morus_Bato_ArtoA.tsv
# Columns: MorusID  BatoID  ArtoAID

# Step 2: add Morus–ArtoB
join -t $'\t' -1 1 -2 1 \
    <(sort -k1,1 quadruplets/tmp_Morus_Bato_ArtoA.tsv) \
    <(sort -k1,1 rbh/Morus_ArtoB.rbh) \
    > quadruplets/Morus_Bato_ArtoA_ArtoB.tsv
# Columns: MorusID  BatoID  ArtoAID  ArtoBID

rm quadruplets/tmp_Morus_Bato_ArtoA.tsv

# Optional: ensure uniqueness (no repeated IDs across rows)
awk '
{
    keep=1
    for (i=1; i<=NF; i++) {
        if (seen[i][$i]++) {
            keep=0
        }
    }
    if (keep) print
}' quadruplets/Morus_Bato_ArtoA_ArtoB.tsv \
  > quadruplets/Morus_Bato_ArtoA_ArtoB.unique.tsv

echo "[`date`] 4-way RBH sets:"
wc -l quadruplets/Morus_Bato_ArtoA_ArtoB.unique.tsv
```

### Extract 4-fold Sites

Extract 4 fold, generate gene trees:

```bash
#!/bin/bash

#SBATCH --job-name=WGD_Quartets_dS_Trees
#SBATCH --time=8-00:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --partition=ceres
#SBATCH --account=coffea_pangenome

# module load miniconda
# source activate beast

JOBS=14
WD=/project/coffea_pangenome/Artocarpus/Comparative_Paper/duplication_orthology/copies_quartet

# Input CDS
MORUS_FA=${WD}/Morus.fa
BATO_FA=${WD}/Bato.fa
ARTOA_FA=${WD}/ArtoA.fa
ARTOB_FA=${WD}/ArtoB.fa

QUADS_TSV="${WD}/quadruplets/Morus_Bato_ArtoA_ArtoB.unique.tsv"   # 4 cols: MorusID BatoID ArtoAID ArtoBID

# Working/output dirs
IN_DIR="${WD}/fastas"             # unaligned 4-way RBH fastas (produced in step A)
OUT_DIR="${WD}/fastas_4fold"      # extracted 4-fold sites (your existing output)
WORK_DIR="${WD}/work_quartet"     # per-gene work
ALIGN_DIR="${WD}/alignments"      # aligned_NT.clean.fasta copied here as *.aln.fa
TREE_WORK="${WD}/tree_work"       # IQ-TREE temp prefix outputs
TREE_DIR="${WD}/trees"            # cleaned final gene trees (*.treefile)
ASTRAL_DIR="${WD}/astral"         # species tree outputs
TREE_NEWICK="(((Bato,ArtoA),ArtoB),Morus);"  # Morus as outgroup

mkdir -p "${IN_DIR}" "${OUT_DIR}" "${WORK_DIR}" "${ALIGN_DIR}" "${TREE_WORK}" "${TREE_DIR}" "${ASTRAL_DIR}"
cd "${WD}"

# Tool checks
for tool in macse degenotate.py seqkit bioawk parallel samtools iqtree java; do
  command -v "$tool" >/dev/null 2>&1 || { echo "ERROR: $tool not found in PATH"; exit 1; }
done

# (A) Build 4-way RBH FASTAs 
echo "[$(date)] Indexing input genomes (if needed)"
for fa in "${MORUS_FA}" "${BATO_FA}" "${ARTOA_FA}" "${ARTOB_FA}"; do
  [ -f "${fa}.fai" ] || samtools faidx "${fa}"
done

# --- Parallel writing ---
echo "[$(date)] Writing 4-way RBH FASTAs to: ${IN_DIR}"

write_quad() {
  local morus_id="$1" bato_id="$2" artoa_id="$3" artob_id="$4"
  local outfa="${IN_DIR}/${morus_id}__${bato_id}__${artoa_id}__${artob_id}.fa"
  echo "  -> ${outfa}"

  {
    # Prefix headers with species so yields Bato/ArtoA/ArtoB/Morus
    samtools faidx "${MORUS_FA}" "${morus_id}" | sed "1s/^>.*/>Morus/"
    samtools faidx "${BATO_FA}"  "${bato_id}"  | sed "1s/^>.*/>Bato/"
    samtools faidx "${ARTOA_FA}" "${artoa_id}" | sed "1s/^>.*/>ArtoA/"
    samtools faidx "${ARTOB_FA}" "${artob_id}" | sed "1s/^>.*/>ArtoB/"
  } > "${outfa}"
}

export -f write_quad
export IN_DIR MORUS_FA BATO_FA ARTOA_FA ARTOB_FA

# Read 4 columns: MorusID BatoID ArtoAID ArtoBID
parallel --colsep '\t' --jobs "${JOBS}" --no-notice \
  write_quad {1} {2} {3} {4} :::: "${QUADS_TSV}"

# Process one fasta
process_fa() {
  local fa="$1"
  local base
  base="$(basename "${fa%.fa}")"

  local wdir="${WORK_DIR}/${base}"
  mkdir -p "${wdir}"
  cd "${wdir}"

  # Validate we have the four species (names Bato/ArtoA/ArtoB/Morus after stripping suffix)
  sed 's/_.*//g' "$fa" > raw.fa
  local required=(Bato ArtoA ArtoB Morus)
  for sp in "${required[@]}"; do
    if ! grep -q "^>${sp}\b" raw.fa; then
      echo "Skipping ${base}: missing sequence for ${sp}" >&2
      return 0
    fi
  done

  # Quartet tree
  printf '%s\n' "${TREE_NEWICK}" > tree.nwk

  # 1) Codon align with MACSE
  macse -prog alignSequences \
    -seq  raw.fa \
    -out_NT aligned_NT.fasta \
    -out_AA aligned_AA.fasta > macse.log 2>&1

  # 2) Clean alignment (remove FS/Stop marks)
  macse -prog exportAlignment \
    -align aligned_NT.fasta \
    -codonForExternalFS --- \
    -codonForFinalStop --- \
    -codonForInternalFS --- \
    -codonForInternalStop --- \
    -charForRemainingFS - \
    -out_NT aligned_NT.clean.fasta \
    -out_AA aligned_AA.clean.fasta

  # Confirm we still have the four samples post-clean
  local nsamp
  nsamp="$(grep -c '^>' aligned_NT.clean.fasta || true)"
  if [[ "$nsamp" -ne 4 ]]; then
    echo "Skipping ${base}: cleaned alignment has ${nsamp} sequences (expected 4)" >&2
    return 0
  fi

  # 3) Degeneracy annotation: extract 4-fold sites
  degenotate.py --overwrite -s aligned_NT.clean.fasta -x 4 -o 4fold > degen.log 2>&1

  # 4) Intersect positions that are 4-fold in ALL FOUR sequences.
  awk '$5 == 4' 4fold/degeneracy-all-sites.bed \
    | cut -f3 \
    | sort -n \
    | uniq -c \
    | awk '$1 == 4 {print $2}' > positions_1based.txt

  # If no positions, emit empty sequences with headers
  if [[ ! -s positions_1based.txt ]]; then
    {
      echo ">Bato";  echo
      echo ">ArtoA"; echo
      echo ">ArtoB"; echo
      echo ">Morus"; echo
    } > 4fold_extracted.fasta
  else
    # 5) Extract those positions from the clean alignment
    pos="$(paste -sd, positions_1based.txt)"
    bioawk -c fastx -v pos="$pos" '
      BEGIN { n = split(pos, P, ",") }
      {
        out = ""
        for (i = 1; i <= n; i++) {
          out = out substr($seq, P[i], 1)
        }
        print ">" $name "\n" out
      }' aligned_NT.clean.fasta > 4fold_extracted.fasta
  fi

  # 6) Ensure order: Bato, ArtoA, ArtoB, Morus
  seqkit grep -r -p Bato  4fold_extracted.fasta \
    | cat - <(seqkit grep -r -p ArtoA 4fold_extracted.fasta) \
            <(seqkit grep -r -p ArtoB 4fold_extracted.fasta) \
            <(seqkit grep -r -p Morus 4fold_extracted.fasta) \
    > ordered.fa

  # 7) Final outputs
  cp ordered.fa "${OUT_DIR}/${base}.fa"

  # (B-pre) Copy the cleaned full alignment (not 4-fold) for tree building later
  cp aligned_NT.clean.fasta "${ALIGN_DIR}/${base}.aln.fa"
}

export -f process_fa
export WD IN_DIR OUT_DIR WORK_DIR ALIGN_DIR TREE_NEWICK

# Parallel
find "${IN_DIR}" -type f -name '*.fa' -print0 \
  | parallel -0 \
      --jobs "${JOBS}" \
      --no-notice \
      --keep-order \
      process_fa {}

# (B) Build gene trees, then ASTRAL species tree
echo "[$(date)] Building IQ-TREE gene trees in parallel..."

build_tree() {
  local aln="$1"
  local base
  base="$(basename "$aln" .aln.fa)"
  local outfile="${TREE_DIR}/${base}.treefile"

  # Skip if already finished
  if [[ -f "$outfile" ]]; then
    echo "Skipping ${base} — already done"
    return 0
  fi

  echo "  -> ${base}"
  iqtree -s "$aln" \
         -m MFP \
         -bb 1000 \
         -nt AUTO \
         -pre "${TREE_WORK}/${base}"

  if [[ ! -f "${TREE_WORK}/${base}.treefile" ]]; then
    echo "WARNING: IQ-TREE failed for ${base}"
  fi
}
export -f build_tree
export TREE_DIR TREE_WORK

# Run on all alignments; keep order of outputs
find "${ALIGN_DIR}" -type f -name '*.aln.fa' -print0 \
  | parallel -0 --jobs "${JOBS}" --no-notice --keep-order build_tree {}

# CONCATENATE ALL GENE TREES FOR ASTRAL
echo "[`date`] Concatenating gene trees..."
cat "${TREE_DIR}"/*.treefile > "${ASTRAL_DIR}/all_gene_trees.tre"

# RUN ASTRAL TO INFER SPECIES TREE
echo "[`date`] Running ASTRAL..."
java -jar ~/symlinks/software/ASTRAL-5.7.1/astral.5.7.1.jar \
    -i "${ASTRAL_DIR}/all_gene_trees.tre" \
    -o "${ASTRAL_DIR}/species_tree.tre" \
    -t 2 \
    2> "${ASTRAL_DIR}/astral.log"

echo "[`date`] ASTRAL species tree written to: ${ASTRAL_DIR}/species_tree.tre"
```

### Count Topos: Gene/Species Tree Discord

Count topologies:

````python
#!/usr/bin/env python3

import argparse
import glob
from collections import Counter
from ete3 import Tree

def classify(tree, outgroup):
    # Root the tree on the user-specified outgroup
    tree.set_outgroup(outgroup)

    # Get all leaf names
    leaves = sorted(tree.get_leaf_names())

    # The ingroup taxa are the three not equal to the outgroup
    ingroup = [t for t in leaves if t != outgroup]

    if len(ingroup) != 3:
        raise ValueError(f"Tree does not contain exactly 3 ingroup taxa: {leaves}")

    a, b, c = ingroup

    # Compute pairwise patristic distances
    d = tree.get_distance
    pairs = {
        (a, b): d(a, b),
        (a, c): d(a, c),
        (b, c): d(b, c)
    }

    # The smallest distance pair defines the topology
    closest = min(pairs, key=pairs.get)

    # Build a canonical topology string
    # Example: Outgroup,(A,(B,C))
    x, y = closest
    z = [t for t in ingroup if t not in closest][0]

    # The topology is: outgroup sister to a clade where the closest pair is nested deepest
    return f"{outgroup},({z},({x},{y}))"

def main():
    parser = argparse.ArgumentParser(
        description="Classify rooted 4‑taxon gene trees into one of three topologies and output counts."
    )
    parser.add_argument(
        "-i", "--input",
        required=True,
        help="Directory containing .treefile gene trees"
    )
    parser.add_argument(
        "-o", "--outgroup",
        required=True,
        help="Name of the outgroup taxon (must appear in every tree)"
    )
    args = parser.parse_args()

    pattern = f"{args.input.rstrip('/')}" + "/*.treefile"
    files = glob.glob(pattern)

    if not files:
        print(f"No .treefile files found in {args.input}")
        # Still write an empty TSV with header for consistency
        with open("topology_counts.tsv", "w") as outfh:
            outfh.write("topology\tcount\n")
        return

    counts = Counter()

    for f in files:
        try:
            t = Tree(f)
            topo = classify(t, args.outgroup)
            counts[topo] += 1
        except Exception as e:
            print(f"Warning: failed to parse {f}: {e}")

    # Screen printing
    print("\nTopology counts:")
    for topo in sorted(counts.keys()):
        print(f"{topo}: {counts[topo]}")

    total = sum(counts.values())
    print(f"\nTotal trees processed: {total}")

    # Write TSV in CWD
    out_path = "topology_counts.tsv"
    with open(out_path, "w") as outfh:
        outfh.write("topology\tcount\n")
        for topo in sorted(counts.keys()):
            outfh.write(f"{topo}\t{counts[topo]}\n")

    print(f"\nWrote TSV: {out_path}")

if __name__ == "__main__":
    main()
````

```
python 03_TopologyCounting.py -i trees/ -o Morus
Topology counts:
Morus,(ArtoB,(ArtoA,Bato)): 7750
Morus,(ArtoA,(ArtoB,Bato)): 788
Morus,(Bato,(ArtoA,ArtoB)): 201

Total trees processed: 8739
```

### Plot Discordance 

Plot:

```R
# install.packages(c("ggplot2", "tibble", "scales", "patchwork", "ape"))  # CRAN
# BiocManager::install("ggtree")  # from Bioconductor
library(ggplot2)
library(tibble)
library(scales)
library(patchwork)
library(ape)
library(ggtree)

# Data
df <- tribble(
  ~Topology, ~Count,
  "Morus,(ArtoB,(ArtoA,Bato))", 7750,
  "Morus,(ArtoA,(ArtoB,Bato))", 788,
  "Morus,(Bato,(ArtoA,ArtoB))", 201
)

df$Topology <- factor(
  df$Topology,
  levels = c(
    "Morus,(ArtoB,(ArtoA,Bato))",
    "Morus,(ArtoA,(ArtoB,Bato))",
    "Morus,(Bato,(ArtoA,ArtoB))"
  )
)

palette <- c("#1f77b4", "#ff7f0e", "#2ca02c")

# Newick trees
t1 <- read.tree(text = "(Morus,(ArtoB,(ArtoA,Bato)));")
t2 <- read.tree(text = "(Morus,(ArtoA,(ArtoB,Bato)));")
t3 <- read.tree(text = "(Morus,(Bato,(ArtoA,ArtoB)));")

# helper compact tree plot 
make_tree_plot <- function(tr, color, title = NULL) {
  ggtree(tr, size = 0.5, color = color) +
    geom_tiplab(size = 3) +
    labs(title = title) +
    theme_void() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 11, face = "bold"),
      plot.margin = margin(4, 4, 4, 4)
    )+
    xlim_tree(0.5)
}

p_t1 <- make_tree_plot(t1, palette[1], "Topology 1")
p_t2 <- make_tree_plot(t2, palette[2], "Topology 2")
p_t3 <- make_tree_plot(t3, palette[3], "Topology 3")

# Bar chart that matches colors
p_bar <- ggplot(df, aes(x = Topology, y = Count, fill = Topology)) +
  geom_col(width = 0.7, color = "grey20") +
  geom_text(aes(label = comma(Count)),
            vjust = -0.35, size = 4.5, fontface = "bold") +
  scale_fill_manual(values = palette) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.10))) +
  labs(x = NULL, y = "Count") +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 11),
    plot.margin = margin(8, 8, 12, 8)
  )

# Stack trees above the bar chart
final_plot <- (p_t1 + p_t2 + p_t3) / p_bar +
  plot_layout(heights = c(0.38, 0.62))

final_plot
ggsave('~/symlinks/comp/figures/20260304_WGD_Topology_Discordance.pdf')

```

### Quartet BEAST

Merge the 4-fold degenerate fasta files and then import them into beauti:

* Gamma model, 4 categories, estimated shape, GTR with estimated frequencies
* Strict clock, log normal default prior
* Yule model, tMRCA prior based on [Williams et al 2017 Out of Borneo](https://academic.oup.com/aob/article/119/4/611/2884288) paper: Morus vs Artocarpus/Batocarpus split: 83.8 74.85-92.65 Ma 
* They place stem of Artocarpus at 40.07 29.8-50.81, with likely split of Batocarpus/Artocarpus around 59.67 55.24-65.03 Ma - although this is based on genetic evidence alone. 
* 50M chains, log every 5k 

Convert the fasta to nex:

```bash
seqkit concat fastas_4fold/*fa > FourFoldDegenerate.fa 2> log
seqret -sequence FourFoldDegenerate.fa -outseq FourFoldDegenerate.nex -osformat nexus
```

Run beast:

```bash
#!/bin/bash
#SBATCH --time=10-00:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=6G
#SBATCH --partition=ceres
#SBATCH --account=coffea_pangenome

RUN=$1
MAX=20

for i in $(seq 1 $MAX); do
  SEED=$RANDOM
  echo "Try $i with seed $SEED"

  rm -f ${RUN}.log ${RUN}.trees ${RUN}.state

  beast -threads 20 -overwrite -beagle_SSE -seed $SEED ${RUN}.xml

  if [ $? -eq 0 ]; then
    treeannotator -b 10 -heights mean ${RUN}.trees ${RUN}.ann
    echo "Success"
    exit 0
  fi
done

echo "Failed after $MAX attempts"
exit 1
```

Plot:

```R
#### Plot BEAST annotated trees 
setwd('/project/coffea_pangenome/Artocarpus/Comparative_Paper/duplication_orthology/copies_quartet/beast')
library(ggtree)
library(phytools)
library(ape)
library(treeio)
library(viridis)
library(ggpubr)
library(RColorBrewer)
library(tidyverse)

#Read in metadata
md = read_tsv('~/artocarpus_comparative_genomics/samples.txt')
md <- data.frame(Accession = c('ArtoA','ArtoB','Bato','Morus'),
                 Species = c('Artocarpus','Artocarpus','Batocarpus','Morus'),
                 Haplotype = c('A','B','Batocarpus','Morus'),
                 col = c('#fc8d62','#83cab3','#8ea1cc','#e78ac3'))

files = list.files('.',paste0('.*ann'))

counter = 0
for (file in files){
  counter = counter +  1 
  iqtree = read.beast(file) 
  gg = ggtree(iqtree,layout='rectangular') %<+% md
  
  #add label for 95% CIs
  lab = gsub('.trees.*','',file)
  heights = gg$data$height_0.95_HPD
  df = as.data.frame(do.call(rbind, heights)) #convert the list to a data frame
  df$node_value = 1:nrow(df) # Add node values as a new column
  colnames(df) = c("value1", "value2", "node")
  df = df[, c("node", "value1", "value2")]
  df = df %>% 
    mutate(
      value1 = if (grepl("mu", file)) value1 / 1e6 else value1,
      value2 = if (grepl("mu", file)) value2 / 1e6 else value2,
      lab = paste0(round(value1,1),' - ',round(value2,1))) %>% 
    select(!c(value1,value2))
  
  leg = md %>% select(Species,Haplotype) %>% unique
  gg$data = left_join(gg$data,df)
  ggp = gg  +
    geom_range(range='height_0.95_HPD', color='red', alpha=.6, size=2) +
    geom_tippoint(aes(fill = Haplotype,shape=Species),size=3)+
    geom_nodelab(aes(label=lab),size=2,vjust=1) +
    ggtitle(lab)+
    #geom_tiplab(size=2)+
    #geom_nodelab(aes(x=branch, label=round(posterior, 2)), vjust=-.5, size=3) +
    scale_fill_manual(values=md$col,breaks=md$Haplotype)+
    scale_shape_manual(values=c(21,22,24))+
    theme(legend.position=c(.1, .8))+
    geom_treescale(x = 5)+
    guides(fill=guide_legend(override.aes=list(shape=21)))+
    theme(legend.position='right')
  ggp
  assign(paste0('p',counter),ggp)
} 

ggarrange(p6,p2,p5,p3,p4,p1,common.legend = TRUE)

pdf('~/symlinks/comp/figures/20260225_BEAST_Divergence_Dating_All.pdf',height=6,width=7.5)
ggarrange(p6,p2,p5,p3,p4,p1,common.legend = TRUE)
dev.off()

```

​	
