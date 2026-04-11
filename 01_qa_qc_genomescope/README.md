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

