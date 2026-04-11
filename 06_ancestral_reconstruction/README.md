# Ancestral Reconstruction

This builds an N=5 comparative dataset (Artocarpus, Batocarpus, Morus, Ficus, Antiaris) by repeat-masking outgroups, lifting over gene annotations, and running OrthoFinder to infer orthogroups and a rooted species tree. Those outputs are then converted into AGORA inputs to reconstruct ancestral gene order/karyotypes (ancestral chromosomes/CARs) and generate karyotype-style plots plus CAR→extant chromosome mapping summaries.

The output will be an karyotype reconstruction:

![agora](/figures/20260410_AGORA_Karyotypes.png)
___

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

