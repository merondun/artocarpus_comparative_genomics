## 10_subgenome_beast/

This workflow repeats the general process from [subgenome_divided_dnds](/09_subgenome_dnds/) except it includes *Ficus carica* as a more distant outgroup. 

Outputs:

Panel a from:

![beast](/figures/panels/04_subgenome_history/20260528_GeneTrees_BEAST_kS.png)



Tracer logs:

Summaries from tracer:

| Run                                            | Chains combined | Posterior ESS | Likelihood ESS | Prior ESS | Tree height mean | Tree height 95% HPD | Clock rate mean | Clock rate 95% HPD     |
| ---------------------------------------------- | --------------- | ------------- | -------------- | --------- | ---------------- | ------------------- | --------------- | ---------------------- |
| 150k Codons                                    | 50M             | 3576.9        | 3471.6         | 8236.5    | 82.9954          | [74.2704, 92.0288]  | 6.38E-04        | [5.7159E-4, 7.0889E-4] |
| Four-fold Degenerate Sites                     | 50M             | 4375.4        | 4109.7         | 9001      | 83.0682          | [74.0576, 91.888]   | 1.76E-03        | [1.5752E-3, 1.956E-3]  |
| Strict: Four-fold Degenerate Sites,  <30% Gaps | 50M             | 3610.6        | 3481.8         | 9001      | 83.012           | [73.7113, 91.492]   | 1.78E-03        | [1.5934E-3, 1.9789E-3] |





___



Take the cds files from: `/project/coffea_pangenome/Artocarpus/Comparative_Paper/annotation/egapx/copies_isoliftoff_longest_transcript_per_gene/cds` 

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
cds_files/Batocarpus.fa  cds_files/HART001B.fa  cds_files/HART058A.fa  cds_files/HART060B.fa  cds_files/HART062A.fa  cds_files/HART063B.fa  cds_files/HART068A.fa  cds_files/N9750A.fa
cds_files/Ficus.fa       cds_files/HART027A.fa  cds_files/HART058B.fa  cds_files/HART061A.fa  cds_files/HART062B.fa  cds_files/HART067A.fa  cds_files/HART068B.fa  cds_files/N9750B.fa
cds_files/HART001A.fa    cds_files/HART027B.fa  cds_files/HART060A.fa  cds_files/HART061B.fa  cds_files/HART063A.fa  cds_files/HART067B.fa  cds_files/Morus.fa

head cds_files/HART061A.fa
>HART061A_000001-R1
ATGATCATGTCTTCAAAAGGGTGTTTAGAGGAGATGGGAATATCTTCAACTAATATCAGT
GATGGTGGGAAAAATTGCTATAGAGGCCATTGGAGACCTGCGGAAGACGAGAAACTCCGA
CAACTCGTCGAACAATACGGTCCTCAGAACTGGAATTTCATCGCCGAGCATCTACAAGGA
AGATCAGGAAAAAGTTGCCGATTGAGATGGTACAACCAACTAGACCCAAACATCAACAAG
AAGCCTTTCACAGAAGAAGAGGAAGAGAGGCTGCTTTCCGCCCACCGGATTTACGGCAAC
AAATGGGCTTACATAGCCAAGTATTTCCAAGGAAGAACCGACAACGCCGTCAAGAACCAT
TACCATGTCGTCATGGCAAGGCGAAAGCGAGAGCGCTTCACGCAGTACCACCatcataat
cataatcataatcatgatcattatcattatcatTATAGCAGGTCTACTTGTACTCCTAAT
CATCTTGATCAAGGTTCTTCTAATGGGAATGTTAATATTCCTCCAAAATTAGGGCTTCTC
(earlgrey722) [justin.merondun@ceres20-compute-22 deeper_subgenome_dnds]$ head cds_files/Ficus.fa 
>Ficus_000001-R1
ATGATGACGGGACACAGCGGGAAGATCACGGAGGAGCAGAAGTCCCAAATATGCAGCTTC
ATCGACTGCAACAAGCTCTCCCACAAGCTCCTCCTTGGCGCCGTCCAAAACCCCCTCATG
CCTCTCCGCTTCGTCGTCCGCGCCATGTTCGCCGACCAGCTTAACACCCGCCGCTGCATC
Atctccgccgccgccaccgccccctccacctcctcccaccaccaccgccgccgccacagc
gacTCTCCTACTCCATCTGCCATGACCCTCGGTGCCCTCCTCCAGCGCGACGCCGCCATC
AGCCAGGCCTCGCAGCTCAAGGCCACCCTCGACGCCACCACCCGCCGCATCCGCAGCCTG
GAGGAGGAGCTCTCCGGCATGAAGAAGCTCCTCATCTTGCAAAACTCCGATCAGGACCGC
CATCGGAGCCTCGTTATGGACCGATCCGCCGGACGATCCGCCAGCTTCCATGTCGGATCG
GAGAATACTAATAAGGTCAAGAAAGAGGACAGATTCTCCGCTTCGTCGGCGAGGTTTCAT
```

1. This takes a reference sample (Ficus), performs reciprocal best‑hit BLAST searches against all other samples, extracts matching CDS sequences and outputs them into a directory named after the Ficus gene in `/genes/` 

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
WD=/project/coffea_pangenome/Artocarpus/Comparative_Paper/deeper_subgenome_dnds
cd ${WD}

# Submit with sample 
B="${1:?usage: $0 <SAMPLE>}"
first=0
[ "$B" = "Batocarpus" ] && first=1
[ "$B" = "Morus" ] && first=1

# dNdS Across tree 
A="Ficus"
Af="cds_files/${A}.fa"

# output dirs
mkdir -p blast db genes

echo "Blasting ${A} against ${B}"
Bf="cds_files/${B}.fa"

# Make blast dbs
if [ ! -f "db/${B}.nhr" ]; then
makeblastdb -in $Bf -dbtype nucl -out db/${B}
fi
if [ ! -f "db/${A}.nhr" ]; then
makeblastdb -in $Af -dbtype nucl -out db/${A}
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
  
  11411 Batocarpus.fa
  10571 Ficus.fa
   9680 HART001A.fa
   8611 HART001B.fa
   9429 HART027A.fa
   8389 HART027B.fa
   9679 HART058A.fa
   8190 HART058B.fa
   9103 HART060A.fa
   7939 HART060B.fa
   9428 HART062A.fa
   8308 HART062B.fa
   9737 HART063A.fa
   8943 HART063B.fa
   9542 HART067A.fa
   8593 HART067B.fa
   9011 HART068A.fa
   7541 HART068B.fa
  10571 Morus.fa
   9085 N9750A.fa
   7861 N9750B.fa
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

module load parallel
WD=/project/coffea_pangenome/Artocarpus/Comparative_Paper/deeper_subgenome_dnds
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

  # require Ficus.fa, Morus.fa, and Batocarpus.fa to exist and be non-empty
  [ -s "${d}/Ficus.fa" ] || exit 0
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

### 4-Fold Sites & Beast

Extract the 4-fold degenerate fasta files:

```bash
#!/bin/bash

#SBATCH --time=3-00:00:00    
#SBATCH --cpus-per-task=48
#SBATCH --mem=64Gb
#SBATCH --partition=ceres
#SBATCH --account=coffea_pangenome

# module load miniconda
# source activate isoseq_ann

JOBS=48
LIST="${1:?usage: $0 <genes_list_file>}"
WD=/project/coffea_pangenome/Artocarpus/Comparative_Paper/deeper_subgenome_dnds
TREE="${WD}/tree.nwk"
GENEDIR="${WD}/raw"
NUM_SAMPS=21

mkdir -p ${WD}/hyphy ${WD}/beast_out ${WD}/beast_out_strict

cd ${WD}

echo "Working on ${LIST}"

# Ensure required tools exist
for tool in macse clipkit parallel; do
  command -v "$tool" >/dev/null 2>&1 || { echo "ERROR: $tool not found in PATH"; exit 1; }
done

SCR=${SLURM_TMPDIR:-/tmp}
export WD TREE RUN GENEDIR SCR NUM_SAMPS

process_gene() {
  local fa="$1"           
  local gene
  gene="$(basename "${fa%.fa}")" 

  # quick skip if not enough sequences (require at least 21)
  seq_count=$(grep -c '^>' "$fa" 2>/dev/null || true)
  if [ "${seq_count:-0}" -lt ${NUM_SAMPS} ]; then
    echo "Skipping ${gene}: only ${seq_count} sequences (need >=${NUM_SAMPS})"
    return 0
  fi

  # Work in a per-gene dir
  (
    # for space limited, replace WD with SCR
    outdir="${WD}/hyphy/${gene}"
    mkdir -p "$outdir"
    cd "$outdir"
    echo "Working on ${gene}"

    # 1) Align with MACSE
    macse -prog alignSequences \
      -seq "${fa}" \
      -out_NT aln_NT.fasta > macse.log 2>&1

    # 2) Clean alignment
    macse -prog exportAlignment \
      -align aln_NT.fasta \
      -codonForExternalFS --- \
      -codonForFinalStop --- \
      -codonForInternalFS --- \
      -codonForInternalStop --- \
      -charForRemainingFS - \
      -out_NT aln_NT.clean.fasta 2>&1

    # 3) Ensure underscores stripped from seq names
    sed 's/_.*//g' aln_NT.clean.fasta > aln_NT.clean.fa

    # 4) Trim gappy codon positions
    clipkit aln_NT.clean.fa -o aln_NT.clipkit.fa --codon -m kpic > clipkit.log 2>&1

    # 4b) Create a "strict" codon alignment that drops sequences with >30% gaps
    seqkit fx2tab aln_NT.clipkit.fa \
    | awk -F'\t' '
        {
          seq = $2
          gsub(/[^-]/, "", gapped)
        }
        {
          total = length($2)
          gaps = gsub(/-/, "", $2)
          if (total > 0 && gaps/total <= 0.30)
              print $1
        }
    ' | seqkit grep -f - aln_NT.clipkit.fa > aln_NT.clipkit.strict.fa

    # Require 21 samples with ATG aligned and extract 
    python /home/justin.merondun/merothon/merothon/scripts/Extract_4Fold.py -i aln_NT.clipkit.fa -m 21 -o extract
    cp extract.4fold.fa "${WD}/beast_out/${gene}.fa"
    
    # and strict...
    python /home/justin.merondun/merothon/merothon/scripts/Extract_4Fold.py -i aln_NT.clipkit.strict.fa -m 21 -o extract_strict \
      || echo "WARN: strict extract failed for ${gene}"
    [ -f extract_strict.4fold.fa ] && cp extract_strict.4fold.fa "${WD}/beast_out_strict/${gene}.fa"
    rm *log aln_NT_AA.fasta aln_NT.clean.fasta aln_NT.fasta aln_NT.clean.fa

  )
}

export -f process_gene

# Run in parallel
parallel --will-cite -j ${JOBS} process_gene :::: "${LIST}"

# concat and convert
seqkit concat beast_out/*fa > Subgenomes_4Fold.fa 2> seqkit.concat.log
seqret -sequence Subgenomes_4Fold.fa -outseq Subgenomes_4Fold.nex -osformat nexus

# concat and convert the STRICT 
seqkit concat beast_out_strict/*fa > Subgenomes_4Fold_Strict.fa 2> seqkit.strict.concat.log
seqret -sequence Subgenomes_4Fold_Strict.fa -outseq Subgenomes_4Fold_Strict.nex -osformat nexus

# also creat a full genes file
find hyphy/ -type f -name "aln_NT.clipkit.strict.fa" -print0 |
while IFS= read -r -d '' f; do
    if [ "$(seqkit stat "$f" | awk 'NR==2{print $4}')" -eq 21 ]; then
        printf '%s\0' "$f"
    fi
done | xargs -0 seqkit concat -o Subgenomes_FullGenes_Strict.fa
seqret -sequence Subgenomes_FullGenes_Strict.fa -outseq Subgenomes_FullGenes_Strict.nex -osformat nexus

# subset 150k codons
seqkit subseq -r 1:450000 Subgenomes_FullGenes_Strict.fa > Subgenomes_FullGenes_Strict_150k.fa
seqret -sequence Subgenomes_FullGenes_Strict_150k.fa -outseq Subgenomes_FullGenes_Strict_150k.nex -osformat nexus
```

Submit:

```
realpath raw/*fa > raw_genes.list
sbatch 03_Extract_4Fold.sh raw_genes.list
```

Count variant and invariant sites in all alignments:

```
#!/usr/bin/env python3
import sys
from itertools import zip_longest

def read_fasta(path):
    seqs = []
    seq = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if seq:
                    seqs.append("".join(seq))
                    seq = []
            else:
                seq.append(line)
        if seq:
            seqs.append("".join(seq))
    return seqs

def count_sites(seqs):
    invariant = 0
    variant = 0
    aln_len = len(seqs[0])
    for col in zip(*seqs):
        bases = set(col) - set("-")   # ignore gaps
        if len(bases) == 1:
            invariant += 1
        else:
            variant += 1
    return invariant, variant

if __name__ == "__main__":
    for fa in sys.argv[1:]:
        seqs = read_fasta(fa)
        inv, var = count_sites(seqs)
        print(fa)
        print("Invariant sites:", inv)
        print("Variant sites:", var)
        print()
```

Count for fastas:

```
python count_invar.py *.fa

Subgenomes_4Fold.fa
Invariant sites: 152127
Variant sites: 100100

Subgenomes_4Fold_Strict.fa
Invariant sites: 128932
Variant sites: 85451

Subgenomes_FullGenes_Strict.fa
Invariant sites: 1645648
Variant sites: 681245

Subgenomes_FullGenes_Strict_150k.fa
Invariant sites: 322415
Variant sites: 127585
```



Merge the 4-fold degenerate fasta files and then import them into beauti:

* Gamma model, 4 categories, estimated shape, GTR with estimated frequencies
* Strict clock, log normal default prior
* Yule model, tMRCA prior based on [Williams et al 2017 Out of Borneo](https://academic.oup.com/aob/article/119/4/611/2884288) paper: Morus vs Artocarpus/Batocarpus split: 83.8 74.85-92.65 Ma, = mean 83.8, sigma = 4.54
* 50M chains, log every 5k 

From the paper: 

> In the mid- to late Cretaceous (83·8 Mya, 74·85–92·65 Mya) the stem node of the tribe Artocarpeae diverged from the rest of the family Moraceae. The biogeographical reconstruction infers a likely origin of the tribe in the Americas. The split between American (*Clarisia* and *Batocarpus*) and Asian Artocarpeae (*Artocarpus*) occurred in the Palaeocene (59·67 Mya, 55·24–65·03 Mya) with a radiation of *Artocarpus* from Borneo in the Eocene to Oligocene (40·07 Mya, 29·8–50·81 Mya). 

ALSO run a BEAST analysis with the extra-filtered 4-fold sites, and one with 150k codons from the full gene alignments. 

### Plot BEAST

Summaries from tracer:

| Run                                            | Chains combined | Posterior ESS | Likelihood ESS | Prior ESS | Tree height mean | Tree height 95% HPD | Clock rate mean | Clock rate 95% HPD     |
| ---------------------------------------------- | --------------- | ------------- | -------------- | --------- | ---------------- | ------------------- | --------------- | ---------------------- |
| 150k Codons                                    | 50M             | 3576.9        | 3471.6         | 8236.5    | 82.9954          | [74.2704, 92.0288]  | 6.38E-04        | [5.7159E-4, 7.0889E-4] |
| Four-fold Degenerate Sites                     | 50M             | 4375.4        | 4109.7         | 9001      | 83.0682          | [74.0576, 91.888]   | 1.76E-03        | [1.5752E-3, 1.956E-3]  |
| Strict: Four-fold Degenerate Sites,  <30% Gaps | 50M             | 3610.6        | 3481.8         | 9001      | 83.012           | [73.7113, 91.492]   | 1.78E-03        | [1.5934E-3, 1.9789E-3] |

```R
#### Plot BEAST annotated trees 
setwd('/project/coffea_pangenome/Artocarpus/Comparative_Paper/deeper_subgenome_dnds/beast')
library(ggtree)
library(phytools)
library(ape)
library(treeio)
library(viridis)
library(ggpubr)
library(RColorBrewer)
library(tidyverse)

# metadata
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
          Accession = c('Ficus','Morus'),
          Group = c('Ficus carica','Morus mongolica')
        ))
hapa <- md1 %>% filter(grepl('A. ',Group)) %>% 
  mutate(Accession = paste0(Accession,'A'),
         Haplotype='A')
hapb <- md1 %>% filter(grepl('A. ',Group)) %>% 
  mutate(Accession = paste0(Accession,'B'),
         Haplotype='B')
ogs <- md1 %>% filter(!grepl('A. ',Group)) %>% mutate(Haplotype=Accession)
md <- rbind(hapa,hapb,ogs)

files = list.files('.',paste0('.*ann'))

counter = 0
for (file in files){
  counter = counter +  1 
  iqtree = read.beast(file) 
  t2 <- drop.tip(iqtree,'Ficus')
  gg = ggtree(t2,layout='rectangular') %<+% md
  
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
  
  leg = md %>% select(Group,Haplotype) %>% unique
  gg$data = left_join(gg$data,df)
  ggp = gg  +
    geom_range(range='height_0.95_HPD', color='red', alpha=.6, size=2) +
    geom_tippoint(aes(fill = Haplotype,shape=Haplotype),size=3)+
    geom_nodelab(aes(label=lab),size=2,vjust=1) +
    ggtitle(lab)+
    #geom_tiplab(size=2)+
    #geom_nodelab(aes(x=branch, label=round(posterior, 2)), vjust=-.5, size=3) +
    scale_fill_manual(values=c('#0072B2','#E69F00','black','white','purple'))+
    scale_shape_manual(values=c(21,22,24,25,23))+
    theme(legend.position=c(.1, .8))+
    geom_treescale(x = 5)+
    guides(fill=guide_legend(override.aes=list(shape=21)))+
    theme(legend.position='right')+
    geom_tiplab(aes(label=gsub('A. ','',Group)),size=1,hjust=-0.25)+xlim(c(-0.5,max(gg@data$x)*1.2))
  ggp
  assign(paste0('p',counter),ggp)
} 

ggarrange(p1,p2,p3,common.legend = TRUE,nrow=1)

pdf('~/symlinks/comp/figures/20260521_BEAST_Divergence_Dating_All.pdf',height=4,width=9.5)
ggarrange(p1,p2,p3,common.legend = TRUE,nrow=1)
dev.off()

```

