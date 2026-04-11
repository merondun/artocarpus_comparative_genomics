## Subgenome-divided dNdS

This workflow splits CDS FASTAs into subgenome (A/B) partitions using chromosome haplotype lists, then identifies orthologous CDS sets via reciprocal best‑hit BLAST using *Morus* as the anchor reference. Per-gene multi-sample CDS alignments are built, filtered, and pruned to matching taxa, and HyPhy (MG94) is run on each gene to estimate branch-specific dN/dS across the tree.

This corresponds to panels B & C & D:

![subgenome_evo](/figures/20260407_SubgenomeEvolution.png)

___

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
md <-  read_tsv('~/symlinks/comp/samples.txt') %>% dplyr::select(ID=Accession,Group) %>% mutate(ID = gsub('_','',ID))

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
labs <- c("×4",
          "×2",
          "equal",
          "×2",
          "×4")

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
md <- read_tsv('~/symlinks/comp/samples.txt') %>% dplyr::select(ID=Accession,Species,Group) %>% mutate(ID = gsub('_','',ID),Group = gsub('A. ','',Group))

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

# gene2go_long <- id %>%
#   filter(!is.na(go_raw), go_raw != "-", str_detect(go_raw, "GO:")) %>%
#   separate_rows(go_raw, sep="\\|") %>%
#   mutate(GO = str_extract(go_raw, "GO:\\d{7}")) %>%
#   filter(!is.na(GO)) %>%
#   distinct(Gene, GO, source_db) %>%
#   count(Gene, GO, name="n_sources") %>%
#   filter(n_sources >= 2) %>%
#   dplyr::select(Gene, GO)

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

sA <- intersect(sA, universe_genes)
sB <- intersect(sB, universe_genes)
uA <- intersect(uA, universe_genes)
uB <- intersect(uB, universe_genes)

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


plot_df <- bind_rows(
  mf_sA %>% mutate(Set = "Shared A"),
  mf_sB %>% mutate(Set = "Shared B"),
  mf_uA %>% mutate(Set = "Unique A"),
  mf_uB %>% mutate(Set = "Unique B")
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

plot_df %>% filter(Set == 'Unique A') %>% dplyr::select(TERM,p_adj,fold_enrich,sig10)
plot_df %>% filter(Set == 'Unique A' & grepl('manno|carb',TERM)) %>% dplyr::select(GO,bg_with_term,cand_with_term,TERM,p_adj,fold_enrich,sig10)
# GO         bg_with_term cand_with_term TERM                         p_adj fold_enrich sig10
# <chr>             <int>          <int> <chr>                        <dbl>       <dbl> <lgl>
#   1 GO:0030246           61             13 carbohydrate binding       0.00957        4.03 TRUE 
# 2 GO:0004559            5              4 alpha-mannosidase activity 0.0171        15.1  TRUE 

# plot
enrich_plot <- plot_df %>% 
  filter(Set == 'Unique A') %>% 
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
  scale_y_discrete(labels = function(x) str_wrap(x, width = 30)) +
  theme_bw(base_size = 8) +
  theme(
    legend.position = "right",
    strip.background = element_rect(fill = "grey95")
  )
enrich_plot

ggsave('~/symlinks/comp/figures/20260406_PurifyingSubgenomeA_Enrichment.pdf',enrich_plot,height=3.5,width=3.5)

# Extract those genes
cands <- gene2go_annot %>% filter(grepl('0030246|0004559',GO)) %>% filter(Gene %in% uA) %>% data.frame
write_tsv(cands,file='TopGenes_Puri_SubgenomeA_20260406.tsv')

```

