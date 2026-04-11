## Quartet Gene Discordance & BEAST

This builds 4-taxon ortholog quartets (Morus–Batocarpus–ArtocarpusA–ArtocarpusB) via reciprocal best-hit BLAST, aligns each quartet, extracts shared 4-fold degenerate sites, and infers gene trees (IQ-TREE) and a quartet species tree (ASTRAL), and summarizes gene-tree/species-tree discordance by topology counting, and concatenates 4-fold sites for BEAST divergence dating.

This corresponds to A & C:

![quartets](/figures/20260304_GeneTrees_BEAST_kS.png)

___

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
md = read_tsv('~/symlinks/comp/samples.txt')
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

