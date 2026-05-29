## Quintets: Gene Discordance

This builds 5-taxon ortholog quintets (Ficus–Morus–Batocarpus–ArtocarpusA–ArtocarpusB) via reciprocal best-hit BLAST, aligns each quintet, and infers gene trees (IQ-TREE), and summarizes gene-tree/species-tree discordance by topology counting.

Primary output:

Panel b from:

![ksrates](/figures/panels/04_subgenome_history/20260528_GeneTrees_BEAST_kS.png)



___



Outputs

- Input quintet: `quintets/Ficus_Morus_Bato_ArtoA_ArtoB.unique.tsv`
- Per-quintetFASTAs: `fastas/*.fa` 
- Clean alignments: `alignments/*.aln.fa`
- Gene trees: `trees/*.treefile` 
- Topology counts: `topology_counts.tsv`

Taking the chr specific CDS files from the annotation section, divide into ArtoA / ArtoB/ Bato / Morus / Ficus  cds fastas. 

Compare the Batocarpus / Artocarpus A / Artocarpus B subgenome maps:

```bash
OD=/project/coffea_pangenome/Artocarpus/Comparative_Paper/duplication_orthology/copies_quartet
awk '{print $1}' Chr_Map.txt | xargs -I {} sh -c 'cat chrs/Artocarpus*{}*.fa' > ${OD}/ArtoA.fa
awk '{print $2}' Chr_Map.txt | xargs -I {} sh -c 'cat chrs/Artocarpus*{}*.fa' > ${OD}/ArtoB.fa
cat chrs/Batocarpus*fa > ${OD}/Bato.fa
cat chrs/Morus*fa > ${OD}/Morus.fa 
cat chrs/Ficus*fa > ${OD}/Ficus.fa 
```

And then to ensure each CDS is unique, ensure the haps for Artocarpus are indicated in CDS headers: 

```bash
sed -i 's/Artocarpus/ArtoA/g' ArtoA.fa
sed -i 's/Artocarpus/ArtoB/g' ArtoB.fa
sed -i 's/Batocarpus/Bato/g' Bato.fa
```

Run initial RBH: 

```bash
#!/bin/bash

#SBATCH --time=1-00:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --partition=ceres
#SBATCH --account=coffea_pangenome

set -euo pipefail

# CONFIG
WD=/project/coffea_pangenome/Artocarpus/Comparative_Paper/duplication_orthology/quintet
cd "${WD}"

THREADS=8
ANCHOR=Ficus

declare -A FASTA
FASTA[Ficus]="${WD}/Ficus.fa"
FASTA[Morus]="${WD}/Morus.fa"
FASTA[Bato]="${WD}/Bato.fa"
FASTA[ArtoA]="${WD}/ArtoA.fa"
FASTA[ArtoB]="${WD}/ArtoB.fa"

TARGETS=(Morus Bato ArtoA ArtoB)

mkdir -p db blast maps rbh quintets fastas work

echo "[`date`] Making BLAST databases..."

for sp in "${ANCHOR}" "${TARGETS[@]}"; do
    if [ ! -f "db/${sp}.nhr" ]; then
        makeblastdb -in "${FASTA[$sp]}" -dbtype nucl -out "db/${sp}"
    fi
done

echo "[`date`] Running BLAST searches..."

for sp in "${TARGETS[@]}"; do

    blastn -query "${FASTA[$ANCHOR]}" -db "db/${sp}" \
        -out "blast/${ANCHOR}_vs_${sp}.tsv" \
        -outfmt "6 qseqid sseqid pident length evalue bitscore" \
        -max_target_seqs 5 -evalue 1e-5 -num_threads "${THREADS}"

    blastn -query "${FASTA[$sp]}" -db "db/${ANCHOR}" \
        -out "blast/${sp}_vs_${ANCHOR}.tsv" \
        -outfmt "6 qseqid sseqid pident length evalue bitscore" \
        -max_target_seqs 5 -evalue 1e-5 -num_threads "${THREADS}"

done

echo "[`date`] Selecting best hits per query..."

best_hit() {
    in=$1
    out=$2

    sort -k1,1 -k6,6nr "${in}" \
        | awk -F'\t' '!seen[$1]++ {print $1"\t"$2}' \
        > "${out}"
}

for sp in "${TARGETS[@]}"; do
    best_hit "blast/${ANCHOR}_vs_${sp}.tsv" "maps/${ANCHOR}2${sp}.best"
    best_hit "blast/${sp}_vs_${ANCHOR}.tsv" "maps/${sp}2${ANCHOR}.best"
done

echo "[`date`] Computing pairwise RBHs with ${ANCHOR} anchor..."

for sp in "${TARGETS[@]}"; do

    join -t $'\t' -1 1 -2 2 \
        <(sort -k1,1 "maps/${ANCHOR}2${sp}.best") \
        <(sort -k2,2 "maps/${sp}2${ANCHOR}.best") \
        | awk -F'\t' '{print $1"\t"$2}' \
        > "rbh/${ANCHOR}_${sp}.rbh"

done

echo "[`date`] Intersecting RBHs to get 5-way 1:1:1:1:1 sets..."

cp "rbh/${ANCHOR}_${TARGETS[0]}.rbh" "quintets/tmp_${ANCHOR}_${TARGETS[0]}.tsv"

current="quintets/tmp_${ANCHOR}_${TARGETS[0]}.tsv"
name="${ANCHOR}_${TARGETS[0]}"

for sp in "${TARGETS[@]:1}"; do

    next="quintets/tmp_${name}_${sp}.tsv"

    join -t $'\t' -1 1 -2 1 \
        <(sort -k1,1 "${current}") \
        <(sort -k1,1 "rbh/${ANCHOR}_${sp}.rbh") \
        > "${next}"

    current="${next}"
    name="${name}_${sp}"

done

cp "${current}" "quintets/${ANCHOR}_Morus_Bato_ArtoA_ArtoB.tsv"

rm quintets/tmp_*.tsv

echo "[`date`] Enforcing no duplicated IDs within species columns..."

awk '
{
    keep=1
    for (i=1; i<=NF; i++) {
        if (seen[i][$i]++) {
            keep=0
        }
    }
    if (keep) print
}' "quintets/${ANCHOR}_Morus_Bato_ArtoA_ArtoB.tsv" \
  > "quintets/${ANCHOR}_Morus_Bato_ArtoA_ArtoB.unique.tsv"

echo "[`date`] 5-way RBH sets:"
wc -l "quintets/${ANCHOR}_Morus_Bato_ArtoA_ArtoB.unique.tsv"

```

Extract those genes and then align, extract 4-fold sites.

```bash
#!/bin/bash

#SBATCH --time=2-00:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=48G
#SBATCH --partition=ceres
#SBATCH --account=coffea_pangenome

# module load miniconda
# source activate beast
module load parallel

JOBS=22
WD=/project/coffea_pangenome/Artocarpus/Comparative_Paper/duplication_orthology/quintet

# Input CDS
FICUS_FA=${WD}/Ficus.fa
MORUS_FA=${WD}/Morus.fa
BATO_FA=${WD}/Bato.fa
ARTOA_FA=${WD}/ArtoA.fa
ARTOB_FA=${WD}/ArtoB.fa

QUIN_TSV="${WD}/quintets/Ficus_Morus_Bato_ArtoA_ArtoB.unique.tsv"   # 5 cols: FicusID MorusID BatoID ArtoAID ArtoBID

# Working/output dirs
RAW_DIR="${WD}/fastas"             # unaligned 5-way RBH fastas
OUT_DIR="${WD}/fastas_4fold"      # extracted 4-fold sites
OUT_DIR_STRICT="${WD}/fastas_4fold_strict"      # extracted 4-fold sites, les than 30% gaps 
WORK_DIR="${WD}/work_quintet"     # per-gene work
ALIGN_DIR="${WD}/alignments"      # aln_NT.clipkit.strict.fa copies here 
TREE_WORK="${WD}/tree_work"       # IQ-TREE temp prefix outputs
TREE_DIR="${WD}/trees"            # cleaned final gene trees (*.treefile)
ASTRAL_DIR="${WD}/astral"         # species tree outputs

mkdir -p "${RAW_DIR}" "${OUT_DIR}" "${OUT_DIR_STRICT}" "${WORK_DIR}" "${ALIGN_DIR}" "${TREE_WORK}" "${TREE_DIR}" "${ASTRAL_DIR}"
cd "${WD}"

# Tool checks
for tool in macse clipkit seqkit bioawk parallel samtools iqtree java; do
  command -v "$tool" >/dev/null 2>&1 || { echo "ERROR: $tool not found in PATH"; exit 1; }
done

# (A) Build 5-way RBH FASTAs 
echo "[$(date)] Indexing input genomes"
for fa in "${FICUS_FA}" "${MORUS_FA}" "${BATO_FA}" "${ARTOA_FA}" "${ARTOB_FA}"; do
  [ -f "${fa}.fai" ] || samtools faidx "${fa}"
done

# --- Parallel writing ---
echo "[$(date)] Writing 5-way RBH FASTAs to: ${RAW_DIR}"

write_quin() {
  local ficus_id="$1" morus_id="$2" bato_id="$3" artoa_id="$4" artob_id="$5"
  local outfa="${RAW_DIR}/${ficus_id}.fa"
  echo "  -> ${outfa}"

  {
    # Prefix headers with species so yields Bato/ArtoA/ArtoB/Morus
    samtools faidx "${FICUS_FA}" "${ficus_id}" | sed "1s/^>.*/>Ficus/"
    samtools faidx "${MORUS_FA}" "${morus_id}" | sed "1s/^>.*/>Morus/"
    samtools faidx "${BATO_FA}"  "${bato_id}"  | sed "1s/^>.*/>Bato/"
    samtools faidx "${ARTOA_FA}" "${artoa_id}" | sed "1s/^>.*/>ArtoA/"
    samtools faidx "${ARTOB_FA}" "${artob_id}" | sed "1s/^>.*/>ArtoB/"
  } > "${outfa}"
}

export -f write_quin
export RAW_DIR FICUS_FA MORUS_FA BATO_FA ARTOA_FA ARTOB_FA WD OUT_DIR OUT_DIR_STRICT

# Read 5 columns: MorusID BatoID ArtoAID ArtoBID
parallel --colsep '\t' --jobs "${JOBS}" --no-notice \
  write_quin {1} {2} {3} {4} {5} :::: "${QUIN_TSV}"

# Process one fasta
process_fa() {
  local fa="$1"
  local base
  base="$(basename "${fa%.fa}")"

  local wdir="${WORK_DIR}/${base}"
  mkdir -p "${wdir}"
  cd "${wdir}"

  # Validate we have the five species (names Bato/ArtoA/ArtoB/Morus/Ficus after stripping suffix)
  sed 's/_.*//g' "$fa" > raw.fa
  local required=(Bato ArtoA ArtoB Morus Ficus)
  for sp in "${required[@]}"; do
    if ! grep -q "^>${sp}\b" raw.fa; then
      echo "Skipping ${base}: missing sequence for ${sp}" >&2
      return 0
    fi
  done

  # 1) Codon align with MACSE
  macse -prog alignSequences \
    -seq  raw.fa \
    -out_NT aligned_NT.fasta > macse.log 2>&1

  # 2) Clean alignment (remove FS/Stop marks)
  macse -prog exportAlignment \
    -align aligned_NT.fasta \
    -codonForExternalFS --- \
    -codonForFinalStop --- \
    -codonForInternalFS --- \
    -codonForInternalStop --- \
    -charForRemainingFS - \
    -out_NT aligned_NT.clean.fasta

  # Confirm we still have the five samples post-clean
  local nsamp
  nsamp="$(grep -c '^>' aligned_NT.clean.fasta || true)"
  if [[ "$nsamp" -ne 5 ]]; then
    echo "Skipping ${base}: cleaned alignment has ${nsamp} sequences (expected 5)" >&2
    return 0
  fi

  # 3) Trim gappy codon positions
  clipkit aligned_NT.clean.fasta -o aln_NT.clipkit.fa --codon -m kpic > clipkit.log 2>&1

  # 4) Create a "strict" codon alignment that drops sequences with >30% gaps
  seqkit fx2tab aln_NT.clipkit.fa \
  | awk -F'\t' '
      {
        total = length($2)
        seq = $2
        gaps = gsub(/-/, "", seq)
        if (total > 0 && gaps / total <= 0.30)
          print $1
      }
    ' \
  | seqkit grep -f - aln_NT.clipkit.fa > aln_NT.clipkit.strict.fa
  
  # Require 5 samples with ATG aligned and extract 
  python /home/justin.merondun/merothon/merothon/scripts/Extract_4Fold.py -i aln_NT.clipkit.fa -m 5 -o extract
  cp extract.4fold.fa "${OUT_DIR}/${base}.fa"
  
  # and strict...
  python /home/justin.merondun/merothon/merothon/scripts/Extract_4Fold.py -i aln_NT.clipkit.strict.fa -m 5 -o extract_strict \
    || echo "WARN: strict extract failed for ${base}"
  [ -f extract_strict.4fold.fa ] && cp extract_strict.4fold.fa "${OUT_DIR_STRICT}/${base}.fa"
  [ -f extract_strict.4fold.fa ] && cp aln_NT.clipkit.strict.fa "${ALIGN_DIR}/${base}.aln.fa"
  rm -f *log raw_AA.fa raw.fa aligned_NT.fasta aligned_NT_AA.fasta
 
}

export -f process_fa
export WD RAW_DIR OUT_DIR OUT_DIR_STRICT WORK_DIR ALIGN_DIR

# Parallel
find "${RAW_DIR}" -type f -name '*.fa' -print0 \
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
         -nt 1 \
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
```

Count topologies:

```python
 
import argparse
import glob
from collections import Counter
from ete3 import Tree
 
 
def canonical_subtree(node):
    """Return a sorted topology string for any rooted subtree."""
    if node.is_leaf():
        return node.name
 
    children = [canonical_subtree(ch) for ch in node.children]
    children = sorted(children)
    return "(" + ",".join(children) + ")"
 
 
def classify(tree, outgroup):
    tree.set_outgroup(outgroup)
 
    leaves = sorted(tree.get_leaf_names())
 
    if outgroup not in leaves:
        raise ValueError(f"Outgroup {outgroup} not found in tree: {leaves}")
 
    ingroup = [t for t in leaves if t != outgroup]
 
    if len(ingroup) != 4:
        raise ValueError(f"Tree does not contain exactly 4 ingroup taxa: {leaves}")
 
    # After rooting on outgroup, the non-outgroup child of the root is the ingroup topology
    root_children = tree.children
    ingroup_nodes = [
        ch for ch in root_children
        if outgroup not in ch.get_leaf_names()
    ]
 
    if len(ingroup_nodes) != 1:
        raise ValueError(f"Could not identify single ingroup clade after rooting: {leaves}")
 
    ingroup_topo = canonical_subtree(ingroup_nodes[0])
 
    return f"{outgroup},{ingroup_topo}"
 
 
def main():
    parser = argparse.ArgumentParser(
        description="Classify rooted 5-taxon gene trees and output topology counts."
    )
    parser.add_argument(
        "-i", "--input",
        required=True,
        help="Directory containing .treefile gene trees"
    )
    parser.add_argument(
        "-o", "--outgroup",
        required=True,
        help="Name of the outgroup taxon, e.g. Ficus"
    )
    parser.add_argument(
        "--output",
        default="topology_counts.tsv",
        help="Output TSV path"
    )
 
    args = parser.parse_args()
 
    pattern = f"{args.input.rstrip('/')}/*.treefile"
    files = glob.glob(pattern)
 
    if not files:
        print(f"No .treefile files found in {args.input}")
        with open(args.output, "w") as outfh:
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
 
    print("\nTopology counts:")
    for topo in sorted(counts.keys()):
        print(f"{topo}: {counts[topo]}")
 
    total = sum(counts.values())
    print(f"\nTotal trees processed: {total}")
 
    with open(args.output, "w") as outfh:
        outfh.write("topology\tcount\n")
        for topo in sorted(counts.keys()):
            outfh.write(f"{topo}\t{counts[topo]}\n")
 
    print(f"\nWrote TSV: {args.output}")
 
 
if __name__ == "__main__":
    main()
```

Summarize:

```
python Count_Topologies.py -i trees/ -o Ficus --output topos
```

```
head topology_counts.tsv
topology        count
Ficus,(((ArtoA,ArtoB),Bato),Morus)      240
Ficus,(((ArtoA,ArtoB),Morus),Bato)      3
Ficus,(((ArtoA,Bato),ArtoB),Morus)      3625
Ficus,(((ArtoA,Bato),Morus),ArtoB)      49
```



### Plot Discordance

```R
setwd('/project/coffea_pangenome/Artocarpus/Comparative_Paper/duplication_orthology/quintet')
library(tidyverse)
library(scales)
library(patchwork)
library(ape)
library(ggtree)

top_n <- 4

df <- read_tsv(
  "topology_counts.tsv",
  col_names = c("Topology", "Count"),
  show_col_types = FALSE
) %>%
  filter(Topology != "topology") %>%
  mutate(Count = as.integer(Count)) %>%
  arrange(desc(Count)) %>%
  slice_head(n = top_n) %>%
  mutate(
    Topology = factor(Topology, levels = Topology),
    Label = paste0("Topology ", row_number())
  )

palette <- scales::hue_pal()(top_n)

make_tree_plot <- function(newick, color, title = NULL) {
  tr <- read.tree(text = paste0("(", newick, ");"))
  gt <- ggtree(tr, size = 0.5, color = color) 
  ggtree(tr, size = 0.5, color = color) +
    geom_tiplab(size = 3) +
    labs(title = title) +
    theme_void() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 11, face = "bold"),
      plot.margin = margin(4, 4, 4, 4)
    ) +
    xlim(0,max(gt@data$x)*.25+max(gt@data$x))
}

tree_plots <- pmap(
  list(df$Topology, palette, df$Label),
  make_tree_plot
)

p_trees <- wrap_plots(tree_plots, nrow = 1)

p_bar <- ggplot(df, aes(x = Label, y = Count, fill = Label)) +
  geom_col(width = 0.7, color = "grey20") +
  geom_text(
    aes(label = comma(Count)),
    vjust = -0.35,
    size = 4.5,
    fontface = "bold"
  ) +
  scale_fill_manual(values = setNames(palette, df$Label)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.10))) +
  labs(x = NULL, y = "Count") +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 11),
    plot.margin = margin(8, 8, 12, 8)
  )

final_plot <- p_trees / p_bar +
  plot_layout(heights = c(0.38, 0.62))

final_plot

ggsave(
  "~/symlinks/comp/figures/20260515_WGD_Topology_Discordance_top4.pdf",
  final_plot,
  width = 6,
  height = 4
)

```

