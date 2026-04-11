## kSRates: Ks distributions

This prepares cleaned CDS FASTAs and the ARtocarpus GFF3, then runs kSRates (in manual mode) to estimate paralog and ortholog Ks distributions across Artocarpus–Batocarpus–Morus and perform rate adjustment before generating plots. 

This corresponds to panel B:

![ksrates](/figures/20260304_GeneTrees_BEAST_kS.png)

___

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

