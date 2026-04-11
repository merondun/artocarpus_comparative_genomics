## Subgenome-specific Expression

This section quantifies RNA-seq/Iso-Seq read expression against the combined HART063 A+B CDS transcriptome with Salmon, maps transcripts to *Morus* orthologs via RBH lookups, and summarizes per-gene subgenome A vs B expression bias across dN/dS-based gene categories.

This corresponds to panel E:

![expressoin](/figures/20260407_SubgenomeEvolution.png)

___

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

