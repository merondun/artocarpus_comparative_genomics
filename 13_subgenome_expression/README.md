## 13_subgenome_expression/

This section quantifies RNA-seq/Iso-Seq read expression against the combined HART063 A+B CDS transcriptome with Salmon, maps transcripts to *Morus* orthologs via RBH lookups, and summarizes per-gene subgenome A vs B expression bias across dN/dS-based gene categories.

This corresponds to panel d:

![expression](/figures/panels/05_subgenome_evolution/20260420_SubgenomeEvolution.png)



____



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
genes <- read_tsv("~/artocarpus_comparative_genomics/09_subgenome_dnds/TopHomeologs_20260414.tsv")

# import dnds for shared copies of all paired genes 
wstats_clean <- fread('~/artocarpus_comparative_genomics/09_subgenome_dnds/Node_dNdS_20260420-RInput.tsv.gz') %>% as_tibble

# read in morus ~ arto/bato gene ID lookups
ids <- read_tsv('~/artocarpus_comparative_genomics/10_subgenome_expression/Morus_Lookup.tsv',col_names = F)
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


### batocarpus
samples_b <- samples %>% filter(prefix != "HART063")
files_b <- samples_b$quant
names(files_b) <- samples_b$sample
txi_b <- tximport(files_b, type = "salmon", txOut = TRUE, ignoreTxVersion = TRUE)

# convert to tpm martrix
bato_tpm <- as_tibble(txi_b$abundance, rownames = "RBH") %>%
  pivot_longer(-RBH, names_to = "sample", values_to = "TPM") %>%
  left_join(samples_b, by = "sample")

# map to hart
bato_tpm2 <- bato_tpm %>%
  left_join(ids, by = "RBH") 

bato_gene <- bato_tpm2 %>%
  filter(!is.na(Morus)) %>%
  group_by(Morus, tissue) %>%
  summarise(TPM_bato = sum(TPM, na.rm = TRUE), .groups = "drop")

# sumarize Arto A/B to Morus
hart_AB <- hart_tpm2 %>%
  filter(!is.na(Morus), !is.na(subgenome)) %>%
  group_by(Morus, tissue, subgenome) %>%
  summarize(TPM = sum(TPM, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = subgenome, values_from = TPM, values_fill = NA) %>%
  mutate(
    total_AB = A + B,
    log2_BA = log2((B + 1) / (A + 1))
  )
write_tsv(hart_AB, "~/artocarpus_comparative_genomics/10_subgenome_expression/hart_AB.tsv")
write_tsv(bato_gene, "~/artocarpus_comparative_genomics/10_subgenome_expression/batocarpus.tsv")
hart_AB <- read_tsv('~/artocarpus_comparative_genomics/10_subgenome_expression/hart_AB.tsv')
bato_gene <- read_tsv('~/artocarpus_comparative_genomics/10_subgenome_expression/batocarpus.tsv')

##### For each selection regime, compare arto/bato expression #####
outgroups <- wstats_clean %>% 
  filter(grepl('HART063|Bato',Branch)) %>% 
  dplyr::select(Gene,Branch,puri) %>% 
  pivot_wider(names_from = Branch,values_from = puri) %>% 
  drop_na(Batocarpus) %>% 
  mutate(purifying_arto = case_when(
    is.na(HART063B) & is.na(HART063A) ~ NA,
    is.na(HART063B) & !is.na(HART063A) ~ HART063A,
    !is.na(HART063B) & is.na(HART063A) ~ HART063B,
    HART063A == TRUE | HART063B == TRUE ~ TRUE,
    HART063A == FALSE & HART063B == FALSE ~ FALSE),
    subset = case_when(
      is.na(HART063B) & is.na(HART063A) ~ NA,
      is.na(HART063B) & !is.na(HART063A) ~ 'A-unique',
      !is.na(HART063B) & is.na(HART063A) ~ 'B-unique',
      !is.na(HART063B) & !is.na(HART063A) ~ 'Homeologs')
  ) %>% 
  drop_na(purifying_arto) %>%
  dplyr::select(Gene,purifying_bato = Batocarpus,purifying_arto,subset) %>% 
  mutate(category = case_when(
    purifying_bato == TRUE & purifying_arto == TRUE ~ 'Both',
    purifying_bato == FALSE & purifying_arto == TRUE ~ 'Arto-only',
    purifying_bato == TRUE & purifying_arto == FALSE ~ 'Bato-only',
    purifying_bato == FALSE & purifying_arto == FALSE ~ 'Neither',
    TRUE ~ NA
  ))

unique_subs <- c("A-unique","Homeologs","B-unique")
unique_cats <- c('Arto-only','Both','Neither','Bato-only')

unique_df <- hart_AB %>%
  dplyr::rename(Gene = Morus) %>%
  mutate(
    subset = case_when(
      is.na(A) ~ 'B-unique',
      is.na(B) ~ 'A-unique',
      !is.na(A) & !is.na(B) ~ 'Homeologs',
      TRUE ~ NA
    ),
    TPM_present = case_when(
      is.na(A) ~ B,
      is.na(B) ~ A,
      !is.na(A) & !is.na(B) ~ total_AB,
      TRUE ~ NA_real_
    )
  ) %>%
  drop_na(TPM_present) %>% 
  dplyr::select(Gene,tissue,subset,TPM_present) %>% 
  left_join(bato_gene %>% dplyr::rename(Gene = Morus), by = c("Gene", "tissue")) %>%
  filter(!is.na(TPM_bato)) %>%
  mutate(log2_Arto_vs_Bato = log2((TPM_present + 1) / (TPM_bato + 1)))

unique_gene <- unique_df %>%
  left_join(outgroups %>% dplyr::select(Gene, subset, category)) %>% 
  drop_na(category) %>% 
  group_by(Gene, subset, category) %>%
  summarise(mean_log2_Arto_vs_Bato = mean(log2_Arto_vs_Bato, na.rm = TRUE), .groups = "drop") %>%
  mutate(subset = factor(subset, levels = unique_subs),
         category = factor(category, levels = unique_cats))
write_tsv(unique_gene,'~/artocarpus_comparative_genomics/10_subgenome_expression/20260420_ExpressionSelectionlog2-ArtocarpusBatocarpus.tsv')

unique_w0 <- unique_gene %>%
  group_by(subset,category) %>%
  summarise(
    n = n(),
    p_value = wilcox.test(mean_log2_Arto_vs_Bato, mu = 0)$p.value,
    .groups = "drop"
  ) %>%
  mutate(
    p_adj = p.adjust(p_value, method = "holm"),
    sig = case_when(
      p_adj <= 0.05  ~ "*",
      TRUE ~ "ns"
    )
  )

unique_sum <- unique_gene %>%
  group_by(subset,category) %>%
  summarise(
    mean = mean(mean_log2_Arto_vs_Bato),
    sd = sd(mean_log2_Arto_vs_Bato),
    n = n(),
    se = sd / sqrt(n),
    conf_low = mean - 1.96 * se,
    conf_high = mean + 1.96 * se,
    .groups = "drop"
  ) %>%
  left_join(unique_w0) %>%
  mutate(y_lab = pmax(conf_high, 0) + 0.05 * diff(range(c(conf_low, conf_high), na.rm = TRUE)))

# subset    category     mean    sd     n     se conf_low conf_high  p_value   p_adj sig   y_lab
# <fct>     <fct>       <dbl> <dbl> <int>  <dbl>    <dbl>     <dbl>    <dbl>   <dbl> <chr> <dbl>
#   1 A-unique  Arto-only -0.0786  1.87    68 0.226   -0.522     0.365  0.246    0.492   ns    0.487
# 2 A-unique  Both      -0.447   1.65    97 0.167   -0.775    -0.119  0.00138  0.0124  *     0.122
# 3 A-unique  Neither   -0.246   1.84   432 0.0888  -0.420    -0.0722 0.000291 0.00320 *     0.122
# 4 A-unique  Bato-only -0.267   1.72   173 0.131   -0.524    -0.0108 0.00687  0.0481  *     0.122
# 5 Homeologs Arto-only  0.383   1.63   164 0.127    0.134     0.633  0.0101   0.0506  ns    0.754
# 6 Homeologs Both       0.182   1.63   425 0.0789   0.0271    0.337  0.127    0.380   ns    0.458
# 7 Homeologs Neither    0.0689  1.91   451 0.0899  -0.107     0.245  0.739    0.739   ns    0.367
# 8 Homeologs Bato-only  0.261   1.67   390 0.0848   0.0945    0.427  0.0151   0.0602  ns    0.549
# 9 B-unique  Arto-only -1.23    2.02    48 0.291   -1.81     -0.663  0.000121 0.00146 *     0.122
# 10 B-unique  Both      -0.489   1.36    67 0.166   -0.814    -0.164  0.00796  0.0481  *     0.122
# 11 B-unique  Neither   -0.364   2.06   263 0.127   -0.613    -0.115  0.00315  0.0252  *     0.122
# 12 B-unique  Bato-only -0.513   1.86   153 0.150   -0.808    -0.218  0.000821 0.00821 *     0.122

p_unique <- ggplot(unique_sum, aes(x = subset, y = mean, fill = category,shape=category)) +
  geom_hline(yintercept = 0, linetype = 2, color = "grey55") +
  geom_errorbar(aes(color = category, ymin = conf_low, ymax = conf_high), width = 0.12, linewidth = 0.7,position=position_dodge(width=0.5)) +
  geom_point(size = 1.8,position=position_dodge(width=0.5)) +
  geom_text(aes(y = y_lab, label = sig), color = "black", size = 2, vjust = 0,position=position_dodge(width=0.5)) +
  scale_fill_manual(values = viridis(4)) +
  scale_color_manual(values = viridis(4)) +
  scale_shape_manual(values=c(21,24,25,22))+
  labs(x = NULL, y = "Mean log2(Arto/Bato expression ratio)") +
  theme_bw(base_size = 10) +
  theme(legend.position = "top", panel.grid.minor = element_blank(),
        axis.text.x = element_text(face = "bold")) +
  coord_flip()

p_unique

ggsave('~/artocarpus_comparative_genomics/figures/20260420_Expression_Bias_Subgenome.pdf',p_unique,height=3.5,width=2.5)

```

