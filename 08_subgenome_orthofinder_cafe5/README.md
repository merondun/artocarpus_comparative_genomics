## Subgenome-divided Orthofinder

This section continues the subgenome-specific evolution analyses, dealing with Artocarpus A/B (HART067), Batocarpus, and Morus.

The results of this section correspond to panel A:

![cafe5](/figures/20260407_SubgenomeEvolution.png)

and the tree from C:

![busco](/figures/20260318_panel_BUSCO_Synteny_SubTreeSimple.png)

___

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
mdraw <-  read_tsv('~/symlinks/comp/samples.txt')
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

### Subgenome CAFE5

Prepare input files:

```R
# Count orthogroups, merge with annotation, plot species tree, plot cafe5 
library(tidyverse)
library(scales)
library(stringr)
library(GO.db)
library(AnnotationDbi)
library(data.table)
library(ggrepel)

ortho_dir="/project/coffea_pangenome/Artocarpus/Comparative_Paper/orthofinder/Subgenome_Divided/OrthoFinder/Results_Feb24/"

# read in orthofinder orthogroup - gene link 
og_map_raw <- read_tsv(
  "/project/coffea_pangenome/Artocarpus/Comparative_Paper/orthofinder/Subgenome_Divided/OrthoFinder/Results_Feb24/Orthogroups/Orthogroups.tsv",
  col_types = cols(.default = col_character()))

# long format: one protein per row with its Orthogroup
og_map_long <- og_map_raw %>%
  pivot_longer(-Orthogroup, names_to = "sample", values_to = "gene_list") %>%
  filter(!is.na(gene_list), gene_list != "") %>%
  separate_rows(gene_list, sep = ",\\s*") %>%
  transmute(Orthogroup, sample, protein_id = gene_list) %>%
  distinct()

# read in functional annotation for each gene 
files <- list.files('/project/coffea_pangenome/Artocarpus/Comparative_Paper/orthofinder/Subgenome_Divided/interproscan/output/',pattern = '*tsv',full.names = TRUE)

df <- list()
for (f in files) {
  id <- gsub('.*/','',gsub('.tsv','',f))
  cat ('Reading in ',id,'\n')
  f0 <- read_tsv(f,col_names=F)
  f1 <- f0 %>% mutate(sample = id) %>% 
    dplyr::select(sample, protein_id = X1, desc = X13, go_raw = X14) 
  df[[id]] <- f1
}
d1 <- rbindlist(df) %>% as_tibble


# parse go ids 
go_long <- d1 %>%
  filter(!is.na(go_raw), go_raw != "-") %>%
  # split multiple GO entries by '|'
  separate_rows(go_raw, sep = "\\|") %>%
  # strip trailing source tags
  mutate(go_id = str_extract(go_raw, "GO:\\d+")) %>%
  filter(!is.na(go_id)) %>%
  distinct(sample, protein_id, go_id, desc) 

# merge 
go_annot <- go_long %>%
  inner_join(og_map_long %>% dplyr::select(Orthogroup, protein_id), by = "protein_id") %>%
  distinct(sample, protein_id, Orthogroup, desc, go_id)

# ALSO, for each orthogroup, I want to add high-level classifications from GO 
go_map <- AnnotationDbi::select(
  GO.db,
  keys = unique(go_long$go_id),
  keytype = "GOID",
  columns = c("GOID", "TERM", "ONTOLOGY")
)

go_mer <- go_annot %>%
  left_join(go_map, by = c("go_id" = "GOID"))

# define categories using GO parent terms 
get_offspring <- function(go_id, ontology = c("BP","MF","CC")) {
  ontology <- match.arg(ontology)
  env <- switch(ontology,
                BP = GOBPOFFSPRING,
                MF = GOMFOFFSPRING,
                CC = GOCCOFFSPRING)
  res <- AnnotationDbi::mget(go_id, envir = env, ifnotfound = NA)
  unique(unlist(res))
}

# Core categories
cats <- tibble::tribble(
  ~label,                   ~go_id,      ~onto,
  "carbohydrate metabolism","GO:0005975","BP",
  "transport",              "GO:0006810","BP",
  "immune/defense",         "GO:0006952","BP",
  "signaling",              "GO:0007165","BP",
  "binding",                "GO:0005488","MF",
  "catalytic activity",     "GO:0003824","MF"
)

# Expand to include offspring sets
cat_offspring <- cats %>%
  rowwise() %>%
  mutate(offs = list(get_offspring(go_id, onto))) %>%
  ungroup()

# add "starch metabolism" explicitly
starch_offs <- get_offspring("GO:0005982", "BP")  # starch metabolic process
cat_offspring <- bind_rows(
  cat_offspring,
  tibble(label = "starch metabolism", go_id = "GO:0005982", onto = "BP", offs = list(starch_offs)))


# Define priority
priority <- c(
  "starch metabolism",
  "carbohydrate metabolism",
  "immune/defense",
  "transport",
  "signaling",
  "binding",
  "catalytic activity"
)

# Flatten offspring lists
cat_map_all <- cat_offspring %>%
  dplyr::select(label, offs) %>%
  unnest(offs) %>%
  dplyr::rename(go_id = offs)

# pre-resolve GO IDs in multiple categories
cat_map_resolved <- cat_map_all %>%
  mutate(label = factor(label, levels = priority, ordered = TRUE)) %>%
  arrange(go_id, label) %>%                  # earlier in priority order first
  distinct(go_id, .keep_all = TRUE) %>%      # keep the highest-priority mapping
  mutate(label = as.character(label))

# Deduplicate protein GO list
prot_go_unique <- go_mer %>%
  distinct(sample, protein_id, desc, Orthogroup, go_id)

# Assign categories, collapse per protein with priority
prot_labels <- prot_go_unique %>%
  left_join(cat_map_resolved, by = "go_id") %>%   # each go_id maps to one category
  group_by(sample,protein_id) %>%
  summarise(
    # pick the first non-NA label
    label = {
      if (all(is.na(label))) "unknown" else {
        # order labels by priority, then take the first
        lab_ord <- factor(label, levels = priority, ordered = TRUE)
        as.character(dplyr::first(label[order(lab_ord)]))
      }
    },
    .groups = "drop"
  ) 

# re-merge these high level labels with protein level descriptions 
full_df <- go_mer %>% 
  inner_join(prot_labels)

merged_df <- full_df %>% 
  dplyr::select(Orthogroup, TERM, label)  %>% 
  group_by(Orthogroup) %>%
  summarize(description = paste(unique(TERM), collapse = "; "),
            category = paste(unique(label), collapse = "; "), .groups = "drop") %>% 
  mutate(d = paste0(category,':',description)) %>% 
  dplyr::select(Orthogroup,d)

# Import gene counts
counts <- read_tsv(paste0(ortho_dir,'/Orthogroups/Orthogroups.GeneCount.tsv'))

cafe_input <- counts %>% 
  left_join(merged_df) %>% 
  dplyr::select(!Total) %>%
  relocate(any_of(c("d", "Orthogroup")), .before = 1) %>% 
  rename_with(~ gsub("_", "", .x)) %>%
  dplyr::select(-HART061A, -HART061B)

# Apply size filter into small and large gene families
species_cols <- colnames(cafe_input)[-(1:2)]
df_small <- cafe_input %>%
  filter(across(all_of(species_cols), ~ . < 100))  # all values < 100 in all species

df_large <- anti_join(cafe_input, df_small)

# Save outputs
write_tsv(df_small, "/project/coffea_pangenome/Artocarpus/Comparative_Paper/cafe5/Subgenome_Divided/counts.filtered_small.tsv")
write_tsv(df_large, "/project/coffea_pangenome/Artocarpus/Comparative_Paper/cafe5/Subgenome_Divided/counts.filtered_large.tsv")
write_tsv(merged_df, "/project/coffea_pangenome/Artocarpus/Comparative_Paper/cafe5/Subgenome_Divided/Orthogroup_Functions_20260304.tsv")
write_tsv(full_df, "/project/coffea_pangenome/Artocarpus/Comparative_Paper/cafe5/Subgenome_Divided/Gene_And_Orthogroup_Functions_20260318.tsv")


```

Run:

```bash
#!/bin/bash

#SBATCH --time=2-00:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=6Gb
#SBATCH --partition=ceres
#SBATCH --account=coffea_pangenome

# module load miniconda
# source activate isoseq_ann

for k in 1 2 3; do

    echo -e "\e[43m~~~~ Running cafe for ${k} ~~~~\e[0m"
    mkdir -p reps/${k}
    cafe5 --infile counts.filtered_small.tsv --tree resolved_ultrametric_binary_tree.nwk --cores 8 -k ${k} -o reps/${k}
    mkdir -p out_reps/${k}
    cafeplotter -i reps/${k}/ -o out_reps/${k} --format pdf

done
```

Save results across replicates:

```bash
find out_reps -type f -name "summary_all_gene_family.pdf" | while read filepath; do   dirname=$(basename "$(dirname "$filepath")");   cp "$filepath" "${dirname}_summary_all_gene_family.pdf"; done

# output lamba and lnl
find reps -type f \( -name "Gamma_results.txt" -o -name "Base_results.txt" \) -exec bash -c 'echo -e "{}\t$(head -n 2 "{}")"' \;
```

Best fit is k=3:

```
reps/1/Base_results.txt Model Base Final Likelihood (-lnL): 174519
Lambda: 0.0031132942119608
reps/2/Gamma_results.txt        Model Gamma Final Likelihood (-lnL): 170239
Lambda: 0.0036744642946421
reps/3/Gamma_results.txt        Model Gamma Final Likelihood (-lnL): 169646
Lambda: 0.0037433186987285
```

Extract top 10 orthogroups:

```
sed '1d' expansions_k2.txt | cut -f1 | head -n 10 > Top10Orthos.list
```

Plot

```R
##### Plot CAFE5 Subgenome #####
setwd('/project/coffea_pangenome/Artocarpus/Comparative_Paper/cafe5/Subgenome_Divided/')
library(tidyverse)
library(ggtree)
library(ape)
library(ggtreeExtra)
library(ggpubr)

# Import metadata
mdraw <-  read_tsv('~/symlinks/comp/samples.txt')
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

counter=0
for (k in c(1,2,3)) { 
  counter=counter+1
  
  dir <- file.path("reps", paste0(k))
  file <- list.files(dir, pattern = "clade_results\\.txt$",full.names = TRUE)
  
  # Import nodes 
  res <- read_tsv(file, comment = "#", col_names = c("Taxon_ID", "Increase", "Decrease"))
  res <- res %>%
    mutate(node_id = as.integer(str_extract(Taxon_ID, "(?<=<)\\d+(?=>)")),
           label = str_remove(Taxon_ID, "<\\d+>"))
  
  # Import tree
  tree <- read.tree('resolved_ultrametric_binary_tree.nwk')
  p <- ggtree(tree)  %<+% md +
    geom_tree()
  
  # Merge and plot 
  tree_df <- p$data
  
  node_data <- left_join(tree_df, res, by = c("node" = "node_id")) 
  tree_plot <- p 
  tree_plot$data <- tree_plot$data %>% mutate(label = ifelse(isTip == TRUE, gsub(' spp.','',gsub('A. ','',Group)), label))
  
  node_data_ordered <- node_data %>%
    dplyr::select(node, x, y, a=Increase, b=Decrease)  # Ensure "Increase" comes first
  bars <- nodebar(node_data_ordered, cols = c("a", "b"), color = c("#1b9e77", "#d95f02"), position='dodge')
  pies <- nodepie(node_data_ordered, cols = c("a", "b"), color = c("#1b9e77", "#d95f02"))
  
  label_df <- node_data %>% mutate(genes = paste0("+", Increase, "/-", Decrease))
  label_coords <- left_join(p$data, label_df) %>%
    drop_na(genes)
  
  # Add pie charts to internal nodes
  max_x <- max(p$data$x)
  cafe_tree <- tree_plot +
    #geom_tippoint(aes(fill = Haplotype, shape = Haplotype), size=2)+
    scale_fill_manual(values=cols)+
    scale_shape_manual(values=c(21,21,4,8))+
    geom_inset(pies,height=0.05,width=0.05) +
    geom_tiplab(offset = 0.06,size=3)+
    geom_text(
      data = label_coords,
      aes(x = x, y = y, label = genes),
      vjust = -0.5, hjust = 0.4, size = 2, fontface = "bold", color = "black"
    ) +
    xlim(c(0,max_x*1.2))+
    ggtitle(paste0('K = ',k))
  cafe_tree
  
  assign(paste0('p',counter),cafe_tree)
  
}

p + geom_nodelab(aes(label=node))+geom_tiplab()+xlim(c(0,max_x*1.2))

ggarrange(p1,p2,p3,nrow=3,ncol=1,common.legend = TRUE)
ggsave('~/symlinks/comp/figures/20260318_Cafe_Tree_Subgenomes.pdf',
       ggarrange(p1,p2,p3,nrow=3,ncol=1,common.legend = TRUE),
       height=20,width=8)

```

**Functional annotation of cafe5 expansion/loss** 

```
setwd("/project/coffea_pangenome/Artocarpus/Comparative_Paper/cafe5/")
# Subgenome-specific expansion cafe5: functional annotation of genes 
library(tidyverse)
library(scales)
library(stringr)
library(GO.db)
library(AnnotationDbi)
library(data.table)
library(ggrepel)
library(igraph)
library(ggraph)
library(RColorBrewer)

# import metadata
mdraw <-  read_tsv('~/symlinks/comp/samples.txt')
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

# After running cafe 5 below, load in the cafeplotter results.tsv and merge with annotation! 
info <- read_tsv('/project/coffea_pangenome/Artocarpus/Comparative_Paper/cafe5/Subgenome_Divided/Orthogroup_Functions_20260304.tsv')
res <- read_tsv('/project/coffea_pangenome/Artocarpus/Comparative_Paper/cafe5/Subgenome_Divided/out_reps/3/result_summary.tsv')
res_top <- res %>% 
  filter(grepl('24$',TaxonID)) %>% 
  filter(Pvalue < 0.01 & abs(Change) >= 3) %>% 
  left_join(info %>% dplyr::rename(FamilyID = Orthogroup))
res_top

# raw counts
ortho_dir="/project/coffea_pangenome/Artocarpus/Comparative_Paper/orthofinder/Subgenome_Divided/OrthoFinder/Results_Feb24/"
counts <- read_tsv(paste0(ortho_dir,'/Orthogroups/Orthogroups.GeneCount.tsv'))
orthos <- res_top %>% drop_na(d) %>% pull(FamilyID)

# ordering
species_order <- c("altilis","camansi","mariannensis","rigidus","odoratissimus","heterophyllus","nitidus","lacucha","dadah")

dat_long <- counts %>%
  filter(Orthogroup %in% 'OG0000015') %>%
  pivot_longer(-Orthogroup, names_to = "sample", values_to = "count") %>%
  filter(!grepl('Total|HART061|Batocarpus|Morus',sample)) %>% 
  mutate(sample_norm = gsub("_", "", sample)) %>%          
  left_join(md, by = c("sample_norm" = "Accession")) %>%
  mutate(Group = gsub('A. ','',Group),
         Group = factor(Group, levels = rev(species_order)))

raw_counts_og15 <- ggplot(dat_long, aes(x = count, y = Group, fill = Haplotype)) +
  geom_col(width = 0.8,position=position_dodge(width=0.8)) +
  facet_wrap(~ Orthogroup, scales = "free_x") +
  labs(x = "Gene count", y = NULL, fill = "Group") +
  scale_y_discrete(labels = function(x) str_wrap(x, 30)) +   # wrap long Group names
  theme_bw(base_size = 11) +
  scale_fill_manual(values=c('#f8766d','#00bf7d'))+
  theme(axis.text.y = element_text(size = 8),
        legend.position = "bottom")
raw_counts_og15
ggsave('~/symlinks/comp/figures/20260326_CountsOrthogroupLoss.pdf',raw_counts_og15,height=3.5,width=2.5)

# Plot network
edges <- res_top %>%
  filter(!is.na(d)) %>%
  separate_rows(d, sep = ";\\s*") %>%
  mutate(term = sub("^[^:]*:\\s*", "", d)) %>%
  dplyr::select(FamilyID, term) %>%
  distinct()

term_names <- na.omit(unique(edges$term))

g <- graph_from_data_frame(edges, directed = FALSE)
V(g)$type <- V(g)$name %in% term_names
V(g)$label <- V(g)$name

network <- ggraph(g, layout = "bipartite") +
  geom_edge_link(alpha = 0.4) +
  geom_node_point(aes(color = factor(type)), size = 3) +    # factor(type) gives a legend
  geom_node_text(aes(label = label), repel = TRUE, size = 3) +
  theme_void()+coord_flip()
network

ggsave('~/symlinks/comp/figures/20260326_NetworkOrthogroupLoss.pdf',network,height=3,width=2.5)
```

