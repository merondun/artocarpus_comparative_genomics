## Subgenome-divided Orthofinder

This section continues the subgenome-specific evolution analyses, dealing with Artocarpus A/B (HART067), Batocarpus, and Morus.

The results of this section correspond to panel B showing missing orthogroups across subgenomes:

![orthogroups](/figures/20260413_Annotation_Counts_Orthogroups.png)



____



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

### Orthogroup missingness

Show orthogroup overlap by sample: 

```R
setwd('/project/coffea_pangenome/Artocarpus/Comparative_Paper/orthofinder/Subgenome_Divided/OrthoFinder/Results_Feb24/Orthogroups')
library(tidyverse)
library(stringr)
library(readr)

gc <- read_tsv("Orthogroups.GeneCount.tsv", show_col_types = FALSE)
md <- read_tsv('~/artocarpus_comparative_genomics/samples.txt') %>% mutate(Accession = gsub('_','',Accession)) %>% dplyr::select(Accession,Group,ord = Accession_Order)
md <- rbind(md, data.frame(Accession = c('Batocarpus','Morus'),Group = c('Batocarpus sp.','Morus mongolica'),ord = c(98,99)))

# First column is usually Orthogroup
og_col <- names(gc)[1]

long <- gc %>%
  pivot_longer(cols = -all_of(og_col), names_to = "Genome", values_to = "n_genes") %>%
  mutate(
    present = n_genes > 0,
    Accession = str_replace(Genome, "_[AB]$", ""),
    Subgenome = str_extract(Genome, "[AB]$") %>% replace_na("NA")
  )

missing_summary <- long %>%
  group_by(Genome, Accession, Subgenome) %>%
  summarise(
    n_orthogroups = n(),
    n_missing = sum(!present),
    pct_missing = 100 * mean(!present),
    .groups = "drop"
  )

# join to your metadata to color by Method / order nicely
missing_summary2 <- missing_summary %>%
  left_join(md %>% select(Accession, ord, Group), by = "Accession") %>% drop_na(Group) %>% 
  mutate(
    # make a nice label: Group (Accession_subgenome)
    ylab = if_else(Subgenome %in% c("A","B"),
                   paste0(Group, " (", Accession, "_", Subgenome, ")"),
                   paste0(Group, " (", Accession, ")")),
    ylab = fct_reorder(ylab, ord, .desc = TRUE)
  )

op <- ggplot(missing_summary2, aes(x = ylab, y = pct_missing, fill = Subgenome)) +
  geom_col(width = 0.8) +
  coord_flip() +
  theme_bw(base_size = 9) +
  scale_fill_manual(values=c('#f8766d','#00bf7d','black'))+
  labs(x = NULL, y = "% orthogroups missing", fill = "Subgenome")

ggsave('~/symlinks/comp/figures/20260413_Orthogroup_Missingness.pdf',op,height=5,width=4.5)

```

