# Artocarpus comparative genomics (assembly → annotation → subgenome evolution)

End-to-end comparative genomics for *Artocarpus* (and outgroups), spanning genome QC/assembly, Iso-Seq annotation, whole-genome alignments, orthology/ancestral reconstruction, and subgenome-aware evolutionary analyses (synteny, dN/dS, expression bias, Ks rate adjustment, and quartet discordance/BEAST dating).

Core sample metadata live in `samples.txt` (with ordering/plot aesthetics) and `Accessions.list` (accession subset used across pipelines).

## Directory map

- [**01_qa_qc_genomescope/**](01_qa_qc_genomescope/) — read/QC and GenomeScope summaries for genome size/heterozygosity context.
- **02_genome_assembly/** — assembly generation and post-processing notes/scripts.
- **03_annotation/** — Iso-Seq prep + eGAPx annotation, longest-isoform extraction, liftoff to other assemblies, standardized CDS/proteomes/GTF/GFF3 formatting.
- **04_whole_genome_alignments/** — nucmer-based WGA, filtering/merging into syntenic blocks, and karyotype-style visualization inputs.
- **05_orthofinder/** — orthogroup inference and species tree building.
- **06_ancestral_reconstruction/** — ancestral chromosome/karyotype reconstructions and lineage fusion/fission summaries.
- **07_subgenome_alignments/** — subgenome partitioning and BUSCO-level synteny comparisons (e.g., chromsyn inputs/plots).
- **08_subgenome_orthofinder_cafe5/** — subgenome-aware orthology and gene family evolution (CAFE5-ready inputs).
- **09_subgenome_dnds/** — subgenome-divided dN/dS along the tree (HyPhy/MG94), candidate scans, and plotting.
- **10_subgenome_expression/** — Salmon quantification + RBH mapping to compare A vs B expression bias across selection categories.
- **11_ksrates/** — kSRates (manual mode) for rate-adjusted Ks distributions and WGD timing context.
- **12_quartet_subgenome_gene_discordance_BEAST/** — Morus–Bato–ArtoA–ArtoB quartets: RBH sets, 4-fold sites, gene trees/ASTRAL, topology discordance, and BEAST divergence dating.

## Qs & Cs

Questions or comments reach out to Justin Merondun heritabilities [@] gmail.com or make an issue here. 