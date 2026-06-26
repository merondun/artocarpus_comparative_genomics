# Artocarpus comparative genomics (assembly → annotation → subgenome evolution)

End-to-end comparative genomics for *Artocarpus* (and outgroups), spanning genome QC/assembly, Iso-Seq annotation, whole-genome alignments, orthology/ancestral reconstruction, and subgenome-aware evolutionary analyses (synteny, dN/dS, expression bias, Ks rate adjustment, and quartet discordance/BEAST dating).

Sample metadata live in `samples.txt` (with ordering/plot aesthetics) and `Accessions.list` (accession subset used across pipelines).

## Directory map

- [**01_qa_qc_genomescope/**](01_qa_qc_genomescope/) — read/QC and GenomeScope summaries for genome size/heterozygosity context.
- [**02_genome_assembly/**](02_genome_assembly/) — assembly generation and post-processing notes/scripts.
- [**03_annotation/**](03_annotation/) — Iso-Seq prep + eGAPx annotation & repeat annotation. Includes longest-isoform extraction, liftoff to other assemblies, standardized CDS/proteomes/GTF/GFF3 formatting, repeat and LTR variation.
- [**04_whole_genome_alignments/**](04_whole_genome_alignments/) — nucmer-based WGA, filtering/merging into syntenic blocks, and karyotype-style visualization inputs.
- [**05_orthofinder/**](05_orthofinder/) — orthogroup inference and species tree building.
- [**06_ancestral_reconstruction/**](06_ancestral_reconstruction/) — ancestral chromosome/karyotype reconstructions and lineage fusion/fission summaries.
- [**07_subgenome_alignments/**](07_subgenome_alignments/) — subgenome partitioning and BUSCO-level synteny comparisons (e.g., chromsyn inputs/plots).
- [**08_subgenome_orthofinder/**](08_subgenome_orthofinder/) — subgenome-aware orthology. 
- [**09_subgenome_dnds/**](09_subgenome_dnds/) — subgenome-divided dN/dS along the tree (HyPhy/MG94), candidate scans, and plotting.
- [**10_subgenome_beast/**](10_subgenome_beast/) — BEAST trees with subgenome-divided RBH hits to Ficus, including Morus + Batocarpus. 
- [**11_ksrates/**](11_ksrates/) kSRates (manual mode) for rate-adjusted Ks distributions and WGD timing context.
- [**12_subgenome_quintet_discordance/**](12_subgenome_quintet_discordance/) — Ficus–Morus–Bato–ArtoA–ArtoB quartets: RBH sets, gene trees, topology discordance.
- [**13_subgenome_expression/**](13_subgenome_expression/) — Salmon quantification + RBH mapping to compare A vs B expression bias across selection categories.

## Qs & Cs

Questions or comments reach out to Justin Merondun heritabilities [@] gmail.com or make an issue here. 
