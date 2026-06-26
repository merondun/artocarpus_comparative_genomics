# 03_annotation/ 

This section deals with both Iso-seq gene annotation on the assemblies, including liftovers, as well as repeat annotation using EarlGrey. 

In the end, can produce these summary figures:

Repeat summaries:

![repeats](/figures/20260410_repeats.png)

Gene annotation: 

![annotation_counts](/figures/20260413_Annotation_Counts_Orthogroups.png)



## Iso-seq Gene Annotation

We have Isoseq data for 2 samples (Artocarpus camansi (6 tissues) and Batocarpus sp. (2 tissues)).

These steps will demultiplex Iso‑Seq reads, convert selected partitions to FASTA, and run eGAPx (with optional short reads just for sensitivity, it was found that Isoseq is sufficient) to generate gene/transcript models. Extract the longest isoform per gene and predict CDS/proteins with TransDecoder, then liftover reference annotations to other assemblies. Finally clean/standardize headers, produce mapping tables, and split CDS/proteomes by chromosome for downstream comparative analyses. 

Primary outputs after this massive codeblock: 

- eGAPx annotations: complete.genomic.gtf and run out/work directories.
- Predicted coding sequences/proteomes: .transdecoder.cds and .transdecoder.pep.
- Longest‑isoform exports: *.longest_transcript_per_gene.gtf and .fa.
- Liftoff transfers: per‑sample GFF3 liftover files.
- Cleaned inputs for downstream: proteomes/, cds/, gtf/ (clean headers) + genes.tsv mappings and chromosome‑split FASTA files.

### Demultiplexing & Input Prep

Split with [pbskera](https://skera.how/) into s-reads, and demultiplex with lima, using example [here](https://skera.how/examples.html). 

This uses the standard mas16 primer file:

```
cat mas16_primers.fasta 
>A
AGCTTACTTGTGAAGA
>B
ACTTGTAAGCTGTCTA
>C
ACTCTGTCAGGTCCGA
>D
ACCTCCTCCTCCAGAA
>E
AACCGGACACACTTAG
>F
AGAGTCCAATTCGCAG
>G
AATCAAGGCTTAACGG
>H
ATGTTGAATCCTAGCG
>I
AGTGCGTTGCGAATTG
>J
AATTGCGTAGTTGGCC
>K
ACACTTGGTCGCAATC
>L
AGTAAGCCTTCGTGTC
>M
ACCTAGATCAGAGCCT
>N
AGGTATGCCGGTTAAG
>O
AAGTCACCGGCACCTT
>P
ATGAAGTGGCTCGAGA
>Q
AGTAGCTGTGTGCA
```

And this isoseq barcode file for demultiplexing:

```
cat Barcodes.fa 
>IsoSeqX_bc01_5p
CTACACGACGCTCTTCCGATCTACTACACGCAATGAAGTCGCAGGGTTGGG
>IsoSeqX_bc02_5p
CTACACGACGCTCTTCCGATCTACTAGTAGCAATGAAGTCGCAGGGTTGGG
>IsoSeqX_bc03_5p
CTACACGACGCTCTTCCGATCTAGTGTACGCAATGAAGTCGCAGGGTTGGG
>IsoSeqX_bc04_5p
CTACACGACGCTCTTCCGATCTATCACTAGCAATGAAGTCGCAGGGTTGGG
>IsoSeqX_bc05_5p
CTACACGACGCTCTTCCGATCTCAGCTGTGCAATGAAGTCGCAGGGTTGGG
>IsoSeqX_bc06_5p
CTACACGACGCTCTTCCGATCTCAGTCACGCAATGAAGTCGCAGGGTTGGG
>IsoSeqX_bc07_5p
CTACACGACGCTCTTCCGATCTCATGTATGCAATGAAGTCGCAGGGTTGGG
>IsoSeqX_bc08_5p
CTACACGACGCTCTTCCGATCTCGTATGTGCAATGAAGTCGCAGGGTTGGG
>IsoSeqX_bc09_5p
CTACACGACGCTCTTCCGATCTGACATGTGCAATGAAGTCGCAGGGTTGGG
>IsoSeqX_bc10_5p
CTACACGACGCTCTTCCGATCTGAGTCTAGCAATGAAGTCGCAGGGTTGGG
>IsoSeqX_bc11_5p
CTACACGACGCTCTTCCGATCTGTAGATAGCAATGAAGTCGCAGGGTTGGG
>IsoSeqX_bc12_5p
CTACACGACGCTCTTCCGATCTGTATGACGCAATGAAGTCGCAGGGTTGGG
>IsoSeqX_3p
AAGCAGTGGTATCAACGCAGAGTAC
```

Run, submit positional for library:

```bash
#!/bin/bash

#SBATCH --time=48:00:00   
#SBATCH --nodes=1  
#SBATCH --cpus-per-task=10
#SBATCH --partition=ceres

module load miniconda
source activate isoseq_ann

#submit e.g. sbatch 00_Skera_Demul.sh m84125_250429_221012_s3.hifi_reads.bcM0002
if [ -z "$1" ]; then
    echo "Error: Library prefix file positional argument is required."
    exit 1
fi

LIBRARY=$1

echo -e "\e[43m~~~~ Demultiplexing isoseq file: ${LIBRARY} ~~~~\e[0m"
skera split ${LIBRARY}.bam mas16_primers.fasta ${LIBRARY}.skera.bam
lima ${LIBRARY}.skera.bam Barcodes.fa ${LIBRARY}.skera.bam --isoseq --overwrite-biosample-names
```

Copy those split and demultiplexed files for annotation:

```bash
#!/bin/bash

#SBATCH --time=2-00:00:00   
#SBATCH --nodes=1  
#SBATCH --ntasks-per-node=16
#SBATCH --partition=ceres
 
samtools fastq -@ 16 m84125_250429_221012_s3.hifi_reads.bcM0001.skera.IsoSeqX_bc01_5p--IsoSeqX_3p.bam | bgzip -c >  /project/coffea_pangenome/Artocarpus/Assemblies/20250101_JustinAssemblies/annotation/isoseq/raw_fastq_nopolyA/N15_23__YL.fastq.gz
samtools fastq -@ 16 m84125_250429_221012_s3.hifi_reads.bcM0001.skera.IsoSeqX_bc02_5p--IsoSeqX_3p.bam | bgzip -c >  /project/coffea_pangenome/Artocarpus/Assemblies/20250101_JustinAssemblies/annotation/isoseq/raw_fastq_nopolyA/N15_23__ML.fastq.gz
samtools fastq -@ 16 m84125_250429_221012_s3.hifi_reads.bcM0001.skera.IsoSeqX_bc03_5p--IsoSeqX_3p.bam | bgzip -c >  /project/coffea_pangenome/Artocarpus/Assemblies/20250101_JustinAssemblies/annotation/isoseq/raw_fastq_nopolyA/HART063__MF.fastq.gz
samtools fastq -@ 16 m84125_250429_221012_s3.hifi_reads.bcM0001.skera.IsoSeqX_bc04_5p--IsoSeqX_3p.bam | bgzip -c >  /project/coffea_pangenome/Artocarpus/Assemblies/20250101_JustinAssemblies/annotation/isoseq/raw_fastq_nopolyA/HART063__FF.fastq.gz
samtools fastq -@ 16 m84125_250429_221012_s3.hifi_reads.bcM0001.skera.IsoSeqX_bc05_5p--IsoSeqX_3p.bam | bgzip -c >  /project/coffea_pangenome/Artocarpus/Assemblies/20250101_JustinAssemblies/annotation/isoseq/raw_fastq_nopolyA/HART063__YL.fastq.gz
samtools fastq -@ 16 m84125_250429_221012_s3.hifi_reads.bcM0001.skera.IsoSeqX_bc06_5p--IsoSeqX_3p.bam | bgzip -c >  /project/coffea_pangenome/Artocarpus/Assemblies/20250101_JustinAssemblies/annotation/isoseq/raw_fastq_nopolyA/HART063__ML.fastq.gz
samtools fastq -@ 16 m84125_250429_221012_s3.hifi_reads.bcM0001.skera.IsoSeqX_bc07_5p--IsoSeqX_3p.bam | bgzip -c >  /project/coffea_pangenome/Artocarpus/Assemblies/20250101_JustinAssemblies/annotation/isoseq/raw_fastq_nopolyA/HART063__FR.fastq.gz
samtools fastq -@ 16 m84125_250429_221012_s3.hifi_reads.bcM0002.skera.IsoSeqX_bc08_5p--IsoSeqX_3p.bam | bgzip -c >  /project/coffea_pangenome/Artocarpus/Assemblies/20250101_JustinAssemblies/annotation/isoseq/raw_fastq_nopolyA/HART063__SD.fastq.gz
```

### eGAPx Annotation

See issues about long read data here: https://github.com/ncbi/egapx/issues/188 

To accomodate this, convert the long read data to fasta:

```bash
#!/bin/bash
#SBATCH --time=1-00:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --partition=ceres
#SBATCH --account=coffea_pangenome

# submit with 
# for f in *.fastq; do sbatch convert_fq_to_fa.sh "$f"; done

FILE="$1"
OUT="${FILE%.fastq}.fasta"

echo "Converting $FILE -> $OUT"
```

After:

```
lt *fasta
-rw-r-----. 1 justin.merondun proj-coffea_pangenome 6.1G Jan 21 18:42 HART063__FF.fasta
-rw-r-----. 1 justin.merondun proj-coffea_pangenome 6.1G Jan 21 18:42 HART063__SD.fasta
-rw-r-----. 1 justin.merondun proj-coffea_pangenome 7.7G Jan 21 18:42 HART063__ML.fasta
-rw-r-----. 1 justin.merondun proj-coffea_pangenome 7.4G Jan 21 18:42 HART063__FR.fasta
-rw-r-----. 1 justin.merondun proj-coffea_pangenome 7.9G Jan 21 18:42 HART063__YL.fasta
-rw-r-----. 1 justin.merondun proj-coffea_pangenome 8.4G Jan 21 18:42 N15_23__YL.fasta
-rw-r-----. 1 justin.merondun proj-coffea_pangenome  12G Jan 21 18:43 N15_23__ML.fasta
-rw-r-----. 1 justin.merondun proj-coffea_pangenome  15G Jan 21 18:43 HART063__MF.fasta
```

I will also try with some available bulk RNAseq from the previous genome:

```bash
total 187G
-rw-r-----. 1 justin.merondun proj-coffea_pangenome 5.3G Jan 20 17:54 external_SRR5997516_1.fastq
-rw-r-----. 1 justin.merondun proj-coffea_pangenome 5.3G Jan 20 17:55 external_SRR5997516_2.fastq
-rw-r-----. 1 justin.merondun proj-coffea_pangenome 5.6G Jan 20 17:55 external_SRR5997517_1.fastq
-rw-r-----. 1 justin.merondun proj-coffea_pangenome 5.6G Jan 20 17:55 external_SRR5997517_2.fastq
-rw-r-----. 1 justin.merondun proj-coffea_pangenome 5.3G Jan 20 17:55 external_SRR5997518_1.fastq
-rw-r-----. 1 justin.merondun proj-coffea_pangenome 5.3G Jan 20 17:55 external_SRR5997518_2.fastq
-rw-r-----. 1 justin.merondun proj-coffea_pangenome 2.2G Jan 20 17:55 external_SRR5997519_1.fastq
-rw-r-----. 1 justin.merondun proj-coffea_pangenome 2.2G Jan 20 17:55 external_SRR5997519_2.fastq
-rw-r-----. 1 justin.merondun proj-coffea_pangenome 2.9G Jan 20 17:55 external_SRR5997520_1.fastq
-rw-r-----. 1 justin.merondun proj-coffea_pangenome 2.9G Jan 20 17:55 external_SRR5997520_2.fastq
-rw-r-----. 1 justin.merondun proj-coffea_pangenome 6.3G Jan 20 17:55 external_SRR5997521_1.fastq
-rw-r-----. 1 justin.merondun proj-coffea_pangenome 6.3G Jan 20 17:55 external_SRR5997521_2.fastq
-rw-r-----. 1 justin.merondun proj-coffea_pangenome 7.4G Jan 20 17:56 external_SRR5997522_1.fastq
-rw-r-----. 1 justin.merondun proj-coffea_pangenome 7.4G Jan 20 17:56 external_SRR5997522_2.fastq
-rw-r-----. 1 justin.merondun proj-coffea_pangenome 5.8G Jan 20 17:56 external_SRR5997523_1.fastq
-rw-r-----. 1 justin.merondun proj-coffea_pangenome 5.8G Jan 20 17:56 external_SRR5997523_2.fastq
-rw-r-----. 1 justin.merondun proj-coffea_pangenome 6.5G Jan 20 17:56 external_SRR5997524_1.fastq
-rw-r-----. 1 justin.merondun proj-coffea_pangenome 6.5G Jan 20 17:56 external_SRR5997524_2.fastq
-rw-r-----. 1 justin.merondun proj-coffea_pangenome 7.8G Jan 20 17:56 external_SRR5997525_1.fastq
-rw-r-----. 1 justin.merondun proj-coffea_pangenome 7.8G Jan 20 17:56 external_SRR5997525_2.fastq
-rw-r-----. 1 justin.merondun proj-coffea_pangenome 7.6G Jan 20 17:56 external_SRR5997526_1.fastq
-rw-r-----. 1 justin.merondun proj-coffea_pangenome 7.6G Jan 20 17:57 external_SRR5997526_2.fastq
-rw-r-----. 1 justin.merondun proj-coffea_pangenome 3.6G Jan 20 17:57 external_SRR5997527_1.fastq
-rw-r-----. 1 justin.merondun proj-coffea_pangenome 3.6G Jan 20 17:57 external_SRR5997527_2.fastq
-rw-r-----. 1 justin.merondun proj-coffea_pangenome 7.6G Jan 20 17:57 external_SRR5997536_1.fastq
-rw-r-----. 1 justin.merondun proj-coffea_pangenome 7.6G Jan 20 17:57 external_SRR5997536_2.fastq
-rw-r-----. 1 justin.merondun proj-coffea_pangenome 8.2G Jan 20 17:57 external_SRR5997537_1.fastq
-rw-r-----. 1 justin.merondun proj-coffea_pangenome 8.2G Jan 20 17:57 external_SRR5997537_2.fastq
-rw-r-----. 1 justin.merondun proj-coffea_pangenome 6.8G Jan 20 17:57 external_SRR5997538_1.fastq
-rw-r-----. 1 justin.merondun proj-coffea_pangenome 6.8G Jan 20 17:57 external_SRR5997538_2.fastq
-rw-r-----. 1 justin.merondun proj-coffea_pangenome 5.0G Jan 20 17:58 external_SRR5997539_1.fastq
-rw-r-----. 1 justin.merondun proj-coffea_pangenome 5.0G Jan 20 17:58 external_SRR5997539_2.fastq
```

Have a long reads file with the paths to the long read fastas:

```bash
cat Hfastas.txt 
HFF     /90daydata/coffea_pangenome/scratch/egapx/isoseq_rename/HART063__FF.fasta
HFR     /90daydata/coffea_pangenome/scratch/egapx/isoseq_rename/HART063__FR.fasta
HMF     /90daydata/coffea_pangenome/scratch/egapx/isoseq_rename/HART063__MF.fasta
HML     /90daydata/coffea_pangenome/scratch/egapx/isoseq_rename/HART063__ML.fasta
HSD     /90daydata/coffea_pangenome/scratch/egapx/isoseq_rename/HART063__SD.fasta
HYL     /90daydata/coffea_pangenome/scratch/egapx/isoseq_rename/HART063__YL.fasta
```

and a short reads file:

```bash
cat Shortreads.txt 
S1      /90daydata/coffea_pangenome/scratch/egapx/rnaseq/external_SRR5997516_1.fastq
S1      /90daydata/coffea_pangenome/scratch/egapx/rnaseq/external_SRR5997516_2.fastq
S2      /90daydata/coffea_pangenome/scratch/egapx/rnaseq/external_SRR5997517_1.fastq
S2      /90daydata/coffea_pangenome/scratch/egapx/rnaseq/external_SRR5997517_2.fastq
S3      /90daydata/coffea_pangenome/scratch/egapx/rnaseq/external_SRR5997518_1.fastq
S3      /90daydata/coffea_pangenome/scratch/egapx/rnaseq/external_SRR5997518_2.fastq
S4      /90daydata/coffea_pangenome/scratch/egapx/rnaseq/external_SRR5997519_1.fastq
S4      /90daydata/coffea_pangenome/scratch/egapx/rnaseq/external_SRR5997519_2.fastq
S5      /90daydata/coffea_pangenome/scratch/egapx/rnaseq/external_SRR5997520_1.fastq
S5      /90daydata/coffea_pangenome/scratch/egapx/rnaseq/external_SRR5997520_2.fastq
S6      /90daydata/coffea_pangenome/scratch/egapx/rnaseq/external_SRR5997521_1.fastq
S6      /90daydata/coffea_pangenome/scratch/egapx/rnaseq/external_SRR5997521_2.fastq
S7      /90daydata/coffea_pangenome/scratch/egapx/rnaseq/external_SRR5997522_1.fastq
S7      /90daydata/coffea_pangenome/scratch/egapx/rnaseq/external_SRR5997522_2.fastq
S8      /90daydata/coffea_pangenome/scratch/egapx/rnaseq/external_SRR5997523_1.fastq
S8      /90daydata/coffea_pangenome/scratch/egapx/rnaseq/external_SRR5997523_2.fastq
S9      /90daydata/coffea_pangenome/scratch/egapx/rnaseq/external_SRR5997524_1.fastq
S9      /90daydata/coffea_pangenome/scratch/egapx/rnaseq/external_SRR5997524_2.fastq
S10     /90daydata/coffea_pangenome/scratch/egapx/rnaseq/external_SRR5997525_1.fastq
S10     /90daydata/coffea_pangenome/scratch/egapx/rnaseq/external_SRR5997525_2.fastq
S11     /90daydata/coffea_pangenome/scratch/egapx/rnaseq/external_SRR5997526_1.fastq
S11     /90daydata/coffea_pangenome/scratch/egapx/rnaseq/external_SRR5997526_2.fastq
S12     /90daydata/coffea_pangenome/scratch/egapx/rnaseq/external_SRR5997527_1.fastq
S12     /90daydata/coffea_pangenome/scratch/egapx/rnaseq/external_SRR5997527_2.fastq
S13     /90daydata/coffea_pangenome/scratch/egapx/rnaseq/external_SRR5997536_1.fastq
S13     /90daydata/coffea_pangenome/scratch/egapx/rnaseq/external_SRR5997536_2.fastq
S14     /90daydata/coffea_pangenome/scratch/egapx/rnaseq/external_SRR5997537_1.fastq
S14     /90daydata/coffea_pangenome/scratch/egapx/rnaseq/external_SRR5997537_2.fastq
S15     /90daydata/coffea_pangenome/scratch/egapx/rnaseq/external_SRR5997538_1.fastq
S15     /90daydata/coffea_pangenome/scratch/egapx/rnaseq/external_SRR5997538_2.fastq
S16     /90daydata/coffea_pangenome/scratch/egapx/rnaseq/external_SRR5997539_1.fastq
S16     /90daydata/coffea_pangenome/scratch/egapx/rnaseq/external_SRR5997539_2.fastq
```

Run egapx:

```bash
#!/bin/bash

#SBATCH --time=2-00:00:00   
#SBATCH --nodes=1  
#SBATCH --cpus-per-task=1
#SBATCH --mem=512Gb
#SBATCH --partition=ceres
#SBATCH --account=coffea_pangenome

module load miniconda
source activate egapx
module load nextflow
module load apptainer

# Variables for genome to annotate, and which isoseq reads to use 
ID=N15_23
ISOdata=ISOfastas
SRdata=None

# Static
RUN="${ID}_${ISOdata}_${SRdata}"
EGAPDIR=/project/coffea_pangenome/Software/Merondun/egapx
WD=/90daydata/coffea_pangenome/scratch/egapx
GENOMES=/project/coffea_pangenome/Artocarpus/Comparative_Paper/assemblies/unmasked
echo "WORKING ON ${RUN}"

# Create yaml for run 
YAML=${WD}/${RUN}.yaml
cat > "$YAML" <<EOF
genome: ${GENOMES}/${ID}.fa
taxid: 241863
long_reads: ${WD}/${ISOdata}.txt
cmsearch:
  enabled: false
trnascan:
  enabled: false
EOF

#long_reads: ${WD}/${ISOdata}.txt
#short_reads: ${WD}/Shortreads.txt
#241863 Bato  194251 Arto altilis 709039 camansi

# Run 
python3 ${EGAPDIR}/ui/egapx.py ${YAML} \
    -c ${EGAPDIR}/egapx_config/ \
    -e slurm \
    -w ${WD}/work/${RUN} -o ${WD}/out/${RUN}
```

Outputs

```bash
./N15_23_Nfastas_None/complete.genomic.gtf 21491
./HART063_Hfastas/complete.genomic.gtf  34139
./HART063_Hfastas_Sreads/complete.genomic.gtf   34914
```

For some analyses (Orthofinder, CAFE5), we only want a single transcript per gene. Use this to extract this, for the references:

```bash
#!/bin/bash

#SBATCH --time=2-00:00:00   
#SBATCH --nodes=1  
#SBATCH --cpus-per-task=6
#SBATCH --mem=32Gb
#SBATCH --partition=ceres
#SBATCH --account=coffea_pangenome

# module load miniconda
# source activate isoseq_ann
WD=/project/coffea_pangenome/Artocarpus/Comparative_Paper/annotation/egapx
mkdir -p ${WD}/only_longest_transcript_per_gene
cd ${WD}/only_longest_transcript_per_gene

for SAMPLE in HART063 N15_23; do 

	TARGET=/project/coffea_pangenome/Artocarpus/Comparative_Paper/assemblies/softmasked/${SAMPLE}.softmasked.fasta
	TDIR=${WD}/${SAMPLE}_IntraspecificISOseq_NoShortRead
	
	agat_sp_keep_longest_isoform.pl --gff ${TDIR}/complete.genomic.gtf -o ${SAMPLE}.longest_transcript_per_gene.gtf
	gffread ${SAMPLE}.longest_transcript_per_gene.gtf -o ${SAMPLE}.gff3 --keep-genes -O -g ${TARGET} -w ${SAMPLE}.fa
	TransDecoder.LongOrfs -t ${SAMPLE}.fa --output_dir .
	TransDecoder.Predict -t ${SAMPLE}.fa --output_dir . --single_best_only

done 
```

### Liftover

Then run liftoff using that gtf onto the other genomes. Run this first once just on 1 sample to generate the database file. I will run using BOTH references for sensitivity, although we will only use liftover from HART063 for all Artocarpus, since Batocarpus is quite diverged. 

```bash
#!/bin/bash

#SBATCH --time=2-00:00:00   
#SBATCH --nodes=1  
#SBATCH --cpus-per-task=12
#SBATCH --mem=128Gb
#SBATCH --partition=ceres
#SBATCH --account=coffea_pangenome

# module load miniconda
# source activate isoseq_ann
t=12
SAMPLE=$1

TARGET=/project/coffea_pangenome/Artocarpus/Comparative_Paper/assemblies/softmasked/${SAMPLE}.softmasked.fasta
WD=/project/coffea_pangenome/Artocarpus/Comparative_Paper/annotation/egapx

for REF in HART063 N15_23; do 

    cd ${WD} 

    rm -rf ${WD}/liftoff_longestiso/${SAMPLE}
    mkdir -p ${WD}/liftoff_longestiso/${SAMPLE} ${WD}/liftoff_longestiso/ref_${REF}
    cd ${WD}/liftoff_longestiso/${SAMPLE}

    REFERENCE=/project/coffea_pangenome/Artocarpus/Comparative_Paper/assemblies/softmasked/${REF}.softmasked.fasta
    TDIR=${WD}/only_longest_transcript_per_gene/
    
    echo "Working on ${SAMPLE} with reference ${REF}"
    liftoff ${TARGET} ${REFERENCE} -g ${TDIR}/${REF}.longest_transcript_per_gene.gtf -o ${WD}/liftoff_longestiso/${SAMPLE}_ref${REF}.gff3 -p ${t}
    cd ${WD}
    rm -rf ${WD}/liftoff_longestiso/${SAMPLE}

done
```

Once the database is made, run for all samples:

```bash
#!/bin/bash

#SBATCH --time=2-00:00:00   
#SBATCH --nodes=1  
#SBATCH --cpus-per-task=12
#SBATCH --mem=128Gb
#SBATCH --partition=ceres
#SBATCH --account=coffea_pangenome

# module load miniconda
# source activate isoseq_ann
t=12
SAMPLE=$1

TARGET=/project/coffea_pangenome/Artocarpus/Comparative_Paper/assemblies/softmasked/${SAMPLE}.softmasked.fasta
WD=/project/coffea_pangenome/Artocarpus/Comparative_Paper/annotation/egapx

for REF in N15_23 HART063; do 

    cd ${WD} 

    rm -rf ${WD}/liftoff_longestiso/${SAMPLE}
    mkdir -p ${WD}/liftoff_longestiso/${SAMPLE} ${WD}/isoliftoff_longest_transcript_per_gene
    cd ${WD}/liftoff_longestiso/${SAMPLE}

    REFERENCE=/project/coffea_pangenome/Artocarpus/Comparative_Paper/assemblies/softmasked/${REF}.softmasked.fasta
    TDIR=${WD}/only_longest_transcript_per_gene/
    
    echo "Working on ${SAMPLE} with reference ${REF}"
    liftoff ${TARGET} ${REFERENCE} -db ${TDIR}/${REF}.longest_transcript_per_gene.gtf_db -o ${WD}/liftoff_longestiso/${SAMPLE}_ref${REF}.gff3 -p ${t} -exclude_partial
    cd ${WD}
    rm -rf ${WD}/liftoff_longestiso/${SAMPLE}
    
    # Extract protein fasta
    cd ${WD}/isoliftoff_longest_transcript_per_gene
    agat_sp_keep_longest_isoform.pl --gff ${WD}/liftoff_longestiso/${SAMPLE}_ref${REF}.gff3 -o ${SAMPLE}_ref${REF}.longest_transcript_per_gene.gtf
    gffread ${SAMPLE}_ref${REF}.longest_transcript_per_gene.gtf -g ${TARGET} -w ${SAMPLE}_ref${REF}.fa
    TransDecoder.LongOrfs -t ${SAMPLE}_ref${REF}.fa --output_dir .
    TransDecoder.Predict -t ${SAMPLE}_ref${REF}.fa --output_dir . --single_best_only
    rm -rf ${SAMPLE}_ref${REF}.fa.transdecoder_dir

done
```

### Formatting for Downstream

First, subset the files we want, and then create a tab delim file with gene info, including chromosome and function.

**For all downstream: I will use HART063 / N15_23 with their isoseq longest isoform annotations**

* **for outgroups (n=3), I will use the N15_23 liftover**
* **for Artocarpus, I will use the HART063 liftover**

```bash
WD=/project/coffea_pangenome/Artocarpus/Comparative_Paper/annotation/egapx/copies_isoliftoff_longest_transcript_per_gene
mkdir -p proteomes cds gtf
# runs.list has a list of all samples EXCEPT HART063 N15_23, because we use the isoseq peptide and not liftover peptide 
for SAMPLE in $(cat run.list); do 
    REF=HART063    
    cp ${SAMPLE}_ref${REF}.fa.transdecoder.pep proteomes/${SAMPLE}.fa;
    cp ${SAMPLE}_ref${REF}.fa.transdecoder.cds cds/${SAMPLE}.fa;
    cp ${SAMPLE}_ref${REF}.longest_transcript_per_gene.gtf gtf/${SAMPLE}.gtf;
done 

# outgroups 
for SAMPLE in Ficus_carica Morus_mongolica Antiaris_toxicaria; do 
	BASE=$(echo ${SAMPLE} | sed 's/_.*//g');
	echo "Copying ${BASE}"
	cp ${SAMPLE}_refN15_23.fa.transdecoder.pep proteomes/${BASE}.fa
	cp ${SAMPLE}_refN15_23.fa.transdecoder.cds cds/${BASE}.fa
	cp ${SAMPLE}_refN15_23.longest_transcript_per_gene.gtf gtf/${BASE}.gtf
done 

# and main targets
for SAMPLE in HART063 N15_23; do 
	cp ../only_longest_transcript_per_gene/${SAMPLE}.fa.transdecoder.pep proteomes/${SAMPLE}.fa
	cp ../only_longest_transcript_per_gene/${SAMPLE}.fa.transdecoder.cds cds/${SAMPLE}.fa
	cp ../only_longest_transcript_per_gene/${SAMPLE}.longest_transcript_per_gene.gtf gtf/${SAMPLE}.gtf	
done 

# and also for N=5 outgroup analysis, rename simply as artocarpus / batocarpus
cp proteomes/HART063.fa proteomes/Artocarpus.fa
cp gtf/HART063.gtf gtf/Artocarpus.gtf
cp cds/HART063.fa cds/Artocarpus.fa
cp proteomes/N15_23.fa proteomes/Batocarpus.fa
cp gtf/N15_23.gtf gtf/Batocarpus.gtf
cp cds/N15_23.fa cds/Batocarpus.fa

# remove underscores since taxon_geneid gets confusing with double underscore 
for dir in gtf cds proteomes; do
  for f in "$dir"/*; do
    mv "$f" "${f//_/}"
  done
done

```

For the N=3 outgroup samples, modify the chr labels for e.g. CM9823XXX becomes Chr01

```bash
cd gtf
for i in Ficus Morus Antiaris; do
    while read old new; do
    	echo "Fixing ${old} chr name to ${new} for ${i}"
        sed -i "s/$old/$new/g" "${i}.gtf"
    done < ../Outgroup_Chr_Map.tsv
done
```

**Format GTF** 

```bash
DIR=/project/coffea_pangenome/Artocarpus/Comparative_Paper/annotation/egapx/copies_isoliftoff_longest_transcript_per_gene/gtf
cd "$DIR"

for f in *.gtf; do
    sample="${f%.gtf}"
    echo "Working on ${sample}"

    awk -v sample="$sample" -F'\t' '
        $3 == "gene" {
            chr = $1
            start = $4
            end = $5
            strand = $7

            # Combine attributes (columns 9..end)
            attr = $9
            for (i = 10; i <= NF; i++) attr = attr " " $i

            # Defaults
            old_gene_id = "NA"
            desc = "NA"
            biotype = "NA"

            # Prefer gene_id, then copy_num_ID, then ID, then locus_tag
            if (match(attr, /gene_id=([^;]+)/, m)) old_gene_id = m[1]
            else if (match(attr, /copy_num_ID=([^;]+)/, m)) old_gene_id = m[1]
            else if (match(attr, /ID=([^;]+)/, m)) old_gene_id = m[1]
            else if (match(attr, /locus_tag=([^;]+)/, m)) old_gene_id = m[1]

            # Clean quotes if any
            gsub(/"/, "", old_gene_id)

            # Description and biotype (optional)
            if (match(attr, /description=([^;]+)/, m)) desc = m[1]
            if (match(attr, /gene_biotype=([^;]+)/, m)) biotype = m[1]

            # Derive gene_num: drop only leading egapxtmp_ prefix, keep copy suffix (_N)
            gene_num = old_gene_id
            sub(/^egapxtmp_/, "", gene_num)

            # First column: sample_geneNumber (ensure uniqueness with copy suffix)
            print sample "_" gene_num "\t" chr "\t" start "\t" end "\t" strand "\t" sample "\t" old_gene_id "\t" biotype "\t" desc
        }
    ' "$f" > "${sample}.genes.tsv"
done
```

Transcripts only:

```bash
for f in *.gtf; do

    sample="${f%.gtf}"
    echo "Working on ${sample}"

    awk -v sample="$sample" -F'\t' '
        $3 == "transcript" {
            chr   = $1; start = $4; end = $5; strand = $7

            # Merge attributes across fields to handle embedded spaces
            attr = $9; for (i = 10; i <= NF; i++) attr = attr " " $i

            # Defaults per record
            tid = old_gene_id = biotype = desc = "NA"

            # Parse key=value pairs exactly
            n = split(attr, a, ";")
            for (i = 1; i <= n; i++) {
                s = a[i]
                gsub(/^[[:space:]]+|[[:space:]]+$/, "", s)
                split(s, kv, "=")
                key = kv[1]; val = kv[2]
                if (key == "transcript_id")        { tid = val; old_gene_id = val }
                else if (key == "transcript_biotype") biotype = val
                else if (key == "product")  { desc = val; gsub(/%2C/, ",", desc) }
            }

            gene_num = tid; sub(/^egapxtmp_/, "", gene_num)

            print sample "_" gene_num "\t" chr "\t" start "\t" end "\t" strand "\t" sample "\t" old_gene_id "\t" biotype "\t" desc
        }
    ' "$f" > "${sample}.genes.tsv"
done
```

Should look like this:

```bash
head gtf/Artocarpus.genes.tsv 
Artocarpus_unassigned_transcript_1      Chr01   66885   69447   +       Artocarpus      unassigned_transcript_1 mRNA    U-box domain-containing protein 28 pseudogene
Artocarpus_004140-R1    Chr01   84232   88381   -       Artocarpus      egapxtmp_004140-R1      mRNA    root UVB sensitive-like protein (Protein of unknown function, DUF647)
Artocarpus_004211-R1    Chr01   88268   92161   -       Artocarpus      egapxtmp_004211-R1      mRNA    pentatricopeptide repeat-containing protein At2g27800, mitochondrial-like isoform 1-T1
Artocarpus_unassigned_transcript_4      Chr01   88276   91925   -       Artocarpus      unassigned_transcript_4 transcript      pentatricopeptide repeat-containing protein At2g27800, mitochondrial-like, transcript variant X7
Artocarpus_029154-R1    Chr01   132092  135486  +       Artocarpus      egapxtmp_029154-R1      mRNA    uridylate kinase PUMPKIN, chloroplastic-like
Artocarpus_014814-R1    Chr01   138641  140098  -       Artocarpus      egapxtmp_014814-R1      mRNA    protein MKS1-like
```

Then, clean up cds and proteins so they match: 

```bash
WD=/project/coffea_pangenome/Artocarpus/Comparative_Paper/annotation/egapx/copies_isoliftoff_longest_transcript_per_gene
cd ${WD}

for dir in cds proteomes; do 
	echo "Cleaning up ${dir}"
	mkdir -p ${WD}/${dir}/original
	cd ${WD}/${dir}
    for fa in *.fa; do
    	cp ${fa} ${WD}/${dir}/original/
        base=$(basename ${fa} .fa)
        echo "Working on ${base}"
        awk -v b=${base} '
            /^>/ {
                id=$1
                sub(/^>/, "", id)
                sub(/^rna-/, "", id)
                sub(/^egapxtmp_/, "", id)

                # remove everything after first "-" OR "."
                sub(/[.].*/, "", id)

                print ">" b "_" id
                next
            }
            { print }
        ' ${fa} > tmp && mv tmp ${fa}
    done

    # ensure there's no duplicates
    grep '^>' *.fa | sort | uniq -d
done 
```

Sanity:

```
head cds/Artocarpus.fa 
>Artocarpus_000001-R1
ATGATCATGTCTTCAAAAGGGTGTTTAGAGGAGATGGGAATATCTTCAACTAATATCAGT
GATGGTGGGAAAAATTGCTATAGAGGCCATTGGAGACCTGCGGAAGACGAGAAACTCCGA
CAACTCGTCGAACAATTCGGTCCTCAGAACTGGAATTTCATCGCCGAGCATCTACAAGGA

head proteomes/Artocarpus.fa 
>Artocarpus_000001-R1
MIMSSKGCLEEMGISSTNISDGGKNCYRGHWRPAEDEKLRQLVEQFGPQNWNFIAEHLQG
RSGKSCRLRWYNQLDPNINKKPFTEEEEERLLSAHRIYGNKWACIAKYFQGRTDNAVKNH
```

Sanity check to ensure the cds headers match the genes.tsv:

```
for fa in *.fa; do
  sample=${fa%.fa}
  tsv="../gtf/${sample}.genes.tsv"
  echo "Checking ${fa} vs ${tsv}"

  awk 'BEGIN{FS="\t"}
    # Read TSV col1 first
    NR==FNR { cnt[$1]++; tsv_total++; next }

    # Read FASTA headers
    /^>/ {
      h = substr($0, 2)
      sub(/[[:space:]]+$/, "", h)
      fasta_total++
      if (!(h in cnt)) { missing[h]=1; nmiss++ }
      else if (cnt[h] != 1) { bad[h]=cnt[h]; nbad++ }
      seen[h]=1
    }

    END {
      # Summary to stderr so the loop can still see failure via exit code
      print "FASTA headers:", fasta_total, " | TSV entries:", tsv_total > "/dev/stderr"

      for (h in missing) print "Missing in TSV:", h > "/dev/stderr"
      for (h in bad)     print "TSV count != 1:", h, "found", bad[h] > "/dev/stderr"
 

      if (nmiss + nbad > 0) exit 1
    }' "$tsv" "$fa" || echo "❌ Problem detected for ${sample}"
done

```

NOW, we need to use that map to divide the CDS into chromosome-specific files:

```bash
WD=/project/coffea_pangenome/Artocarpus/Comparative_Paper/annotation/egapx/copies_isoliftoff_longest_transcript_per_gene/cds
cd ${WD}

mkdir -p chrs
for fa in *.fa; do
    sample=${fa%.fa}
    echo "Splitting ${sample}"
    awk -v s="$sample" '
    NR==FNR { map[$1]=$2; next }
    /^>/ {
        id = substr($1, 2)
        chr = map[id]
        file = "chrs/" s "_" chr ".fa"
    }
    { if (chr != "") print > file }
    ' "${sample}.genes.tsv" "$fa"
done

```

Do the same for proteins:

```bash
WD=/project/coffea_pangenome/Artocarpus/Comparative_Paper/annotation/egapx/copies_isoliftoff_longest_transcript_per_gene/proteomes
cd ${WD}

mkdir -p chrs
for fa in *.fa; do
    sample=${fa%.fa}
    echo "Splitting ${sample}"
    awk -v s="$sample" '
    NR==FNR { map[$1]=$2; next }
    /^>/ {
        id = substr($1, 2)
        chr = map[id]
        file = "chrs/" s "_" chr ".fa"
    }
    { if (chr != "") print > file }
    ' "${sample}.genes.tsv" "$fa"
done
```

### Final Counts

Take the `.gtf` for each species and extract counts:

```R
setwd('/project/coffea_pangenome/Artocarpus/Comparative_Paper/annotation/egapx/copies_isoliftoff_longest_transcript_per_gene/gtf')
library(tidyverse)
library(data.table)
library(stringr)

read_gtf_like <- function(path) {
  dt <- data.table::fread(
    path, sep = "\t", skip = 3, header = FALSE, quote = "", fill = TRUE,
    data.table = TRUE, showProgress = FALSE
  )
  dt <- dt[!grepl("^#", V1)]
  if (ncol(dt) < 9) stop("File seems malformed (<9 columns): ", path)
  
  setnames(dt, c("seqid","source","type","start","end","score","strand","phase","attributes"))
  dt[, start := as.integer(start)]
  dt[, end   := as.integer(end)]
  dt
}

get_attr <- function(x, key) {
  str_match(x, paste0("(^|;)", key, "=([^;]+)"))[,3]
}

split_on_comma <- function(df, col) {
  df %>%
    mutate("{col}" := str_split(.data[[col]], ",")) %>%
    unnest(.data[[col]])
}

summarize_one_file <- function(path) {
  g <- read_gtf_like(path) %>%
    as_tibble() %>%
    mutate(
      ID     = get_attr(attributes, "ID"),
      Parent = get_attr(attributes, "Parent"),
      gene_biotype = get_attr(attributes, "gene_biotype"),
      transcript_biotype = get_attr(attributes, "transcript_biotype"),
      pseudo = get_attr(attributes, "pseudo")
    )
  
  genes <- g %>%
    filter(type == "gene", !is.na(ID)) %>%
    transmute(
      gene_id = ID,
      gene_biotype = gene_biotype,
      pseudo = pseudo,
      # normalize / coarsen biotypes for reporting
      gene_biotype2 = case_when(
        is.na(gene_biotype) | gene_biotype == "" ~ "unknown",
        gene_biotype %in% c("pseudogene", "transcribed_pseudogene") ~ "pseudogene",
        gene_biotype == "protein_coding" ~ "protein_coding",
        gene_biotype == "lncRNA" ~ "lncRNA",
        TRUE ~ gene_biotype
      )
    )
  
  tx <- g %>%
    filter(type %in% c("transcript","mRNA"), !is.na(ID)) %>%
    transmute(
      tx_id = ID,
      gene_id = Parent,
      tx_start = start,
      tx_end = end
    ) %>%
    split_on_comma("gene_id")
  
  exons <- g %>%
    filter(type == "exon", !is.na(Parent)) %>%
    transmute(
      tx_id = Parent,
      exon_len = (end - start + 1L)
    ) %>%
    split_on_comma("tx_id")
  
  cds <- g %>%
    filter(type == "CDS", !is.na(Parent)) %>%
    transmute(
      tx_id = Parent,
      cds_len = (end - start + 1L)
    ) %>%
    split_on_comma("tx_id")
  
  ex_by_tx <- exons %>%
    group_by(tx_id) %>%
    summarise(
      exon_count = n(),
      exon_bp = sum(exon_len),
      .groups = "drop"
    )
  
  cds_by_tx <- cds %>%
    group_by(tx_id) %>%
    summarise(
      cds_bp = sum(cds_len),
      .groups = "drop"
    )
  
  tx2 <- tx %>%
    left_join(ex_by_tx, by = "tx_id") %>%
    left_join(cds_by_tx, by = "tx_id") %>%
    mutate(
      exon_count = replace_na(exon_count, 0L),
      exon_bp    = replace_na(exon_bp, 0L),
      cds_bp     = replace_na(cds_bp, 0L),
      tx_len_bp = if_else(exon_bp > 0L, exon_bp, (tx_end - tx_start + 1L)),
      is_multiexon = exon_count >= 2L
    )
  
  # representative transcript per gene (longest)
  gene_rep <- tx2 %>%
    filter(!is.na(gene_id)) %>%
    group_by(gene_id) %>%
    slice_max(order_by = tx_len_bp, n = 1, with_ties = FALSE) %>%
    ungroup()
  
  gene_rep2 <- gene_rep %>%
    left_join(genes, by = "gene_id") %>%
    mutate(
      is_protein_coding = gene_biotype2 == "protein_coding",
      is_lncRNA         = gene_biotype2 == "lncRNA",
      is_pseudogene     = gene_biotype2 == "pseudogene"
    )
  
  tibble(
    Accession = str_replace(basename(path), "\\.gtf$", ""),
    file = basename(path),
    n_genes = n_distinct(genes$gene_id),
    n_transcripts = n_distinct(tx2$tx_id),
    
    # gene biotype counts (using the normalized gene_biotype2)
    n_protein_coding_genes = sum(genes$gene_biotype2 == "protein_coding", na.rm = TRUE),
    n_pseudogenes          = sum(genes$gene_biotype2 == "pseudogene", na.rm = TRUE),  # includes transcribed_pseudogene
    n_lncRNA_genes         = sum(genes$gene_biotype2 == "lncRNA", na.rm = TRUE),
    n_other_or_unknown_genes = sum(!genes$gene_biotype2 %in% c("protein_coding","pseudogene","lncRNA"), na.rm = TRUE),
    
    # structural stats (using representative transcript per gene)
    mean_tx_len_bp = mean(gene_rep2$tx_len_bp, na.rm = TRUE),
    mean_cds_len_bp = mean(gene_rep2$cds_bp[gene_rep2$is_protein_coding], na.rm = TRUE),
    mean_exons_per_gene = mean(gene_rep2$exon_count, na.rm = TRUE),
    pct_multiexon_genes = 100 * mean(gene_rep2$is_multiexon, na.rm = TRUE)
  )
}

gtfs <- list.files(pattern = "\\.gtf$", full.names = TRUE)

metrics <- purrr::map_dfr(gtfs, summarize_one_file) %>%
  mutate(
    Method = case_when(
      Accession %in% c("N1523","HART063") ~ "Iso-Seq annotation",
      TRUE ~ "Liftover"
    )
  )

metrics

#removes dups
mets <- metrics %>% filter(!grepl('Arto|Bato',Accession))
readr::write_tsv(metrics, "~/artocarpus_comparative_genomics/03_annotation/gene_annotation/AnnotationQC_fromGTF.tsv")

# What leads to the discrepancy between n_genes != pcgs + psuedo?
g <- data.table::fread("HART063.gtf", sep="\t", skip=3, header=FALSE, fill=TRUE, data.table=FALSE)
colnames(g)[1:9] <- c("seqid","source","type","start","end","score","strand","phase","attr")
genes <- g[g$type=="gene", ]
biotype <- stringr::str_match(genes$attr, "(^|;)gene_biotype=([^;]+)")[,3]
sort(table(biotype), decreasing=TRUE)[1:10]
# biotype
# protein_coding             pseudogene                 lncRNA transcribed_pseudogene                   <NA>                   <NA>                   <NA> 
#   27475                   5956                    703                      5                                                                      
# <NA>                   <NA>                   <NA> 

mets <- read_tsv("~/artocarpus_comparative_genomics/03_annotation/gene_annotation/AnnotationQC_fromGTF.tsv")
md <- read_tsv('~/artocarpus_comparative_genomics/samples.txt') %>% mutate(Accession = gsub('_','',Accession))
df <- md %>%
  left_join(mets, by = "Accession") %>%
  mutate(
    Species_short = gsub('Artocarpus','A.',Group),
    ylab = paste0(Species_short, " (", Accession, ")"),
    ylab = fct_reorder(ylab, Accession_Order,.desc = TRUE)
  )

plot_df <- df %>%
  transmute(
    ylab, Method,
    protein_coding = n_protein_coding_genes,
    pseudogene = n_pseudogenes,
    lncRNA = n_lncRNA_genes,
    other = n_other_or_unknown_genes
  ) %>%
  pivot_longer(cols = c(protein_coding, pseudogene, lncRNA, other),
               names_to = "biotype", values_to = "n") %>%
  mutate(
    biotype = factor(biotype, levels = c("protein_coding", "pseudogene", "lncRNA", "other"))
  )
lab_df <- plot_df %>%
  filter(biotype == "protein_coding") %>%
  group_by(ylab,Method) %>%
  summarise(n = sum(n), .groups = "drop") %>%
  mutate(ypos = n / 2)

p <- ggplot(plot_df, aes(x = ylab, y = n, fill = biotype)) +
  geom_col(width = 0.8) +
  geom_text(
    data = lab_df,
    aes(x = ylab, y = ypos, label = scales::comma(n)),
    inherit.aes = FALSE,
    color = "white",
    size = 2.8
  ) +
  coord_flip() +
  scale_fill_manual(values = c(
    protein_coding = "#377eb8",
    pseudogene = "#e41a1c",
    lncRNA = "#4daf4a",
    other = "grey70"
  )) +
  facet_grid(Method ~ ., scales = "free_y", space = "free_y") +
  theme_bw(base_size = 9) +
  theme(
    legend.position = "right",
    axis.title = element_blank()
  ) +
  labs(fill = NULL, y = "Genes")

p
ggsave('~/artocarpus_comparative_genomics/figures/20260413_Annotation_Counts.pdf',p,height=5,width=4.5)

mets %>% select(matches('Acc|^n_ge|n_transcripts|protein'))
# # A tibble: 16 × 4
# Accession  n_genes n_transcripts n_protein_coding_genes
# <chr>        <dbl>         <dbl>                  <dbl>
#   1 Antiaris     12693         12736                  11337
# 2 Artocarpus   34139         35714                  27475
# 3 Batocarpus   21491         22508                  18374
# 4 Ficus        13259         13303                  11886
# 5 HART001      33194         33590                  26816
# 6 HART027      31338         31679                  25707
# 7 HART058      31600         31955                  25898
# 8 HART060      29803         30109                  24609
# 9 HART061      31240         31600                  25590
# 10 HART062      31398         31760                  25705
# 11 HART063      34139         35714                  27475
# 12 HART067      32903         33300                  26610
# 13 HART068      29716         30031                  24543
# 14 Morus        15116         15175                  13625
# 15 N1523        21491         22508                  18374
# 16 N9750        29577         29890                  24477

```

## Repeat Annotation

Annotate each genome with earlgrey v6.3.0

```bash
#!/bin/bash

#SBATCH --time=14-00:00:00   
#SBATCH --nodes=1  
#SBATCH --cpus-per-task=20
#SBATCH --mem=112Gb
#SBATCH --partition=ceres
#SBATCH --account=coffea_pangenome

t=20

#module load miniconda
#source activate earlgrey
#source activate eg6

SAMPLE=$1
WD=/project/coffea_pangenome/Artocarpus/Assemblies/20250101_JustinAssemblies/primary_asm

mkdir -p ${WD}/EarlGrey_SampleLibrary
cd ${WD}/EarlGrey_SampleLibrary

FASTA="${WD}/${SAMPLE}.fa" 

if [ ! -s ${WD}/EarlGrey_SampleLibrary/${SAMPLE}_EarlGrey/${SAMPLE}_summaryFiles/${SAMPLE}.softmasked.fasta ] && [ -s ${FASTA} ]; then
    echo -e "\e[43m~~~~ Starting repeat annotation for ${SAMPLE} ~~~~\e[0m"
    # Run earlgrey with eudicotyledons repeatmasker search time, generating soft-masked genome and run helitrons. 
    earlGrey -r eudicotyledons -d yes -e yes -t ${t} -g ${FASTA} -s ${SAMPLE} -o ${WD}/EarlGrey_SampleLibrary

else
    echo -e "\e[42m~~~~ Skipping repeat annotation for ${SAMPLE}, already exists ~~~~\e[0m"
fi 
```

And copy:

```bash
DIR=/90daydata/coffea_pangenome/scratch/repeats
REP=/project/coffea_pangenome/Artocarpus/Comparative_Paper/repeats/inhouse_genomes_only_compiled/

for SAMPLE in $(cat CompSamples.list) ; do 
cp ${DIR}/${SAMPLE}_EarlGrey/${SAMPLE}_summaryFiles/* ${REP}/
done 
```

Add accession to each output:

```bash
for SAMPLE in $(ls *.familyLevelCount.txt | sed 's/.familyLevelCount.txt//g'); do 
echo "${SAMPLE}"
awk -v s=${SAMPLE} '{OFS="\t"}{print $0, s}' ${SAMPLE}.familyLevelCount.txt > ${SAMPLE}.families.out
awk -v s=${SAMPLE} '{OFS="\t"}{print $0, s}' ${SAMPLE}.highLevelCount.txt > ${SAMPLE}.summary.out
awk -v s=${SAMPLE} '{OFS="\t"}{print $0, s}' ${SAMPLE}_divergence_summary_table.tsv > ${SAMPLE}.divergence.out
done 

mergem *families.out > Repeat_Families.txt
mergem *summary.out > Repeat_Summaries.txt
mergem *divergence.out > Divergence_Summaries.txt
```



### Repeat variation & plots

- Summarize high‑level TE composition per assembly (stacked proportions from EarlGrey / RepeatMasker).
- PGLS: `log(GSize) ~ log(Repeats)` (phylogeny‑corrected using species tree; λ estimated by ML).
- PGLS: `young_frac_LTR ~ log(GSize)` (young_frac = bp‑weighted fraction of LTR bp below chosen Kimura cutoff).

```R
setwd('/project/coffea_pangenome/Artocarpus/Comparative_Paper/repeats/')
library(tidyverse)
library(RColorBrewer)
library(ggpubr)
library(meRo) #devtools::install_github('merondun/meRo')
library(vegan)
library(broom)
library(ggrepel)
library(caper)

# Add metadata information
md <- read_tsv('~/artocarpus_comparative_genomics/samples.txt')

##### High Level #####
high_level <- read_tsv('~/artocarpus_comparative_genomics/03_annotation/repeat_annotation/Repeat_Summaries.txt') 
names(high_level) <- c('Classification','Coverage','Count','Proportion','Gen','Distinct_Classifications','Accession')
non_repeat <- high_level %>% 
  group_by(Accession) %>% 
  summarize(Proportion = 1-(sum(Proportion)),
            Classification = "Non Repeat",
            Coverage=NA,Count=NA,Gen=NA,Distinct_Classifications=NA)
full_level <- rbind(high_level,non_repeat) %>% filter(!grepl('Anti|Ficus|Morus',Accession))
fl <- left_join(full_level,md)
ord <- md %>% dplyr::select(Accession,Group,Accession_Order) %>% arrange(Accession_Order)
grpord <- md %>% dplyr::select(Group,Accession_Order) %>% arrange(Accession_Order) %>% dplyr::select(Group) %>% distinct
fl$Accession <- factor(fl$Accession,levels=rev(ord$Accession))
fl$Group <- factor(fl$Group,levels=grpord$Group)
fl$Classification <- factor(fl$Classification,levels=c('Non Repeat','Unclassified','Other (Simple Repeat, Microsatellite, RNA)','DNA','Penelope','Rolling Circle','LTR','LINE','SINE'))
cols <- fl %>% dplyr::select(Classification) %>% distinct %>% mutate(Color = brewer.pal(9,'Paired'))
gcols <- md %>% dplyr::select(Group,Color,Shape) %>% distinct

# ltr labs
fl_labels <- fl %>%
  filter(Classification == "LTR") %>%
  mutate(
    label = paste0(round(Proportion * 100, 1), "%"),
    text_color = ifelse(Proportion > 0.08, "black", "black")
  )


# Plot landscape 
all <- fl %>% 
  mutate(Coverage = Coverage / 1e6) %>% 
  pivot_longer(c(Proportion)) %>%
  filter( !(name == 'Distinct_Classifications' & (Classification == 'Unclassified' | Classification == "Other (Simple Repeat, Microsatellite, RNA)")) ) %>% 
  ggplot(aes(y = Accession, x = value, fill = Classification)) +
  geom_bar(stat = 'identity', position = position_stack()) +
  # add LTR percent labels
  geom_text(
    data = fl_labels,
    aes(y = Accession, x = Proportion, label = label),
    position = position_stack(vjust = 0.5),
    color = fl_labels$text_color,
    size = 2.5
  ) +
  theme_bw() +
  facet_grid(Group ~ name, scales = 'free', space = 'free_y') +
  scale_fill_manual(values = cols$Color, breaks = cols$Classification) +
  theme(strip.text.y = element_text(angle = 0)) +
  ylab('') + xlab('Distinct Classifications') +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 3))

all
ggsave('~/symlinks/comp/figures/20260410_RepeatsHighLevelSummary.pdf',
       all,dpi=300,height=4,width=6.5)

# fl %>% dplyr::select(Accession,Group,Classification,Proportion) %>% pivot_wider(names_from = Classification,values_from=Proportion)
# # A tibble: 11 × 11
# Accession Group               DNA    LINE   LTR `Other (Simple Repeat, Microsatellite, RNA)` `Rolling Circle`      SINE Unclassified Penelope `Non Repeat`
# <fct>     <fct>             <dbl>   <dbl> <dbl>                                        <dbl>            <dbl>     <dbl>        <dbl>    <dbl>        <dbl>
#   1 HART001   A. altilis       0.0385 0.00807 0.414                                       0.0199          0.00531   3.81e-4        0.200 NA              0.314
# 2 HART027   A. heterophyllus 0.0230 0.00449 0.466                                       0.0190          0.00320  NA              0.204  1.24e-4        0.280
# 3 HART058   A. rigidus       0.0342 0.00776 0.388                                       0.0198          0.00349   2.53e-4        0.232 NA              0.315
# 4 HART060   A. dadah         0.0386 0.00721 0.336                                       0.0322          0.00192  NA              0.257  1.71e-3        0.325
# 5 HART061   A. odoratissimus 0.0378 0.00347 0.394                                       0.0184          0.00681   2.42e-5        0.226 NA              0.313
# 6 HART062   A. odoratissimus 0.0381 0.00374 0.403                                       0.0175          0.00462  NA              0.221 NA              0.312
# 7 HART063   A. camansi       0.0467 0.00601 0.382                                       0.0194          0.00299  NA              0.208 NA              0.335
# 8 HART067   A. mariannensis  0.0418 0.00676 0.378                                       0.0222          0.00574   3.21e-4        0.209 NA              0.336
# 9 HART068   A. nitidus       0.0375 0.00495 0.364                                       0.0203          0.00438   4.30e-6        0.222  3.55e-4        0.346
# 10 N15_23    Batocarpus sp.   0.0593 0.00617 0.360                                       0.0379          0.00573   1.62e-4        0.214 NA              0.317
# 11 N97_50    A. lacucha       0.0437 0.00337 0.346                                       0.0186          0.00117   6.78e-4        0.230  6.51e-4        0.356
write.csv(fl %>% dplyr::select(Accession,Accession_Order,Group,Classification,Coverage,Count,Proportion,Distinct_Classifications),'~/artocarpus_comparative_genomics/03_annotation/repeat_annotation/Repeat_proportions_coverage_summarized.csv',quote = F,row.names = F)

#### 	PGLS ####
# Significance genome size ~ LTRs 
rep_df <- fl %>% filter(!grepl('Non Repeat',Classification)) %>% 
  dplyr::rename(GSize = Genome_Size_Assembly_Mb,
                phylo_order = Accession_Order) %>% 
  group_by(Accession,Group,phylo_order,GSize) %>% 
  summarize(Repeats = sum(Coverage)/1e6)
rep_df

# Spearman correlation
cor_res <- cor.test(
  ~ Repeats + GSize,
  data = rep_df,
  method = "spearman"
) %>% tidy()
cor_res
# # A tibble: 1 × 5
# estimate statistic p.value method                          alternative
# <dbl>     <dbl>   <dbl> <chr>                           <chr>      
#   1    0.991      2.00       0 Spearman's rank correlation rho two.sided  

# pgls
pgls_in <- rep_df %>% ungroup %>% dplyr::select(Accession,Repeats,GSize) %>% mutate( Accession = gsub("_","", Accession) ) %>% as.data.frame
nwk <- read.tree('~/artocarpus_comparative_genomics/05_orthofinder/SpeciesTree_rooted_node_labels.txt')
nwk$node.label <- NULL
comp <- comparative.data(phy = nwk, data = pgls_in, names.col = Accession, vcv = TRUE, na.omit = FALSE )

pgls_model <- pgls( log(GSize) ~ log(Repeats) , data = comp, lambda = "ML" ) 
summary(pgls_model)
# Call:
#   pgls(formula = log(GSize) ~ log(Repeats), data = comp, lambda = "ML")
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.18570 -0.06371  0.07070  0.14834  0.25227 
# 
# Branch length transformations:
#   
#   kappa  [Fix]  : 1.000
# lambda [ ML]  : 0.951
# lower bound : 0.000, p = 0.054709
# upper bound : 1.000, p = 0.59073
# 95.0% CI   : (NA, NA)
# delta  [Fix]  : 1.000
# 
# Coefficients:
#   Estimate Std. Error t value  Pr(>|t|)    
# (Intercept)  1.271753   0.336977   3.774  0.004389 ** 
#   log(Repeats) 0.854150   0.055247  15.461 8.672e-08 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.1539 on 9 degrees of freedom
# Multiple R-squared: 0.9637,	Adjusted R-squared: 0.9597 
# F-statistic:   239 on 1 and 9 DF,  p-value: 8.672e-08


pgls_line <- data.frame(
  Repeats = pgls_in$Repeats,
  fitted = fitted(pgls_model)
)

profile <- pgls.profile(pgls_model) 
plot(profile)

# add fitted values aligned to rows used by the model
pdat <- comp$data %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "Accession") %>%
  mutate(
    logG = log(GSize),
    logR = log(Repeats),
    fitted_logG = as.numeric(fitted(pgls_model))
  ) %>% 
  left_join(., 
            md %>% dplyr::select(Accession,Group,Color,Shape) %>% mutate(Accession = gsub('_','',Accession)))

# for a line, sort by x
pdat_line <- pdat %>% arrange(logR)

s <- summary(pgls_model)
beta <- s$coefficients["log(Repeats)", "Estimate"]
se   <- s$coefficients["log(Repeats)", "Std. Error"]
pval <- s$coefficients["log(Repeats)", "Pr(>|t|)"]
lam  <- pgls_model$param["lambda"]

label_pgls <- sprintf("PGLS: b = %.3f ± %.3f\np = %.2e\nl = %.2f", beta, se, pval, lam)

p_pgls <- ggplot(pdat, aes(x = logR, y = logG, fill = Group, shape = Group, label =Group)) +
  geom_point(size = 2) +
  geom_text_repel(size = 2, max.overlaps = Inf) +
  geom_line(data = pdat_line, aes(x = logR, y = fitted_logG), inherit.aes = FALSE,
            color = "blue", linewidth = 0.8) +
  theme_bw(base_size = 8) +
  scale_shape_manual(values=md$Shape,breaks=md$Group)+
  scale_fill_manual(values=md$Color,breaks=md$Group)+
  annotate("text", x = min(pdat$logR)+0.01, y = max(pdat$logG)-0.05, label = label_pgls,
           hjust = 0, vjust = 0, size = 1.2) +
  labs(
    x = "log(Repeat masked bp)",
    y = "log(Genome size (bp))"
  ) +
  theme(legend.position = "none")

p_pgls

ggsave('~/symlinks/comp/figures/20260410_RepeatsHighLevelSummary-Repeats-GenomeSize_PGLS.pdf',
       p_pgls,dpi=300,height=3,width=2.75)

###### Divergence Summaries ######
t <- read_tsv('~/artocarpus_comparative_genomics/03_annotation/repeat_annotation/Divergence_Summaries.txt')
t <- t %>% dplyr::rename(Accession = HART001)
tm <- left_join(t,md)
tm$Accession <- factor(tm$Accession,levels=ord$Accession)
tm$Group <- factor(tm$Group,levels=rev(grpord$Group))

dom <- c("HART001","HART027")

# helper build div_traits + run pgls for a given cutoff since we run few models sensitivity
run_div_pgls <- function(young_cut_pct, tm, rep_df, nwk) {
  young_cut <- young_cut_pct / 100
  
  div_traits <- tm %>%
    filter(subclass %in% c("LTR","LINE","DNA")) %>%
    mutate(domest = if_else(as.character(Accession) %in% dom, "domesticated", "wild")) %>%
    group_by(Accession, domest, subclass) %>%
    summarise(
      te_mb = sum(total_bp, na.rm = TRUE) / 1e6,
      young_mb = sum(total_bp[mean_div <= young_cut], na.rm = TRUE) / 1e6,
      young_frac = young_mb / te_mb,
      wmean_div = weighted.mean(mean_div, w = total_bp, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    pivot_wider(
      names_from = subclass,
      values_from = c(te_mb, young_mb, young_frac, wmean_div),
      names_sep = "_"
    ) %>%
    mutate(Accession = gsub("_","", Accession)) %>%
    left_join(
      rep_df %>%
        dplyr::select(Accession, GSize, Repeats) %>%
        mutate(Accession = gsub("_","", Accession)),
      by = "Accession"
    ) %>%
    as.data.frame()
  
  comp_div <- comparative.data(
    phy = nwk,
    data = div_traits,
    names.col = "Accession",
    vcv = TRUE,
    na.omit = FALSE
  )
  
  # ensure domest is a factor
  comp_div$data$domest <- factor(comp_div$data$domest, levels = c("wild","domesticated"))
  
  m_dom <- pgls(wmean_div_LTR ~ domest, data = comp_div, lambda = "ML")
  m_div <- pgls(young_frac_LTR ~ log(GSize), data = comp_div, lambda = "ML")
  
  list(div_traits = div_traits, comp_div = comp_div, m_dom = m_dom, m_div = m_div)
}

res10 <- run_div_pgls(young_cut_pct = 10, tm = tm, rep_df = rep_df, nwk = nwk)
m_div <- res10$m_div

summary(res10$m_div)

# Call:
#   pgls(formula = young_frac_LTR ~ log(GSize), data = comp_div, 
#        lambda = "ML")
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.94304 -0.54048 -0.04824  0.37944  0.83333 
# 
# Branch length transformations:
#   
#   kappa  [Fix]  : 1.000
# lambda [ ML]  : 1.000
# lower bound : 0.000, p = 0.019737
# upper bound : 1.000, p = 1    
# 95.0% CI   : (0.393, NA)
# delta  [Fix]  : 1.000
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)  
# (Intercept) -4.11747    1.64560 -2.5021  0.03375 *
#   log(GSize)   0.66154    0.25404  2.6041  0.02855 *
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.6447 on 9 degrees of freedom
# Multiple R-squared: 0.4297,	Adjusted R-squared: 0.3663 
# F-statistic: 6.781 on 1 and 9 DF,  p-value: 0.02855 

pdat_young <- res10$comp_div$data %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "Accession") %>%
  mutate(
    logG = log(GSize),
    young = young_frac_LTR,
    fitted_young = as.numeric(fitted(m_div))
  )

pdat_young_line <- pdat_young %>% arrange(logG)

# label (β, p, λ)
s <- summary(m_div)
beta <- s$coefficients["log(GSize)", "Estimate"]
se   <- s$coefficients["log(GSize)", "Std. Error"]
pval <- s$coefficients["log(GSize)", "Pr(>|t|)"]
lam  <- m_div$param["lambda"]

label_pgls_young <- sprintf(
  "PGLS: b = %.3f ± %.3f\np = %.2e\nl = %.2f\nyoung cutoff = %d%%",
  beta, se, pval, lam, 10
)

# Plot 
p_young <- ggplot(
  pdat_young,
  aes(x = logG, y = young, fill = Group, shape = Group, label = Accession)
) +
  geom_point(size = 2) +
  geom_text_repel(size = 2, max.overlaps = Inf) +
  geom_line(
    data = pdat_young_line,
    aes(x = logG, y = fitted_young),
    inherit.aes = FALSE,
    color = "blue",
    linewidth = 0.8
  ) +
  theme_bw(base_size = 8) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 2)) +  
  scale_shape_manual(values=md$Shape,breaks=md$Group)+
  scale_fill_manual(values=md$Color,breaks=md$Group)+
  annotate(
    "text",
    x = min(pdat_young$logG, na.rm = TRUE) + 0.01,
    y = max(pdat_young$young, na.rm = TRUE) - 0.02,
    label = label_pgls_young,
    hjust = 0, vjust = 1, size = 1.2
  ) +
  labs(
    x = "log(Genome size (bp))",
    y = "LTR young fraction (bp-weighted)"
  ) +
  theme(legend.position = "none")

p_young

ggsave('~/symlinks/comp/figures/20260410_LTRYoungFrac_vs_GSize_PGLS.pdf',
       p_young, dpi=300, height=3, width=2.75)

# domesticates?
summary(res10$m_dom)

# Call:
#   pgls(formula = wmean_div_LTR ~ domest, data = comp_div, lambda = "ML")
# 
# Residuals:
#   Min        1Q    Median        3Q       Max 
# -0.109326 -0.050600  0.007363  0.046943  0.071145 
# 
# Branch length transformations:
#   
#   kappa  [Fix]  : 1.000
# lambda [ ML]  : 1.000
# lower bound : 0.000, p = 0.0017301
# upper bound : 1.000, p = 1    
# 95.0% CI   : (0.776, NA)
# delta  [Fix]  : 1.000
# 
# Coefficients:
#   Estimate Std. Error t value  Pr(>|t|)    
# (Intercept)         0.1559046  0.0086090 18.1095 2.176e-08 ***
#   domestdomesticated -0.0128897  0.0051162 -2.5194    0.0328 *  
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.06965 on 9 degrees of freedom
# Multiple R-squared: 0.4136,	Adjusted R-squared: 0.3484 
# F-statistic: 6.347 on 1 and 9 DF,  p-value: 0.0328

## Sensitivity loop over young_cut_pct values (5, 10, 15)
sens_vals <- c(5, 10, 15)
sens_out <- list()

for (sens in sens_vals) {
  
  res <- run_div_pgls(young_cut_pct = sens, tm = tm, rep_df = rep_df, nwk = nwk)
  
  # extract stats: genome size effect on young_frac_LTR
  sm_div <- summary(res$m_div)
  gs_beta <- sm_div$coefficients["log(GSize)", "Estimate"]
  gs_se   <- sm_div$coefficients["log(GSize)", "Std. Error"]
  gs_p    <- sm_div$coefficients["log(GSize)", "Pr(>|t|)"]
  gs_lam  <- res$m_div$param["lambda"]
  
  sens_out[[paste0(as.character(sens), "_2")]] <- tibble(
    young_cut_pct = sens,
    
    model = "young_frac_LTR ~ log(GSize)",
    beta = gs_beta,
    se = gs_se,
    p = gs_p,
    lambda = gs_lam
  )
}

sens_tbl <- bind_rows(sens_out) %>%
  arrange(model, young_cut_pct)

sens_tbl
# # A tibble: 3 × 6
# young_cut_pct model                        beta    se      p lambda
# <dbl> <chr>                       <dbl> <dbl>  <dbl>  <dbl>
#   1             5 young_frac_LTR ~ log(GSize) 0.374 0.132 0.0196      1
# 2            10 young_frac_LTR ~ log(GSize) 0.662 0.254 0.0285      1
# 3            15 young_frac_LTR ~ log(GSize) 0.388 0.224 0.117       1


```

