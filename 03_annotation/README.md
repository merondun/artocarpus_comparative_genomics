# Iso-seq Gene Annotation

We have Isoseq data for 2 samples (Artocarpus camansi (6 tissues) and Batocarpus sp. (2 tissues)).

These stpes will demultiplex Iso‑Seq reads, convert selected partitions to FASTA, and run eGAPx (with optional short reads just for sensitivity, it was found that Isoseq is sufficient) to generate gene/transcript models. Extract the longest isoform per gene and predict CDS/proteins with TransDecoder, then liftover reference annotations to other assemblies. Finally clean/standardize headers, produce mapping tables, and split CDS/proteomes by chromosome for downstream comparative analyses. 

In the end, can produce these summary figures:

![annotation_counts](/figures/20260413_Annotation_Counts_Orthogroups.png)



Primary outputs after this massive codeblock: 

- eGAPx annotations: complete.genomic.gtf and run out/work directories.
- Predicted coding sequences/proteomes: .transdecoder.cds and .transdecoder.pep.
- Longest‑isoform exports: *.longest_transcript_per_gene.gtf and .fa.
- Liftoff transfers: per‑sample GFF3 liftover files.
- Cleaned inputs for downstream: proteomes/, cds/, gtf/ (clean headers) + genes.tsv mappings and chromosome‑split FASTA files.

___

## Demultiplexing & Input Prep

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

## eGAPx Annotation

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

Sanity checks: count genes/transcripts/cds/proteins from the files:

```bash
awk 'BEGIN{genes=0; transcripts=0; cds=0} !/^#/{
  if($3=="gene") genes++;
  if($3=="transcript" || $3=="mRNA") transcripts++;
  if($3=="CDS") cds++;
  if(match($0,/protein_id "[^"]+"/)){
    pid[substr($0,RSTART+12,RLENGTH-13)]=1;
  }
} END{
  print "Genes:",genes,"Transcripts:",transcripts,"CDS:",cds,"Unique_Proteins:",length(pid)
}' complete.genomic.gtf
Genes: 34139 Transcripts: 55865 CDS: 304032 Unique_Proteins: 45455

```

## Liftover

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

## Formatting for Downstream

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

## Final Counts

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
readr::write_tsv(metrics, "AnnotationQC_fromGTF.tsv")

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

md <- read_tsv('~/symlinks/comp/samples.txt') %>% mutate(Accession = gsub('_','',Accession))
df <- md %>%
  left_join(mets, by = "Accession") %>%
  mutate(
    Species_short = gsub('Artocarpus','A.',Group),
    ylab = paste0(Species_short, " (", Accession, ")"),
    ylab = fct_reorder(ylab, `Accession Order`,.desc = TRUE)
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
ggsave('~/symlinks/comp/figures/20260413_Annotation_Counts.pdf',p,height=5,width=4.5)

```

And after running orthofinder later... return here to show orthogroup overlap by sample: 

```R
setwd('/project/coffea_pangenome/Artocarpus/Comparative_Paper/orthofinder/Subgenome_Divided/OrthoFinder/Results_Feb24/Orthogroups')
library(tidyverse)
library(stringr)
library(readr)

gc <- read_tsv("Orthogroups.GeneCount.tsv", show_col_types = FALSE)
md <- read_tsv('~/symlinks/comp/samples.txt') %>% mutate(Accession = gsub('_','',Accession)) %>% dplyr::select(Accession,Group,ord = `Accession Order`)
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

