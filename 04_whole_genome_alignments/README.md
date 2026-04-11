# Whole Genome Alignments

Per‑pair whole‑genome alignments are generated with nucmer, filtered by length/identity and merged into syntenic blocks (within a MAXGAP), then converted to BED/link files (.simple / .cols.simple). These files plus a layout/chromosome table are used by jcvi.graphics.karyotype to render karyotype/link plots that visualize intra‑ and inter‑chromosomal synteny across samples.

The output of this section will be a WGA: 

![wga](/figures/20260409_WGA_N11.png)

___

Run the alignment:

```bash
#!/bin/bash

#SBATCH --time=2-00:00:00   
#SBATCH --nodes=1  
#SBATCH --cpus-per-task=4
#SBATCH --mem=16Gb
#SBATCH --partition=ceres
#SBATCH --account=coffea_pangenome

# Check if the correct number of arguments is provided
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <first_pair> <second_pair>"
    exit 1
fi

#Submit as cat Pairs.list  | xargs -L 1 sbatch 01_WGA.sh
#module load miniconda
#source activate wga #for mummer
#source activate merothon #for last

# This comes from improvements on the parakeet synteny scripts from https://github.com/lurebgi/monkParakeet
# Use a file with $SAMPLE_1 \t $SAMPLE_2 within pairs.list, submit
target=$1  #reference
query=$2 #target
GENOME_DIR=/project/coffea_pangenome/Artocarpus/Comparative_Paper/wga
echo "Working on $target aligned to $query"

mkdir -p 1kb_merge1mb
cd 1kb_merge1mb
MAXGAP=1000000
MINLEN=2500
MINIDENT=95

# generating nucmer alignments
# nucmer -l 50 -b 400  -t 10 ${GENOME_DIR}/${target}.pri.chr.fa  ${GENOME_DIR}/${query}.pri.chr.fa -p ${target}-${query}
# delta-filter -1 -l ${MINLEN} ${target}-${query}.delta > ${target}-${query}.delta.filt
# show-coords -H -c -l -o -r -T ${target}-${query}.delta.filt > ${target}-${query}.delta.filt.coords

# Filter 
awk -v OFS="\t" -v minlen=$MINLEN -v minid=$MINIDENT '
{
  tchr=$12; tstart=$1; tend=$2
  qchr=$13; qstart=$3; qend=$4
  tlen=$5; qlen=$6
  if (tlen<minlen || qlen<minlen || $7<=minid) next

  torient = (tstart < tend ? "+" : "-")
  qorient = (qstart < qend ? "+" : "-")
  strand  = (torient == qorient ? "+" : "-")

  if (qorient == "-") { tmp=qstart; qstart=qend; qend=tmp }
  if (torient == "-") { tmp=tstart; tstart=tend; tend=tmp }

  print tchr, tstart, tend, qchr, qstart, qend, strand
}' ${target}-${query}.delta.filt.coords | sort -k1,1 -k2,2n -k4,4 -k5,5n > ${target}-${query}.delta.filtered.sorted.tsv

# merge any links that are within MAXGAP on the same chr, for both ref and query 
awk -v OFS="\t" -v maxgap=$MAXGAP '
function flush() {
    if (rt != "") {
        print rt, rs, re, qt, qs, qe, str
    }
}
NR==1 {
    rt=$1; rs=$2; re=$3;
    qt=$4; qs=$5; qe=$6;
    str=$7;
    next
}
{
    if ($1==rt && $4==qt && $7==str && $2 - re <= maxgap && $5 - qe <= maxgap) {
        # extend block
        re = ($3 > re ? $3 : re)
        qe = ($6 > qe ? $6 : qe)
    } else {
        flush()
        rt=$1; rs=$2; re=$3;
        qt=$4; qs=$5; qe=$6;
        str=$7;
    }
}
END { flush() }
' ${target}-${query}.delta.filtered.sorted.tsv > ${target}-${query}.merged.tsv

# Reference and query BEDs
> ${target}-${query}_ref.bed
> ${target}-${query}_qry.bed
> ${target}-${query}.simple

awk -v OFS="\t" -v prefix="${target}-${query}" '
{
    # reference anchors
    print $1,$2,$2+1,$1"__"$2,0,"+" >> (prefix "_ref.bed")
    print $1,$3,$3+1,$1"__"$3,0,"+" >> (prefix "_ref.bed")
    # query anchors
    print $4,$5,$5+1,$4"__"$5,0,"+" >> (prefix "_qry.bed")
    print $4,$6,$6+1,$4"__"$6,0,"+" >> (prefix "_qry.bed")
}
' ${target}-${query}.merged.tsv

# Simple link file
awk -v OFS="\t" -v prefix="${target}-${query}" '
{
    print $1"__"$2,$1"__"$3,$4"__"$5,$4"__"$6,($3-$2),$7
}
' ${target}-${query}.merged.tsv > ${target}-${query}.simple

# Color inter-chromosomal hits
awk '{
  split($1,a,"__"); split($3,b,"__");
  if(a[1] == b[1]) print $0;
  else print "g*"$0;
}' ${target}-${query}.simple > ${target}-${query}.cols.simple

raw=$(cat ${target}-${query}.delta.filt.coords | wc -l)
filtered=$(cat ${target}-${query}.delta.filtered.sorted.tsv | wc -l)
merged=$(cat ${target}-${query}.simple | wc -l)
echo "${raw} raw links; ${filtered} filtered links; ${merged} merged blocks"
```

Submit with:

```
sed '1d' Pairs.list  | xargs -L 1 sbatch 02_Mummer.sh
```

Combine the outputs: 

Create a layout file:

```
cat chr_layout.txt 
#       y,      xstart, xend,   rotation,       color,  label,  va,     bed
        0.92,   0.2,    0.95,   0,      ,       HART001,        top,    HART001.bed
        0.83,   0.2,    0.95,   0,      ,       HART067,        top,    HART067.bed
        0.75,   0.2,    0.95,   0,      ,       HART063,        top,    HART063.bed
        0.67,   0.2,    0.95,   0,      ,       HART061,        top,    HART061.bed
        0.58,   0.2,    0.95,   0,      ,       HART062,        top,    HART062.bed
        0.5,    0.2,    0.95,   0,      ,       HART058,        top,    HART058.bed
        0.42,   0.2,    0.95,   0,      ,       HART027,        top,    HART027.bed
        0.33,   0.2,    0.95,   0,      ,       HART060,        top,    HART060.bed
        0.25,   0.2,    0.95,   0,      ,       N97_50, top,    N97_50.bed
        0.17,   0.2,    0.95,   0,      ,       HART068,        top,    HART068.bed
        0.08,   0.2,    0.95,   0,      ,       N15_23, bottom, N15_23.bed
#       edges
e,      0,      1,      HART001-HART067.cols.simple
e,      1,      2,      HART067-HART063.cols.simple
e,      2,      3,      HART063-HART061.cols.simple
e,      3,      4,      HART061-HART062.cols.simple
e,      4,      5,      HART062-HART058.cols.simple
e,      5,      6,      HART058-HART027.cols.simple
e,      6,      7,      HART027-HART060.cols.simple
e,      7,      8,      HART060-N97_50.cols.simple
e,      8,      9,      N97_50-HART068.cols.simple
e,      9,      10,     HART068-N15_23.cols.simple
```

And the chrs.txt file:

```
cat chrs.txt 
Chr01,Chr02,Chr03,Chr04,Chr05,Chr06,Chr07,Chr08,Chr09,Chr10,Chr11,Chr12,Chr13,Chr14,Chr15,Chr16,Chr17,Chr18,Chr19,Chr20,Chr21,Chr22,Chr23,Chr24,Chr25,Chr26,Chr27,Chr28
Chr01,Chr02,Chr03,Chr04,Chr05,Chr06,Chr07,Chr08,Chr09,Chr10,Chr11,Chr12,Chr13,Chr14,Chr15,Chr16,Chr17,Chr18,Chr19,Chr20,Chr21,Chr22,Chr23,Chr24,Chr25,Chr26,Chr27,Chr28
Chr01,Chr02,Chr03,Chr04,Chr05,Chr06,Chr07,Chr08,Chr09,Chr10,Chr11,Chr12,Chr13,Chr14,Chr15,Chr16,Chr17,Chr18,Chr19,Chr20,Chr21,Chr22,Chr23,Chr24,Chr25,Chr26,Chr27,Chr28
Chr01,Chr02,Chr03,Chr04,Chr05,Chr06,Chr07,Chr08,Chr09,Chr10,Chr11,Chr12,Chr13,Chr14,Chr15,Chr16,Chr17,Chr18,Chr19,Chr20,Chr21,Chr22,Chr23,Chr24,Chr25,Chr26,Chr27,Chr28
Chr01,Chr02,Chr03,Chr04,Chr05,Chr06,Chr07,Chr08,Chr09,Chr10,Chr11,Chr12,Chr13,Chr14,Chr15,Chr16,Chr17,Chr18,Chr19,Chr20,Chr21,Chr22,Chr23,Chr24,Chr25,Chr26,Chr27,Chr28
Chr01,Chr02,Chr03,Chr04,Chr05,Chr06,Chr07,Chr08,Chr09,Chr10,Chr11,Chr12,Chr13,Chr14,Chr15,Chr16,Chr17,Chr18,Chr19,Chr20,Chr21,Chr22,Chr23,Chr24,Chr25,Chr26,Chr27,Chr28
Chr01,Chr02,Chr03,Chr04,Chr05,Chr06,Chr07,Chr08,Chr09,Chr10,Chr11,Chr12,Chr13,Chr14,Chr15,Chr16,Chr17,Chr18,Chr19,Chr20,Chr21,Chr22,Chr23,Chr24,Chr25,Chr26,Chr27,Chr28
Chr01,Chr02,Chr03,Chr04,Chr05,Chr06,Chr07,Chr08,Chr09,Chr10,Chr11,Chr12,Chr13,Chr14,Chr15,Chr16,Chr17,Chr18,Chr19,Chr20,Chr21,Chr22,Chr23,Chr24,Chr25,Chr26,Chr27,Chr28
Chr01,Chr02,Chr03,Chr04,Chr05,Chr06,Chr07,Chr08,Chr09,Chr10,Chr11,Chr12,Chr13,Chr14,Chr15,Chr16,Chr17,Chr18,Chr19,Chr20,Chr21,Chr22,Chr23,Chr24,Chr25,Chr26,Chr27,Chr28
Chr01,Chr02,Chr03,Chr04,Chr05,Chr06,Chr07,Chr08,Chr09,Chr10,Chr11,Chr12,Chr13,Chr14,Chr15,Chr16,Chr17,Chr18,Chr19,Chr20,Chr21,Chr22,Chr23,Chr24,Chr25,Chr26,Chr27,Chr28
Chr01,Chr04,Chr06,Chr08,Chr10,Chr11,Chr13,Chr15,Chr18,Chr19,Chr22,Chr24,Chr26,Chr28
```

And create karyoplot:

```bash
for i in $(cat Samples.list ); do cat ${i}.bed*-* > ${i}.bed; done

# And add a color (e.g. g=green) for inter-chromosomal translocations
for i in $(ls *simple | sed 's/.simple//g'); do 
awk '{
  split($1,a,"__"); split($3,b,"__");
  if(a[1] == b[1]) print $0;
  else print "g*"$0;
}' ${i}.simple > ${i}.cols.simple
done 

python -m jcvi.graphics.karyotype chrs.txt chr_layout.txt
```

