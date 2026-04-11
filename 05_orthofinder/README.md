# Orthofinder

Run orthofinder on our n=10 Artocarpus accessions and our Batocarpus sample, using the protein data from above. 

This is pretty much just to get the initial species tree among our samples. 

```bash
#!/bin/bash

#SBATCH --time=2-00:00:00   
#SBATCH --nodes=1  
#SBATCH --ntasks-per-node=20
#SBATCH --mem=128Gb
#SBATCH --partition=ceres
#SBATCH --account=coffea_pangenome

#module load miniconda
#source activate orthofinder

WD=/project/coffea_pangenome/Artocarpus/Comparative_Paper/orthofinder

t=20
RUN=$1
cd ${WD} 

orthofinder -f ${RUN} -a ${t}

# Afterwards, clean up for long term storage since this generates thousdands of files 
#find proteomes/${REF} -type d \( -name "MultipleSequenceAlignments" -o -name "Gene_Trees" -o -name "WorkingDirectory" -o -name "Orthogroup_Sequences" -o -name "Orthologues" -o -name "Resolved_Gene_Trees" \) -exec rm -rf {} +
```

