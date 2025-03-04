#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH -D . # set working directory to .
#SBATCH -p bioseq # partition to use
#SBATCH --time=24:00:00 # maximum walltime
#SBATCH -A Research_Project-BioProduction # research project to submit under
#SBATCH --nodes=1 # specify number of nodesq pq
#SBATCH --exclusive
#SBATCH --job-name=cr_count


# multiome specific reference
GRCh38=/data2/databases/10x/refdata-cellranger-arc-GRCh38-2020-A-2.0.0

GRCh38=/lustre/projects/Research_Project-BioProduction/references/refdata-cellranger-arc-GRCh38-2020-A-2.0.0

# for bioseqfs02:
CRARC=/data1/software/cellranger-arc-2.0.2/cellranger-arc

# for ISCA
CRARC=/lustre/projects/Research_Project-BioProduction/software/cellranger-arc-2.0.2/cellranger-arc

sample=11062_NP16_239_PFC
echo ${sample}
$CRARC count \
   --id ${sample} \
   --reference=${GRCh38} \
   --libraries=libraries.csv
