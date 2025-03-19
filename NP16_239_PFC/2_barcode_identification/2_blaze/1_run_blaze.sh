#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH -D . # set working directory to .
#SBATCH -p mrcq # submit to the parallel queue
#SBATCH --time=20:00:00 # maximum walltime for the job
#SBATCH -A Research_Project-MRC148213 # research project to submit under
#SBATCH --nodes=1 # specify number of nodes
#SBATCH --ntasks-per-node=16 # specify number of processors per node
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=sl693@exeter.ac.uk # email address
#SBATCH --output=1_run_blaze.o
#SBATCH --error=1_run_blaze.e

module load Miniconda2
source activate nanopore

#cd /lustre/projects/Research_Project-MRC148213/lsl693/software
#wget -O cellranger-arc-2.0.2.tar.gz "https://cf.10xgenomics.com/releases/cell-arc/cellranger-arc-2.0.2.tar.gz?Expires=1741930379&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA&Signature=FsvZ-pdW2ciMj1ID9H8mGzMlKAu2cndFMg5MBK~PHNfqT7-20yeZ9VqkoWnr67Kx4np9IY2MC374eySGow1aOcriQHBgd8dF3zyRFnJ51malsJY0bTqsaA1XsQXAWOjn83Lj2JDM-0P1FxiMqrmb5LDMj48-eI2fzBykCppFBQxnYoKbxlWS53DFqAPZvzJaaOlGeYXyBES3IgWW7Kqjo9h5-OBzQ7unsP96jDnA2nElonBloW-yTzN-Ppz9fgNHnx7T0ebM3hTstcEd7EufhOzBMKzffzxvxbEfi6kwRxHv7Sfb3RZw16M9idBm~RDENKlvPtcsBtNx7022nNPaNQ__"

# Human reference (GRCh38) - 2024-A
cd 
wget "https://cf.10xgenomics.com/supp/cell-arc/refdata-cellranger-arc-GRCh38-2024-A.tar.gz"
tar -xvzf refdata-cellranger-arc-GRCh38-2024-A.tar.gz

# raw short-read FASTQ
SRFASTQ=/lustre/projects/Research_Project-MRC190311/SingleCell/BMIAMP/1_raw/project_11062/11062/V0311-gex-demux/11062

# raw long-read FASTQ basecalled
LRFASTQ=/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/Gina/basecalled_sequence

# run BLAZE
# number of cells from short-read sequencing data
cd /lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/Gina/analysis/6_blaze/1_blaze
whitelist=/lustre/projects/Research_Project-MRC148213/lsl693/software/cellranger-arc-2.0.2/lib/python/cellranger/barcodes/737K-arc-v1.txt 
blaze --force-cells 1378 --output-prefix NP16_239_PFC --threads=12 ${LRFASTQ} --full-bc-whitelist ${whitelist} &> blaze.log