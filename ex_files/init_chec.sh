#!/bin/bash

#SBATCH --job-name=Med14_30s     		### Job Name
#SBATCH --output=M30s.out         ### File in which to store job output
#SBATCH --error=M30s.err          ### File in which to store job error messages
#SBATCH --time=0-01:00:00       	### Wall clock time limit in Days-HH:MM:SS
#SBATCH --nodes=1               	### Number of nodes needed for the job
#SBATCH --ntasks-per-node=1     	### Number of tasks to be launched per Node
#SBATCH --cpus-per-task=8			### Number of CPU's used to complete submitted tasks within each requested node
#SBATCH --mem=12GB				### Amount of preallocated memory for submitted job
#SBATCH --account=mcknightlab     ### Account used for job submission

#This script initiates a set of actions converting fastq -> bam + bam.bai
#Can be run as sbatch or in chunks to do only what you need
#Before running, edit file path at line 17 to be your data directory and copy edited file to new directory

cd /projects/mcknightlab/obanks/analysis/SpLiT_ChEC/Med14/30s/SCAR_test/

mkdir -p genome && cd genome

wget ftp://ftp.ensembl.org/pub/release-100/fasta/saccharomyces_cerevisiae/dna/Saccharomyces_cerevisiae.R64-1-1.dna_sm.toplevel.fa.gz

gunzip Saccharomyces_cerevisiae.R64-1-1.dna_sm.toplevel.fa.gz

wget ftp://ftp.ensembl.org/pub/release-100/gtf/saccharomyces_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.100.gtf.gz

gunzip Saccharomyces_cerevisiae.R64-1-1.100.gtf.gz

cd /projects/mcknightlab/obanks/analysis/common/

module load singularity

if [ ! -e ./scar_v* ]; then
	singularity pull library://banksorion/default/scar_v0.1:sha256.12d36f4f6c72aeb2bc7ffe1c17079ac64eaeea4e5d5e117bc0e1b115bc51da4e

else
	break

fi

cd /projects/mcknightlab/obanks/analysis/SpLiT_ChEC/Med14/30s/SCAR_test/

singularity exec -eCB `pwd` -H `pwd` /projects/mcknightlab/obanks/analysis/common/scar_v* Rscript chec_process.r

exit
