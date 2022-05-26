#!/usr/bin/bash

#SBATCH --job-name=SCAR	    		 	### Job Name
#SBATCH --output=scripts/logs/run.out           ### File in which to store job output DO NOT EDIT
#SBATCH --error=scripts/logs/run.err            ### File in which to store job error messages DO NOT EDIT
#SBATCH --time=0-01:00:00           		### Wall clock time limit in Days-HH:MM:SS
#SBATCH --nodes=1                   		### Number of nodes needed for the job
#SBATCH --ntasks-per-node=1         		### Number of tasks to be launched per Node
#SBATCH --cpus-per-task=8	    		### Number of CPU's used to complete submitted tasks within each requested node
#SBATCH --mem=12GB		    		### Amount of preallocated memory for submitted job
#SBATCH --account=your_account       		### Account used for job submission

## Be sure to edit the above to fit your HPC job submisssion framework
## For example directory structure, see the README at https://github.com/Bankso/SCAR or https://github.com/Bankso/mock_home
## Usage: sbatch run.sh 1 2 3 4 5 6
## Requires 6 inputs at the command line:
	## For start bash script:
		## First input - Home directory that holds sample directories and will be target for sif, genome downloads
		## Second input - Path to options.R, determines settings used during processing

	## For Singularity environment, used in Rscript:
		## Third input - Name of tag to be used on end of plotting directories and for logs, should be unique for each run 
		## Fourth input - Sample file directory path, relative to home (i.e., for $home_dir/protein; protein=sample_dir=$2)
		## Fifth input - path to a BED formatted file to be used for plotting
		## Sixth input - fragment set to be used in processing small, large, full, or all. Bounds are set in options.R
		

home_dir=$1  
process_script=$2
job_name=$3
sample_dir=$4
samples_file=$4/samples.txt
plot_regions=$5
frag_set=$6

cd $home_dir

#mkdir -p genome && cd genome

#wget ftp://ftp.ensembl.org/pub/release-100/fasta/saccharomyces_cerevisiae/dna/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa.gz

#gunzip Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa.gz

#wget ftp://ftp.ensembl.org/pub/release-100/gff3/saccharomyces_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.100.gff3.gz

#gunzip Saccharomyces_cerevisiae.R64-1-1.100.gff3.gz

#cd $home_dir

module load singularity

if [ ! -e ./scar* ]; then
	
	singularity cache clean -f
	singularity pull --arch amd64 library://banksorion/default/scar_software:1.5_latest
fi

# Start singularity-based processing run 
singularity exec --env sample_dir=$sample_dir,samples_file=$samples_file,plot_regions=$plot_regions,job_name=$job_name,frag_set=$frag_set -eCB $home_dir:/opt/conda/process -H /opt/conda/process *.sif Rscript $process_script

# Move sbatch logs to sample dir for storage
mkdir -p $sample_dir/logs/ && mv scripts/logs/run.out $sample_dir/logs/$job_name.out && mv scripts/logs/run.err $sample_dir/logs/$job_name.err

# Clean-up from processing
rm -r ./temp/

exit
