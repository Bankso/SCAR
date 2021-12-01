# SpLiT-ChEC Analysis with R (SCAR)
## Version 1.3
### Why does this exist?

This repository houses SCAR, a NGS data analysis pipeline based on ZentTools from gzentner and rpolicastro. 

The goal of this project is to optimize analysis of low-background, whole-genome sequencing datasets derived from the in-development factor mapping method SpLiT-ChEC.

- A Singularity container is used to bundle software and simplify reproducibility  
- The code is primarily in R and is installed as an R package within the virtual environment generated by Singularity  
- Scripts to interface with SCAR can be found in the 'ex_files' directory and you can find out how to use them below!

### How to use this software

It is recommended that the analysis outlined here is performed on a HPC, considering many of the calculations are quite intensive. Run submission can be carried out via the job management system for you cluster. However; I have had some success running several of the post-alignment functions on a single-processor desktop with sufficient RAM.

Singularity must be available for this software to function. A singularity image file (SIF) allows you to create a new "clone" environment for each run, which improves consistency and reproducibility.

Instructions to install Singularity can be found [here](https://sylabs.io/guides/3.5/user-guide/quick_start.html#quick-installation-steps).  

**Note:** You probably can't install this (or anything else) on your HPC if you are a non-admin user, so contact your admin if it's not available to you and you want it! For the UO Talapas cluster, the availability of modules can be determined with the command
```
module avail
```
### Processing data with SCAR
**Note/Disclaimer:** This is my own example working directory, meaning it might not be right for you, but the organization described here influences the run, process, and sample file organization, so this system should be used unless you are willing to redesign the processing files.

To begin, generate a clean directory. This will serve as a quasi HOME directory for data processing.
```
mkdir -p /path/to/dir
```
Within this directory (let's call it new_home), make three new directories:
```
cd /path/to/new_home
mkdir -p scripts/ samples/ controls/
```
Each of these directories should contain a specific set of files for processing: 

**FASTQ :**  
Raw sequencing data output files from ChEC/SpLiT-ChEC, ChIP, etc. The pipeline is designed for paired-end reads, but can run single-end too!  
Place FASTQ files directly into 'sample' or 'control' folder, depending on if the file is from a sample or control 

**samples.txt :**  
A tab-separated file or data frame with sample names and paths to files  
(Template in ex_files)  
Place a **single** samples.txt file in the 'samples' directory, nowhere else.

**run_start.sh :**  
Bash script to submit jobs via slurm, manage external files, and run Singularity  
(Template in ex_files)  
Place directly in scripts directory

**process_options.R :**  
Specifies the functions and inputs that will be passed to SCAR. Change the inputs to fit your analysis, but if you change the code  
format or structure, you are on your own if it goes weird!  
(Template in ex_files)  
Place directly in scripts directory

Once all the files are in place, a run can be initiated with the following command format:
```
sbatch scripts/run_start.sh $PWD scripts/process_options.R Run_name sample_dir sample_dir/samples.txt
```
The analysis types flagged as 'TRUE' in process_settings.R will be performed, at the moment these include:  
FASTQ -> BAM (samples and controls, paired or unpaired reads)  
  
BAM -> bigwig/bedgraph (samples and controls via bamCompare, bamCoverage, bigwigCompare)
  
BAM/bedgraph -> peaks (requires both sample and control bams or bedgraphs, via SEACR)  
  
Visualization options with deeptools will be added later!  
  
Following the conclusion of each processing step, the resultant file will be dumped to the outdir defined in the options script. By default, this is some variant of sample/new_dir/out_file, where new_dir is a term corresponding to the outfile type (e.g., aligned/, coverage/, etc.)

*More info to be added*
