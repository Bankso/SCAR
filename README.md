# SpLiT-ChEC Analysis with R (SCAR)
## Version 1.5
### Overview

This repository houses a next-generation sequencing (NGS) data analysis package, SpLiT-ChEC Analysis with R (SCAR). SCAR was designed to function as a reproducible bioinformatics pipeline with easily accessible runtime settings and self-documenting directories/filenames. These decisions have been made to simplify and optimize analysis for data derived from novel NGS-coupled techniques and organization of processing outputs into coherent directory structures. The set of processes available is being expanded and optimized as we design methods to analyze SpLiT-ChEC data from DNA-associating proteins in *S. cerevisiae*. 
The framework of SCAR is derived from ZentTools, created by gzentner and rpolicastro (https://github.com/rpolicastro/ZentTools)

While SCAR has been designed for use with DNA-centric NGS techniques, the pipeline is highly amenable to modification and can easily accomodate RNA-centric analysis tools. 

Current feature set and runtime options for each process (noted in parentheses):
- Retrieve and index a reference genome (primarily use SacCer3, can be adapted to other databases and organisms)
- FASTA/Q quality control analysis with fastQC
- Alignment to a reference genome using bowtie2 (options to filter input fragment sizes, change alignment mode, apply MAPQ filtering)
- Coverage calculation and normalization with deeptools2 (options for centering/not centering coverage, normalization types, bin size)
- Calculating enrichment over background with deeptools2 (options for function used in comparison)  
- Peak finding using MACS (https://github.com/macs3-project/MACS) (options for broad/standard cutoffs)
- Matrix building, heatmap and profile plotting for peaks and other input BED files using deeptools2 (options for input coverage file, input BED regions, plotting region bounds, number of kmeans clusters to calculate/plot, data sorting type and direction)

In-progress features (can be found here https://github.com/Bankso/SEAPE):
- Tools written in R, python, and Bash/shell to simplify statistical analysis, counting processes, comparison to external information sources, and pre-processing/plotting of outputs
- Implementing clustering/machine learning algorithms with R, Python and/or Tensorflow to robustly identify sub-groups in complex NGS datasets and generate reliable models for identifying features/comparison to experimental results

### How to use this software

It is recommended that SCAR is used on your local high-performance cluster or a cloud-based service, considering how long it takes to run on most single-CPU computers. I have had some success running some post-alignment functions on a single-processor desktop with sufficient RAM, but alignments generally take far too long to be done regularly on standard hardware.

A Singularity container is used to bundle software and simplify reproducibility. As such, Singularity must be available for this software to function. Instructions to install Singularity can be found [here](https://sylabs.io/guides/3.5/user-guide/quick_start.html#quick-installation-steps). 

For the UO Talapas cluster, using slurm job handling, the availability of modules can be determined with the command:
```
module avail
```
**Note:** You probably can't install Singularity (or anything else) on your HPC if you are a non-admin user, so contact your admin to get help if it isn't available.

The Singularity recipe file in this repository is used to build a signed and verified Singularity container, hosted on sylabs.io, which can be downloaded with the following command: 
```
singularity pull --arch amd64 library://banksorion/default/scar_software:1.5_latest
```
The container allows processing to be carried out independently of local software by incorporating all the programs and packages necessary for all SCAR processes to run. Each run, the SIF as accessed as a naive computing environment to carry out the desired tasks, most of which can be easily toggled at run time. This also means directory structure is important to pay attention to. 
I've outlined a basic directory structure below, but have also put together a repository that can be copied to your local environment and used in a "drag-and-drop" format (https://github.com/Bankso/mock_home/). From the start, I intended for SCAR to automate most organization decisions when handling processing output so users can focus more on what questions to ask next and less on where to fit everything in a file directory. 

### Processing data with SCAR

To begin, generate a clean directory. This will serve as a quasi HOME directory for data processing.
```
mkdir -p /path/to/dir
```
Within this directory (let's call it new_home), make a few new directories:
```
cd /path/to/new_home
mkdir -p scripts/process scripts/start samples/ controls/
```
These are the types of files used when initiating a processing run using SCAR:

**FASTA/Q :**  
Raw sequencing data output files from ChEC/SpLiT-ChEC, ChIP, etc. The pipeline is designed for paired-end reads, but can analyze single-end sequencing results as well 
Place FASTQ files directly into 'sample' or 'control' folder, depending on if the file is from a sample or control 

**samples.txt :**  
A tab-separated file or data frame with sample and/or control file names + paths to FASTA/Q files 
Place a **single** samples.txt file in the 'samples' directory, nowhere else.

**run.sh :**  
Bash script to submit jobs to a HPC, collect and manage external files, and download SIF + initiate a Singularity environment with the necessary variables and settings
Place in 'scripts/start' directory

**options.R :**  
Specifies the functions and inputs that will be passed to SCAR. I've made many useful settings easily changeable to simplify optimization of data analysis across a wide set of input paramaters. This seems to be a generally useful format for process optimization and exploration of data features.
Place 'scripts/process' directory

Once all the files are in place, a run can be initiated with the following command format (using slurm job handling):
```
sbatch scripts/start/run.sh $path_to_new_home scripts/process/options.R $name_for_output $path_to_sample_files $path_to_BED_file $fragment_range
``` 
Following the conclusion of each processing step, the resultant directories and files will be dumped to the sample directory. By default, this is some variant of sample/new_dir/out_file, where new_dir is a term corresponding to the outfile type (e.g., aligned/, coverage/, etc.)

SCAR was designed to be primarily self-documenting, so most options will be recorded in outputs or file logs. After every run, file logs are moved to a 'logs' directory in the sample directory. To make sure the log isn't overwritten with a new run, either modify the name of the log to a unique identifier or be sure to use a a new file output tag for each run with otherwise identical settings.

*More info to be added*
