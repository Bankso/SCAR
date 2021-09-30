# SpLiT-ChEC Analysis with R (SCAR)

### Why does this exist?

This repository houses SCAR, a pipeline based on ZentTools from gzentner and rpolicastro. 

The goal of this project is to optimize analysis of low-background whole-genome sequencing datasets derived from the in-development factor mapping method SpLiT-ChEC.

- A Singularity container is used to bundle software and simplify reproducibility
- The code is primarily in R and is installed as an R package within the virtual environment generated by Singularity
- Scripts to interface with SCAR can be found in the 'ex_files' directory and you can find out how to use them below!

### How to use this software

#### On a HPC (strongly recommended)

It is recommended that you use SCAR on a HPC. Run submission can be carried out via slurm or similar.

Singularity must be available for this software to function. A singularity image file (SIF) creates an identical environment for each run, allowing you to consistently analyze data in an identical environment.

Instructions to install Singularity can be found [here](https://sylabs.io/guides/3.5/user-guide/quick_start.html#quick-installation-steps).  

**Note:** You probably can't install this  (or anything else) on your HPC if you are a non-admin user, so contact your admin if it's not available to you and you want it! For the UO Talapas cluster, the availability of modules can be determined with the command
```
module avail

```

### Processing data with SCAR

To begin, generate a clean directory to work within and include the files detailed below

#### File types required:

**FASTQ**  
Raw sequencing data output files for samples and controls (if desired)

**samples.txt**  
A tab-separated sheet or data frame with sample names and file names to be analyzed

**init.sh**  
Bash script to submit jobs via slurm, acquire genome files, and run Singularity

**process.R**  
The "heavy-lifter" script, which contains all the commands to process the data by interacting with the SCAR package

The processing steps are outlined below (to be completed expanded on later...)

#### Collect yeast genome assembly and annotation

#### Initiate Singularity

#### Create the ZentTools object

#### Perform FASTQ quality qontrol with FastQC

#### Align reads to genome with Tophat2

#### Make tracks with Deeptools

#### Peak calling with SEACR and MACS2





