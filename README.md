# FLAQ-SC2-Meta
FL BPHL's SARS-CoV-2 (SC2) analysis pipeline for Illumina paired-end, whole-genome tiled-amplicon data enriched from environmental metagenomic samples (e.g., wastewater). 

## About
FLAQ-SC2-Meta was developed to analyze Illumina paired-end, whole-genome tiled-amplicon data (i.e., [ARTIC protocol](https://artic.network/ncov-2019)) from wastewater samples. The pipeline generates variant files along with reports including read/mapping quality metrics and estimated Pango lineage abundances. The current version will run only on [HiPerGator](https://www.rc.ufl.edu/about/hipergator/)(HPG) using local Singularity containers for each pipeline process.

## Dependencies
- Python3
- Singularity/Apptainer
- [Freyja](https://github.com/andersen-lab/Freyja)
- Git

To load python3 into your current environment on HiPerGator, either use `module load python` to get the lastest version of python or activate your base conda environment. For more information on how to set up your base conda environment on HPG, see the [HiPerGator Analysis Reference Guide](https://github.com/StaPH-B/southeast-region/tree/master/hipergator)).

Singularity/Apptainer will be loaded as a module during your job execution on HPG using the sbatch job script in this repository. 

Freyja can be installed via [conda](https://anaconda.org/bioconda/freyja). Activate your freyja conda environment with `conda activate freyja`.

Git is already installed in your HPG environment upon login.

## Usage

For first time use, clone this repository to a directory in blue on HPG, such as in /blue/bphl-\<state\>/\<user\>/repos/bphl-molecular/.
```
cd /blue/bphl-<state>/<user>/repos/bphl-molecular/
git clone https://github.com/BPHL-Molecular/flaq_sc2_meta.git
```
For future use, update any changes to your local repository on HPG by navigating to the flaq_sc2_meta repository and pulling any changes.
```
cd flaq_sc2_meta/
git pull
```
To run the FLAQ-SC2-Meta pipeline, copy all files from the flaq_sc2_meta local repository to your analysis folder. Make an input directory and copy your fastqs.
```
mkdir <analysis_dir>
cd <analysis_dir>
cp /blue/bphl-<state>/<user>/repos/bphl-molecular/flaq_sc2_meta/* .
mkdir fastqs/
cp /path/to/fastqs/*.fastq.gz fastqs/
```
Rename your fastq files to the following format: sample_1.fastq.gz, sample_2.fastq.gz. See below for a helpful command to rename your R1 and R2 files.
```
cd fastqs/
for i in *_R1_001.fastq.gz; do mv -- "$i" "${i%[PATTERN to REMOVE]}_1.fastq.gz"; done
for i in *_R2_001.fastq.gz; do mv -- "$i" "${i%[PATTERN to REMOVE]}_2.fastq.gz"; done
```
Edit your sbatch job submission script to include your email to receive an email notification upon job END or FAIL. Replace ENTER EMAIL in `#SBATCH --mail-user=ENTER EMAIL` with your email address. Make sure there is no space between = and your email address. Edit additional sbatch parameters as needed to run your job succesfully, such as the length of time the job will run.

Submit your job.
```
sbatch sbatch_flaq_sc2_meta.sh
```

## Main processes
- [Fastqc](https://github.com/s-andrews/FastQC)
- [Trimmomatic](https://github.com/usadellab/Trimmomatic)
- [BBDuk](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/bbduk-guide/)
- [BWA](https://github.com/lh3/bwa)
- [Samtools](https://github.com/samtools/samtools)
- [iVar](https://github.com/andersen-lab/ivar)
- [Freyja](https://github.com/andersen-lab/Freyja)

## Primary outputs

Outputs from each process for each individual sample can be found in a sample-specific subdirectory within the FLAQ-SC2-Meta analysis directory. Report.txt contains the main summary report with read/mapping quality metrics. Additional details can be found in the report outputs from each process such as variant files (.variant.tsv) and Freyja reports.
```
analysis_dir/
|__ <date>_flaq_run/
     |__ report.txt
     |__ sample1/
     |__ sample2/
```

## Developed by:
[@SESchmedes](https://www.github.com/SESchmedes)<br />

Please email bphl16BioInformatics@flhealth.gov for questions and support.
