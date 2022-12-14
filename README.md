## Stenglein lab taxonomic assessment pipeline


This is a nextflow implementation of the pipeline used in the [Stenglein lab](http://www.stengleinlab.org) to taxonomically classify sequences in NGS datasets.

It is mainly designed to identify virus sequences in a metagenomic dataset but it also performs a general taxonomic classification of sequences.

A [previous bash-based version of this pipeline](https://github.com/stenglein-lab/taxonomy_pipeline/releases/tag/2022_3_30_bash) has been reporpted in a number of [publications](https://www.stengleinlab.org/papers/) from our lab. 

## How to run the pipeline


#### Step 1. Clone the repository:

First you need to get the pipeline code.  You do this using git:
```
git clone https://github.com/stenglein-lab/taxonomy_pipeline.git 2022_3_analysis
```
This will download a copy of the pipeline files and put them in the directory `2022_3_analysis` but this is just an example name: you can name the directory anything you want. 
 
#### Step 2. Move your fastq files to the new directory 

The pipeline expects your fastq files to be in a subdirectory named input/fastq.  The next step will be to move your fastq files there.

```
# transfer fastq files from an example location to the pipeline directory
mv /home/mdstengl/raw_data/*.fastq.gz 2022_3_analysis/input/fastq
```

The fastq files should contain Illumina reads, but can be single or paired-end, or a mix, or cogzip compressed (end in .gz) or not, or a mix of compressed and uncompressed.  It would be a best practice to keep your fastq files gzipped to take up less space on server storage.  

Note that the pipeline looks for files with names that match the pattern `*_R[12]_*.fastq*`.  You can change this pattern using the argument `--fastq_pattern` as input to the nextflow command.

#### Step 3. Setup host filtering

The pipeline optionally removes host-derived reads because generally they are not what we are interested in and removing these reads makes downstream taxonomic classification steps go faster.  To perform host filtering, you will need a bowtie2 index of a host reference genome on the same server as you are running the pipeline.  [This tutorial section](https://github.com/stenglein-lab/taxonomy_pipeline/blob/main/docs/tutorial.md#section_genome) describes how to download a new host genome and build a bowtie2 index if there is not one already available.

To specify which bowtie2 index will be used for your datasets you should setup a file that maps dataset names to host genomes.  An [example of this file is here](map://github.com/stenglein-lab/taxonomy_pipeline/blob/nextflow/input/host_mapping.txt).  This example contains two columns separated by a tab.  The first column contains a pattern that will be searched for in fastq file names.  The second columns contains the location of a bowtie2 index.

```
# view an example host filtering map file
$ cat input/host_mapping.txt
_M_	/home/databases/fly/Dmel_genome_transcriptome
_F_	/home/databases/fly/Dmel_genome_transcriptome
Hela	/home/databases/human/GCRh38
HeLa	/home/databases/human/GCRh38
```

In this example, any datasets whose fastq file names contain the text `_M_` or `_F_` will be host-filtered using the bowtie2 index located at `/home/databases/fly/Dmel_genome_transcriptome` and any datasets whose fastq file names contain the text `Hela` or `HeLa` will be filtered using the GCRh38 human reference genome.  Note that bowtie2 indexes are actually made up of files and the paths listed in the example above are the prefix of all those files (they are the path that will be passed to the bowtie2 `-x` command line argument).

```
# bowtie2 indexes are made up of 6 files 
$ ls /home/databases/human/GCRh38.*
/home/databases/human/GCRh38.1.bt2
/home/databases/human/GCRh38.2.bt2
/home/databases/human/GCRh38.3.bt2
/home/databases/human/GCRh38.4.bt2
/home/databases/human/GCRh38.rev.1.bt2
/home/databases/human/GCRh38.rev.2.bt2
```

You will need to setup this host_mapping file to specify which host reference genomes (bowtie2 indexes) should be used to remove host reads from your datasets.  You can edit the provided `input/host_mapping.txt` file or you can create a new file and specify its location using the `--host_map_file` command line argument when you run nextflow.

If you fail to specify a host mapping file the pipeline will warn you that no host filtering will be performed.  If you do specify host filtering, you don't necessarily need to perform filtering for all datasets: Any dataset whose fastq file name doesn't match one of the patterns in the host map file will not be host filtered.

#### Step 4. Actually run the pipeline

Now you actually need to run the pipeline.  An example command to run it:

```
# assuming you are in the pipeline main directory (2022_3_analyis in the example above)
nextflow run main.nf -resume -profile singularity,local --host_map_file input/host_mapping.txt
```

This assumes you are working on a server with nextflow and conda installed.  See the Dependencies below for more information on these dependencies.  

Because the pipeline will likely run for a while (depending on dataset sizes), you will want to run it via a screen session.  See [this tutorial section](https://github.com/stenglein-lab/taxonomy_pipeline/blob/main/docs/tutorial.md#section_screen) for more information on using screen.

## Pipeline output 

The pipeline puts output in a `results` directory.

- Virus sequences: in the `results/virus_sequences` directory
- Taxonomic tabulation: in the `tallies` directory.  These are tab-delimited tables of the taxa observed in each dataset, the # of mapping reads, average percent identity to database sequences, etc.
- Contigs: in the `contigs` directory.  Contigs from assembly of reads remaining after host filtering.
- QC reports: `initial_qc_report.html` and `post_trim_qc_report.html
- Plots of # and fraction reads remaining after filtering steps: `filtering_plots.pdf`

TODO - flesh out this section

## Dependencies

### Nextflow

This pipeline is implemented in Nextflow.  To run the pipeline you will need to be working on a computer that has nextflow installed. Installation instructions can be [found here](https://www.nextflow.io/docs/latest/getstarted.html#installation).  To test if you have nextflow installed, run:

```
nextflow -version
```

### Other software dependencies

This pipeline deals with other software dependencies in 2 possible ways: by using Singularity containers or by using an all-in-one conda environment.  You are free to use either one of these when running the pipeline.

#### Singularity containers

The pipeline can use singularity containers to run programs like bowtie2 and blast.  To use these containers you must be running the pipeline on a computer that has [singularity](https://sylabs.io/singularity) [installed](https://sylabs.io/guides/latest/admin-guide/installation.html).  To test if singularity is installed, run:

```
singularity --version
```

To run with singularity containers include the option `-profile singularity` in the nextflow command line, for instance:

```
nextflow run main.nf -profile singularity
```

Singularity containers will be automatically downloaded and stored in a directory named `singularity_cacheDir` in your home directory.  They will only be downloaded once. 

#### Conda environment

The pipeline can also use an all-in-one conda environment.  This requires conda to be installed on your computer.  To test if conda is installed, run:

```
conda -V
```

The conda environment is defined in [this yaml file](./conda/taxonomy_conda_environment.yaml) and will be automatically created if you run the pipeline using the conda profile.  To run the pipeline with conda, include `-profile conda` in the nextflow command line, for instance:

```
nextflow run main.nf -profile conda
```

The conda environment will be created in a directory in your home directory named `conda_cacheDir` and will only be created once.

You should specify either `-profile conda` or `-profile singularity` or the pipeline will output an error message and halt.

#### Custom scripts

The pipeline also uses custom R and perl scripts, which are located in the [scripts](./scripts) directory of this repository.


### Database dependencies

This pipeline uses two databases for taxonomic classification:

(1) The [NCBI nt nucleotide database](https://ftp.ncbi.nlm.nih.gov/blast/db/) for use in blastn searches.  

(2) A [diamond](https://github.com/bbuchfink/diamond) database created from the NCBI nr protein sequence database.

(3) A local version of the [NCBI Taxonomy](https://ftp.ncbi.nih.gov/pub/taxonomy/) database.

These databases need to be installed locally on a computer or server to run this pipeline.  This requires ~1.2 Tb of disk space (as of early 2022).  A script to download these databases from NCBI [is included](./scripts/download_and_process_sequence_databases) in this repository.  These databases take up a lot of space, so before doing this make sure that these databases are not already downloaded.

#### Database locations

##### nt database

The default location of the blast nt databases is: `/home/databases/nr_nt/nt`. This value will be passed to the -db option when running blastn. This default value can be overridden using the `--nt_blast_db` parameter, for instance:

```
nextflow run main.nf -profile conda --nt_blast_db /path/to/my/local/db/nt
```
In this example, there should be a file named `/path/to/my/local/db/nt.nal` in addition to the other nt database files.

##### diaomond nr database

The default location of the diamond database is: `/home/databases/nr_nt/nr.dmnd`.  This path will be passed to the --db option when running diamond.  This default value can be overridden using the `--nt_diamond_db` parameter.

##### Taxonomy db

TODO - fill out this section

## Assumptions about input data

1. This pipeline is designed to work with Illumina data.  

## Testing

TODO - fill out this section

## Tutorial

There is [a tutorial](./docs/tutorial.md) that describes the main steps of this pipeline.  This tutorial refers to the older, pre-nextflow bash version of this pipeline, but the steps are similar and so the information in the tutorial remains applicable to understanding how the pipeline works in principle.


