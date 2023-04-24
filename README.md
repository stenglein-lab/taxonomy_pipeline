## Stenglein lab taxonomic assessment pipeline


This is a nextflow implementation of the pipeline used in the [Stenglein lab](http://www.stengleinlab.org) to taxonomically classify sequences in NGS datasets.

It is mainly designed to identify virus sequences in a metagenomic dataset but it also performs a general taxonomic classification of sequences.

A [previous bash-based version of this pipeline](https://github.com/stenglein-lab/taxonomy_pipeline/releases/tag/2022_3_30_bash) has been reporpted in a number of [publications](https://www.stengleinlab.org/papers/) from our lab. 

## How to run the pipeline

This pipeline is run using [nextflow]()

1. [The path to the directory containing your fastq files (`--fastq_dir`)](#Fastq directory)
2. [The path to a file defining how host filtering will be performed (`--host_map_file`)](#Host filtering)
3. [Whether you will use singularity or conda to provide required software tools (`-profile singularity` or `profile conda`).](#Dependencies)

Here's a basic example of a command line to run the pipeline.

```
nextflow run stenglein-lab/taxonomy_pipeline -resume -profile singularity --host_map_file host_mapping.txt --fastq_dir /path/to/directory/containing/fastq/
```
 
#### Fastq directory

You will need to provide the location of a directory containing fastq files using the `--fastq_dir` command line argument to nextflow.   These need not be in the same directory as you are running the pipeline.  For instance:

```
nextflow run stenglein-lab/taxonomy_pipeline --fastq_dir /path/to/fastq/containing/directory/ --host_map_file /path/to/host_map_file.txt
```

The fastq files should contain Illumina reads, but can be single or paired-end, or a mix, or gzip compressed (end in .gz) or not, or a mix of compressed and uncompressed.  It would be a best practice to keep your fastq files gzipped to take up less space on server storage.  

Note that the pipeline looks for files with names that match the pattern `*_R[12]_*.fastq*`.  You can change this pattern using the argument `--fastq_pattern` as input to the nextflow command.  For instance:

```
nextflow run stenglein-lab/taxonomy_pipeline --fastq_dir /path/to/fastq/containing/directory/ --fastq_pattern "*_R1*fastq.gz"
```

The above pattern will match fastq files with R1 in their names, so would only use single-end data even if R2 files were present.

#### Step 2. Setup host filtering

The pipeline optionally removes host-derived reads because generally they are not what we are interested in and removing these reads makes downstream taxonomic classification steps go faster.  To perform host filtering, you will need a bowtie2 index of a host reference genome.  [This tutorial section](https://github.com/stenglein-lab/taxonomy_pipeline/blob/main/docs/tutorial.md#section_genome) describes how to download a new host genome and build a bowtie2 index if there is not one already available.

Host mapping is described in a file whose location is provided to the pipeline using the `--host_map_file` command line argument.  For instance:  

nextflow run stenglein-lab/taxonomy_pipeline --fastq_dir /path/to/fastq/containing/directory/ --host_map_file /path/to/host_map_file.txt

An [example of this file is here](map://github.com/stenglein-lab/taxonomy_pipeline/blob/nextflow/input/host_mapping.txt).  This file contains two columns separated by a tab.  The first column contains a pattern that will be searched for in fastq file names.  The second columns contains the location of a bowtie2 index.

```
# view an example host filtering map file
$ cat input/host_mapping.txt
_M_	/home/databases/fly/Dmel_genome_transcriptome
_F_	/home/databases/fly/Dmel_genome_transcriptome
Hela	/home/databases/human/GCRh38
HeLa	/home/databases/human/GCRh38
```

In this example, any datasets whose fastq file names contain the text `_M_` or `_F_` will be host-filtered using the bowtie2 index located at `/home/databases/fly/Dmel_genome_transcriptome` and any datasets whose fastq file names contain the text `Hela` or `HeLa` will be filtered using the GCRh38 human reference genome.  Note that bowtie2 indexes are actually made up of files and the paths listed in the example above are the *prefix* of all those files (they are the path that will be passed to the bowtie2 `-x` command line argument).

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

If you fail to specify a host mapping file the pipeline will warn you that no host filtering will be performed.  Note that you don't necessarily need to perform filtering: Any dataset whose fastq file name doesn't match one of the patterns in the host map file will not be host filtered.

#### Step 3. Actually run the pipeline

Now you actually need to run the pipeline.  An example command to run it:

```
# assuming you are in the pipeline main directory (2022_3_analyis in the example above)
nextflow run stenglein-lab/taxonomy_pipeline -resume -profile singularity,local --host_map_file input/host_mapping.txt
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

This pipeline deals with other software dependencies in 2 possible ways: by using Singularity containers or by using an all-in-one conda environment.  You are free to use either one of these when running the pipeline, but singularity is preferred over conda.

#### Singularity containers

The pipeline can use singularity containers to run programs like bowtie2 and blast.  To use these containers you must be running the pipeline on a computer that has [singularity](https://sylabs.io/singularity) [installed](https://sylabs.io/guides/latest/admin-guide/installation.html).  To test if singularity is installed, run:

```
singularity --version
```

To run with singularity containers include the option `-profile singularity` in the nextflow command line, for instance:

```
nextflow run stenglein-lab/taxonomy_pipeline -profile singularity
```

Singularity containers will be automatically downloaded and stored in a directory named `singularity_cacheDir` in your home directory.  They will only be downloaded once. 

#### Conda environment

The pipeline can also use an all-in-one conda environment.  This requires conda to be installed on your computer.  To test if conda is installed, run:

```
conda -V
```

The conda environment is defined in [this yaml file](./conda/taxonomy_conda_environment.yaml) and will be automatically created if you run the pipeline using the conda profile.  To run the pipeline with conda, include `-profile conda` in the nextflow command line, for instance:

```
nextflow run stenglein-lab/taxonomy_pipeline -profile conda
```

The conda environment will be created in a directory in your home directory named `conda_cacheDir` and will only be created once.

You must specify either `-profile conda` or `-profile singularity` or the pipeline will output an error message and halt.

#### Custom scripts

The pipeline also uses custom R and perl scripts, which are located in the [scripts](./scripts) directory of this repository.


### Sequence database dependencies

This pipeline uses two databases for taxonomic classification:

(1) The [NCBI nt nucleotide database](https://ftp.ncbi.nlm.nih.gov/blast/db/) for use in blastn searches.  [This script](./scripts/download_and_process_sequence_databases) will download the NCBI nt BLAST database (and create the diamond database described next).

(2) A [diamond](https://github.com/bbuchfink/diamond) database created from the NCBI nr protein sequence database.  This is downloaded and created as part of [this same script](./scripts/download_and_process_sequence_databases)

These databases need to be installed locally on a computer or server to run this pipeline.  This requires ~1.2 Tb of disk space (as of early 2022).  These databases take up a lot of space, so before doing this make sure that these databases are not already available locally.

#### Database locations

##### nt database

The default location of the blast nt databases is: `/home/databases/nr_nt/nt`. This value will be passed to the -db option when running blastn. This default value can be overridden using the `--nt_blast_db` parameter, for instance:

```
nextflow run stenglein-lab/taxonomy_pipeline -profile singularity --nt_blast_db /path/to/my/local/db/nt
```
In this example, there should be a file named `/path/to/my/local/db/nt.nal` in addition to the other nt database files.

##### diaomond nr database

The default location of the diamond database is: `/home/databases/nr_nt/nr.dmnd`.  This path will be passed to the --db option when running diamond.  This default value can be overridden using the `--nt_diamond_db` parameter.

##### Taxonomy db

The NCBI taxonomy databases are downloaded automatically as part of the pipeline and stored in the work directory created by nextflow.  These files use ~1.5 Gb of disk storage.

## Assumptions about input data

1. This pipeline is designed to work with Illumina data.  

## Testing

The pipeline is provided with small test fastq files that can be used to test that the pipeline is working properly.  To run these test datasets, run with `profile test`, for instance:

```
nextflow run stenglein-lab/taxonomy_pipeline -profile test,singularity
```

## Tutorial

There is [a tutorial](./docs/tutorial.md) that describes the main steps of this pipeline.  This tutorial refers to the older, pre-nextflow bash version of this pipeline, but the steps are similar and so the information in the tutorial remains applicable to understanding how the pipeline works in principle.


