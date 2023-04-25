## Stenglein lab taxonomic assessment pipeline

This is a nextflow implementation of the pipeline used in the [Stenglein lab](http://www.stengleinlab.org) to taxonomically classify sequences in NGS datasets.

It is mainly designed to identify virus sequences in a metagenomic dataset but it also performs a general taxonomic classification of sequences.

A [previous bash-based version of this pipeline](https://github.com/stenglein-lab/taxonomy_pipeline/releases/tag/2022_3_30_bash) has been reported in a number of [publications](https://www.stengleinlab.org/papers/) from our lab. 

## How to run the pipeline

This pipeline is implemented in [nextflow](https://www.nextflow.io) and [requires nextflow](#Dependencies) to be run.

The pipeline requires several main inputs:

1. [The path to the directory containing your fastq files (`--fastq_dir`)](#FASTQ_directory)
2. [The path to a file defining how host filtering will be performed (`--host_map_file`)](#Host_filtering)
3. [Local copies of NCBI nucleotide and protein sequence databases](#Sequence_databases)
4. [Singularity or conda to provide required software tools (`-profile singularity` or `-profile conda`).](#Other_software_dependencies)


**Here's an example of a command line to run the pipeline:**

```
nextflow run stenglein-lab/taxonomy_pipeline -resume -profile singularity --fastq_dir /path/to/directory/containing/fastq/ --host_map_file /path/to/host_map_file.txt
```

Nextflow will [automatically download](https://www.nextflow.io/docs/latest/sharing.html) the [pipeline code](https://github.com/stenglein-lab/taxonomy_pipeline)  from github and run using the provided parameter values. `-resume` will [resume pipeline execution](https://www.nextflow.io/docs/latest/cli.html#run) if it stopped for some reason (or if you, for instance, added additional fastq files to the fastq-containing directory).  `-profile singularity` tells nextflow to use Singularity to [handle software dependencies](#Other_software_dependencies).  `--fastq_dir` specifies the location of a [directory containing input fastq](#FASTQ_directory).  `--host_map_file` specifies the location of a file that defines [how host reads](#Host_filtering) should be filtered.

 
#### FASTQ_directory

You must provide the location of a directory containing fastq files using the `--fastq_dir` command line argument to nextflow as in the above example.   This directory need not be a sub-directory of the directory from which you are running the pipeline.  

The fastq files should contain Illumina reads, but can be single or paired-end, or a mix, or gzip compressed (end in .gz) or not, or a mix of compressed and uncompressed.  It is best practice to keep your fastq files gzipped to take up less storage space.

Note that the pipeline looks for files with names that match the pattern `*_R[12]_*.fastq*`.  You can change this pattern using the argument **`--fastq_pattern`** as input to the nextflow command.  For instance:

```
nextflow run stenglein-lab/taxonomy_pipeline -resume -profile singularity --fastq_dir /path/to/directory/containing/fastq/ --host_map_file /path/to/host_map_file.txt --fastq_pattern "*R1*.fastq*"
```

In this example, the pattern will match fastq files with R1 in their names, so would only use single-end data even if R2 files were present.

#### Host_filtering

The pipeline optionally removes host-derived reads from datasets because generally host reads are not what we are interested in and removing these reads makes downstream taxonomic classification steps go **much** faster.  To perform host filtering, you will need one or more bowtie2 indexes of host reference genomes.  [This tutorial section](https://github.com/stenglein-lab/taxonomy_pipeline/blob/main/docs/tutorial.md#section_genome) describes how to download a new host genome and build a bowtie2 index.

The location of bowtie2 indexes are provided in a file specified by the `--host_map_file` command line argument as in the above example.

An [example of this host mapping file is here](map://github.com/stenglein-lab/taxonomy_pipeline/blob/nextflow/input/host_mapping.txt).  This file must contain two columns separated by a tab.  The first column contains a pattern that will be searched for in fastq file names.  The second columns contains the location of a bowtie2 index.

```
# an example host filtering map file
$ cat host_mapping.txt
_M_	/home/databases/fly/Dmel_genome_transcriptome
_F_	/home/databases/fly/Dmel_genome_transcriptome
Hela	/home/databases/human/GCRh38
HeLa	/home/databases/human/GCRh38
```

In this example, any datasets whose fastq file names contain the text `_M_` or `_F_` will be host-filtered using the bowtie2 index located at `/home/databases/fly/Dmel_genome_transcriptome` and any datasets whose fastq file names contain the text `Hela` or `HeLa` will be filtered using the GCRh38 human reference genome.  Note that bowtie2 indexes are actually made up of 6 files and the paths listed in the example above are the prefix of all those files (this prefix is passed to the bowtie2 `-x` command line argument).

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

If you fail to specify a host mapping file the pipeline will warn you that no host filtering will be performed.  Note that you don't necessarily need to perform filtering: Any dataset whose fastq file name doesn't match one of the patterns in the host map file will not be host filtered.  But beware that not host filtering when you can will cause the pipeline to run much more slowly than if you removed host reads.

#### Running with screen

Because the pipeline will likely run for a while (depending on dataset sizes), you will want to run it via a screen session.  See [this tutorial section](https://github.com/stenglein-lab/taxonomy_pipeline/blob/main/docs/tutorial.md#section_screen) for more information on using screen.

## Pipeline output 

The pipeline puts output in a `results` directory.

- Virus sequences: in the `results/virus_sequences` directory
- Taxonomic tabulation: in the `tallies` directory.  These are tab-delimited tables of the taxa observed in each dataset, the # of mapping reads, average percent identity to database sequences, etc.
- Host-filtered fastq: in the `host_filtered_fastq` directory.  This directory contains reads after quality and (optional) host-filtering.
- Contigs: in the `contigs` directory.  Contigs from assembly of reads remaining after host filtering.
- Remapping of reads to putative virus sequences: in the `virus_remapping` directory.  The pipeline automatically remaps host-filtered reads to putative virus contigs.  SAM files from this remapping are in this directory as well as files containing mapping stats.
- QC reports: `initial_qc_report.html` and `post_trim_qc_report.html`
- Plots of # and fraction reads remaining after filtering steps: `filtering_plots.pdf`

## Dependencies

### Nextflow

This pipeline is implemented in Nextflow.  To run the pipeline you will need to be working on a computer that has nextflow installed. Installation instructions can be [found here](https://www.nextflow.io/docs/latest/getstarted.html#installation).  To test if you have nextflow installed, run:

```
nextflow -version
```

**This pipeline currently requires Nextflow version >= 22.10.4 and < 23.**  Nextflow version 23 will not work because it [requires DSL2 syntax](https://www.nextflow.io/blog/2022/evolution-of-nextflow-runtime.html), and this pipeline is still in DSL1.  A [future enhancement](https://github.com/stenglein-lab/taxonomy_pipeline/issues/2) will convert it to DSL2.  

you can download Nextflow version 22.10.8 from [this link](https://github.com/nextflow-io/nextflow/releases/download/v22.10.8/nextflow)

#### Updating cached versions of the pipeline

Nextflow will download and cache the pipeline code in a subdirectory of your home directory: 

```
# See the downloaded pipeline code:
$ ls ~/.nextflow/assets/stenglein-lab/taxonomy_pipeline/
bin/    conf/        docs/   LICENSE  make_clean*      README.md      run_test*  server/
conda/  containers/  input/  main.nf  nextflow.config  run_pipeline*  scripts/   test/
```

If the pipeline has been updated and you want to ensure that you are running the latest code you can run:

```
nextflow drop stenglein-lab/taxonomy_pipeline
```

To [force nextflow](https://www.nextflow.io/docs/latest/sharing.html#deleting-a-downloaded-project) to delete the downloaded pipeline code and re-download the latest version.  You can also run specific versions of the pipeline, as [described here in more detail](https://www.nextflow.io/docs/latest/sharing.html#handling-revisions).

### Other_software_dependencies

This pipeline deals with other software dependencies in 2 possible ways: by using Singularity containers or by using an all-in-one conda environment.  You are free to use either one of these when running the pipeline, but singularity is preferred over conda.

#### Singularity containers

The pipeline can use singularity containers to run programs like bowtie2 and blast.  To use these containers you must be running the pipeline on a computer that has [singularity](https://sylabs.io/singularity) [installed](https://sylabs.io/guides/latest/admin-guide/installation.html).  To test if singularity is installed, run:

```
singularity --version
```

There is no prescribed minimum version of Singularity, but older version (>1-2 years old) are likely to have problems.

To run with singularity containers include the option `-profile singularity` in the nextflow command line, for instance:

```
nextflow run stenglein-lab/taxonomy_pipeline -profile singularity ...
```

Singularity containers will be automatically downloaded and stored in a directory named `singularity_cacheDir` in your home directory.  They will only be downloaded once. 

#### Conda environment

The pipeline can also use an all-in-one conda environment.  This requires conda to be installed on your computer.  To test if conda is installed, run:

```
conda -V
```

The conda environment is defined in [this yaml file](./conda/taxonomy_conda_environment.yaml) and will be automatically created if you run the pipeline using the conda profile.  To run the pipeline with conda, include `-profile conda` in the nextflow command line, for instance:

```
nextflow run stenglein-lab/taxonomy_pipeline -profile conda ...
```

The conda environment will be created in a directory in your home directory named `conda_cacheDir` and will only be created once.

You must specify either `-profile conda` or `-profile singularity` or the pipeline will output an error message and halt.

#### Custom scripts

The pipeline also uses custom R and perl scripts, which are located in the [scripts](./scripts) directory of this repository.

### Sequence_databases

This pipeline uses two databases for taxonomic classification.  These must exist locally.

(1) The [NCBI nt nucleotide database](https://ftp.ncbi.nlm.nih.gov/blast/db/) for use in BLASTN searches.  

(2) A [diamond](https://github.com/bbuchfink/diamond) database created from the NCBI nr protein sequence database.  

[This script](./scripts/download_and_process_sequence_databases) will download both databases.  

These databases must be present locally on a computer or server to run this pipeline.  This requires ~1.2 Tb of disk space (as of early 2022).  These databases take up a lot of space, so before doing this make sure that these databases are not already available locally.

#### Database locations

##### nt database

The default location of the blast nt databases is: `/home/databases/nr_nt/` and the default name of the database is `nt`. These values will be used to set the `-db` option when running blastn. These default value can be overridden using the `local_nt_database_dir` and `local_nt_database_name` parameters, for instance:

```
nextflow run stenglein-lab/taxonomy_pipeline -resume -profile singularity --fastq_dir /path/to/directory/containing/fastq/ --host_map_file /path/to/host_map_file.txt --local_nt_database_dir /path/to/my/local/db/
```
In this example, there should be a file named `/path/to/my/local/db/nt.nal` in addition to the other nt database files.

##### diaomond nr database

The default location of the diamond database is: `/home/databases/nr_nt/nr.dmnd`.  This path will be passed to the `--db` option when running diamond.  This path is specififed by the `local_diamond_database_dir` and `local_diamond_database_name` parameters, which can be overridden as described above for the nt BLASTN database. 

##### Taxonomy db

The NCBI taxonomy databases are downloaded automatically as part of the pipeline and stored in the work directory created by nextflow.  These files use ~1.5 Gb of disk storage.

## Assumptions about input data

1. This pipeline is designed to work with Illumina data.  
2. Input fastq can be single-end or paired-end, or a mix.
3. Input fastq can be compressed (.gz) or not, but you might as well keep your input compressed.

## Testing

The pipeline is provided with small test fastq files that can be used to test that the pipeline is working properly.  To run these test datasets, run with `profile test`, for instance:

```
nextflow run stenglein-lab/taxonomy_pipeline -profile test,singularity
```

## Pipeline options

A number of pipeline defaults are specified in [nextflow.config](./nextflow.config). These options can be [overridden in a variety of ways](https://www.nextflow.io/docs/latest/config.html).  Overridable options include:

- `outdir`: the main output directory. Default: results.
- `fastq_pattern`: the pattern to match for finding input fastq files.  Default: "*_R[12]*.fastq*"
- Output sub-directories: there a number of these, for instance `contigs_out_dir` (default: $outdir/contigs), and they can all be overridden.
- `minimum_contig_length`: the minimum contig length: contigs shorter than this will be discarded.  Default: 200 bases.
- `classify_singletons`: should non-assembling reads (singletons) be taxonomically classified?  Default: no.
- Resource allocation:  it's possible to specify the maximium available number of processers and memory. A good way to use this is to use a [custom profile](https://www.nextflow.io/docs/latest/config.html#config-profiles).  For example, see profile named local in [nextflow.config](./nextflow.config). 

## Tutorial

There is [a tutorial](./docs/tutorial.md) that describes the main steps of this pipeline.  This tutorial refers to the older, pre-nextflow [bash version](https://github.com/stenglein-lab/taxonomy_pipeline/releases/tag/2022_3_30_bash) of this pipeline, but the steps are similar and so the information in the tutorial remains applicable to understanding how the pipeline works in principle.


