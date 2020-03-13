## Stenglein lab taxonomic assessment pipeline

This is the pipeline used in the [Stenglein lab](http://www.stengleinlab.org) to taxonomically classify sequences in NGS datasets.

This pipeline has been reported in a number of [publications](https://www.stengleinlab.org/papers/) from our lab.

This pipeline has a number of dependencies and is not necessarily easily portable.  The dependencies include:

- [fastqc](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
- [bowtie2](hhttp://bowtie-bio.sourceforge.net/bowtie2/manual.shtml)
- [cutadapt](https://cutadapt.readthedocs.io/en/stable/)
- [blast](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
- [diamond](https://github.com/bbuchfink/diamond)
- [spades](http://cab.spbu.ru/software/spades/)

As well as a number of scripts that we've written, which can be found [here](https://github.com/stenglein-lab/stenglein_lab_scripts)

The pipeline also expects local installations of the NCBI nt/nr databases, as well as the NCBI taxonomy database for accession->taxid mapping.

#### Dependency installation 

We setup our linux bioinformatics servers on which this pipeline runs using [this script](https://github.com/stenglein-lab/stenglein_lab_scripts/blob/master/setup_server.sh)

We download and setup NCBI databases using these scripts:

- https://github.com/stenglein-lab/stenglein_lab_scripts/blob/master/fetch_NCBI_Taxonomy_db
- https://github.com/stenglein-lab/stenglein_lab_scripts/blob/master/fetch_nr_nt
- https://github.com/stenglein-lab/stenglein_lab_scripts/blob/master/build_indexes

The [run_pipeline_one_sample](./run_pipeline_one_sample) script is the main entry point to the pipeline.

The main steaps in the pipeline are:

- quality assessment w/ FastQC
- trimming or removal of low quality sequences and collapse of non-unique read pairs [run_preprocessing_pipeline_one_sample](./run_preprocessing_pipeline_one_sample)
- quality assessment w/ FastQC
- filtering of host sequences via the script specified by $host_filtering_script in [run_pipeline_one_sample](./run_pipeline_one_sample)
- de novo assembly of remaining reads and taxonomic assement of contigs, first by nt to nt alignments (blastn) and then by translated-nt to protein alignments (diamond)  [contig_based_taxonomic_assessment](./contig_based_taxonomic_assessment)
- taxonomic assement of non contig-mapping reads, first by nt to nt alignments (blastn) and then by translated-nt to protein alignments (diamond)  [contig_based_taxonomic_assessment](./contig_based_taxonomic_assessment)



We are working to make this pipeline more portable and welcome your feedback in the meantime.

This script can be run in parallel on multiple datasets using a tool like [gnu parallel](https://www.gnu.org/software/parallel/).  We usually run it via our own [simple scheduler](https://github.com/stenglein-lab/stenglein_lab_scripts/blob/master/simple_scheduler)
