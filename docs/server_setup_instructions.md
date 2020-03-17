## Server setup instructions for the Stenglein lab taxonomic assessment pipeline

This is the pipeline used in the [Stenglein lab](http://www.stengleinlab.org) to taxonomically classify sequences in NGS datasets.

The goal of this document is to describe how to setup a computer to run this pipeline.  


### Conda installation

The dependencies for this pipeline (see [below](#section_dependencies)) can be installed via [conda](https://docs.conda.io/en/latest/).  The [conda recipe file for 64 bit linux is here](https://github.com/stenglein-lab/taxonomy_pipeline/blob/master/server/taxo_recipe.yml).  After downlading this recipe file, create a conda environment with this command:
```
conda create --name taxonomy --file taxo_recipe.yaml
```

#### <a name="section_dependencies"></a> Dependencies.

This pipeline has several main dependencies (included in the conda environment described above), including: 

- [fastqc](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
- [bowtie2](hhttp://bowtie-bio.sourceforge.net/bowtie2/manual.shtml)
- [cutadapt](https://cutadapt.readthedocs.io/en/stable/)
- [blast](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
- [diamond](https://github.com/bbuchfink/diamond)
- [spades](http://cab.spbu.ru/software/spades/)

As well as some perl modules and a the scripts in this repository.

The pipeline also expects local installations of the NCBI nt/nr databases, as well as the NCBI taxonomy database for accession->taxid mapping.

### Installing NCBI databases

The pipeline needs databases of nucleotide (nt) and protein (nr) sequences from [NCBI](https://www.ncbi.nlm.nih.gov/).  These databases are pretty big.  Expect that they will take up something like ~1Tb of disk space (and they continue to grow) . 

We download and setup NCBI databases using the scripts in this repository in the server directory.  The main scripts to run are:

#### Download and process the NCBI taxonomy database
- https://github.com/stenglein-lab/taxonomy_pipeline/blob/master/server/fetch_NCBI_Taxonomy_db:

#### Download and process the NCBI nt and nr databases
- https://github.com/stenglein-lab/taxonomy_pipeline/blob/master/server/fetch_nr_nt
