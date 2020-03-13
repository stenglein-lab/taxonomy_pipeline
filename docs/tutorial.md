## Stenglein lab taxonomy pipeline tutorial

This tutorial describes how to run the [Stenglein lab's](http://www.stengleinlab.org) [taxonomy pipeline](http://github.com/stenglein-lab/taxonomy_pipeline) to perform metagenomic classification.

This tutorial is tailored for users of the on the aidlngs01 server but could be useful to users who wish to use the pipeline on another server.   See [setup instructions](./setup_instructions.md). 


### Pipeline overview

The goal of the pipeline is to taxonomically classify sequences in shotgun sequencing datasets.  The overall strategy is similar to that used by the [SURPI pipeline](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4079973/) and the [ID-Seq](https://github.com/chanzuckerberg/idseq-web) platform.  (This is because all are derivatives of a generic strategy developed in the [DeRisi lab](http://derisilab.ucsf.edu/)).  


### Logging in and getting started

First, you'll need to login to the `aidlngs01` server.  On a MacOS computer, open the Terminal app and enter:

```
ssh <your_csu_eid>@aidlngs01.cvmbs.colostate.edu
```

Enter your CSU EID password when prompted.

Note that if you are off the CSU campus, you will have [VPN in](https://www.acns.colostate.edu/security/#pulse-connect) to be able to login to the server.  If you are on a Windows computer, you can use [MobaXterm](https://mobaxterm.mobatek.net/download-home-edition.html).  

When you first login, you will be asked a question about an `"ECDSA key fingerprint"`.  Answer yes to continue.  You should only see this the first time you login.  You'll also see some information about a conda environment installing on your first login, which may run for a few minutes. 

If you're analyzing your own dataset, you'll have to transfer the sequence data (fastq files) to the server.  For the purposes of this tutorial, we will analyze an existing dataset.  This dataset is shotgun 1x150 RNA sequencing data from wild caught [Drosophila spp.](http://obbard.bio.ed.ac.uk/photos.html) flies.  (Thanks to Reyes Murrieta for catching flies and Emily Fitzmeyer for preparing the sequencing library).  

Let's start getting things into place.  You should be in your home directory on the aidlngs01 server, so your terminal should look like this:

```
mdstengl@aidlngs01:~$ 
```

First, let's make a directory for our analysis and move (change) into that directory:
```
mkdir taxo_tutorial
cd taxo_tutorial
```

Now, let's copy this example dataset to the directory we are in (our present working directory = .)
```
cp /data/training/taxonomy_pipeline/brazil_dros_pool_R1.fastq .
```

You should see this file if you type the command `ls -lh` and you should see that its size is about 1.4Gb.  

**Question:** How many reads are in this dataset? (Hint: the `wc -l` command outputs the # of lines in a file and each read uses 4 lines in typical [fastq format](https://en.wikipedia.org/wiki/FASTQ_format))


#### Getting the pipeline

The taxonomy pipeline is hosted at github, [here](https://github.com/stenglein-lab/taxonomy_pipeline).   If you visit that page, you will see the scripts and other files that make up the pipeline.  You can download them all at once using this command:

```
git clone https://github.com/stenglein-lab/taxonomy_pipeline.git
```

An advantage of using git clone to obtain the pipeline is that it effectively timestamps the version of it you are using.  You should see a new directory named taxonomy_pipeline when you type the command `ls`.

Run the following command to get the pipeline files into your present working directory:

```
cp taxonomy_pipeline/bin/* .
```



###

The main default pipeline steps are:

1. [Initial QC of reads.](#section1)
2. [Filtering of low quality and adapter sequences.](#section2)
3. [Post-filtering QC check.](#section3)
4. [Collapsing of duplicate reads (likely PCR duplicates).](#section4)
5. [Filtering of host-derived reads.](#section5)
6. [Assembly of remaining host-filtered reads.](#section6)
7. [BLASTN search of contigs against the NCBI nucleotide (nt) database.](#section7)
8. [Taxonomic assignment of contigs based on nucleotide-level BLASTN alignments and tabulation of results.](#section8)
9. [Extraction of putative virus contigs.](#section9)
10. [Diamond (BLASTX) search of NCBI protein (nr) database.](#section10)
11. [Taxonomic assignment of contigs based on protein-level diamond alignments and tabulation of results.](#section11)
12. [Extraction of putative virus contigs.](#section12)
13. [Repeat of steps 7-12 for any reads that didn't assemble into contigs (singletons).](#section13)


Let's go through these step by step:

### <a name="section1"></a> 1. Initial QC of reads. 

The first thing the pipeline does is 

### <a name="section2"></a> 2. Filtering of low quality and adapter sequences
### <a name="section3"></a> 3. Post-filtering QC check
### <a name="section4"></a> 4. Collapsing of duplicate reads (likely PCR duplicates).
### <a name="section5"></a> 5. Filtering of host-derived reads.
### <a name="section6"></a> 6. Assembly of remaining host-filtered reads.
### <a name="section7"></a> 7. BLASTN search of contigs against the NCBI nucleotide (nt) database.
### <a name="section8"></a> 8. Taxonomic assignment of contigs based on nucleotide-level BLASTN alignments and tabulation of results.
### <a name="section9"></a> 9. Extraction of putative virus contigs
### <a name="section10"></a> 10. Diamond (BLASTX) search of NCBI protein (nr) database.
### <a name="section11"></a> 11. Taxonomic assignment of contigs based on protein-level diamond alignments and tabulation of results.
### <a name="section12"></a> 12. Extraction of putative virus contigs
### <a name="section13"></a> 13. Repeat of steps 7-12 for any reads that didn't assemble into contigs (singletons).





### Filtering a different host species
### Using your own data
