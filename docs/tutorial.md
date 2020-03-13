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

Run the following command to copy the pipeline files into your present working directory:
```
cp taxonomy_pipeline/bin/* .
```

You should see a bunch of scripts now in the present working directory that constitute the pipeline if you type `ls`.

#### To run the pipeline on this example dataset, you can run these commands:
```
conda activate taxonomy 
./run_pipeline_single_end brazil_dros_pool
```
The conda command will activate a conda environment that contains software like BLAST that the pipeline uses.  (The file that used to create this environment is [here](../server/taxo_recipe.yaml).

Note that the input to the pipeline is the name of the dataset, and the pipeline script expects a file whose name is `<dataset>_R1.fastq`.  For a paired-end analysis, the equivalent command would be `./run_pipeline <dataset_name>` and the pipeline would expect that two files exist: `<dataset_name>_R1.fastq` and `<dataset_name>_R2.fastq`

The pipeline will take a while to run, maybe 30 minutes for a dataset this size on this server.   

Running the pipeline this way doesn't really help you understand what's happening, so let's consider the steps individually.

### The main pipeline steps are:

1. [Initial QC of reads.](#section1)
2. [Filtering of low quality and adapter sequences and duplicate collapsing.](#section2)
3. [Post-filtering QC check.](#section3)
4. [Filtering of host-derived reads.](#section4)
5. [Assembly of remaining host-filtered reads.](#section5)
6. [BLASTN search of contigs against the NCBI nucleotide (nt) database.](#section6)
7. [Taxonomic assignment of contigs based on nucleotide-level BLASTN alignments and tabulation of results.](#section7)
8. [Extraction of putative virus contigs.](#section8)
9. [Diamond (BLASTX) search of NCBI protein (nr) database.](#section9)
10. [Taxonomic assignment of contigs based on protein-level diamond alignments and tabulation of results.](#section10)
11. [Extraction of putative virus contigs.](#section11)
12. [Repeat of steps 7-12 for any reads that didn't assemble into contigs (singletons).](#section12)

Let's consider these steps one at a time:

### <a name="section1"></a> 1. Initial QC of reads. 

The first thing the pipeline does is use [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) to do a preliminary QC on the input sequence data.  The output of the fastqc command will be an HTML file in a new `fastqc_pre` directory.  

Running this command will recapitulate what the pipeline does:
```
mkdir -p fastqc_pre 
fastqc -o fastqc_pre brazil_dros_pool_R1.fastq
```

In the pipeline, this command is run as one command in the larger [run_pipeline](../bin/run_pipeline) bash script.  See if you can find where it's run in that script by looking at [the file](./bin/run_pipeline) on github. 

After fastqc completes, you should see a new html file if you run this command:
```
ls fastqc_pre
```

[Transfer this file](#section_transfer) to your own computer and open it in a web browser.  Here are some things to look for in the fastq report:

- Basic Statistics: this will indicate how many reads there are in the dataset, what their lengths are, etc.
- Per base sequence quality: an indication of the average basecall quality scores as a function of read length.  Hopefully much of your data will be in the green, but in any case, the pipeline trims off low quality parts off the ends of reads.
- Per base sequence content: this the average percent of each of the 4 bases at each position in the reads.  These lines should be flat (base composition should be equal across reads in shotgun data) and the overall percentages should match the average GC content of the data.
- Sequence duplication levels:  this will give you an indication of how many duplicated sequences are in your dataset.  The more cycles of PCR you do during library prep, the more duplicated sequences you will have (which is why it's good to minimize PCR during library prep).  Duplicated sequences don't really add independent information, so the pipeline removes them.
- Over-represented sequences: sometimes these are highly duplicated sequences and sometimes they are adapter or rRNA sequences.  If you see over-represented sequences, you can blast them at the [NCBI website](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastSearch) to try to get an idea of what they are.
- Adapter content: if your read length is longer than the insert size in your libraries, you will see adapter sequences in your dataset.  You definitely want to trim off adapter sequences because they derive from the sample you are sequencing, and the pipeline does this.

### <a name="section2"></a> 2. Filtering of low quality and adapter sequences and duplicate collapsing

In this step, the pipeline uses [cutadapt](https://cutadapt.readthedocs.io/en/stable/) to remove low quality and adapter bases from reads.  It then uses [cd-hit](https://github.com/weizhongli/cdhit) to collapse duplicate reads.

These commands are part of the [run_preprocessing_pipeline_single_end script](../bin/run_preprocessing_pipeline_single_end), which is called from the main run_pipeline script.  For paired-end data, the script is [run_preprocessing_pipeline](../bin/run_preprocessing_pipeline).

You can reproduce these commands manually by running:

```
 cutadapt -a AGATCGGAAGAGC -g GCTCTTCCGATCT -a AGATGTGTATAAGAGACAG -g CTGTCTCTTATACACATCT -q 30,30 --minimum-length 80 -u 1 -o brazil_dros_pool_R1_f.fastq brazil_dros_pool_R1.fastq
```

Breaking down this command (see the excellent [cutadapt documentation](https://cutadapt.readthedocs.io/en/stable/) for more info):

- -a AGATCGGAAGAGC, etc.:  These are Illumina adapter sequences to be trimmed
- q30,30: a quality score threshhold for trimming low quality bases 
- --minimum-length 80: any read that is shorter than this length after trimming will be discarded
- -u 1: always trim the last 3' base of the read, which is usually low quality . 
- -o: Name of the output file.  Here `<dataset_name>_R1_f.fastq`  The `_f` indicates the reads have now been trimmed.

Note that you are completely free to change these parameters.  Feel free to play around with them and re-run the command.

Cutadapt will run for a couple of minutes.  When it's done it will output statistics about how much it trimmed.  Have a look at that output.  Does it make sense?  Did it trim adapters?  Low quality bases?  Where any reads too short?

Now let's run [cd-hit-dup](https://github.com/weizhongli/cdhit/blob/master/doc/cdhit-user-guide.wiki)

```
cd-hit-dup -i brazil_dros_pool_R1_f.fastq -o brazil_dros_pool_R1_fu.fastq -e 2 -u 50
```

The command syntax here means that cd-hit-dup will collapse any reads that have 2 or fewer differences in their first 50 bases, likely PCR duplicates.

You should now see 2 new files in your directory when you run `ls`: `brazil_dros_pool_R1_f.fastq` and `brazil_dros_pool_R1_fu.fastq`

How big are these files compared to the original pre-trimming file?  How many reads are in each? 

### <a name="section3"></a> 3. Post-filtering QC check

The pipeline next runs fastqc again to double-check that trimming and duplicate collapsing actually worked.  Run fastqc again like this:

```
mkdir -p fastqc_post
fastqc brazil_dros_pool_R1_fu.fastq -o fastqc_post
```

There should be a new html file in the new fastqc_post directory.  Transfer it to your computer and look at it.  Did trimming and collapsing steps have the desired impact?

### <a name="section4"></a> 4. Filtering of host-derived reads.

The next step in the pipeline is to remove host-derived reads. The host-derived reads might be interesting, but the pipeline's primary purpose is to taxonomically classify non-host reads.

The [run_pipeline_single_end script](../bin/run_pipeline_single_end) is by default to filter out Drosophila (fly) reads.  This is the appropriate type of host filtering for our example dataset but obviously wouldn't always be the appropriate choice.  See [below](#section_different_host) for information about how to run this filtering for other hosts.

The pipeline uses [bowtie2]() to map the reads to the [Drosophila melanogaster reference genome]().  Although these flies are not necessarily D. melanogaster, they are closely enough related that using this reference geneome will work well.  

The main pipeline script calls [filter__fly_reads_single_end](../bin/filter_fly_reads_single_end) to accomplish this.  See if you can find both where the main pipeline calls this script and have a look at this host filtering script itself.

Let's recreate that host filtering by running this command.
``` 
```

### <a name="section5"></a> 5. Assembly of remaining host-filtered reads.
### <a name="section6"></a> 6. BLASTN search of contigs against the NCBI nucleotide (nt) database.
### <a name="section7"></a> 7. Taxonomic assignment of contigs based on nucleotide-level BLASTN alignments and tabulation of results.
### <a name="section8"></a> 8. Extraction of putative virus contigs
### <a name="section9"></a> 9. Diamond (BLASTX) search of NCBI protein (nr) database.
### <a name="section10"></a> 10. Taxonomic assignment of contigs based on protein-level diamond alignments and tabulation of results.
### <a name="section11"></a> 11. Extraction of putative virus contigs
### <a name="section12"></a> 12. Repeat of steps 7-12 for any reads that didn't assemble into contigs (singletons).





### <a name="section_different_host"></a>Filtering a different host species
### <a name="section_own_data"></a>Using your own data
### <a name="section_transfer"></a>Transferring files
### <a name="section_simple_scheduler"></a>Running the pipeline on multiple datasets
