## Stenglein lab taxonomy pipeline tutorial

This tutorial describes how to run the [Stenglein lab's](http://www.stengleinlab.org) [taxonomy pipeline](http://github.com/stenglein-lab/taxonomy_pipeline) to perform metagenomic classification.

This tutorial is tailored for users on the aidlngs01 server but could be useful to users who wish to use the pipeline on another server.   See [setup instructions](./setup_instructions.md). 


### Pipeline overview

The goal of the pipeline is to taxonomically classify sequences in shotgun sequencing datasets.  The overall strategy is similar to that used by the [SURPI pipeline](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4079973/) and the [ID-Seq](https://github.com/chanzuckerberg/idseq-web) platform.  (This is because all are more or less descendents of the [DeRisi lab](http://derisilab.ucsf.edu/)).  


### Logging in and getting started

First, you'll need to login to the `aidlngs01` server.  On a MacOS computer, open the Terminal app and enter:

```
ssh <your_csu_eid>@aidlngs01.cvmbs.colostate.edu
```

Enter your CSU EID password when prompted.

Note that if you are off the CSU campus, you will have [VPN in](https://www.acns.colostate.edu/security/#pulse-connect) to be able to login to the server.  If you are on a Windows computer, you can use [MobaXterm](https://mobaxterm.mobatek.net/download-home-edition.html).  

When you first login, you will be asked a question about an `"ECDSA key fingerprint"`.  Answer yes to continue.  You should only see this the first time you login.  You'll also see some information about a conda environment installing on your first login, which may run for a few minutes. 

If you're analyzing your own dataset, you'll have to transfer the sequence data (fastq files) to the server.  For this tutorial, we will analyze an existing dataset.  This dataset is shotgun 1x150 RNA sequencing data from wild caught [Drosophila spp.](http://obbard.bio.ed.ac.uk/photos.html) flies.  (Thanks to Reyes Murrieta for catching flies and Emily Fitzmeyer for preparing the sequencing library and generating the data).  

Let's start getting things into place.  You should be in your home directory on the aidlngs01 server, so your terminal should look something like this:

```
mdstengl@aidlngs01:~$ 
```

First, let's make a directory for our analysis and move (change) into that directory:
```
# make a directory
mkdir taxo_tutorial

# move (change) into it
cd taxo_tutorial
```

If you disconnect from the server and log back in, you'll want to change into that directory to continue to the tutorial.


Now, let's copy this example dataset to the directory we are in (our present working directory = .)
```
# copy a fastq file into the present directory
cp /data/training/taxonomy_pipeline/dros_pool_R1.fastq .
```

You should see this file if you type the command `ls -lh` and you should see that its size is about 1.4Gb.  

**Question:** How many reads are in this dataset? (Hint: the `wc -l` command outputs the # of lines in a file. Each fastq read uses 4 lines in typical [fastq format](https://en.wikipedia.org/wiki/FASTQ_format))


#### Getting the pipeline

The taxonomy pipeline is hosted [at github, here](https://github.com/stenglein-lab/taxonomy_pipeline).   If you visit that page, you will see the scripts and other files that make up the pipeline.  You can download them all at once using this command:

```
git clone https://github.com/stenglein-lab/taxonomy_pipeline.git
```

An advantage of using git clone to obtain the pipeline is that it effectively timestamps the version you are using.  You should see a new directory named taxonomy_pipeline when you type the command `ls`.

Run the following command to copy the pipeline scripts into your present working directory:
```
cp taxonomy_pipeline/bin/* .
```

You should see a bunch of scripts now in the present working directory that constitute the pipeline if you type `ls`.

#### To run the pipeline on this example dataset, you could run these commands:

#### But don't do this yet!
```
# activate the conda environment needed for the pipeline
conda activate taxonomy 

# run the entire pipeline on one dataset
# don't actually run this yet!
./run_pipeline_single_end dros_pool
```

The conda command will activate a conda environment that contains software like BLAST that the pipeline uses.  (The file that used to create this environment is [here](../server/taxo_recipe.yaml).

Note that the input to the pipeline is the name of the dataset, and the pipeline script expects a file whose name is `<dataset>_R1.fastq`.  For a paired-end analysis, the equivalent command would be `./run_pipeline <dataset_name>` and the pipeline would expect that two files exist: `<dataset_name>_R1.fastq` and `<dataset_name>_R2.fastq`

The pipeline will take a while to run, maybe 30 minutes for a dataset this size on this server.   

**Running the pipeline this way doesn't really help you understand what's happening, so let's consider the steps individually.**

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

The first thing the pipeline does is use [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) to do a preliminary QC on the input sequence data.  The output of the fastqc command will be an HTML file.

Note that you need to have the conda taxonomy environment activated.  If your command line prompt looks like: `(taxonomy) mdstengl@aidlngs01:~$ ` you are good to go.  If not, you need to run:

```
# activate the taxonomy conda environment
conda activate taxonomy
```

Running the following commands will recapitulate what the pipeline does:
```
# use FastQC to analyze the data in one fastq file
fastqc -o fastqc_pre dros_pool_R1.fastq
```

In the pipeline, this command is run as one command in the larger [run_pipeline](../bin/run_pipeline) bash script.  See if you can find where it's run in that script by looking at [the file](../bin/run_pipeline) on github. 

After fastqc completes, you should see a new html file if you run the `ls` command.

[Transfer this file](#section_transfer) to your own computer and open it in a web browser.  Here are some things to look for in the fastq report:

- **Basic Statistics:** this will indicate how many reads there are in the dataset, what their lengths are, etc.
- **Per base sequence quality:** an indication of the average basecall quality scores as a function of read length.  Hopefully much of your data will be in the green, but in any case, the pipeline trims off low quality parts off of reads.
- **Per base sequence content:** this the average percent of each of the 4 bases at each position in the reads.  These lines should be flat (base composition should be equal across reads in shotgun data) and the overall percentages should match the average GC content of the data.
- **Sequence duplication levels:**  this will give you an indication of how many duplicated sequences are in your dataset.  The more cycles of PCR you do during library prep, the more duplicated sequences you will have (which is why it's good to minimize PCR during library prep).  Duplicated sequences don't really add independent information, so the pipeline removes them.
- **Over-represented sequences:** sometimes these are highly duplicated sequences and sometimes they are adapter or rRNA sequences.  If you see over-represented sequences, you can blast them at the [NCBI website](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastSearch) to try to get an idea of what they are.
- **Adapter content:** if your read length is longer than the insert size in your libraries, you will see adapter sequences in your dataset.  You definitely want to trim off adapter sequences because they don't derive from the sample you are sequencing, and the pipeline does this.

### <a name="section2"></a> 2. Filtering of low quality and adapter sequences and duplicate collapsing

In this step, the pipeline uses [cutadapt](https://cutadapt.readthedocs.io/en/stable/) to remove low quality and adapter bases from reads.  It then uses [cd-hit](https://github.com/weizhongli/cdhit) to collapse duplicate reads.

These commands are part of the [run_preprocessing_pipeline_single_end script](../bin/run_preprocessing_pipeline_single_end), which is called from the main run_pipeline script.  For paired-end data, the script is [run_preprocessing_pipeline](../bin/run_preprocessing_pipeline).

You can reproduce these commands manually by running:

```
 # use cutadapt to trim low quality bases and adapter-derived bases
 cutadapt -a AGATCGGAAGAGC -g GCTCTTCCGATCT -a AGATGTGTATAAGAGACAG -g CTGTCTCTTATACACATCT -q 30,30 --minimum-length 80 -u 1 -o dros_pool_R1_f.fastq dros_pool_R1.fastq
```

Breaking down this command (see the excellent [cutadapt documentation](https://cutadapt.readthedocs.io/en/stable/) for more info):

- **-a AGATCGGAAGAGC, etc.:**  These are Illumina adapter sequences to be trimmed
- **-q30,30:** a quality score threshhold for trimming low quality bases 
- **--minimum-length 80:** any read that is shorter than this length after trimming will be discarded
- **-u 1:** always trim the last 3' base of the read, which is usually low quality. 
- **-o:** Name of the output file.  Here `<dataset_name>_R1_f.fastq`  The `_f` indicates the reads have now been trimmed.

Note that you are completely free to change these parameters for your own analyses (this comment applies to all steps of the pipeline).  Feel free to play around with them and re-run the command.

Cutadapt will run for a couple of minutes.  When it's done it will output statistics about how much it trimmed.  Have a look at that output.  Does it make sense?  Did it trim adapters?  Low quality bases?  Where any reads too short?  What fraction of the data was removed?

Now let's run [cd-hit-dup](https://github.com/weizhongli/cdhit/blob/master/doc/cdhit-user-guide.wiki)

```
# collapse duplicate reads with cd-hit
cd-hit-dup -i dros_pool_R1_f.fastq -o dros_pool_R1_fu.fastq -e 2 -u 50
```

The command syntax here means that cd-hit-dup will collapse any reads that have 2 or fewer differences in their first 50 bases, likely PCR duplicates.

You should now see 2 new files in your directory when you run `ls`: `dros_pool_R1_f.fastq` and `dros_pool_R1_fu.fastq`

How big are these files compared to the original pre-trimming file?  How many reads are in each? 

### <a name="section3"></a> 3. Post-filtering QC check

The pipeline next runs fastqc again to double-check that trimming and duplicate collapsing actually worked.  Run fastqc again like this:

```
# use FastQC to visualize effects of trimming
fastqc dros_pool_R1_fu.fastq -o fastqc_post
```

There should be a new html file in the new fastqc_post directory.  Transfer it to your computer and look at it.  Did trimming and collapsing steps have the desired impact?

### <a name="section4"></a> 4. Filtering of host-derived reads.

The next step in the pipeline is to remove host-derived reads. The host-derived reads might be interesting, but the pipeline's primary purpose is to taxonomically classify non-host reads.

The [run_pipeline_single_end script](../bin/run_pipeline_single_end) by default filters out Drosophila (fly) reads.  This is the appropriate type of host filtering for our example dataset but obviously wouldn't be the appropriate choice for most datasets.  See [below](#section_different_host) for information about how to filtering reads from other hosts.

The pipeline uses [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml) to map the reads to the [Drosophila melanogaster reference genome](https://www.ncbi.nlm.nih.gov/genome/?term=txid7227[Organism:noexp]).  Although these wild drosophilid flies are not necessarily _D. melanogaster_, they are closely enough related that using this reference geneome will work well (we determined this empirically).  

The main pipeline script calls [filter_fly_reads_single_end](../bin/filter_fly_reads_single_end) to accomplish this.  See if you can find both where the main pipeline calls this script and have a look at this host filtering script itself.

We can recreate that host filtering by running this command.
``` 
# map reads to the D. melanogaster genome to remove host reads
bowtie2 -x /home/databases/fly/fly_genome --local -q -U dros_pool_R1_fu.fastq --sensitive --score-min C,60,0 --time --un dros_pool_R1_fuh.fastq -p 12 > /dev/null 2> dros_pool_R1_fu.fastq.fly_genome_bt.log
```

This is a long command!  Let's break it down:

- **-x /home/databases/fly/fly_genome:** this specify the reference genome to which bowtie will map reads.  See [below](#section_different_host) for info on using other genomes.
- **--local:** this tells bowtie to not force the whole read to map.  
- **-q:** the input is fastq
- **-U dros_pool_R1_fu.fastq:** the name of the file containing (U)npaired reads that will be mapped.
- **--sensitive:** this tells bowtie to map with increased sensitivity.  This causes it to run slower than it could but will be more likely to map reads.
- **--score-min C,60,0:** this requires a minimum mapping score of 60: about 30 bp of perfect alignment (though mismatches are allowed).  
- **--un dros_pool_R1_fuh.fastq:** this tells bowtie to output unmapped reads (non-host reads) to this `_fuh.fastq` file.  **These are the reads we will analyze in subsequent steps.**
- **-p 12:** use 12 threads (~CPUs) to run faster.  
- **>/dev/null:** don't actually store the mapping information.  If you cared more about the host reads, you could change this to store the information. 

See the [bowtie2 manual](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml) for more info about all these settings.

The main output that we care about for this command will be the `dros_pool_R1_fuh.fastq` reads that didn't map to the fly genome. These are (presumably) the non-host reads that we will analyze in the rest of the pipeline.  How many reads are in this file? In other words, how many reads remain after host filtering?  What fraction of the original dataset is that?

This command will also output a log file: `dros_pool_R1_fu.fastq.fly_genome_bt.log`. Have a look at it. What percentage of reads mapped to the host genome?  


### <a name="section5"></a> 5. Assembly of remaining host-filtered reads.

The next 6 steps are all done by a single script, [contig_based_taxonomic_assessment](../bin/contig_based_taxonomic_assessment).  The commands run here are a bit more than it is convenient to run manually, so let's run the script as is, and talk about the various steps and their output.  

Let's run the command.  
```
# assemble non-host reads and taxonomically assign contigs
./contig_based_taxonomic_assessment dros_pool
```

This is going to be the slowest part of the pipeline and may take 20-30 min to complete (or longer).  

The first thing you should see is the output from the [spades assembler](http://cab.spbu.ru/software/spades/). Assemblers like Spades stitch together overlapping short reads into longer contigs.  This process should take a couple minutes for the number of reads remaining in the `_fuh.fastq` file.   Spades has verbose output so you will see a lot of information cascade down your screen.

After assembly is complete, you will see a new file in the directory named `dros_pool_spade_contigs.fa`. This is a fasta file containing the assembled contigs. 

The pipeline then uses bowtie to map host filtered reads (those in the `fuh.fastq` file) to the contigs.  The goals of this are:

1. To determine how many reads were collapsed into each contig.  When the pipeline taxonomically assigns the contigs, it weights them by the number of reads that contributed to each contig.
2. To identify reads that didn't assemble into contigs (assembly requires a certain amount of coverage (overlap), so not all reads will assemble).  These "singleton" reads can also be taxonomically assigned.

Once spade has completed, look for the `dros_pool_spade_contigs.fa` file, and have a look at it: `less dros_pool_spade_contigs.fa`.  How long are the longest contigs?  Can you tell how much coverage the different contigs have?

### <a name="section6"></a> 6. BLASTN search of contigs against the NCBI nucleotide (nt) database.

Once contigs have been created, it is time to try to taxonomically categorize them.  The pipeline first uses BLASTN to identify existing sequences in the NCBI nucleotide (nt) database that have a certain amount of similarity with each contig.  The output file created by the pipeline will be named `dros_pool_spade_contigs.fa.bn_nt`.  Note that this file will exist with a size of 0 bytes before it is populated with results (BLASTN creates the file immediately before it begins writing results to the files). Use the `ls -l` command to look for this file and note its file size.

Can you find the line in the [contig_based_taxonomic_assessment](../bin/contig_based_taxonomic_assessment) script where blastn is run?

After blastn completes, look at the blast output file (`dros_pool_spade_contigs.fa.bn_nt`) by running the command `less dros_pool_spade_contigs.fa.bn_nt`.  Can you interpret the output?   (Hint: [this page](http://www.metagenomics.wiki/tools/blast/blastn-output-format-6) describes the column in the blast output). Did the first contig (named `NODE_1_...`) map at a nucleotide level to a nt database sequence?  What is the NCBI accession of that sequence?  What is that sequence (you can check [here](https://www.ncbi.nlm.nih.gov/genbank/)).  



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
### <a name="section_tally"></a>Running custom tabulations
### <a name="section_distribute"></a>Getting sequences for various taxa
