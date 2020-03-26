# Stenglein lab taxonomy pipeline tutorial

This tutorial describes how to run the [Stenglein lab's](http://www.stengleinlab.org) [taxonomy pipeline](http://github.com/stenglein-lab/taxonomy_pipeline) to perform metagenomic classification.

This tutorial is tailored for users on the aidlngs01 server but could be useful to users who wish to use the pipeline on another server.   See [setup instructions](./server_setup_instructions.md) for how to get a computer setup to run the pipeline.  


## Pipeline overview

The goal of the pipeline is to taxonomically classify sequences in shotgun sequencing datasets.  The overall strategy is similar to that used by the [SURPI pipeline](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4079973/) and the [ID-Seq](https://github.com/chanzuckerberg/idseq-web) platform.  

The pipeline is very well suited for identifying virus and virus-like sequences in metagenomic datasets, but can be used for other purposes as well.

This pipeline has been reported in a number of [publications](https://www.stengleinlab.org/papers/) from our lab.


## Table of contents

- [Logging in and getting started on the aidlngs01 server](#section_login)
- [Setting up to run the tutorial](#section_setup)
- [Getting the pipeline](#section_get_pipeline)
- [Running the pipeline](#section_run_pipeline)
- [Main pipeline steps explained](#section_steps)
  - [Initial QC of reads.](#section1)
  - [Filtering of low quality and adapter sequences and duplicate collapsing.](#section2)
  - [Post-filtering QC check.](#section3)
  - [Filtering host reads.](#section4)
  - [Assembly of remaining host-filtered reads.](#section5)
  - [BLASTN search of contigs against the NCBI nucleotide (nt) database.](#section6)
  - [Taxonomic assignment of contigs based on nucleotide-level BLASTN alignments and tabulation of results.](#section7)
  - [Extraction of putative virus contigs.](#section8)
  - [Diamond (BLASTX) search of NCBI protein (nr) database.](#section9)
  - [Taxonomic assignment of contigs based on protein-level diamond alignments and tabulation of results.](#section10)
  - [Extraction of putative virus contigs.](#section11)
- Other topics
  - [Filtering a different host species](#section_different_host)
  - [Downloading a new host genome](#section_genome)
  - [Running the pipeline on multiple datasets](#section_simple_scheduler)
  - [Validating putative hits](#section_validation)
  - [Merging the results from multiple datasets](#section_matrix)
  - [Using the screen utility to avoid dropped connections](#section_screen)


### <a name="section_login"></a>Logging in and getting started on the aidlngs01 server

This section describes how to login to the aidlngs01 server, hosted at [CSU](http://www.colostate.edu).  This will not be that useful for folks not at CSU.

First, you'll need to login to the `aidlngs01` server.  On a MacOS computer, open the Terminal app and enter:

```
ssh <your_csu_eid>@aidlngs01.cvmbs.colostate.edu
```

Enter your CSU EID password when prompted.

Note that if you are off the CSU campus, you will have [VPN in](https://www.acns.colostate.edu/security/#pulse-connect) to be able to login to the server.  If you are on a Windows computer, you can use [MobaXterm](https://mobaxterm.mobatek.net/download-home-edition.html).  

When you first login, you will be asked a question about an `"ECDSA key fingerprint"`.  Answer yes to continue.  You should only see this the first time you login.  You'll also see some information about a conda environment installing on your first login, which may run for a few minutes. 

If you're analyzing your own dataset, you'll have to transfer the sequence data (fastq files) to the server.  For this tutorial, we will analyze an existing dataset.  This dataset is shotgun 1x150 RNA sequencing data from wild caught [Drosophila spp.](http://obbard.bio.ed.ac.uk/photos.html) flies.  (Thanks to Reyes Murrieta for catching the flies and Emily Fitzmeyer for preparing the sequencing library and generating the data).  

### <a name="section_setup"></a>Setting up to run the tutorial

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

If you disconnect from the server and log back in, you'll want to change into that directory to continue the tutorial.


Now, let's copy this example dataset to the directory we are in (our present working directory = .)
```
# copy a fastq file into the present directory
cp /data/training/taxonomy_pipeline/dros_pool_R1.fastq .
```

You should see this file if you type the command `ls -lh` and you should see that its size is about 1.4Gb.  

**Question:** How many reads are in this dataset? (Hint: the `wc -l` command outputs the # of lines in a file, so you can run `wc -l dros_pool_R1.fastq` to count the # of lines in that file. Each fastq read uses 4 lines in typical [fastq format](https://en.wikipedia.org/wiki/FASTQ_format))


#### <a name="section_get_pipeline"></a>Getting the pipeline

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

#### <a name="section_run_pipeline"></a> To run the pipeline on this example dataset, you could run these commands:

##### But don't do this yet!
```
# activate the conda environment needed for the pipeline
conda activate taxonomy 

# run the entire pipeline on one dataset
# ** don't actually run this yet! **
./run_pipeline_single_end dros_pool
```

The `conda activate` command will activate a conda environment that contains software like BLAST that the pipeline uses.  (The file that used to create this environment is [here](../server/taxo_recipe.yaml).

Note that the input to the pipeline is the name of the dataset, and the pipeline expects a file whose name is `<dataset>_R1.fastq` to exist in your present directory.  

For a paired-end analysis, the equivalent command would be `./run_pipeline <dataset_name>` and the pipeline would expect that two files exist: `<dataset_name>_R1.fastq` and `<dataset_name>_R2.fastq`

The pipeline would take a while to run, maybe 30 minutes for a dataset this size on this server.   

**Running the pipeline this way doesn't really help you understand what's happening, so let's consider the steps individually.**   We will run the commands individually rather than all at once via the main pipeline script.  

### <a name="section_steps"></a>The main pipeline steps are:

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

Every time you log out then log back in the server, remember to activate that taxonomy conda environment.  (You will also have to activate this environment again if you start a [screen session](#section_screen).

Running the following commands will recapitulate how the pipeline runs fastqc:
```
# use FastQC to analyze the data in one fastq file
fastqc dros_pool_R1.fastq
```

In the pipeline, this command is run from the main [run_pipeline_single_end](../bin/run_pipeline_single_end) bash script.  See if you can find where it's run in that script by looking at [the file](../bin/run_pipeline_single_end) on github. 

After fastqc completes, you should see a new html file if you run the `ls` command.

[Transfer this file](#section_transfer) to your own computer and open it in a web browser.  Here are some things to look for in the fastq report:

- **Basic Statistics:** this will indicate how many reads there are in the dataset, what their lengths are, etc.
- **Per base sequence quality:** an indication of the average basecall quality scores as a function of read length.  Hopefully much of your data will be in the green, but in any case, the pipeline trims off low quality parts off of reads.
- **Per base sequence content:** this the average percent of each of the 4 bases at each position in the reads.  These lines should be more or less flat (base composition should be equal across reads in shotgun data) and the overall percentages should match the average GC content of the data.
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
fastqc dros_pool_R1_fu.fastq 
```

There should be a new html file in your present directory.  Transfer it to your computer and look at it.  Did trimming and collapsing steps have the desired impact?

### <a name="section4"></a> 4. Filtering host reads.

The next step in the pipeline is to remove host-derived reads. The host-derived reads might well be interesting, but the pipeline's primary purpose is to taxonomically classify non-host reads.

The [run_pipeline_single_end script](../bin/run_pipeline_single_end) by default filters out _Drosophila_ (fly) reads.  This is the appropriate type of host filtering for our example dataset but obviously wouldn't be the appropriate choice for most datasets.  See [below](#section_different_host) for information about how to filtering reads from other hosts.

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

The main output that we care about for this command will be the `dros_pool_R1_fuh.fastq` file, which contains reads that _did not_ map to the fly genome. These are (presumably) the non-host reads that we will analyze in the rest of the pipeline.  

How many reads are in this file? In other words, how many reads remain after host filtering?  What fraction of the original dataset is that?  Could this set of non-mapping reads contain host-derived reads?  What are some of the reasons that host reads might not have been filtered by this step? 

This command will also output a log file: `dros_pool_R1_fu.fastq.fly_genome_bt.log`. Have a look at it. What percentage of reads mapped to the host genome?  


### <a name="section5"></a> 5. Assembly of remaining host-filtered reads.

The next 6 steps are all done by a single script, [contig_based_taxonomic_assessment_single_end](../bin/contig_based_taxonomic_assessment_single_end).  The commands run here are a bit more than it is convenient to run manually, so let's run the script as is, and walk through the steps and their output.  

Let's run the command.  
```
# assemble non-host reads and taxonomically assign contigs
./contig_based_taxonomic_assessment_single_end dros_pool
```
The first thing you should see is the output from the [spades assembler](http://cab.spbu.ru/software/spades/). Assemblers like Spades stitch together overlapping short reads into longer contigs.  This process should take a couple minutes for the number of reads remaining in the `_fuh.fastq` file.   Spades has verbose output so you will see a lot of information cascade down your screen.

Once spades has completed, look for the `dros_pool_contigs_singletons.fa` file, and have a look at it: `less dros_pool_contigs_singletons.fa`.  How long are the longest contigs?  Can you tell how much coverage the different contigs have?

The pipeline then uses bowtie to map host filtered reads (those in the `fuh.fastq` file) to the contigs.  The goals of this are:

1. To determine how many reads were collapsed into each contig.  When the pipeline taxonomically assigns the contigs, it weights them by the number of reads that contributed to each contig.
2. To identify reads that didn't assemble into contigs (assembly requires a certain amount of coverage (overlap), so not all reads will assemble).  These "singleton" reads can also be taxonomically assigned.  In fact, the pipeline concatenates the contigs with all the singtoln reads into a new merged file named `dros_pool_contigs_singletons.fa`.

The output of this step is a file named `dros_pool_contig_weights.txt`.  Have a look at the first 20 lines by running `head -20 dros_pool_contig_weights.txt`.  How many reads mapped to the top 20 contigs? 


### <a name="section6"></a> 6. BLASTN search of contigs and singletons against the NCBI nucleotide (nt) database.

Once contigs have been created, it is time to try to taxonomically categorize them and the non-assembling singleton reads (below, I refer to contigs, but really both singletons and contigs are being analyzed).  The pipeline first uses [BLASTN](https://en.wikipedia.org/wiki/BLAST_(biotechnology)) to identify existing sequences in the NCBI nucleotide (nt) database that exceen a certain similarity threshold for each contig.  The output file created by the pipeline will be named `dros_pool_contigs_singletons.fa.bn_nt`.  Note that this file will exist with a size of 0 bytes before it is populated with results (BLASTN creates the file immediately before it begins writing results to the files). Use the `ls -l` command to look for this file and note its file size.

Can you find the line in the [contig_based_taxonomic_assessment_single_end](../bin/contig_based_taxonomic_assessment_single_end) script where blastn is run?  What is the minimum E-value threshhold used? 

After blastn completes, look at the blast output file (`dros_pool_contigs_singletons.fa.bn_nt`) by running the command `less dros_pool_contigs_singletons.fa.bn_nt`.  Can you interpret the output?   (Hint: [this page](http://www.metagenomics.wiki/tools/blast/blastn-output-format-6) describes the column in the blast output). 

Did the first contig (named `NODE_1_...`) map at a nucleotide level to a nt database sequence?  What is the NCBI accession of that sequence?  What is that sequence? (You can check [here](https://www.ncbi.nlm.nih.gov/genbank/) by pasting in the accession).  

### <a name="section7"></a> 7. Taxonomic assignment of contigs based on nucleotide-level BLASTN alignments and tabulation of results.

You could go through the contigs one at a time, but that's not very practical.  The next step of the pipeline taxonomically categorizes the blast "hits" and tabulates the results in a couple different formats.  The script that does this is called [tally_blast_hits](../bin/tally_blast_hits).  Can you find where this is called in the [contig_based_taxonomic_assessment_single_end](../bin/contig_based_taxonomic_assessment_single_end) script? Note that it is called multiple times in slightly different ways.

The [tally_blast_hits](../bin/tally_blast_hits) script does a couple things.  First it goes through the blast results and maps the database sequences to their taxa.  Note that this depends on annotation in the NCBI database. For instance, visit [this database sequence](https://www.ncbi.nlm.nih.gov/nuccore/827047338).  What species is this sequence assigned to?  Can you see where the entire taxonomic lineage of this species is noted?  Note also that sometimes this annotation is incorrect.  

If a contig produces equally high scoring blast alignments to multiple database sequences, then it will assign this contig to the "lowest common ancestor" (LCA) of all the hits.  This means that sometimes contigs will be assigned at a higher taxonomic level (e.g. at the genus or family level).  

The [tally_blast_hits](../bin/tally_blast_hits) script also tabulates the hits and outputs tables of the number of hits to each taxa, the average percent identity of these hits, etc.  i

Have a look at the output of tally_blast_hits by running `less dros_pool_contigs_singletons.fa.bn_nt.tally`.  The output is sorted by the number of reads assigned to each taxon.  (Remember that this number is weighted by the number of reads that mapped to contigs).  

What was the most abunundant non-host taxon?   How many reads were assigned to it?  What was the median percent identity of the blast alignments for this taxon?    

This tally file is tab-delimited, so you can open it in programs like Excel or R for further analysis.  

The [tally_blast_hits](../bin/tally_blast_hits) script is also aware of the entire taxonomy heirarchy, so it can keep track of things like how many reads mapped to viruses, bacteria, etc.  The output file `dros_pool_contigs_singletons.fa.bn_nt.tab_tree_tally` has this information.  Look at this file by running `less dros_pool_contigs_singletons.fa.bn_nt.tab_tree_tally`.  

How many reads mapped to viruses?  How many to bacteria?  Note that the lines of this file are indented at the beginning to reflect the taxonomy heirarchy, meaning it's *not* suitable for being opened in Excel or R (these initial indents can be removed - see below).  

The [tally_blast_hits](../bin/tally_blast_hits) script is highly configurable.  It is run a couple different ways by the pipeline, but you can run it many different ways to suit your own analysis needs.  For instance, try running the command these way:

```
# run the script by itself to see all the options for running it
./tally_blast_hits 


# tally hits at the genus level
./tally_blast_hits -lca -r genus -w dros_pool_contig_weights.txt dros_pool_contigs_singletons.fa.bn_nt


# tally only for viruses (NCBI taxid 10239) 
# see: https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?mode=Info&id=10239&lvl=3&lin=f&keep=1&srchmode=1&unlock
./tally_blast_hits -lca -it 10239 -w dros_pool_contig_weights.txt dros_pool_contigs_singletons.fa.bn_nt


# same as .tab_tree_tally file, but with no indents at beginning of lines
./tally_blast_hits -lca -w dros_pool_contig_weights.txt -t dros_pool_contigs_singletons.fa.bn_nt > dros_pool_contigs_singletons.fa.bn_nt.tree_tally
```

Play around with this script to produce different types of output.  If you want to capture the output in a new file, you can use the bash [output redirection operator](https://thoughtbot.com/blog/input-output-redirection-in-the-shell) `>` as in the last example above.  

### <a name="section8"></a> 8. Extraction of putative virus contigs

The pipeline by default also outputs putative virus-mapping contigs and singleton reads into separate fasta files.  

The command responsible for doing this is named [distribute_fasta_by_blast_taxid](../bin/distribute_fasta_by_blast_taxid).  The pipeline uses this script to create one fasta file for each virus taxon identified.  

Look in your directory to identify these files.  For instance, you should see a file named `dros_pool_contigs_singletons.fa_1654579_Galbut_virus.fa`. Output the contents of this file by running `cat dros_pool_contigs_singletons.fa_1654579_Galbut_virus.fa`.  

How many galbut virus-mapping contigs were there?  How high are the coverage levels of these? 

Like [tally_blast_hits](../bin/tally_blast_hits), the [distribute_fasta_by_blast_taxid](../bin/distribute_fasta_by_blast_taxid) script is configurable.  You can run it multiple ways to get the sequences of contigs that were assigned at any taxonomic level.  Here are some examples:

```
# run the script by itself to see usage information
./distribute_fasta_by_blast_taxid

# output all bacteria-mapping (bacteria = taxid 2) reads into a single file:
./distribute_fasta_by_blast_taxid -t 2 dros_pool_contigs_singletons.fa.bn_nt
```

**Exercise:** you should see in your .tally file that the pipeline identified Wolbachia-mapping contigs in this example dataset.  Can you create a file containing all the Wolbachia-mapping contigs? 

### <a name="section9"></a> 9. Diamond (BLASTX) search of NCBI protein (nr) database.
### <a name="section10"></a> 10. Taxonomic assignment of contigs based on protein-level diamond alignments and tabulation of results.
### <a name="section11"></a> 11. Extraction of putative virus contigs

The next three steps are also run as part of [contig_based_taxonomic_assessment_single_end](../bin/contig_based_taxonomic_assessment_single_end).  

After taxonomically categorizing contigs by nucleotide-level similarity, the pipeline attempts to classify the remaining non-assigned contigs and singletons using protein-level similarity.  The pipeline uses the [diamond aligner](http://www.diamondsearch.org/index.php) to do this. Diamond is essentially equivalent to BLASTX, but runs faster.  Diamond identifies open reading frames in the contigs and unassembled reads, translates these in silico, then compares the resulting protein sequences to databases of protein sequences.  Comparisons at a protein level have the ability to identify more distant homologies (**Q:** why is this so?).  

The contigs and singletons that didn't produce a nucleotide-level similarity alignment are in the `dros_pool_contigs_singletons_n.fa`  

The [contig_based_taxonomic_assessment_single_end](../bin/contig_based_taxonomic_assessment_single_end) script runs diamond using this file as input.  Can you identify the line in that script where diamond is run?  

After running diamond the pipeline performs taxonomic assignment, result tabulation, and outputting of virus-mapping contigs and reads as above.  

Have a look at the tally file for the diamond assignments: `less dros_pool_contigs_singletons_n.fa.dmd_nr.tally`.  What taxa were identified this way?  How many reads did they account for?  What was the median % identity of the alignments?  

The pipeline also outputs putative viral contigs and singletons.  For example, you should see a file named `dros_pool_contigs_singletons_n.fa_33724_Nilaparvata_lugens_reovirus.fa`.   This file contains putative reovirus contigs that were identified by protein-level similarity to a reovirus from a brown planthopper insect (`Nilaparvata lugens`).  Note that since these were identified by protein level alignments, this means that these new sequences represent a 'new' virus: one for which closely related database sequences don't exist.  Cool, huh?

Note that some contigs or reads are not assigned either by nucleotide or protein level similarity.  These contigs are in a file named `dros_pool_contigs_singletons_nn.fa`.  Have a look at this file.  There may be interesting sequences in there that just weren't assigned but could be recognizable by other methods.


## Other considerations

### <a name="section_different_host"></a>Filtering a different host species

If you are going to be running this pipeline on your own datasets, it is likely that you will want to filter reads from a different host genome.  It's pretty straightforward to modify the pipeline to do this.

The host filtering in the examples above occurs via the script `filter_fly_reads_single_end`.  You need to do 2 things to change filtering to a different host:

1. Make a new host filtering script.
2. Modify the pipeline to use this new script.

#### 1. Making a new host filtering script

The easiest way to make a new host filtering script is to copy an existing one.  Say, for example, that you wanted to filter reads from `Aedes aegypti` mosquitoes.  First, you would copy the host filtering script:

```
# copy the host filtering script
cp filter_fly_reads_single_end filter_aedes_reads_single_end
```

Then, you would want to edit this script. You could use a text editor such as nano (`nano filter_aedes_reads_single_end`).  Or you could edit it remotely using a text editor like [BBedit](https://www.barebones.com/products/bbedit/).  

The lines that you will want to change in the original `filter_fly_reads_single_end` script are these lines:
```
# Location of bowtie2 index of genome
btindex=/home/databases/fly/fly_genome

# output suffix of files created
output_suffix=fly_genome
```

You would want to change them so they read as follows.  
```
# Location of bowtie2 index of genome
btindex=/home/mdstengl/databases/mosquito/aedes_aegypti

# output suffix of files created
output_suffix=aedes_genome
```

Note that in this example assumes the existence of a bowtie index in the directory `/home/mdstengl/databases/mosquito/`.  In real life, you would actually have needed to download the new host genome and build a bowtie2 index for it: [see below](#section_genome)

If there is no genome for your host of interest, often using a relatively closely related species will work sufficiently well for the purposes of host filtering.  

#### 2. Changing the main pipeline script to filter the appropriate host

Once you've created a new host filtering script, you'll need to change the main pipeline file to use it.  The relevant line in [run_pipeline_single_end](../bin/run_pipeline_single_end) is:
```
host_filtering_script=filter_fly_reads_single_end
```

In this example, you would want to change this by editing the file to:
```
host_filtering_script=filter_aedes_reads_single_end
```

Note that the comments in the [script](../bin/run_pipeline_single_end) describe a strategy to use different hosts for different datasets.

##### Host filtering paired-end data:

For paired end datasets, you'd want to change the [run_pipeline](../bin/run_pipeline) script.  [filter_human_reads](../bin/filter_human_reads) is an example of a paired-end host filtering script.

#### If you don't want to do host filtering

If you don't want to do host filtering, you can create a bash script like this that doesn't actually do any filtering between the `_fu.fastq` and `_fuh.fastq` files.
```
#!/bin/bash
id=$1
cp ${id}_R1_fu.fastq ${id}_R1_fuh.fastq
```

Create this script, name it something like `dummy_host_filtering` and be sure to give it executable permissions using the chmod command:
```
chmod +x dummy_host_filtering
```
Then modify the appropriate line in the `run_pipeline_single_end` script as above for a different host

### <a name="section_genome"></a>Downloading a new host genome

In the example [above](#section_different_host), the example assumed the existince of a bowtie index named `aedes_aegypti` in the directory `/home/mdstengl/databases/mosquito/`.  

These example commands would allow you to create a create such a database:
```
# change to your user's home directory
cd

# make a new directory to hold the genome and bowtie index
mkdir -p databases/mosquito/

# change to that new directory
cd databases/mosquito/

# download the fasta file containing the Aedes aegypti genome
curl -OL https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/204/515/GCF_002204515.2_AaegL5.0/GCF_002204515.2_AaegL5.0_genomic.fna.gz

# unzip the compressed fasta file
gunzip GCF_002204515.2_AaegL5.0_genomic.fna.gz

# make a bowtie2 index named aedes_aegypti from this fasta file
bowtie2-build GCF_002204515.2_AaegL5.0_genomic.fna aedes_aegypti
```

Note that an excellent way to see if there is a genome for your species of interest (or a related species) is via the NCBI taxonomy page for that species.  For instance, here is the page for the [Aedes aegypti mosquito](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?mode=Info&id=7159&lvl=3&lin=f&keep=1&srchmode=1&unlock).  On the table on the upper right side of that page, there is a link to the page in the NCBI Genome database for this species.  If you click that link, you will be taken to [this page](https://www.ncbi.nlm.nih.gov/genome/?term=txid7159[Organism:exp]), where there is an FTP link for the fasta file for the reference genome for the species.   (There are also things like the transcriptome, annotation, etc.)

Note also that the bowtie2-build step may take a while to run.  [See below](#section_screen) for an explanation of how to use the screen utility to avoid dropping a connection and interruping a long-running process like this.   


### <a name="section_simple_scheduler"></a>Running the pipeline on multiple datasets

It is useful to be able to run the pipeline on a number of datasets.  There are a number of ways to do this.  For instance, you could use the gnu [parallel command](https://www.gnu.org/software/parallel/).  We are working on implementing a [nextflow](https://www.nextflow.io/) version of the pipeline, which would also accomplish this. 

In the meantime, in the Stenglein lab, we often take advantage of utility called [simple_scheduler](https://github.com/stenglein-lab/stenglein_lab_scripts/blob/master/simple_scheduler) that accomplishes this.  

To run the pipeline in parallel, you could use simple scheduler as follows:
```
# run the pipeline in 4-way parallel on the datasets named in dataset_names.txt
simple_scheduler -a 4 ./run_pipeline_single_end `cat dataset_names.txt` 
```

In this example, you would create a text file named `dataset_names.txt` that contained the names of your datasets (one per line).  For example:
```
dros_pool_1
dros_pool_2
dros_pool_3
etc
```

**Note that it would not be advisable to run the pipeline in parallel on more than ~4-6 datasets at once.  Some of the steps use a fair amount of memory (RAM), and you want to avoid bogging down the server.**

### <a name="section_validation"></a>Validating putative hits

You should not blindly trust the results of the pipeline (or any bioinformatics analysis for that matter).  Here are some relatively simple things you can do to validate certain aspects of the results.

1. Is a taxonomic assignment correct?

If the pipeline indicates that a contig was assigned to a particular organism, one simple thing to do is to copy the contig sequence and blast it on the [NCBI website](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastSearch).  If the contig was assigned at a protein level, I will often do this both as a blastn search and a blastx search.

2. Is a contig assembled correctly?

Say you've pulled out a putative virus-derived contig. One of the first steps we will do to validate it is to re-map reads from the dataset back to the contig.  For example, say you wanted to validate the assembly of the galbut virus-derived contigs from the examples above.  You could run:
```
# make a new bowtie2 index from the putative galbut virus contigs
bowtie2-build dros_pool_contigs_singletons.fa_1654579_Galbut_virus.fa galbut_index

# use bowtie2 to map reads back:
bowtie2 -x galbut_index -q -U dros_pool_R1_fuh.fastq --local --score-min C,120,1 --no-unal --threads 12 -S dros_pool_R1_fuh.fastq.bt_galbut.sam
```

In this example, the sam file created by bowtie2 describes alignment of reads to the putative galbut virus contigs.  This sam file can be imported into a tool like Geneious to visualize the mapped reads and look for assembly errors.

### <a name="section_matrix"></a>Merging the results from multiple datasets

The pipeline contains a script named [make_taxa_matrix](https://github.com/stenglein-lab/taxonomy_pipeline/blob/master/bin/make_taxa_matrix) that can be used to combine the results from multiple datasets.  This script takes as input the .tally files described above and will output a matrix (table), where rows are datasets and columns are taxa.  For example:
```
# run the script by itself for usage information
make_taxa_matrix

# merge the results from nucleotide-level taxa tabulation for multiple datasets
make_taxa_matrix *.bn_nt.tally > taxa_matrix.txt
```
The resulting tab-delimited file `taxa_matrix.txt` can be imported into programs like Excel or R for further analysis.


### <a name="section_screen"></a> Using the screen utility to avoid dropped connections

If you are running a process that is going to take a long time to complete, it is useful to run it via the [gnu screen utility](https://www.gnu.org/software/screen/).  Screen allow you to create a persistent session on a server that will not terminate even if you disconnect from the server.  This avoids the artificial termination of a long-running process before it completes.  

Some useful screen command examples:

```
# start a screen session
screen 

# start a named screen session
screen -S taxonomy

# detach from a session
ctrl-a d  

# stop (terminate) a session
exit  

# list your screen sessions
screen -ls 

# re-attach to a detached session (if there is only one)
screen -r  

# re-attach to a particular detached session, of more than one exists 
screen -r session_id  

# for instance:
screen -r taxonomy
```

### <a name="section_transfer"></a> Transferring files to and from servers

You will sometimes need to transfer files to and from a server.  There are a number of ways to do this:

- There are some GUI-based tools that allow you to upload and download files.  Some examples include:
  - [Cyberduck](https://cyberduck.io/)
  - [Filezila](https://filezilla-project.org/)
- It may be possible to connect remotely to a server using a samba interface or something similar
- [sftp](https://www.digitalocean.com/community/tutorials/how-to-use-sftp-to-securely-transfer-files-with-a-remote-server) is a command line tool and protocol for transferring files between servers.  It will be worth your time to learn how to use sftp.  Especially for transferring large files (such as large .fastq files), sftp can be quite useful. 
