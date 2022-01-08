#!/usr/bin/env nextflow

/*
    Stenglein lab metagenomic classification pipeline 
    implemented in nextflow

    March 9, 2021

    Mark Stenglein
*/


// TODO: command line options and parameter checking

// output usage message
params.help = false
params.h = false

params.input_dir = "$baseDir/input/"
params.fastq_dir = "${params.input_dir}/fastq/"
params.outdir = "$baseDir/results"                                                       

params.initial_fastqc_dir = "${params.outdir}/initial_fastqc/" 
params.post_trim_fastqc_dir = "${params.outdir}/post_trim_fastqc/" 
params.host_filtered_out_dir = "${params.outdir}/host_filtered_fastq/"                        
params.contigs_out_dir = "${params.outdir}/contigs/"                        
params.blast_out_dir = "${params.outdir}/blast_output/"                        
params.tally_out_dir = "${params.outdir}/tallies/"                        
params.virus_seq_out_dir = "${params.outdir}/virus_sequences/"                        
params.counts_out_dir = "${params.outdir}/fastq_counts/"                        
params.fastq_out_dir = "${params.outdir}/trimmed_fastq/"                        
params.bam_out_dir = "${params.outdir}/bam/"                                    

// ------------------
// Trimming settings
// ------------------
params.always_trim_5p_bases = "0" 
params.always_trim_3p_bases = "1" 
params.post_trim_min_length = "30" 

// --------------------
// Host cell filtering
// --------------------
// Define one of the 2 following parameters:
// 
// 1. A 2-column tab-delimited file with:
//    - the first column defining dataset IDs or patterns that will
//      match dataset IDs
//    - the second column will be the path of a bowtie index that will be
//      used to filter out host reads
//  
//    This enables different filtering for different datasets

// params.host_bt_index_map_file = "${params.input_dir}/host_index_map.txt"
params.host_bt_index_map_file = ""

// 2. The path to a bowtie index that will be used to filter host reads
//    for all datasets
// 
// params.host_bt_index = "/home/databases/fly/combined_fly_index"
params.host_bt_index = ""

// min bowtie alignment score to be considered a host-derived read
params.host_bt_min_score = "60"

// where are R scripts found...
params.R_bindir="${baseDir}/scripts"
params.scripts_bindir="${baseDir}/scripts"

// cd-hit-dup cutoff for collapsing reads with >= this much fractional similarity
params.mismatches_allowed = "2"


// Blast e-value cutoffs                                                        
params.max_blast_nt_evalue = "1e-10"  
params.max_blasx_nr_evalue = "1e-3"  

params.blast_db_dir = "/home/databases/nr_nt/"
params.nt_blast_db = "${params.blast_db_dir}/nt"
params.nr_blast_db = "${params.blast_db_dir}/nr"
params.nr_diamond_db = "${params.blast_db_dir}/nr.dmnd"

// TODO: command line arg processing and validating 

// TODO: check that appropriate refseq files exist (host bt index, nr/nt, etc.)

// TODO: handle fastq.gz compressed files 

/*
 These fastq files represent the main input to this workflow
*/
Channel
  .fromFilePairs("${params.fastq_dir}/*_R{1,2}*.fastq*", size: -1, checkIfExists: true, maxDepth: 1)
  .into {samples_ch_qc; samples_ch_trim; samples_ch_count; host_setup_ch}


// define usage output message
// TODO:
def helpMessage() {
  log.info """
        Usage:
        The typical command for running the pipeline is as follows:
        nextflow run main.nf --query QUERY.fasta --dbDir "blastDatabaseDirectory" --dbName "blastPrefixName"

        Mandatory arguments:
        One of: 
        -- host_bt_index_map_file      The path to a file that contains host filtering
                                       info for all datasets or subsets of datasets.   

                                       This enables different filtering for different datasets
          
                                       This file should be a 2-column tab-delimited file 
                                       with:

                                       - the first column defining dataset IDs or patterns that will
                                         match dataset IDs

                                       - the second column will be the path of a bowtie index that will be
                                         used to filter out host reads
                                       

        -- host bt_index               The path to a bowtie index that will be used 
                                       to filter host reads for all datasets

        Optional arguments:
        --outdir                       Output directory into which to place final results files 

        --help                         This usage statement.
        --h                            This usage statement.
        """
}

if (params.help || params.h) {
    helpMessage()
    exit 0
}


/* 
  Check input parameters 
*/
def check_params () {

  // must specify one and only one of these 2 host mapping 
  if (!params.host_bt_index_map_file && !params.host_bt_index){
    log.info """
      Error: you may specify one of these two parameters:
        1. host_bt_index_map_file ($params.host_bt_index_map_file) and 
        2. host_bt_index ($params.host_bt_index) 
    """
    helpMessage()
  }

  if (params.host_bt_index_map_file && params.host_bt_index){
    log.info """
      Error: you may specify only one of these two parameters:
        1. host_bt_index_map_file ($params.host_bt_index_map_file) and 
        2. host_bt_index ($params.host_bt_index) 
    """
    helpMessage()
  }

  if (params.host_bt_index_map_file) {
    host_map_file = file(params.host_bt_index_map_file)
    if (!host_map_file.exists()) {  
      log.info  """
        Error: host_bt_index_map_file ($params.host_bt_index_map_file) does not exist 
      """
      helpMessage() 
    }
  }
}

host_map = [:]

/* 
  This function defines bowtie indexes to be used to do host filtering

  This can be one index for all datasets or different indexes for 
  different datasets.
*/
def setup_host_indexes () {
   
  // Store the bowtie index to be used for each dataset ID or dataset ID pattern
  if (params.host_bt_index_map_file) {

    host_map_file = file(params.host_bt_index_map_file)

    // read through the file
    allLines  = host_map_file.readLines()
    for( line : allLines ) {
        // split by tabs
        fields = line.tokenize()
        host_map.put(fields[0], fields[1])
    }
  }
  else {
    // because of the way pattern matching works below
    // this empty string will match all sample IDS
    // and host will be set to params.host_bt_index
    host_map.put("" , params.host_bt_index)
  }
}

/*
  Get the bowtie index for host filtering for one dataset
*/
def get_dataset_host_index (sample_id) {

  def index = ""

  // first, check for a direct match
  if (host_map.containsKey(sample_id)) {
       index = host_map.get(sample_id)
  }

  // if that doesn't work, check for a pattern match
  // using groovy's find operator to do pattern matching
  // see: https://groovy-lang.org/operators.html#_find_operator
  if (!index) {
    host_map.each {
      key, value -> 
      // Pattern matching 
      // TODO: if pattern ($key) is empty, this will evaluate to true
      //       is this the expected behavior?
      if (sample_id =~ /$key/) {
         index = value
      }
    }
  }

  // no index found: no host filtering will be performed 
  if (!index) {
    log.info "Info: no host-filtering index found for $sample_id.  No host filtering will be performed\n"
    index = "dummy"
  }

  return (index)
}

// check parameters
check_params()

// setup indexes for host filtering
setup_host_indexes()


/*
   Setup some initial indexes and dictionaries needed by downstream processes.
   Only do this once at beginning.
*/
process setup_indexes {

  output:
  val("indexes_complete") into post_index_setup_ch

  script:
  """

  """
}

/*
 Run fastqc on input fastq 
*/
process initial_qc {
  label 'lowmem_non_threaded'                                                                

  input:
  tuple val(sample_id), path(initial_fastq) from samples_ch_qc

  output:
  val(sample_id) into post_initial_qc_ch
  // TODO: count

  script:
  """
  mkdir -p  ${params.initial_fastqc_dir} 
  fastqc -o ${params.initial_fastqc_dir} $initial_fastq 
  """
}

/*
 Count initial fastq
*/
process initial_fastq_count {
  label 'lowmem_non_threaded'                                                                
  publishDir "${params.counts_out_dir}", mode:'link'

  input:
  tuple val(sample_id), path(initial_fastq) from samples_ch_count

  output:
  path("${sample_id}_initial_count.txt") into post_count_initial_ch

  shell:
  // only count the first file because for paired-read data both files
  // will have the same # of reads.

  // for an explanation of the xargs command used for arithmetic in a pipe, see: 
  // https://askubuntu.com/questions/1203063/how-can-i-pipe-the-result-of-the-wc-command-into-an-arithmetic-expansion
  '''
  zcat -f !{initial_fastq[0]} | wc -l | xargs bash -c 'echo $(($0 / 4))' | awk '{print "!{sample_id}" "\tinitial\t" $1}' > "!{sample_id}_initial_count.txt"
  '''
}

/*
 Use multiqc to merge initial fastqc reports
*/
process initial_multiqc {
  publishDir "${params.outdir}", mode:'link'

  input:
  val(all_sample_ids) from post_initial_qc_ch.collect()

  output: 
  path("initial_qc_report.html")

  script:
  """
  multiqc -n "initial_qc_report.html" -m fastqc ${params.initial_fastqc_dir}
  """
}

/*
 Use cutadapt to trim off adapters and low quality bases
*/
process trim_adapters_and_low_quality {
  label 'lowmem_non_threaded'                                                                

  input:
  tuple val(sample_id), path(initial_fastq) from samples_ch_trim

  output:
  tuple val(sample_id), path("*_f.fastq") optional true into post_trim_ch
  tuple val(sample_id), path("*_f.fastq") optional true into post_trim_count_ch

  // TODO: put adapters as a param (?)
  script:

  // this handles paired-end data, in which case must specify a paired output file
  def paired_output   = initial_fastq[1] ? "-p ${sample_id}_R2_f.fastq" : ""
  def paired_adapters = initial_fastq[1] ? "-A AGATCGGAAGAGC -G GCTCTTCCGATCT -A AGATGTGTATAAGAGACAG -G CTGTCTCTTATACACATCT" : ""
  // TODO: don't trim this much for non-amplicon data!
  def paired_trimming = initial_fastq[1] ? "-U $params.always_trim_5p_bases -U -${params.always_trim_3p_bases}" : ""

  """
  cutadapt \
   -a AGATCGGAAGAGC -g GCTCTTCCGATCT -a AGATGTGTATAAGAGACAG -g CTGTCTCTTATACACATCT \
   $paired_adapters \
   -q 30,30 \
   --minimum-length ${params.post_trim_min_length} \
   -u ${params.always_trim_5p_bases} \
   -u -${params.always_trim_3p_bases} \
   $paired_trimming \
   -o ${sample_id}_R1_f.fastq \
   $paired_output \
   $initial_fastq 
  """

}

/*
 Count post-trimming fastq
*/
process trimmed_fastq_count {
  label 'lowmem_non_threaded'                                                                
  publishDir "${params.counts_out_dir}", mode:'link'

  input:
  tuple val(sample_id), path(trimmed_fastq) from post_trim_count_ch

  output:
  path("${sample_id}_trimmed_count.txt") into post_count_trim_ch

  shell:
  // only count the first file because for paired-read data both files
  // will have the same # of reads.
  '''
  zcat -f !{trimmed_fastq[0]} | wc -l | xargs bash -c 'echo $(($0 / 4))' |  awk '{print "!{sample_id}" "\tpost_trimming\t" $1}' > "!{sample_id}_trimmed_count.txt"
  '''
}

/*
 Collapse duplicate reads - likely PCR duplicates: duplicated in any case
*/

process collapse_duplicate_reads {
  label 'lowmem_threaded'                                                                

  input:
  // the filter{size()} functionality here checks if fastq is empty, 
  // which causes cd-hit-dup to choke
  // see: https://stackoverflow.com/questions/47401518/nextflow-is-an-input-file-empty
  tuple val(sample_id), path(input_fastq) from post_trim_ch.filter{ it[1].first().size() > 0}

  output:
  tuple val(sample_id), path("*_fu.fastq") optional true into post_collapse_count_ch
  tuple val(sample_id), path("*_fu.fastq") optional true into post_collapse_qc_ch
  tuple val(sample_id), path("*_fu.fastq") optional true into host_filtering_ch

  script:

  // this handles paired-end data, in which case must specify a paired output file
  def r1              = input_fastq[1] 
  def prefix_param    = input_fastq[1] ? "-u 30" : "-u 50"
  def paired_input    = input_fastq[1] ? "-i2 $r1" : ""
  def paired_output   = input_fastq[1] ? "-o2 ${sample_id}_R2_fu.fastq" : ""

  """
  cd-hit-dup \
   -e $params.mismatches_allowed \
   $prefix_param \
   -i ${input_fastq[0]} \
   $paired_input \
   -o ${sample_id}_R1_fu.fastq \
   $paired_output \
  """
}

/*
 Count post-collapse fastq
*/
process collapsed_fastq_count {
  label 'lowmem_non_threaded'                                                                
  publishDir "${params.counts_out_dir}", mode:'link'

  input:
  tuple val(sample_id), path(collapsed_fastq) from post_collapse_count_ch

  output:
  path("${sample_id}_collapsed_count.txt") into post_count_collapse_ch

  shell:
  // only count the first file because for paired-read data both files
  // will have the same # of reads.
  '''
  zcat -f !{collapsed_fastq[0]} | wc -l | xargs bash -c 'echo $(($0 / 4))' |  awk '{print "!{sample_id}" "\tpost_collapse\t" $1}' > "!{sample_id}_collapsed_count.txt"
  '''
}

/*
 Use fastqc to do QC on post-collapsed fastq
*/
process post_collapse_qc {
  label 'lowmem_non_threaded'

  input:
  tuple val(sample_id), path(input_fastq) from post_collapse_qc_ch

  output:
  val(sample_id) into post_collapse_multiqc_ch

  script:

  """
  mkdir -p  ${params.post_trim_fastqc_dir} 
  fastqc -o ${params.post_trim_fastqc_dir} $input_fastq
  """
}

/*
 Use multiqc to merge post-trimming reports
*/
process post_preprocess_multiqc {
  publishDir "${params.outdir}", mode:'link'

  input:
  val(all_sample_ids) from post_collapse_multiqc_ch.collect()

  output:
  path("post_trim_qc_report.html")

  """
  multiqc -n "post_trim_qc_report.html" -m fastqc -m cutadapt ${params.post_trim_fastqc_dir}
  """
}


/*
  Use bowtie2 to remove host-derived reads
*/
process host_filtering {
  label 'lowmem_threaded'                                                                
  publishDir "${params.host_filtered_out_dir}", mode:'link'

  input:
  tuple val(sample_id), path(input_fastq) from host_filtering_ch

  output:
  tuple val(sample_id), path("*_fuh.fastq") optional true into post_host_filter_count_ch
  tuple val(sample_id), path("*_fuh.fastq") optional true into post_host_filter_ch

  script:

  def bt_index = get_dataset_host_index(sample_id)

  // handle single-end or paired-end inputs
  def r1 = input_fastq[0] 
  def r2 = input_fastq[1] ? input_fastq[1] : ""
  def bowtie_file_input  = input_fastq[1] ? "-1 $r1 -2 $r2" : "-U $r1"
  def bowtie_file_output = input_fastq[1] ? "--un-conc ${sample_id}_R%_fuh.fastq" : "--un ${sample_id}_R1_fuh.fastq"

  if (bt_index == "dummy") {
    // this is the case where we don't have a host index specified 
    // for this sample.  Don't do host filtering.  Instead, just
    // use the ln (link) command to effectively copy the input fastq file
    // to the output file

    // handle paired end data case
    def r2_dummy_cmd = input_fastq[1]? "ln $r2 ${sample_id}_R2_fuh.fastq" : ""

    """
    ln $r1 ${sample_id}_R1_fuh.fastq
    $r2_dummy_cmd
    """
  }
  else {
    // case where there is a bowtie index specified for host filtering
    """
    bowtie2 \
    -x "${bt_index}" \
    --local \
    -q \
    $bowtie_file_input \
    --sensitive \
    --score-min "C,${params.host_bt_min_score},0" \
    -p ${task.cpus} \
    $bowtie_file_output 2> "${sample_id}.host_filtering_bt.log" > /dev/null 
    """
  }
}


/*
 Count post-host-filtering fastq
*/
process host_filtered_fastq_count {
  label 'lowmem_non_threaded'                                                                
  publishDir "${params.counts_out_dir}", mode:'link'

  input:
  tuple val(sample_id), path(filtered_fastq) from post_host_filter_count_ch

  output:
  path("${sample_id}_host_filtered_count.txt") into post_count_host_ch

  shell:
  // only count the first file because for paired-read data both files
  // will have the same # of reads.
  '''
  zcat -f !{filtered_fastq[0]} | wc -l | xargs bash -c 'echo $(($0 / 4))' |  awk '{print "!{sample_id}" "\tpost_host_filtered\t" $1}' > "!{sample_id}_host_filtered_count.txt"
  '''
}

/* 
  Collect all the fastq counts from various steps of processing and create a plot / excel summary
*/
process tabulate_fastq_counts {
  publishDir "${params.outdir}", mode: 'link'

  input:
  path(all_count_files) from post_count_initial_ch.concat(post_count_collapse_ch, post_count_trim_ch, post_count_host_ch).collect()

  output:
  path ("all_read_counts.txt") 
  path ("filtering_plots.pdf") 

  script:
  """
  Rscript ${params.R_bindir}/process_fastq_counts.R ${params.R_bindir} ${all_count_files}
  """
}

/*
  Assemble reads remaining after host filtering
*/
process assemble_remaining_reads {
  publishDir "${params.contigs_out_dir}", mode:'link'
  label 'highmem'

  // Spades will fail if, for instance, there are very few reads in the post-host-filtered datasets
  // this is expected for some kinds of datasets, such as water negative control datasets
  // we don't want the whole pipeline to stop in this case, so set errorStrategy to ignore
  // in this case, the pipeline 
  // see: https://www.nextflow.io/docs/latest/process.html#errorstrategy
  errorStrategy 'ignore'


  input:
  tuple val(sample_id), path(input_fastq) from post_host_filter_ch

  output:
  tuple val(sample_id), path("${sample_id}_contigs.fa"), path(input_fastq) into post_assembly_ch

  script:

  // handle single-end or paired-end inputs
  def r1 = input_fastq[0] 
  def r2 = input_fastq[1] ? input_fastq[1] : ""
  def spades_input  =      input_fastq[1] ? "-1 $r1 -2 $r2" : "-s $r1"
  // def bowtie_file_output = input_fastq[1] ? "--un-conc ${sample_id}_R%_fuh.fastq" : "--un ${sample_id}_R1_fuh.fastq"
  """
  # run spades
  spades.py -o ${sample_id}.spades ${spades_input} -t ${task.cpus} -m ${task.memory_gb}
  
  # consolidate output
  # this forces it into 1-line format
  # also remove contigs shorter than 150 bases
  seqtk seq -A ${sample_id}.spades/contigs.fasta | ${params.scripts_bindir}/filter_fasta_by_size > ${sample_id}_contigs.fa
  """
}

process quantify_read_mapping_to_contigs {
  label 'lowmem_threaded'
  publishDir "${params.contigs_out_dir}", mode:'link'

  input:
  tuple val(sample_id), path(contigs), path(input_fastq) from post_assembly_ch

  output:
  tuple val(sample_id), path("${sample_id}_contig_weights.txt"), path("${sample_id}_contigs_and_singletons.fasta") into post_contigs_singletons_ch

  script:
  def r1 = input_fastq[0] 
  """
  bowtie2-build $contigs contig_bt_index

  # C,120,1 makes min score 120 -> ~corresponds to 100% identity over ~60 bases
  # TODO: make score configurable?
  bowtie2 -x contig_bt_index --local --score-min C,120,1 -q -U $r1 -p ${task.cpus} -S contig_mapping_sam 

  # count the # of reads that hit each contig 
  # tally_sam_subjects script is just a series of piped bas commands: could move it here. 
  ${params.scripts_bindir}/tally_sam_subjects contig_mapping_sam > ${sample_id}_contig_weights.txt

  # create file of non-aligning aka singleton reads
  samtools view -f 4 contig_mapping_sam | samtools fasta > ${sample_id}_non_contig_mapping_reads.fasta

  # concatenate contigs and singleton
  cat ${sample_id}_non_contig_mapping_reads.fasta $contigs > ${sample_id}_contigs_and_singletons.fasta
  """
}

process blastn_contigs_and_singletons {
  label 'highmem'
  publishDir "${params.blast_out_dir}", mode:'link'

  input:
  tuple val(sample_id), path(contig_weights), path(contigs_and_singletons) from post_contigs_singletons_ch

  output:
  tuple val(sample_id), path(contig_weights), path("${sample_id}_contigs_and_singletons_n.fasta") into post_blastn_blastx_ch
  path("${sample_id}_contigs_and_singletons_n.fasta") into post_blastn_blastx_merge_ch
  tuple val(sample_id), path(contig_weights), path("${contigs_and_singletons}.bn_nt") into post_blastn_tally_ch
  tuple val(sample_id), path(contigs_and_singletons), path("${contigs_and_singletons}.bn_nt") into post_blastn_distribute_ch

  script:
  def blastn_columns = "qaccver saccver pident length mismatch gaps qstart qend sstart send evalue bitscore staxid ssciname scomname sblastname sskingdom"
  /*                                                                              
            staxid means Subject Taxonomy ID                                    
          ssciname means Subject Scientific Name                                
          scomname means Subject Common Name                                    
        sblastname means Subject Blast Name                                     
         sskingdom means Subject Super Kingdom                                  
  */                                                                              

  // TODO: run blast with different # of threads?  Would be good to benchmark first...
  // this post suggests past 4 threads you get diminishing returns
  // note that individual blastn processes here will use up to ~70Gb of RAM, possibly more
                                                                                
  """                                                                           
  # run megablast
  # using db nt here implies that you have a 

  # have to set this environmental variable so blast can be taxonomically aware 
  # poorly documented feature of command line blast 
  export "BLASTDB=${params.blast_db_dir}"

  # run the megablast 
  # blastn -num_threads $task.cpus -db ${params.nt_blast_db} -task megablast -evalue ${params.max_blast_nt_evalue} -query $contigs_and_singletons -outfmt "6 $blastn_columns" -out ${contigs_and_singletons}.bn_nt.no_header
  blastn -num_threads 4 -db ${params.nt_blast_db} -task megablast -evalue ${params.max_blast_nt_evalue} -query $contigs_and_singletons -outfmt "6 $blastn_columns" -out ${contigs_and_singletons}.bn_nt.no_header

  # prepend blast output with the column names so we don't have to manually name them later
  # and replace spaces with tabs                                                    
  echo $blastn_columns | perl -p -e 's/ /\t/g' > blast_header                   

  # concatenate header line and actual data
  # TODO: prepend with sample ID and type of blast (?)
  cat blast_header ${contigs_and_singletons}.bn_nt.no_header > ${contigs_and_singletons}.bn_nt

  # get rid of unneeded temporary file
  rm ${contigs_and_singletons}.bn_nt.no_header 

  # pull out remaining sequences 
  ${params.scripts_bindir}/fasta_from_blast -r ${contigs_and_singletons}.bn_nt > ${sample_id}_contigs_and_singletons_n.fasta
  """                                                                           
}

process tally_blastn_results {
  label 'lowmem_nonthreaded'
  publishDir "${params.tally_out_dir}", mode:'link'

  input:
  tuple val(sample_id), path(contig_weights), path(blast_out) from post_blastn_tally_ch

  output:
  tuple val(sample_id), path("*.tally") into post_tally_ch

  script:
  // TODO: move tally_blast_hits logic into an R script?
  // TODO: tally_blast_hits shouldn't need to do accession->taxid lookups anymore
  """
  ${params.scripts_bindir}/tally_blast_hits -lca -w $contig_weights $blast_out > ${blast_out}.tally
  ${params.scripts_bindir}/tally_blast_hits -lca -w $contig_weights -t -ti ${blast_out} > ${blast_out}.tab_tree.tally
  """
}

process distribute_blastn_results {
  label 'lowmem_nonthreaded'
  publishDir "${params.virus_seq_out_dir}", mode:'link'

  input:
  tuple val(sample_id), path(contigs_and_singletons), path(blast_out) from post_blastn_distribute_ch

  output:
  tuple val(sample_id), path("*.fasta_*") optional true

  script:
  // TODO: move distribute logic into an R script?
  """
  # pull out virus-derived reads - this will create a file for each viral taxon
  ${params.scripts_bindir}/distribute_fasta_by_blast_taxid -v $contigs_and_singletons $blast_out
  """
}

// merge fasta for faster combined blastx search 

post_blastn_blastx_merge_ch
      .collectFile(name: 'merged_contigs_and_singletons.fasta')
      .set{merged_blastx_ch}


process blastx_merged_contigs_and_singletons {
  label 'highmem'

  input:
  path(merged_contigs_and_singletons) from merged_blastx_ch

  output:
  path("${merged_contigs_and_singletons}.bx_nr") into merged_blastx_results_ch

  script:
  def blastx_columns = "qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids sscinames sskingdoms"
                                                                                
  """                                                                           
  # run diamond to do a blastx-style search

  # have to set this environmental variable so blast can be taxonomically aware 
  # poorly documented feature of command line blast 
  export "BLASTDB=${params.blast_db_dir}"

  diamond blastx --threads ${task.cpus} -c 1 -b 8 --db ${params.nr_diamond_db} --query $merged_contigs_and_singletons --more-sensitive --outfmt 6 $blastx_columns --evalue ${params.max_blasx_nr_evalue} --out ${merged_contigs_and_singletons}.bx_nr
  """
}

merged_blastx_results_ch
  .combine(post_blastn_blastx_ch)
  .set{merged_blastx_with_input_ch}

process split_merged_blastx_results {
  label 'highmem'
  publishDir "${params.blast_out_dir}", mode:'link'

  input:
  tuple path(merged_blastx_results), val(sample_id), path(contig_weights), path(contigs_and_singletons) from merged_blastx_with_input_ch

  output:
  tuple val(sample_id), path(contig_weights), path("${contigs_and_singletons}.bx_nr") into post_blastx_tally_ch
  tuple val(sample_id), path(contigs_and_singletons), path("${contigs_and_singletons}.bx_nr") into post_blastx_distribute_ch
  tuple val(sample_id), path("${sample_id}_contigs_and_singletons_nn.fasta") into post_blastx_unassigned_ch

  script:
  def blastx_columns = "qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids sscinames sskingdoms"
                                                                                
  """                                                                           
  # split out blastx hits for the contigs and singletons in this sample's file
  # grep -e '^>' $contigs_and_singletons | tr -d '>' > fasta_headers
  # pull out matching blastx results 
  # grep -f fasta_headers $merged_blastx_results > ${contigs_and_singletons}.bx_nr.no_header

  blast_from_fasta -f $contigs_and_singletons $merged_blastx_results > ${contigs_and_singletons}.bx_nr.no_header

  # prepend blast output with the column names so we don't have to manually name them later
  # and replace spaces with tabs                                                    
  echo $blastx_columns | perl -p -e 's/ /\t/g' > blast_header                   

  # concatenate header line and actual data
  # TODO: prepend with sample ID and type of blast (?)
  cat blast_header ${contigs_and_singletons}.bx_nr.no_header > ${contigs_and_singletons}.bx_nr

  # get rid of unneeded temporary file
  rm ${contigs_and_singletons}.bx_nr.no_header 

  # pull out remaining sequences 
  ${params.scripts_bindir}/fasta_from_blast -r ${contigs_and_singletons}.bx_nr > ${sample_id}_contigs_and_singletons_nn.fasta
  """                                                                           
}

process tally_blastx_results {
  label 'lowmem_nonthreaded'
  publishDir "${params.tally_out_dir}", mode:'link'

  input:
  tuple val(sample_id), path(contig_weights), path(blast_out) from post_blastx_tally_ch

  output:
  tuple val(sample_id), path("*.tally") 

  script:
  // TODO: move tally_blast_hits logic into an R script?
  // TODO: tally_blast_hits shouldn't need to do accession->taxid lookups anymore
  """
  ${params.scripts_bindir}/tally_blast_hits -lca -w $contig_weights $blast_out > ${blast_out}.tally
  ${params.scripts_bindir}/tally_blast_hits -lca -w $contig_weights -t -ti ${blast_out} > ${blast_out}.tab_tree.tally
  """
}

process distribute_blastx_results {
  label 'lowmem_nonthreaded'
  publishDir "${params.virus_seq_out_dir}", mode:'link'

  input:
  tuple val(sample_id), path(contigs_and_singletons), path(blast_out) from post_blastx_distribute_ch

  output:
  tuple val(sample_id), path("*.fasta_*") optional true

  script:
  // TODO: move distribute logic into an R script?
  """
  # pull out virus-derived reads - this will create a file for each viral taxon
  ${params.scripts_bindir}/distribute_fasta_by_blast_taxid -v $contigs_and_singletons $blast_out
  """
}
