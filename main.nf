#!/usr/bin/env nextflow

/*
    Stenglein lab metagenomic classification pipeline 
    implemented in nextflow

    Feb 25, 2022 

    Mark Stenglein
*/


// output usage message
params.help = false
params.h = false

// TODO: command line arg processing and validating 

// TODO: check that appropriate refseq files exist (host bt index, nr/nt, etc.)

/*
 These fastq files represent the main input to this workflow

 See here for info on glob pattern matching:
    https://docs.oracle.com/javase/tutorial/essential/io/fileOps.html#glob
*/
Channel
  .fromFilePairs("${params.fastq_dir}/${params.fastq_pattern}", size: -1, checkIfExists: true, maxDepth: 1)
  // .into {samples_ch_qc; samples_ch_trim; samples_ch_count; host_setup_ch}
  .into {initial_fastq_subsample_ch; initial_fastq_ch}


// define usage output message
// TODO:
def helpMessage() {
  log.info """
        Usage:
        The typical command for running the pipeline is as follows:
        nextflow run main.nf --

        Mandatory arguments:
        One of: 
        --host_map_file                The path to a file that contains host filtering
                                       info for all datasets or subsets of datasets.   

                                       This enables different filtering for different datasets
          
                                       This file should be a 2-column tab-delimited file 
                                       with:

                                       - the first column defining dataset IDs or patterns that will
                                         match dataset IDs

                                       - the second column will be the path of a bowtie index that will be
                                         used to filter out host reads
                                       

        --host bt_index                The path to a bowtie index that will be used 
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
  if (!params.host_map_file && !params.host_bt_index){
    log.info """
      Warning: you have not specified either of these two parameters.  
               No host filtering will be performed.

        1. host_map_file (--host_map_file) or 
        2. host_bt_index (--host_bt_index) 
     
    """
    helpMessage()
  }

  if (params.host_map_file && params.host_bt_index){
    log.info """
      Error: you may specify only one of these two parameters:

        1. host_map_file (--host_map_file) and 
        2. host_bt_index (--host_bt_index) 
    """
    helpMessage()
    // Stop execution
    System.exit(1)
  }

  if (params.host_map_file) {
    host_map_file = file(params.host_map_file)
    if (!host_map_file.exists()) {  
      log.info  """
        Error: host_map_file ($params.host_map_file) does not exist 
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
  if (params.host_map_file) {

    host_map_file = file(params.host_map_file)

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
process subsample_input {
  label 'lowmem'

  // singularity info for this process                                          
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
      container "https://depot.galaxyproject.org/singularity/seqtk:1.3--h5bf99c6_3"
  } else {                                                                      
      container "quay.io/biocontainers/seqtk:1.3--h5bf99c6_3"
  }   

  when: 
  params.subsample_fraction

  input:
  tuple val (sample_id), path(input_fastq) from initial_fastq_subsample_ch

  output:
  tuple val (sample_id), path("*.ss.fastq*") into subsampled_fastq_ch

  script:
  def r1 = input_fastq[0]
  def r1_ss = r1.name.replaceAll("fastq", "ss.fastq")
  def r2 = input_fastq[1] 
  def r2_ss = input_fastq[1] ? r2.name.replaceAll("fastq", "ss.fastq") : ""
  def r1_command = "seqtk sample $r1 ${params.subsample_fraction} | gzip > $r1_ss" 
  def r2_command = input_fastq[1] ? "seqtk sample $r2 ${params.subsample_fraction} | gzip > $r2_ss" : ""
  """
  $r1_command
  $r2_command
  """
}

/*
 Fork fastq input channel
 */
if (params.subsample_fraction) {
  subsampled_fastq_ch.into {samples_ch_qc; samples_ch_trim; samples_ch_count; host_setup_ch}
} else {
  initial_fastq_ch.into {samples_ch_qc; samples_ch_trim; samples_ch_count; host_setup_ch}
}


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

  // singularity info for this process                                          
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
      container "https://depot.galaxyproject.org/singularity/fastqc:0.11.9--0"  
  } else {                                                                      
      container "quay.io/biocontainers/fastqc:0.11.9--0"                        
  }     

  input:
  tuple val(sample_id), path(initial_fastq) from samples_ch_qc

  output:
  // path("*.html") into initial_qc_fastqc_html_ch
  path("*.zip") into initial_qc_fastqc_zip_ch

  script:
  """
  fastqc $initial_fastq
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

  // singularity info for this process                                          
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
      container "https://depot.galaxyproject.org/singularity/multiqc:1.11--pyhdfd78af_0"
  } else {                                                                      
      container "quay.io/biocontainers/multiqc:1.11--pyhdfd78af_0"              
  }   

  input:
  path("fastqc_zips/*") from initial_qc_fastqc_zip_ch.collect()

  output: 
  path("initial_qc_report.html")

  script:
  """
  multiqc -n "initial_qc_report.html" -m fastqc fastqc_zips
  """
}

/*
 Use cutadapt to trim off adapters and low quality bases
*/
process trim_adapters_and_low_quality {
  label 'lowmem_non_threaded'                                                                

  // singularity info for this process                                          
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
    container "https://depot.galaxyproject.org/singularity/cutadapt:3.5--py39h38f01e4_0"
  } else {                                                                      
    container "quay.io/biocontainers/cutadapt:3.5--py39h38f01e4_0"              
  }        

  input:
  tuple val(sample_id), path(initial_fastq) from samples_ch_trim

  output:
  tuple val(sample_id), path("*_f.fastq") optional true into post_trim_ch
  tuple val(sample_id), path("*_f.fastq") optional true into post_trim_count_ch

  // TODO: parameterize adapter sequences
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

  // singularity info for this process                                          
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
    container "https://depot.galaxyproject.org/singularity/cd-hit-auxtools:4.8.1--h7d875b9_1"
  } else {                                                                      
    container "quay.io/biocontainers/cd-hit-auxtools:4.8.1--h7d875b9_1"
  }        

  input:
  // the filter{size()} functionality here checks if fastq is empty, 
  // which causes cd-hit-dup to choke
  // see: https://stackoverflow.com/questions/47401518/nextflow-is-an-input-file-empty  
  // this means that fastq that are empty at this stage will just stop going through pipeline 
  // TODO: there shouldn't be an asterisk (spread operator) in this filter step ?  
  tuple val(sample_id), path(input_fastq) from post_trim_ch.filter{ it[1]*.getAt(0).size() > 0}

  output:
  tuple val(sample_id), path("*_fu.fastq") optional true into post_collapse_count_ch
  tuple val(sample_id), path("*_fu.fastq") optional true into post_collapse_qc_ch
  tuple val(sample_id), path("*_fu.fastq") optional true into host_filtering_ch

  script:

  // this handles paired-end data, in which case must specify a paired output file
  def r1 = input_fastq[0]
  def r2 = input_fastq[1] 
  def prefix_param    = input_fastq[1] ? "-u 30" : "-u 30"
  def paired_input    = input_fastq[1] ? "-i2 $r2" : ""
  def paired_output   = input_fastq[1] ? "-o2 ${sample_id}_R2_fu.fastq" : ""

  """
  cd-hit-dup \
   -e $params.mismatches_allowed \
   $prefix_param \
   -i $r1 \
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

  // singularity info for this process                                          
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
      container "https://depot.galaxyproject.org/singularity/fastqc:0.11.9--0"  
  } else {                                                                      
      container "quay.io/biocontainers/fastqc:0.11.9--0"                        
  }     

  input:
  tuple val(sample_id), path(input_fastq) from post_collapse_qc_ch

  output:
  // path("*.html") into post_trim_fastqc_html_ch
  path("*.zip") into post_trim_fastqc_zip_ch

  script:

  """
  fastqc $input_fastq
  """
}

/*
 Use multiqc to merge post-trimming reports
*/
process post_preprocess_multiqc {
  publishDir "${params.outdir}", mode:'link'

  // singularity info for this process                                          
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
      container "https://depot.galaxyproject.org/singularity/multiqc:1.11--pyhdfd78af_0"
  } else {                                                                      
      container "quay.io/biocontainers/multiqc:1.11--pyhdfd78af_0"              
  }   

  input:
  path("fastqc_zips/*") from post_trim_fastqc_zip_ch.collect()

  output:
  path("post_trim_qc_report.html")

  """
  multiqc -n "post_trim_qc_report.html" -m fastqc -m cutadapt fastqc_zips
  """
}

/*
  Use bowtie2 to remove host-derived reads
*/
process host_filtering {
  label 'lowmem_threaded'                                                                
  publishDir "${params.host_filtered_out_dir}", mode:'link'

  // singularity info for this process                                          
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
      container "https://depot.galaxyproject.org/singularity/bowtie2:2.4.5--py39ha4319a6_1"
  } else {                                                                      
      container "quay.io/biocontainers/bowtie2:2.4.5--py39ha4319a6_1"
  }   

  input:
  tuple val(sample_id), path(input_fastq) from host_filtering_ch

  output:
  tuple val(sample_id), path("*_fuh.fastq") optional true into post_host_filter_ch
  tuple val(sample_id), path("*_fuh.fastq") optional true into post_host_filter_count_ch
  tuple val(sample_id), path("*_fuh.fastq") optional true into post_host_filter_remap_ch


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

  // singularity info for this process                                          
  if (workflow.containerEngine == 'singularity') {
      container "library://stenglein-lab/r_taxonomy_tools/r_taxonomy_tools:1.0.0"                 
  }          

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
  label 'highmem'

  // singularity info for this process                                          
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
      container "https://depot.galaxyproject.org/singularity/spades:3.15.4--h95f258a_0"
  } else {                                                                      
      container "quay.io/biocontainers/spades:3.15.4--h95f258a_0"
  }   

  // Spades will fail if, for instance, there are very few reads in the post-host-filtered datasets
  // this is expected for some kinds of datasets, such as water negative control datasets
  // we don't want the whole pipeline to stop in this case, so set errorStrategy to ignore
  // in this case, the pipeline 
  // see: https://www.nextflow.io/docs/latest/process.html#errorstrategy
  errorStrategy 'ignore'

  input:
  tuple val(sample_id), path(input_fastq) from post_host_filter_ch

  output:
  tuple val(sample_id), path("${sample_id}.spades"), path(input_fastq) into post_assembly_ch

  script:

  // handle single-end or paired-end inputs
  def r1 = input_fastq[0] 
  def r2 = input_fastq[1] ? input_fastq[1] : ""
  def spades_input  =      input_fastq[1] ? "-1 $r1 -2 $r2" : "-s $r1"
  // def bowtie_file_output = input_fastq[1] ? "--un-conc ${sample_id}_R%_fuh.fastq" : "--un ${sample_id}_R1_fuh.fastq"
  """
  # run spades
  spades.py -o ${sample_id}.spades ${spades_input} -t ${task.cpus} 
  """
}

/*
  Get contigs from assembly, converting into 1-line fasta format
*/
process retrieve_contigs {
  publishDir "${params.contigs_out_dir}", mode:'link', pattern:"*.fa"
  label 'lowmem'

  // singularity info for this process                                          
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
      container "https://depot.galaxyproject.org/singularity/seqtk:1.3--h5bf99c6_3"
  } else {                                                                      
      container "quay.io/biocontainers/seqtk:1.3--h5bf99c6_3"
  }   

  input:
  tuple val(sample_id), path("${sample_id}.spades"), path(input_fastq) from post_assembly_ch

  output:
  tuple val(sample_id), path("${sample_id}_contigs.fa"), path(input_fastq) into post_contigs_ch

  script:
  """
  # consolidate assembly output
  # this forces contigs fasta into 1-line format
  # also remove contigs shorter than minimum_contig_length bases
  seqtk seq -A -L ${params.minimum_contig_length} ${sample_id}.spades/contigs.fasta > ${sample_id}_contigs.fa
  """
}

/*
  Map reads to contigs so that contigs can be weighted by the # of 
  reads they account for.
*/
process quantify_read_mapping_to_contigs {
  label 'lowmem_threaded'
  // publishDir "${params.contigs_out_dir}", mode:'link'

  // singularity info for this process                                          
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
      container "https://depot.galaxyproject.org/singularity/bowtie2:2.4.5--py39ha4319a6_1"
  } else {                                                                      
      container "quay.io/biocontainers/bowtie2:2.4.5--py39ha4319a6_1"
  }   

  input:
  tuple val(sample_id), path(contigs), path(input_fastq) from post_contigs_ch

  output:
  tuple val(sample_id), path("${sample_id}_contig_weights.txt"), path(contigs), path(contig_mapping_sam) into post_contigs_weight_ch

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
  """
}

/*
  Optionally merge contigs and singletons
*/
process merge_contigs_and_singletons {
  label 'lowmem_threaded'
  // publishDir "${params.contigs_out_dir}", mode:'link'

  // singularity info for this process                                          
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
      container "https://depot.galaxyproject.org/singularity/samtools:1.17--h00cdaf9_0"
  } else {                                                                      
      container "quay.io/biocontainers/samtools:1.17--h00cdaf9_0"
  }     

  input:
  tuple val(sample_id), path(contig_weights), path(contigs), path(contig_mapping_sam) from post_contigs_weight_ch

  output:
  tuple val(sample_id), path(contig_weights), path("${sample_id}_contigs_and_singletons.fasta") into post_contigs_singletons_ch

  script:
  // classify singletons and contigs - slower but more thorough
  if (params.classify_singletons)
  """
  # create file of non-aligning aka singleton reads
  samtools view -f 4 contig_mapping_sam | samtools fasta > ${sample_id}_non_contig_mapping_reads.fasta

  # concatenate contigs and singleton
  cat ${sample_id}_non_contig_mapping_reads.fasta $contigs > ${sample_id}_contigs_and_singletons.fasta
  """

  // don't classify singletons - only classify contigs
  else
  """
  # simply create a link of the contigs file named contigs_and_singletons
  # TODO: fix this misleading naming scheme in case when not considering 
  #       singletons
  ln $contigs ${sample_id}_contigs_and_singletons.fasta
  """
}

/*
  Create a channel containing the path to the directory containing a local copy of the blast nt
  database.  This is so that this directory can be staged (by link) in the work directory.
 */
blast_db_ch   = Channel.fromPath( params.local_nt_database_dir, type: 'dir', checkIfExists: true)

/*
  This process confirms that the local NT database is valid 
 */

process check_local_blast_database {
  label 'process_low'

  // singularity info for this process
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
      container "https://depot.galaxyproject.org/singularity/blast:2.12.0--pl5262h3289130_0"
  } else {
      container "quay.io/biocontainers/blast:2.12.0--pl5262h3289130_0"
  }

  input:
  path(local_nt_database_dir) from blast_db_ch

  output:
  // tuple path("${local_nt_database_dir}"), val("${params.local_nt_database_name}") into post_blast_db_check_ch
  path(local_nt_database_dir) into post_blast_db_check_ch


  script:                                                                       
  def local_nt_database = "${local_nt_database_dir}/${params.local_nt_database_name}"

  """
  # check BLAST database
  # check for expected .nal file: if not present, output a helpful warning message
  if [ ! -f "${local_nt_database}.nal" ]                                 
  then                                                                          
    echo "ERROR: it does not appear that ${local_nt_database} is a valid BLAST database."
  fi   
  
  # check validity of database with blastdbcmd.  If not valid, this will error 
  # and stop pipeline.
  blastdbcmd -db ${local_nt_database} -info 

  """
}

/*
  Create a channel containing the path to the directory containing a local copy of the diamond
  database.  This is so that this directory can be staged (by link) in the work directory.
 */
diamond_db_ch = Channel.fromPath( params.local_diamond_database_dir, type: 'dir', checkIfExists: true)

/*
  This process confirms that the local diamond database is valid 
 */

process check_local_diamond_database {
  label 'process_low'

  // singularity info for this process                                          
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
      container "https://depot.galaxyproject.org/singularity/diamond:2.0.14--hdcc8f71_0"
  } else {                                                                      
      container "quay.io/biocontainers/diamond:2.0.14--hdcc8f71_0"
  }                                                                             

  input:
  path(local_diamond_database_dir) from diamond_db_ch

  output:
  path(local_diamond_database_dir) into post_diamond_db_check_ch

  script:                                                                       
  def local_dmd_database = "${local_diamond_database_dir}/${params.local_diamond_database_name}"

  """
  # check Diamond database
  # check for expected file: if not present, output a helpful warning message
  if [ ! -f "${local_dmd_database}" ]                                 
  then                                                                          
    echo "ERROR: it does not appear that Diamond database ${local_dmd_database} exists."
  fi   
  
  # check validity of database with diamond dbinfo.  If not valid, this will error 
  # and stop pipeline.
  diamond dbinfo --db ${local_dmd_database}
  """
}

/*
  Create a channel containing the path to an (optional) local copy of the 
  NCBI taxonomy databases
 */
ncbi_tax_ch = Channel.empty()

// if this path was provided as a parameter, then create a channel
// from this path and set a boolean to true to indicate it's an existing
// directory
if (params.ncbi_tax_dir) {
   ncbi_tax_ch = Channel.fromPath( params.ncbi_tax_dir )
                         .map { path -> [ path , true ] }  
} else {
   // if this path was *not* provided as a parameter, then create a channel
   // from a bogus path "ncbi_tax_dir" and set a boolean to false 
   // to indicate it *doesn't* refer to an existing directory
   ncbi_tax_ch = Channel.fromPath( "ncbi_tax_dir" )
                         .map { path -> [ path , false ] }  
}

/* 
  Download NCBI Taxonomy database files if necessary
 */
process download_ncbi_taxonomy_if_necessary {
  label 'process_low'

  // singularity info for this process
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
      container "https://depot.galaxyproject.org/singularity/curl:7.80.0"
  } else {
      container "quay.io/biocontainers/curl:7.80.0"
  }

  input:
  tuple path(ncbi_tax_dir), val(existing_db) from ncbi_tax_ch

  output:
  tuple path (ncbi_tax_dir), val(existing_db) into post_ncbi_tax_download_ch

  script:
  // if a local directory containing the ncbi taxonomy database  is specified, 
  // check that it exists and contains the expected sqlite db file
  existing_db ? 
  """
    # check that the directory exists
    if [ ! -d "${ncbi_tax_dir}" ] ; then
      echo "ERROR: BLAST taxonomy directory ${ncbi_tax_dir} (--ncbi_tax_dir) does not exist."
      exit 1
    fi 
    # check that appropriate files exist
    if [ ! -f "${ncbi_tax_dir}/${params.ncbi_tax_db}" ] ; then
      echo "ERROR: required NCBI taxonomy database file ${params.ncbi_tax_db} (--ncbi_tax_db) not in directory ${ncbi_tax_dir} (--ncbi_tax_dir)."
      exit 1
    fi 
  """ :
  // if tax db doesn't already exist : download the necessary files and keep track of directory 
  // log.info("Downloading the NCBI Taxonomy database files.")
  """
     rm $ncbi_tax_dir
     mkdir $ncbi_tax_dir
     ${params.scripts_bindir}/download_taxonomy_databases 
     mv *.dmp $ncbi_tax_dir
  """
}

/* 
  Pre-process NCBI Taxonomy database files 
 */
process preprocess_ncbi_taxonomy_files {
  label 'process_low'

  // singularity info for this process
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
      container "https://depot.galaxyproject.org/singularity/perl:5.26.2"
  } else {
      container "quay.io/biocontainers/perl:5.26.2"
  }

  input:
  tuple path (ncbi_tax_dir), val(existing_db) from post_ncbi_tax_download_ch

  output:
  tuple path ("prepped_ncbi_tax_dir"), val(existing_db) into post_ncbi_tax_prep_ch

  script:
  existing_db ? 
  // if the path to a local copy of the ncbi taxonomy database is specified, simply link to the directory
  """
     echo "nothing to do"   
     ln -s $ncbi_tax_dir prepped_ncbi_tax_dir
  """ :
  // deal with the newly downloaded *.dmp files 
  // pull out the info we need: this will be dumped into a sqlite db in the
  // next downstream process
  """
     rm -rf prepped_ncbi_tax_dir
     mkdir -p prepped_ncbi_tax_dir
     ${params.scripts_bindir}/preprocess_dmp_files $ncbi_tax_dir prepped_ncbi_tax_dir ${params.scripts_bindir}
  """
}

/* 
  Setup NCBI Taxonomy databases: use existing installation if available
  or download and create new one
 */
process create_ncbi_tax_sqlite_db {
  label 'process_low'
  publishDir "${params.tax_db_out_dir}", mode:'link', pattern:"*.sqlite3"

  // singularity info for this process
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
      container "https://depot.galaxyproject.org/singularity/sqlite:3.33.0"
  } else {
      container "quay.io/biocontainers/sqlite:3.33.0"
  }

  input:
  tuple path (ncbi_tax_dir), val(existing_db) from post_ncbi_tax_prep_ch

  output:
  path ("${params.ncbi_tax_db}") into post_ncbi_tax_check_ch

  script:
  // if a local blast_tax_dir is specified, check that it contains the expected files
  // and that these are valid SQLite databases
  existing_db ? 
  """
    # validate that a SQLite database is minimally valid
    # see: https://stackoverflow.com/questions/3888529/how-to-tell-if-sqlite-database-file-is-valid-or-not
    # This sqlite3 command will exit with 0 status (success) if everything looks OK
    # and will exit with non-zero exit status otherwise
    sqlite3 -batch ${ncbi_tax_dir}/${params.ncbi_tax_db} <<"EOF"
    pragma schema_version;
    EOF

    # create a link in the current directory to the existing sqlite db
    ln -s "$ncbi_tax_dir/${params.ncbi_tax_db}" ${params.ncbi_tax_db} 

  """ :
  // create a sqlite database containing NCBI taxonomy files if it doesn't already exist
  """
     ${params.scripts_bindir}/create_sqlite_taxonomy_databases $ncbi_tax_dir $params.ncbi_tax_db 
  """
}

post_ncbi_tax_check_ch.first().into{tax_db_ch_1; tax_db_ch_2; tax_db_ch_3; tax_db_ch_4}

/*
  Create a channel containing the path to an (optional) local copy of the 
  NCBI blast/taxonomy db to make BLAST taxonomically aware
 */
blast_tax_ch = Channel.empty()

// if this path was provided as a parameter, then create a channel
// from this path and set a boolean to true to indicate it's an existing
// directory
if (params.blast_tax_dir) {
   blast_tax_ch = Channel.fromPath( params.blast_tax_dir )
                         .map { path -> [ path , true ] }  
} else {
   // if this path was *not* provided as a parameter, then create a channel
   // from a bogus path and set a boolean to false 
   // to indicate it *doesn't* refer to an existing directory
   blast_tax_ch = Channel.fromPath( "does_not_exist" )
                         .map { path -> [ path , false ] }  
}

/* 
  Confirm that local blast will be taxonomically aware
 */
process check_blast_tax {
  label 'process_low'
  publishDir "${params.tax_db_out_dir}", mode:'link'

  // singularity info for this process
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
      container "https://depot.galaxyproject.org/singularity/curl:7.80.0"
  } else {
      container "quay.io/biocontainers/curl:7.80.0"
  }

  input:
  tuple path(blast_tax_dir), val(existing_db) from blast_tax_ch

  output:
  path ("checked_blast_tax_dir") into post_blast_tax_check_ch
  path ("checked_blast_tax_dir") into post_blast_tax_check_dmd_ch

  script:
  // if a local blast_tax_dir is specified, check that it contains the expected files
  existing_db ? 
  """
    # check that the directory exists
    if [ ! -d "${blast_tax_dir}" ] ; then
      echo "ERROR: BLAST taxonomy directory ${blast_tax_dir} (--blast_tax_dir) does not exist."
      exit 1
    fi 
    # check that appropriate files exist
    if [ ! -f "${blast_tax_dir}/taxdb.btd" ] ; then
      echo "ERROR: required BLAST taxonomy file taxdb.btd not in directory ${blast_tax_dir} (--blast_tax_dir)."
      exit 1
    fi 
    if [ ! -f "${blast_tax_dir}/taxdb.bti" ] ; then
      echo "ERROR: required BLAST taxonomy file taxdb.bti not in directory ${blast_tax_dir} (--blast_tax_dir)."
      exit 1
    fi 

    # create a link to the actual dir and output that link
    ln $blast_tax_dir checked_blast_tax_dir                                    
  
  """ :
  // if tax db doesn't already exist : download the necessary files and keep track of directory 
  """
    # make a new local directory to contain the files
    # first remove broken link to non-existent file
    rm $blast_tax_dir
    # make a new directory, into which we'll put the blast tax files
    mkdir checked_blast_tax_dir
    # download taxdb files
    curl -OL ${params.blast_tax_url}
    # unpack archive
    tar xvf taxdb.tar.gz
    # move files to blast_tax_dir
    mv taxdb.??? checked_blast_tax_dir
    # get rid of archive
    rm taxdb.tar.gz
  """
}

/*
  Use blastn to identify closest related database sequences
*/
process blastn_contigs_and_singletons {
  label 'highmem'
  // publishDir "${params.blast_out_dir}", mode:'link'

  // singularity info for this process                                          
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
      container "https://depot.galaxyproject.org/singularity/blast:2.12.0--pl5262h3289130_0"
  } else {                                                                      
      container "quay.io/biocontainers/blast:2.12.0--pl5262h3289130_0"          
  }                                                                             

  input:
  tuple val(sample_id), path(contig_weights), path(contigs_and_singletons) from post_contigs_singletons_ch
  // the collect operators for the next 2 steps are required to convert the paths to value channels so they
  // repeat for each input.  
  path(blast_tax_dir), stageAs: 'local_blast_tax_dir' from post_blast_tax_check_ch.collect()
  path(local_nt_database_dir), stageAs: 'local_nt_db_dir'  from post_blast_db_check_ch.collect()

  output:
  tuple path("${contigs_and_singletons}.bn_nt"), val(sample_id), path(contig_weights), path(contigs_and_singletons) into post_blastn_ch

  script:
  def local_nt_database = "${local_nt_database_dir}/${params.local_nt_database_name}"

  def blastn_columns = "qaccver saccver pident length mismatch gaps qstart qend sstart send evalue bitscore staxid ssciname scomname sblastname sskingdom"
  /*                                                                              
            staxid means Subject Taxonomy ID                                    
          ssciname means Subject Scientific Name                                
          scomname means Subject Common Name                                    
        sblastname means Subject Blast Name                                     
         sskingdom means Subject Super Kingdom                                  
  */                                                                              

  // TODO: run blast with different # of threads?  Would be good to benchmark first...
  // note that individual blastn processes here will use up to ~70Gb of RAM, possibly more
                                                                                
  """                                                                           
  # run megablast
  # using db nt here implies that you have a local copy of this db installed
  # TODO: remote blast option (may not work with taxids, etc.)

  # have to set this environmental variable so blast can be taxonomically aware 
  # poorly documented feature of command line blast 
  export BLASTDB="$blast_tax_dir"

  # run the megablast 
  blastn -num_threads $task.cpus -db $local_nt_database -task megablast -evalue ${params.max_blast_nt_evalue} -query $contigs_and_singletons -outfmt "6 $blastn_columns" -out ${contigs_and_singletons}.bn_nt
  """                                                                           
}

/*
  This splits the merged blastn results into results for each sample

  The split is based on the query IDs in the blast results and whether
  they match read or contig IDs present in the input contigs and singletons fasta
  for individual samples.  Uses blast_from_fasta script.
*/
process process_blastn_output {
  label 'highmem'
  publishDir "${params.blast_out_dir}", mode:'link'

  // singularity info for this process                                          
  if (workflow.containerEngine == 'singularity') {
      container "https://depot.galaxyproject.org/singularity/ubuntu:20.04"
  }          

  input:
  tuple path(blastn_results), val(sample_id), path(contig_weights), path(contigs_and_singletons) from post_blastn_ch

  output:
  tuple val(sample_id), path(contig_weights), path("${sample_id}_contigs_and_singletons_n.fasta") into post_blastn_blastx_ch
  path("${sample_id}_contigs_and_singletons_n.fasta") into post_blastn_blastx_merge_ch
  tuple val(sample_id), path(contig_weights), path("${contigs_and_singletons}.bn_nt") into post_blastn_tally_ch
  tuple val(sample_id), path(contigs_and_singletons), path("${contigs_and_singletons}.bn_nt") into post_blastn_distribute_ch

  script:
  // TODO: this is defined twice: possible issue
  def blastn_columns = "qaccver saccver pident length mismatch gaps qstart qend sstart send evalue bitscore staxid ssciname scomname sblastname sskingdom"
  /*                                                                              
            staxid means Subject Taxonomy ID                                    
          ssciname means Subject Scientific Name                                
          scomname means Subject Common Name                                    
        sblastname means Subject Blast Name                                     
         sskingdom means Subject Super Kingdom                                  
  */                                                                              
                                                                                
  """                                                                           
  mv $blastn_results ${blastn_results}.no_header

  # TODO: move this code block up to the blast process to avoid defining blastn_columns twice
  # prepend blast output with the column names so we don't have to manually name them later
  # and replace spaces with tabs                                                    
  echo $blastn_columns | perl -p -e 's/ /\t/g' > blast_header                   

  # concatenate header line and actual data
  cat blast_header ${blastn_results}.no_header > ${blastn_results}

  # get rid of unneeded temporary file
  rm ${blastn_results}.no_header 

  # pull out remaining sequences 
  ${params.scripts_bindir}/fasta_from_blast -r ${blastn_results} > ${sample_id}_contigs_and_singletons_n.fasta
  """                                                                           
}

/*
  Use the tally_blast_hits script to tally the number of reads mapping to different taxa
*/
process tally_blastn_results {
  label 'lowmem_nonthreaded'
  publishDir "${params.tally_out_dir}", mode:'link'

  // singularity info for this process                                          
  if (workflow.containerEngine == 'singularity') {
      container "library://stenglein-lab/perl_taxonomy_tools/perl_taxonomy_tools:1.0.0"
  }          

  input:
  tuple val(sample_id), path(contig_weights), path(blast_out) from post_blastn_tally_ch
  path(tax_db) from tax_db_ch_1

  output:
  path("*bn_nt.tally") into blastn_tally_file_ch
  tuple val ("genus"),       path("*.genus.tally")        into genus_tally_ch
  tuple val ("family"),      path("*.family.tally")       into family_tally_ch
  tuple val ("order"),       path("*.order.tally")        into order_tally_ch
  tuple val ("class"),       path("*.class.tally")        into class_tally_ch
  tuple val ("kingdom"),     path("*.kingdom.tally")      into kingdom_tally_ch
  tuple val ("superkingdom"),path("*.superkingdom.tally") into superkingdom_tally_ch

  script:

  // TODO: move tally_blast_hits logic into an R script?
  // TODO: tally_blast_hits shouldn't need to do accession->taxid lookups if taxid are in blast results
  // ${params.scripts_bindir}/tally_blast_hits -ntd $local_tax_db_dir/${params.ncbi_tax_db} -lca -w $contig_weights $blast_out > ${blast_out}.tally
  // ${params.scripts_bindir}/tally_blast_hits -ntd $local_tax_db_dir/${params.ncbi_tax_db} -lca -w $contig_weights -t -ti ${blast_out} > ${blast_out}.tab_tree.tally
  """
  ${params.scripts_bindir}/tally_blast_hits -ntd $tax_db -lca -w $contig_weights $blast_out > ${blast_out}.tally
  ${params.scripts_bindir}/tally_blast_hits -ntd $tax_db -lca -w $contig_weights -t -ti ${blast_out} > ${blast_out}.tab_tree.tally
  ${params.scripts_bindir}/tally_blast_hits -ntd $tax_db -lca -w $contig_weights -t -r genus ${blast_out} > ${blast_out}.genus.tally
  ${params.scripts_bindir}/tally_blast_hits -ntd $tax_db -lca -w $contig_weights -t -r family ${blast_out} > ${blast_out}.family.tally
  ${params.scripts_bindir}/tally_blast_hits -ntd $tax_db -lca -w $contig_weights -t -r order ${blast_out} > ${blast_out}.order.tally
  ${params.scripts_bindir}/tally_blast_hits -ntd $tax_db -lca -w $contig_weights -t -r class ${blast_out} > ${blast_out}.class.tally
  ${params.scripts_bindir}/tally_blast_hits -ntd $tax_db -lca -w $contig_weights -t -r kingdom ${blast_out} > ${blast_out}.kingdom.tally
  ${params.scripts_bindir}/tally_blast_hits -ntd $tax_db -lca -w $contig_weights -t -r superkingdom ${blast_out} > ${blast_out}.superkingdom.tally
  """
}

/*
   Use the distribute_fasta_by_blast_taxid script to pull out putative virus sequences
*/
process distribute_blastn_results {
  label 'lowmem_nonthreaded'
  publishDir "${params.virus_seq_out_dir}", mode:'link'

  // singularity info for this process                                          
  if (workflow.containerEngine == 'singularity') {
      container "library://stenglein-lab/perl_taxonomy_tools/perl_taxonomy_tools:1.0.0"
  }          

  input:
  tuple val(sample_id), path(contigs_and_singletons), path(blast_out) from post_blastn_distribute_ch
  path(tax_db) from tax_db_ch_2
                                                                                
  output:                                                                       
  tuple val(sample_id), path("*.fasta_*") optional true into virus_fasta_ch

  script:
  """
  # pull out virus-derived reads - this will create a file for each viral taxon
  ${params.scripts_bindir}/distribute_fasta_by_blast_taxid -ntd $tax_db -v $contigs_and_singletons $blast_out
  """
}

/*
  This creates a new channel that merges all of the contigs and singletons
  that remain unassigned after blastn-based assignments
  into a single merged file that will be the input to a single diamond (blastx) call
  
  Run diamond as a merged single search because it takes so long to load
  the database into memory and because each individual process
  uses so much memory.  So do it all as one then split up results afterwards.
*/
post_blastn_blastx_merge_ch
      .collectFile(name: 'merged_contigs_and_singletons.fasta')
      .set{merged_blastx_ch}

/*
  Run diamond to attempt to assign contigs (and singletons) remaining
  after blastn assignment
*/
process blastx_merged_contigs_and_singletons {
  label 'highmem'

  // singularity info for this process                                          
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
      container "https://depot.galaxyproject.org/singularity/diamond:2.0.14--hdcc8f71_0"
  } else {                                                                      
      container "quay.io/biocontainers/diamond:2.0.14--hdcc8f71_0"
  }                                                                             

  input:
  // the filter here will prevent empty inputs from running
  path(merged_contigs_and_singletons) from merged_blastx_ch.filter{ it.size() > 0 }  
  path(local_diamond_database_dir) from post_diamond_db_check_ch.collect()
  path(blast_tax_dir), stageAs: 'local_blast_tax_dir' from post_blast_tax_check_dmd_ch.collect()

  output:
  path("${merged_contigs_and_singletons}.bx_nr") into merged_blastx_results_ch

  script:

  def diamond_database = "${local_diamond_database_dir}/${params.local_diamond_database_name}"
  def blastx_columns = "qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids sscinames sskingdoms"
                                                                                
  """                                                                           
  # run diamond to do a blastx-style search

  # have to set this environmental variable so diamond can be taxonomically aware 
  export BLASTDB="$blast_tax_dir"

  diamond blastx --threads ${task.cpus} -c 1 -b 8 --db ${diamond_database} --query $merged_contigs_and_singletons --more-sensitive --outfmt 6 $blastx_columns --evalue ${params.max_blasx_nr_evalue} --out ${merged_contigs_and_singletons}.bx_nr
  """
}

/*
  This combines the diamond output file with the individual fasta from each
  dataset that were mushed together for a single diamond search
*/
merged_blastx_results_ch
  .combine(post_blastn_blastx_ch)
  .set{merged_blastx_with_input_ch}


/*
  Split up diamond results into results for each dataset.   Uses blast_from_fasta script.
*/
process split_merged_blastx_results {
  label 'highmem'
  publishDir "${params.blast_out_dir}", mode:'link'

  // singularity info for this process                                          
  if (workflow.containerEngine == 'singularity') {
      container "https://depot.galaxyproject.org/singularity/ubuntu:20.04"
  }          

  input:
  tuple path(merged_blastx_results), val(sample_id), path(contig_weights), path(contigs_and_singletons) from merged_blastx_with_input_ch

  output:
  tuple val(sample_id), path(contig_weights), path("${contigs_and_singletons}.bx_nr") into post_blastx_tally_ch
  tuple val(sample_id), path(contigs_and_singletons), path("${contigs_and_singletons}.bx_nr") into post_blastx_distribute_ch
  tuple val(sample_id), path("${sample_id}_contigs_and_singletons_nn.fasta") into post_blastx_unassigned_ch

  script:
  def blastx_columns = "qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids sscinames sskingdoms"
                                                                                
  """                                                                           
  # split out blastx hits for the contigs and singletons in *this* sample's file
  ${params.scripts_bindir}/blast_from_fasta -f $contigs_and_singletons $merged_blastx_results > ${contigs_and_singletons}.bx_nr.no_header

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

/*
  Use the tally_blast_hits script to tally the number of reads mapping to different taxa
*/
process tally_blastx_results {
  label 'lowmem_nonthreaded'
  publishDir "${params.tally_out_dir}", mode:'link'

  // singularity info for this process                                          
  if (workflow.containerEngine == 'singularity') {
      // container "https://depot.galaxyproject.org/singularity/ubuntu:20.04"
      container "library://stenglein-lab/perl_taxonomy_tools/perl_taxonomy_tools:1.0.0"
  }          

  input:
  tuple val(sample_id), path(contig_weights), path(blast_out) from post_blastx_tally_ch
  path(tax_db) from tax_db_ch_3

  output:
  path("*bx_nr.tally") into blastx_tally_file_ch

  script:
  // TODO: move tally_blast_hits logic into an R script?
  // TODO: tally_blast_hits shouldn't need to do accession->taxid lookups anymore
  """
  ${params.scripts_bindir}/tally_blast_hits -ntd $tax_db -lca -w $contig_weights $blast_out > ${blast_out}.tally
  ${params.scripts_bindir}/tally_blast_hits -ntd $tax_db -lca -w $contig_weights -t -ti ${blast_out} > ${blast_out}.tab_tree.tally
  """
}

/*
   Use the distribute_fasta_by_blast_taxid script to pull out putative virus sequences
*/
process distribute_blastx_results {
  label 'lowmem_nonthreaded'
  publishDir "${params.virus_seq_out_dir}", mode:'link'

  // singularity info for this process                                          
  if (workflow.containerEngine == 'singularity') {
      container "library://stenglein-lab/perl_taxonomy_tools/perl_taxonomy_tools:1.0.0"
  }          

  input:
  tuple val(sample_id), path(contigs_and_singletons), path(blast_out) from post_blastx_distribute_ch
  path(tax_db) from tax_db_ch_4

  output:
  tuple val(sample_id), path("*.fasta_*") optional true

  script:
  // TODO: move distribute logic into an R script?
  """
  # pull out virus-derived reads - this will create a file for each viral taxon
  ${params.scripts_bindir}/distribute_fasta_by_blast_taxid -ntd $tax_db -v $contigs_and_singletons $blast_out
  """
}

/* 
 * We need to deal with the fact that each sample can produce multiple virus
 * fasta files.  This will create a channel that passes sample ID, one fasta
 * and the input fastq to the remap process below
 * see: https://www.nextflow.io/docs/latest/process.html#input-repeaters-each
 *      https://www.nextflow.io/docs/latest/operator.html#operator-combine
 *      https://github.com/nextflow-io/nextflow/issues/440
 */
virus_remap_ch = virus_fasta_ch.transpose().combine(post_host_filter_remap_ch, by: 0)

/*
   Remap host-filtered reads to putative virus contigs
*/
process remap_reads_to_virus_seqs {
  label 'lowmem_threaded'
  publishDir "${params.virus_remap_out_dir}", mode:'link'

  // singularity info for this process
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
      container "https://depot.galaxyproject.org/singularity/bowtie2:2.4.5--py39ha4319a6_1"
  } else {
      container "quay.io/biocontainers/bowtie2:2.4.5--py39ha4319a6_1"
  }

  input:
  tuple val(sample_id), path(virus_fasta), path(input_fastq) from virus_remap_ch

  output:
  tuple val(sample_id), path("*.sam") into virus_sam_ch

  script:
  // handle single-end or paired-end inputs
  def r1 = input_fastq[0]
  def r2 = input_fastq[1] ? input_fastq[1] : ""
  def bowtie_file_input  = input_fastq[1] ? "-1 $r1 -2 $r2" : "-U $r1"

  """
  bowtie2-build $virus_fasta "${virus_fasta}_index"

  bowtie2 \
  -x "${virus_fasta}_index" \
  --local \
  -q \
  --no-unal \
  $bowtie_file_input \
  --sensitive \
  --score-min "C,${params.host_bt_min_score},0" \
  -p ${task.cpus} \
  -S "${virus_fasta}.sam" \
  2> "${sample_id}.host_filtering_bt.log"
  """
}

/*
   Calculate virus remapping coverage stats 
*/
process virus_remapping_coverage_stats {
  label 'lowmem_nonthreaded'
  publishDir "${params.virus_remap_out_dir}", mode:'link'

  // singularity info for this process                                          
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
      container "https://depot.galaxyproject.org/singularity/samtools:1.17--h00cdaf9_0"
  } else {                                                                      
      container "quay.io/biocontainers/samtools:1.17--h00cdaf9_0"
  }     

  input:
  tuple val(sample_id), path(sam) from virus_sam_ch

  output:
  tuple val(sample_id), path("*.txt") 

  // TODO: version with prepended output
  // samtools depth !{bam} | awk '{ print !{sample_id} "\t" $0; }' > !{depth}
  // samtools coverage !{bam} | awk '{ print !{sample_id} "\t" $0; }' > !{cov}

  script:
  def bam   = sam.getName().replaceAll('.sam$', '.bam')
  def depth = sam.getName().replaceAll('.sam$', '.depth.txt')
  def cov   = sam.getName().replaceAll('.sam$', '.cov.txt')
  """
  # create sorted bam from sam
  samtools sort ${sam} > ${bam}

  # create depth file with a prepended column containing sample ID
  samtools depth ${bam} > ${depth}

  # create coverage stats file with a prepended column containing sample ID
  samtools coverage ${bam}  > ${cov}

  """
}

/*
   Output a matrix of # of virus-mapping reads
*/
process  output_virus_mapping_matrix {
  label 'lowmem_nonthreaded'
  publishDir "${params.virus_matrix_out_dir}", mode:'link'

  // singularity info for this process
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
      container "https://depot.galaxyproject.org/singularity/perl:5.26.2"
  } else {
      container "quay.io/biocontainers/perl:5.26.2"
  }

  input:
  path(tally_files) from blastn_tally_file_ch.mix(blastx_tally_file_ch).collect()

  output:
  path("*.txt") 

  script:
  """
  ${params.scripts_bindir}/make_taxa_matrix -v -c ${params.min_matrix_reads} $tally_files > virus_matrix.txt
  """
}

all_tally_ch = genus_tally_ch.mix(family_tally_ch, order_tally_ch, class_tally_ch, kingdom_tally_ch, superkingdom_tally_ch)
   .groupTuple()

/*
   Output matrices of # of reads mapping at different taxonomic levels
*/
process  output_taxa_matrices {
  label 'lowmem_nonthreaded'
  publishDir "${params.virus_matrix_out_dir}", mode:'link'

  // singularity info for this process
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
      container "https://depot.galaxyproject.org/singularity/perl:5.26.2"
  } else {
      container "quay.io/biocontainers/perl:5.26.2"
  }

  input:
  tuple val(rank), path(tally_files) from all_tally_ch

  output:
  path("*.txt") 

  script:
  """
  ${params.scripts_bindir}/make_taxa_matrix -r -c ${params.min_matrix_reads} $tally_files > ${rank}_taxa_matrix.txt
  """
}


