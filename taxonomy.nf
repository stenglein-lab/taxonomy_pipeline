#!/usr/bin/env nextflow

/*
    Stenglein lab metagenomic classification pipeline 
    implemented in nextflow

    March 9, 2021

    Mark Stenglein
*/


// TODO: command line options and parameter checking


params.fastq_dir = "$baseDir/fastq/"
params.outdir = "$baseDir/results"                                                       
params.initial_fastqc_dir = "${params.outdir}/initial_fastqc/" 
params.post_trim_fastqc_dir = "${params.outdir}/post_trim_fastqc/" 
params.counts_out_dir = "${params.outdir}/fastq_counts/"                        
params.fastq_out_dir = "${params.outdir}/trimmed_fastq/"                        
params.bam_out_dir = "${params.outdir}/bam/"                                    

// ------------------
// Trimming settings
// ------------------
params.always_trim_5p_bases = "0" 
params.always_trim_3p_bases = "1" 
params.post_trim_min_length = "60" 

// --------------------
// Host cell filtering
// --------------------
// Vero cell: African green monkey genome for host filtering
params.host_bt_index = "/home/databases/fly/combined_fly_index"
params.host_bt_suffix = "fly_genome"
params.host_bt_min_score = "80"


// where are R scripts found...
params.R_bindir="${baseDir}/scripts"
params.scripts_bindir="${baseDir}/scripts"

// cd-hit-dup cutoff for collapsing reads with >= this much fractional similarity
params.mismatches_allowed = "2"


// TODO: command line arg processing and validating 

// TODO: check that appropriate refseq files exist (host bt index, nr/nt, etc.)

// TODO: handle fastq.gz compressed files 

/*
 These fastq files represent the main input to this workflow
*/
Channel
    .fromFilePairs("${params.fastq_dir}/*_R{1,2}*.fastq*", size: -1, checkIfExists: true, maxDepth: 1)
    .into {samples_ch_qc; samples_ch_trim; samples_ch_count}


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
  cat !{initial_fastq[0]} | wc -l | xargs bash -c 'echo $(($0 / 4))' | awk '{print "!{sample_id}" "\tinitial\t" $1}' > "!{sample_id}_initial_count.txt"
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
  cat !{trimmed_fastq[0]} | wc -l | xargs bash -c 'echo $(($0 / 4))' |  awk '{print "!{sample_id}" "\tpost_trimming\t" $1}' > "!{sample_id}_trimmed_count.txt"
  '''
}

/*
 Collapse duplicate reads - likely PCR duplicates: duplicated in any case
*/

process collapse_duplicate_reads {
  label 'lowmem_threaded'                                                                

  input:
  tuple val(sample_id), path(input_fastq) from post_trim_ch

  output:
  tuple val(sample_id), path("*_fu.fastq") optional true into post_collapse_count_ch
  tuple val(sample_id), path("*_fu.fastq") optional true into post_collapse_qc_ch
  tuple val(sample_id), path("*_fu.fastq") optional true into post_collapse_ch

  script:

  // this handles paired-end data, in which case must specify a paired output file
  def prefix_param    = input_fastq[1] ? "-u 30" : "-u 50"
  def paired_input    = input_fastq[1] ? "-i2 input_fastq[1]" : ""
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
  cat !{collapsed_fastq[0]} | wc -l | xargs bash -c 'echo $(($0 / 4))' |  awk '{print "!{sample_id}" "\tpost_collapse\t" $1}' > "!{sample_id}_collapsed_count.txt"
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
 Use multiqc to merge post-trimming fastq reports
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
// TODO: make host filtering optional
// TODO: switch to bwa for host filtering too?  The reason to use bowtie2 in
//       instead of bwa is that it has convenient built-in options for 
//       outputting unmapped reads (the --un or --un-conc) options, whereas
//       for bwa you have to go through additional steps with samtools / bedtools 
process host_filtering {
  label 'lowmem_threaded'                                                                

  input:
  tuple val(sample_id), path(input_fastq) from post_collapse_ch
  val("indexes_complete") from post_index_setup_ch

  output:
  // tuple val(sample_id), path("*_fuh.fastq") optional true into post_host_ch_dvg
  tuple val(sample_id), path("*_fuh.fastq") optional true into post_host_ch_count


  // TODO: multiqc analysis of bowtie output (host filtering) (?)

  script:

  // handle single-end or paired-end inputs
  def r1 = input_fastq[0] 
  def r2 = input_fastq[1] ? input_fastq[1] : ""
  def bowtie_file_input  = input_fastq[1] ? "-1 $r1 -2 $r2" : "-U $r1"
  def bowtie_file_output = input_fastq[1] ? "--un-conc ${sample_id}_R%_fuh.fastq" : "--un ${sample_id}_R1_fuh.fastq"

  """
  bowtie2 \
  -x "${params.host_bt_index}" \
  --local \
  -q \
  $bowtie_file_input \
  --sensitive \
  --score-min "C,${params.host_bt_min_score},0" \
  -p ${task.cpus} \
  $bowtie_file_output 2> "${sample_id}.host_filtering_bt.log" > /dev/null 
  """
}


/*
 Count post-host-filtering fastq
*/
process host_filtered_fastq_count {
  label 'lowmem_non_threaded'                                                                
  publishDir "${params.counts_out_dir}", mode:'link'

  input:
  tuple val(sample_id), path(filtered_fastq) from post_host_ch_count

  output:
  path("${sample_id}_host_filtered_count.txt") into post_count_host_ch

  shell:
  // only count the first file because for paired-read data both files
  // will have the same # of reads.
  '''
  cat !{filtered_fastq[0]} | wc -l | xargs bash -c 'echo $(($0 / 4))' |  awk '{print "!{sample_id}" "\tpost_host_filtered\t" $1}' > "!{sample_id}_host_filtered_count.txt"
  '''
}

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
