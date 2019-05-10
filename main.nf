#!/usr/bin/env nextflow

VERSION="0.2"

log.info "===================================================================="
log.info "GATK4 Best Practice Nextflow Pipeline (v${VERSION})                        "
log.info "===================================================================="

params.help = ""
if (params.help) {
  log.info " "
  log.info "USAGE: "
  log.info " "
  log.info "nextflow run oliverSI/GATK4_Best_Practice --fastq1 read_R1.fastq.gz --fastq2 read_R2.fastq.gz"
  log.info " "
  log.info "Mandatory arguments:"
  log.info "    --fastq1        FILE               Fastq(.gz) file for read1"
  log.info "    --fastq2        FILE               Fastq(.gz) file for read2"
  log.info " "
  log.info "Optional arguments:"
  log.info "    --outdir        DIR                Output directory(default: ./Results)"
  log.info "    --samplename    STRING             Sample name(dafault: fastq1 basename)"
  log.info "    --rg            STRING             Read group tag(dafault: fastq1 basename)"
  log.info " "
  log.info "===================================================================="
  exit 1
}

// Validate inputs
params.fasta = params.genome ? params.genomes[ params.genome ].fasta ?: false : false
if (params.fasta) {
    Channel.fromPath(params.fasta)
           .ifEmpty { exit 1, "fasta annotation file not found: ${params.fasta}" }
           .into { fasta_scatter_intervals; fasta_bwa; fasta_baserecalibrator; fasta_haplotypecaller; fasta_genotypegvcfs; fasta_variantrecalibrator_snps; fasta_variantrecalibrator_tranches; fasta_variant_eval; fasta_structural_variantcaller }
}
params.fai = params.genome ? params.genomes[ params.genome ].fai ?: false : false
if (params.fai) {
    Channel.fromPath(params.fai)
           .ifEmpty { exit 1, "fai annotation file not found: ${params.fai}" }
           .into { fai_scatter_intervals; fai_bwa; fai_baserecalibrator; fai_haplotypecaller; fai_genotypegvcfs; fai_variantrecalibrator_snps; fai_variantrecalibrator_tranches; fai_variant_eval; fai_structural_variantcaller }
}
params.dict = params.genome ? params.genomes[ params.genome ].dict ?: false : false
if (params.dict) {
    Channel.fromPath(params.dict)
           .ifEmpty { exit 1, "dict annotation file not found: ${params.dict}" }
           .into { dict_scatter_intervals; dict_bwa; dict_baserecalibrator; dict_haplotypecaller; dict_genotypegvcfs; dict_variantrecalibrator_snps; dict_variantrecalibrator_tranches; dict_variant_eval }
}
params.dbsnp_gz = params.genome ? params.genomes[ params.genome ].dbsnp_gz ?: false : false
if (params.dbsnp_gz) {
    Channel.fromPath(params.dbsnp_gz)
           .ifEmpty { exit 1, "dbsnp annotation file not found: ${params.dbsnp_gz}" }
           .set { dbsnp_gz}
}
params.dbsnp_idx_gz = params.genome ? params.genomes[ params.genome ].dbsnp_idx_gz ?: false : false
if (params.dbsnp_idx_gz) {
    Channel.fromPath(params.dbsnp_idx_gz)
           .ifEmpty { exit 1, "dbsnp_idx_gz annotation file not found: ${params.dbsnp_idx_gz}" }
           .set { dbsnp_idx_gz}
}
params.golden_indel_gz = params.genome ? params.genomes[ params.genome ].golden_indel_gz ?: false : false
if (params.golden_indel_gz) {
    Channel.fromPath(params.golden_indel_gz)
           .ifEmpty { exit 1, "golden_indel_gz annotation file not found: ${params.golden_indel_gz}" }
           .set { golden_indel_gz }
}
params.golden_indel_idx_gz = params.genome ? params.genomes[ params.genome ].golden_indel_idx_gz ?: false : false
if (params.golden_indel_idx_gz) {
    Channel.fromPath(params.golden_indel_idx_gz)
           .ifEmpty { exit 1, "golden_indel_idx_gz annotation file not found: ${params.golden_indel_idx_gz}" }
           .set { golden_indel_idx_gz }
}

params.bwa_index_amb = params.genome ? params.genomes[ params.genome ].bwa_index_amb ?: false : false
if (params.bwa_index_amb) {
    Channel.fromPath(params.bwa_index_amb)
           .ifEmpty { exit 1, "bwa_index_amb annotation file not found: ${params.bwa_index_amb}" }
           .set { bwa_index_amb }
}
params.bwa_index_ann = params.genome ? params.genomes[ params.genome ].bwa_index_ann ?: false : false
if (params.bwa_index_ann) {
    Channel.fromPath(params.bwa_index_ann)
           .ifEmpty { exit 1, "bwa_index_ann annotation file not found: ${params.bwa_index_ann}" }
           .set { bwa_index_ann }
}
params.bwa_index_bwt = params.genome ? params.genomes[ params.genome ].bwa_index_bwt ?: false : false
if (params.bwa_index_bwt) {
    Channel.fromPath(params.bwa_index_bwt)
           .ifEmpty { exit 1, "bwa_index_bwt annotation file not found: ${params.bwa_index_bwt}" }
           .set { bwa_index_bwt }
}
params.bwa_index_pac = params.genome ? params.genomes[ params.genome ].bwa_index_pac ?: false : false
if (params.bwa_index_pac) {
    Channel.fromPath(params.bwa_index_pac)
           .ifEmpty { exit 1, "bwa_index_pac annotation file not found: ${params.bwa_index_pac}" }
           .set { bwa_index_pac }
}
params.bwa_index_sa = params.genome ? params.genomes[ params.genome ].bwa_index_sa ?: false : false
if (params.bwa_index_sa) {
    Channel.fromPath(params.bwa_index_sa)
           .ifEmpty { exit 1, "bwa_index_sa annotation file not found: ${params.bwa_index_sa}" }
           .set { bwa_index_sa }
}
if (params.intervals) {
    Channel.fromPath(params.intervals)
           .ifEmpty { exit 1, "Interval list file for HaplotypeCaller not found: ${params.intervals}" }
           .set { intervals_input }
}
if (params.bai) {
    Channel.fromPath(params.bai)
           .ifEmpty { exit 1, "BAM index file not found: ${params.bai}" }
           .set { bai }
}

// set threadmem equal to total memory divided by number of threads
int threads    = Runtime.getRuntime().availableProcessors()
threadmem      = (((Runtime.getRuntime().maxMemory() * 4) / threads) as nextflow.util.MemoryUnit)

// More memory for RunBam processes
threadmem_more = 4 * threadmem

// Added soft-coded method but hard-coded value of cpu-sage percentage for StructuralVariantCallers process
// ToDo: Expose the hard-coded value as parameter if needed in the future for user to allocate resources at will

// Declaring percentage of total cpus (aka 'threads' var) to be allocated to StructuralVariantCallers process
cpu_percentage = 1

// Multiplying & converting java.math.BigDecimal object to java.lang.Integer
// Check object type with 'my_object.getClass()' method
// More info here: https://www.geeksforgeeks.org/bigdecimal-intvalue-method-in-java/
cpus_to_use_StructVarCall    = (cpu_percentage * threads).intValue()

/*
 * Create a channel for input read files
 * Dump can be used for debugging purposes, e.g. using the -dump-channels operator on run
 */
if (params.reads) {
  Channel
      .fromPath(params.reads)
      .map { file -> tuple(file.baseName, file) }
      .ifEmpty { exit 1, "Cannot find any reads matching: ${reads}\nNB: Path needs to be enclosed in quotes!\nIf this is single-end data, please specify --singleEnd on the command line." }
      .into { reads_fastqc; reads_files }
  reads_files
    .combine(fasta_bwa)
    .dump(tag:'input')
    .set { reads_bwa }
} else if (params.reads_folder){
  reads="${params.reads_folder}/${params.reads_prefix}_{1,2}.${params.reads_extension}"
  Channel
      .fromFilePairs(reads, size: 2)
      .ifEmpty { exit 1, "Cannot find any reads matching: ${reads}\nNB: Path needs to be enclosed in quotes!\nIf this is single-end data, please specify --singleEnd on the command line." }
      .into { reads_fastqc; reads_files}
  reads_files
    .combine(fasta_bwa)
    .dump(tag:'input')
    .set { reads_bwa }
} else if (params.bam) {
  Channel.fromPath(params.bam)
         .map { file -> tuple(file.baseName, file) }
         .ifEmpty { exit 1, "BAM file not found: ${params.bam}" }
         .set { bam_bqsr }
} else {
  exit 1, "Please specify either --reads singleEnd.fastq, --reads_folder pairedReads or --bam myfile.bam"
}

scattered_intervals_ref = fasta_scatter_intervals.merge(fai_scatter_intervals, dict_scatter_intervals)

process ScatterIntervals {
    tag "$fasta"
    container 'broadinstitute/gatk:latest'

    input:
    set file(fasta), file(fai), file(dict) from scattered_intervals_ref

    output:
    file("skip_Ns.interval_list") into scattered_intervals

    script:
    """
    gatk ScatterIntervalsByNs \
    -O skip_Ns.interval_list \
    -R ${fasta} \
    --MAX_TO_MERGE 1000000 \
    --OUTPUT_TYPE ACGT \
    --VERBOSITY ERROR
    """
}

process SkipIntervals {
    tag "$intervals"
    container 'broadinstitute/gatk:latest'

    input:
    file(intervals) from scattered_intervals

    output:
    file("skip_Ns.bed") into skipped_intervals

    script:
    """
    grep -v @ ${intervals}
    #awk '{print \$1":"\$2"-"\$3}'
    awk 'BEGIN { OFS = "" }{ print \$1,":",\$2,"-",\$3 }' ${intervals} > intervals.interval_list

    gatk IntervalListToBed \
    -I ${intervals} \
    -O skip_Ns.bed
    """
}

process MakeIntervals {
    tag "$intervals"
    container 'broadinstitute/gatk:latest'
    publishDir "${params.outdir}", mode: 'copy'

    input:
    file intervals from intervals_input
    file bed from skipped_intervals

    output:
    file("ACGTmers_interval_size.interval_list") into expanded_intervals

    script:
    """
    cut -f 1-3 $bed > ACGTmers.per-interval.bed
    awk '{print \$1, \$2, \$3, \$3-\$2}' $bed | sort -k4nr | cut -f 1-3 > ACGTmers_interval_size.bed
    awk '{print \$1":"\$2+1"-"\$3}' ACGTmers.per-interval.bed > ACGTmers.per-interval.txt

    grep -v @ $intervals | sed 's/:/\\t/' | sed 's/-/\\t/' | cut -f 1-3 > target_temp

    bedtools intersect -a ACGTmers.per-interval.bed -b target_temp -wa | awk '{print \$1":"\$2+1"-"\$3}' | sort -u > ACGTmers_interval_target_intersect.txt

    cp ACGTmers_interval_size.bed ACGTmers_interval_size_temp.txt

    awk 'FNR==NR{a[\$1] = \$1;next}{if (\$1 in a) print \$1}' ACGTmers_interval_size_temp.txt ACGTmers_interval_target_intersect.txt > ACGTmers.per-interval_temp.txt
    awk 'FNR==NR{a[\$1] = \$1;next}{if (\$1 in a) print \$1}' ACGTmers_interval_target_intersect.txt ACGTmers.per-interval.txt > ACGTmers_interval_size.interval_list

    mv ACGTmers.per-interval_temp.txt ACGTmers.per-interval.txt
    """
}

expanded_intervals.into { intervals_file; split_intervals}
    
split_intervals
        .splitText()
        .map { it -> it.trim() }
        .set { interval }

if (!params.bam) {
  
  process fastqc {
    tag "$name"
    publishDir "${params.outdir}/fastqc", mode: 'copy',
        saveAs: {filename -> filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename"}
    container 'flowcraft/fastqc:0.11.7-1'

    input:
    set val(name), file(reads) from reads_fastqc

    output:
    file "*_fastqc.{zip,html}" into fastqc_results

    when: !params.skip_multiqc

    script:
    """
    fastqc -q $reads
    """
  }
  
  process gunzip_dbsnp {
    tag "$dbsnp_gz"

    input:
    file dbsnp_gz from dbsnp_gz
    file dbsnp_idx_gz from dbsnp_idx_gz

    output:
    file "*.vcf" into dbsnp, dbsnp_variantrecalibrator_snps, dbsnp_variantrecalibrator_indels
    file "*.vcf.idx" into dbsnp_idx, dbsnp_idx_variantrecalibrator_snps, dbsnp_idx_variantrecalibrator_indels

    script:
    if ( "${dbsnp_gz}".endsWith(".gz") ) {
     """
     gunzip -d --force $dbsnp_gz
     gunzip -d --force $dbsnp_idx_gz
     """
   } else {
     """
     cp $dbsnp_gz dbsnp.vcf
     cp $dbsnp_idx_gz dbsnp.vcf.idx
     """
   }
  }

  process gunzip_golden_indel {
    tag "$golden_indel_gz"

    input:
    file golden_indel_gz from golden_indel_gz
    file golden_indel_idx_gz from golden_indel_idx_gz

    output:
    file "*.vcf" into golden_indel, golden_indel_variantrecalibrator_indels
    file "*.vcf.idx" into golden_indel_idx, golden_indel_idx_variantrecalibrator_indels

    script:
    if ( "${golden_indel_gz}".endsWith(".gz") ) {
     """
     gunzip -d --force $golden_indel_gz
     gunzip -d --force $golden_indel_idx_gz
     """
    } else {
     """
     cp $golden_indel_gz golden_indel.vcf
     cp $golden_indel_idx_gz golden_indel.vcf.idx
     """
   }
  }

  bwa_index = bwa_index_amb.merge(bwa_index_ann, bwa_index_bwt, bwa_index_pac, bwa_index_sa)
  bwa = reads_bwa.combine(bwa_index)

process BWA {
  tag "$reads"
  container 'kathrinklee/bwa:latest'

  input:
  set val(name), file(reads), file(fasta), file(amb), file(ann), file(bwt), file(pac), file(sa) from bwa

  output:
  set val(name), file("${name}.sam") into sam

  """
  bwa mem -M -R '@RG\\tID:${name}\\tSM:${name}\\tPL:Illumina' $fasta $reads > ${name}.sam
  """
  }


process BWA_sort {
  tag "$sam"
  container 'lifebitai/samtools:latest'

  input:
  set val(name), file(sam) from sam

  output:
  set val(name), file("${name}-sorted.bam") into bam_sort, bam_sort_qc

  """
  samtools sort -o ${name}-sorted.bam -O BAM $sam
  """

}

process RunBamQCmapped {
    tag "$bam"
    container 'maxulysse/sarek:latest'
    memory threadmem_more 
    cpus 4

    input:
    set val(name), file(bam) from bam_sort_qc

    output:
    file("${name}") into bamQCmappedReport

    when: !params.skip_multiqc

    script:
    """
    qualimap \
    bamqc \
    -bam ${bam} \
    --paint-chromosome-limits \
    --genome-gc-distr HUMAN \
    -nt ${task.cpus} \
    -skip-duplicated \
    --skip-dup-mode 0 \
    -outdir ${name} \
    -outformat HTML
    """
}

process MarkDuplicates {
  tag "$bam_sort"
  container 'broadinstitute/gatk:latest'

  input:
  set val(name), file(bam_sort) from bam_sort

  output:
  set val(name), file("${name}_MarkDup.bam") into bam_markdup_baserecalibrator, bam_markdup_applybqsr
  file "metrics.txt" into markdup_multiqc

  """
  gatk MarkDuplicates -I $bam_sort -M metrics.txt -O ${name}_MarkDup.bam
  """

}

baserecalibrator_index = fasta_baserecalibrator.merge(fai_baserecalibrator, dict_baserecalibrator, dbsnp, dbsnp_idx, golden_indel, golden_indel_idx)
baserecalibrator = bam_markdup_baserecalibrator.combine(baserecalibrator_index)

process BaseRecalibrator {
  tag "$bam_markdup"
  container 'broadinstitute/gatk:latest'

  input:
  set val(name), file(bam_markdup), file(fasta), file(fai), file(dict), file(dbsnp), file(dbsnp_idx), file(golden_indel), file(golden_indel_idx) from baserecalibrator

  output:
  set val(name), file("${name}_recal_data.table") into baserecalibrator_table
  file ("*data.table") into baseRecalibratorReport

  """
  gatk BaseRecalibrator \
  -I $bam_markdup \
  --known-sites $dbsnp \
  --known-sites $golden_indel \
  -O ${name}_recal_data.table \
  -R $fasta
  """
}

applybqsr = baserecalibrator_table.join(bam_markdup_applybqsr)

process ApplyBQSR {
  tag "$baserecalibrator_table"
  container 'broadinstitute/gatk:latest'

  input:
  set val(name), file(baserecalibrator_table), file(bam_markdup) from applybqsr

  output:
  set val(name), file("${name}_bqsr.bam") into bam_bqsr

  script:
  """
  gatk ApplyBQSR -I $bam_markdup -bqsr $baserecalibrator_table -O ${name}_bqsr.bam
  """
}
}

if (!params.bai){
  process IndexBam {
  tag "$bam"
  publishDir "${params.outdir}/Bam", mode: 'copy'
  container 'lifebitai/samtools:latest'

  input:
  set val(name), file(bam) from bam_bqsr

  output:
  set val(name), file("ready/${bam}"), file("ready/${bam}.bai") into indexed_bam_bqsr, indexed_bam_qc, indexed_bam_structural_variantcaller

  script:
  """
  mkdir ready
  [[ `samtools view -H ${bam} | grep '@RG' | wc -l`   > 0 ]] && { mv $bam ready;}|| { picard AddOrReplaceReadGroups \
  I=${bam} \
  O=ready/${bam} \
  RGID=${params.rgid} \
  RGLB=${params.rglb} \
  RGPL=${params.rgpl} \
  RGPU=${params.rgpu} \
  RGSM=${params.rgsm};}
  cd ready ;samtools index ${bam};
  """
  }
} else {
  bam_bqsr.merge(bai).into { indexed_bam_bqsr; indexed_bam_qc; indexed_bam_structural_variantcaller }
}


process RunBamQCrecalibrated {
    tag "$bam"
    container 'maxulysse/sarek:latest'
    memory threadmem_more
    cpus 4

    input:
    set val(name), file(bam), file(bai) from indexed_bam_qc

    output:
    file("${name}_recalibrated") into bamQCrecalibratedReport

    when: !params.skip_multiqc

    script:
    """
    qualimap \
    bamqc \
    -bam ${bam} \
    --paint-chromosome-limits \
    --genome-gc-distr HUMAN \
    -nt ${task.cpus} \
    --java-mem-size=${task.memory.toGiga()}G \
    -skip-duplicated \
    --skip-dup-mode 0 \
    -outdir ${name}_recalibrated \
    -outformat HTML
    """
}

haplotypecaller_index = fasta_haplotypecaller.merge(fai_haplotypecaller, dict_haplotypecaller, indexed_bam_bqsr, intervals_file)
haplotypecaller = interval.combine(haplotypecaller_index)

process HaplotypeCaller {
  tag "$interval"
  container 'broadinstitute/gatk:latest'

  memory threadmem

  input:
  set val(interval), file(fasta), file(fai), file(dict), val(sample), file(bam_bqsr), file(bai), file(intervals_file) from haplotypecaller

  output:
  file("${sample}.g.vcf") into haplotypecaller_gvcf
  file("${sample}.g.vcf.idx") into index
  val(sample) into sample

  script:
  int mem = (Runtime.getRuntime().totalMemory()) >> 30
  """
  gatk HaplotypeCaller \
    --java-options -Xmx${task.memory.toMega()}M \
    -R $fasta \
    -O ${sample}.g.vcf \
    -I $bam_bqsr \
    -ERC GVCF \
    -L $interval \
    -isr INTERSECTION \
    --native-pair-hmm-threads 1 \
    -ip 100 \
    --max-alternate-alleles 3 \
    -contamination 0 \
    -L $intervals_file \
    --QUIET
  """
}

process MergeVCFs {
  tag "${name[0]}.g.vcf"
  publishDir "${params.outdir}", mode: 'copy'
  container 'broadinstitute/gatk:latest'

  input:
  file ('*.g.vcf') from haplotypecaller_gvcf.collect()
  file ('*.g.vcf.idx') from index.collect()
  val name from sample.collect()

  output:
  set val("${name[0]}"), file("${name[0]}.g.vcf"), file("${name[0]}.g.vcf.idx") into vcf_bcftools, vcf_variant_eval

  script:
  """
  ## make list of input variant files
  for vcf in \$(ls *vcf); do
    echo \$vcf >> input_variant_files.list
  done

  gatk MergeVcfs \
  --INPUT= input_variant_files.list \
  --OUTPUT= ${name[0]}.g.vcf
  """
}

// Adding structural variant callers  with parliament2
// Input data: -- fasta --fai --bam --bai
// channel with inpu files: fasta, fai from params in the beginning, bam and bai from IndexBam process 

input_structural_variantcaller =  indexed_bam_structural_variantcaller.merge(fasta_structural_variantcaller, fai_structural_variantcaller)

process StructuralVariantCallers {
  tag "$bam"
  container 'lifebitai/parliament2:latest'
  publishDir "${params.outdir}/parliament2", mode: 'copy'

  // Allocate cpus to be utilised in this process
  cpus cpus_to_use_StructVarCall

  input:
  set val(name), file(bam), file(bai), file(fasta), file(fai) from input_structural_variantcaller

  output:
  file("*") into output_structural_variantcaller
  
  when: !params.skip_structural_variants

  script:
  // TODO: --filter_short_contigs (include when using real data) --svviz_only_validated_candidates (both filterings to reduce computations)
  """
  nf_work_dir=\$(pwd)
  cp ${fasta} ref.fa
  cp ${fai} ref.fa.fai
  cp ${bam} input.bam
  cp ${bai} input.bai
  gzip ref.fa

  mv * /home/dnanexus/in
  cd /home/dnanexus

  parliament2.py \
    --bam input.bam \
    --bai input.bai \
    --fai ref.fa.fai \
    --ref_genome ref.fa.gz \
    --prefix ${name} \
    --delly_deletion \
    --delly_insertion \
    --delly_inversion \
    --delly_duplication \
    --breakseq \
    --breakdancer \
    --manta \
    --lumpy \
    --cnvnator \
    --genotype \
    --svviz 


  mv /home/dnanexus/out/* \$nf_work_dir
  """
}


process bcftools{
  tag "$vcf"

  container 'lifebitai/bcftools:latest'

  input:
  set val(name), file(vcf), file(index) from vcf_bcftools
  output:
  file("*") into bcftools_multiqc

  when: !params.skip_multiqc

  script:
  """
  bcftools stats $vcf > bcfstats.txt
  """
}

variant_eval_ref = fasta_variant_eval.merge(fai_variant_eval, dict_variant_eval)
variant_eval = vcf_variant_eval.combine(variant_eval_ref)

process VariantEval {
    tag "$vcf"
    container 'broadinstitute/gatk:latest'

    input:
    set val(name), file(vcf), file(index), file(fasta), file(fai), file(dict) from variant_eval

    output:
    file("${name}.eval.grp") into variantEvalReport

    when: !params.skip_multiqc

    script:
    // TODO: add dbsnp & gold standard
    """
    touch ${name}.eval.grp
    gatk VariantEval \
    -R ${fasta} \
    --eval:${name} $vcf \
    -O ${name}.eval.grp
    """
}


if (!params.skip_fastqc && !params.skip_multiqc) {

  if (!params.bam) {
      fastqc_multiqc = fastqc_results.collect().ifEmpty([])
      multiqc_data = markdup_multiqc.merge(fastqc_multiqc, baseRecalibratorReport, variantEvalReport, bamQCmappedReport, bamQCrecalibratedReport)
      multiqc = bcftools_multiqc.combine(multiqc_data)
  } else {
    multiqc = bcftools_multiqc.combine(variantEvalReport.merge(bamQCrecalibratedReport))
  }

  process multiqc {
    tag "multiqc_report.html"

    publishDir "${params.outdir}/MultiQC", mode: 'copy'
    container 'ewels/multiqc:v1.7'

    input:
    file multiqc from multiqc

    output:
    file("*") into viz

    when: !params.skip_multiqc

    script:
    """
    multiqc . -m fastqc -m qualimap -m picard -m gatk -m bcftools 
    """
  }
}

