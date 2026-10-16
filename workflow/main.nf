#!/usr/bin/env nextflow
/*
 * Upstream RNA-seq workflow: raw FASTQ -> QC -> trimming -> HISAT2 alignment
 * -> QualiMap -> HTSeq-count, with a MultiQC report.
 * Same steps and parameters as scripts/01-07, run in parallel per sample.
 *
 *   nextflow run workflow/main.nf -profile conda
 *   nextflow run workflow/main.nf -profile conda --samples config/samples.tsv --outdir results_nf
 */
nextflow.enable.dsl = 2

process FASTQC {
    tag "${sample_id}:${stage}"
    publishDir "${params.outdir}/fastqc_${stage}", mode: 'copy'

    input:
    tuple val(sample_id), path(reads), val(stage)

    output:
    path "*_fastqc.{zip,html}", emit: reports

    script:
    """
    fastqc --threads ${task.cpus} ${reads}
    """

    stub:
    """
    touch ${sample_id}_${stage}_fastqc.zip ${sample_id}_${stage}_fastqc.html
    """
}

process TRIM_GALORE {
    tag "${sample_id}"
    publishDir "${params.outdir}/trimmed", mode: 'copy', pattern: '*_trimming_report.txt'

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}_trimmed.fq.gz"), emit: reads
    path "*_trimming_report.txt",                             emit: report

    script:
    """
    trim_galore --cores ${task.cpus} --basename ${sample_id} ${reads}
    """

    stub:
    """
    touch ${sample_id}_trimmed.fq.gz ${reads}_trimming_report.txt
    """
}

process HISAT2_ALIGN {
    tag "${sample_id}"
    publishDir "${params.outdir}/alignments", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)
    path index_files

    output:
    tuple val(sample_id), path("${sample_id}.sorted.bam"), path("${sample_id}.sorted.bam.bai"), emit: bam
    path "${sample_id}.hisat2.log",                                                             emit: log

    script:
    """
    hisat2 -p ${task.cpus} -k ${params.hisat2_k} -x ${params.hisat2_index_name} -U ${reads} \\
        2> ${sample_id}.hisat2.log \\
        | samtools sort -@ ${task.cpus} -o ${sample_id}.sorted.bam -
    samtools index ${sample_id}.sorted.bam
    """

    stub:
    """
    touch ${sample_id}.sorted.bam ${sample_id}.sorted.bam.bai ${sample_id}.hisat2.log
    """
}

process QUALIMAP_RNASEQ {
    tag "${sample_id}"
    publishDir "${params.outdir}/qualimap", mode: 'copy'

    input:
    tuple val(sample_id), path(bam), path(bai)
    path gtf

    output:
    path "${sample_id}", emit: results

    script:
    """
    qualimap rnaseq -bam ${bam} -gtf ${gtf} --java-mem-size=${params.qualimap_mem} -outdir ${sample_id}
    """

    stub:
    """
    mkdir ${sample_id}
    """
}

process HTSEQ_COUNT {
    tag "${sample_id}"
    publishDir "${params.outdir}/counts", mode: 'copy'

    input:
    tuple val(sample_id), path(bam), path(bai)
    path gtf

    output:
    path "${sample_id}.counts.tsv", emit: counts

    script:
    """
    htseq-count --format bam --order pos --stranded ${params.stranded} --idattr gene_id \\
        ${bam} ${gtf} > ${sample_id}.counts.tsv
    """

    stub:
    """
    touch ${sample_id}.counts.tsv
    """
}

process COUNT_SUMMARY {
    publishDir "${params.outdir}/counts", mode: 'copy'

    input:
    path counts
    path script

    output:
    path "htseq_assignment_summary.*"

    script:
    """
    python ${script} ${counts} --out htseq_assignment_summary.png --table htseq_assignment_summary.tsv
    """

    stub:
    """
    touch htseq_assignment_summary.png htseq_assignment_summary.tsv
    """
}

process MULTIQC {
    publishDir "${params.outdir}/multiqc", mode: 'copy'

    input:
    path '*'

    output:
    path "multiqc_report.html"

    script:
    """
    multiqc .
    """

    stub:
    """
    touch multiqc_report.html
    """
}

workflow {
    samples = Channel
        .fromPath(params.samples)
        .splitCsv(header: true, sep: '\t')
        .map { row -> tuple(row.run_accession, file("${params.reads_dir}/${row.run_accession}.fastq.gz", checkIfExists: true)) }

    gtf         = file(params.gtf, checkIfExists: true)
    index_files = Channel.fromPath("${params.hisat2_index_dir}/${params.hisat2_index_name}*.ht2").collect()

    trimmed = TRIM_GALORE(samples)

    // FastQC on raw and trimmed reads (a process can only be invoked once per workflow)
    fastqc = FASTQC(
        samples.map { id, reads -> tuple(id, reads, 'raw') }
            .mix(trimmed.reads.map { id, reads -> tuple(id, reads, 'trimmed') })
    )

    aligned = HISAT2_ALIGN(trimmed.reads, index_files)

    QUALIMAP_RNASEQ(aligned.bam, gtf)
    counts = HTSEQ_COUNT(aligned.bam, gtf)
    COUNT_SUMMARY(counts.counts.collect(), file("${projectDir}/../scripts/07_plot_count_summary.py"))

    MULTIQC(
        fastqc.reports
            .mix(trimmed.report, aligned.log, counts.counts)
            .collect()
    )
}
