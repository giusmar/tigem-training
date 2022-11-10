ch_fasta = Channel.fromPath( "/workspace/tigem-training/rnaseq/data/genome/sacCer3.fasta" )
ch_fastq = Channel.fromPath( "/workspace/tigem-training/rnaseq/data/fastq/*")
                    .map { tuple( it.name.split('_')[0], it ) }

process fastqc {
    label 'fastqc'
    tag 'FASTQC'
    publishDir "$params.outdir" , mode: 'copy',
    saveAs: {filename ->
            if (filename.endsWith("zip"))      "fastqc/zip/$filename"    
        else if (filename.endsWith("html"))     "fastqc/html/$filename"
        else null            
    }

    input:
    tuple val(sample_id), path(read)

    output:
    tuple val(sample_id), path("${sample_id}_fastqc.{zip,html}"), emit: qc

    script:
    """
    ln -s ${read} ${sample_id}.fastq.gz
    fastqc ${sample_id}.fastq.gz
    """
}

process align {
    label 'alignment'
    tag 'BWA'
    publishDir "$params.outdir" , mode: 'copy',
    saveAs: {filename ->
            if (filename.indexOf("sam") > 0)     "bwa/align/$filename"
        else null            
    }

    input:
    tuple val(sample_id), path(read)
    path(fasta)

    output:
    tuple val(sample_id), path("${sample_id}_aligned.sam"), emit: aligned

    script:
    """
    bwa index ${fasta}
    bwa mem ${fasta} ${read} -o ${sample_id}_aligned.sam
    """
}

process report {
    label 'report'
    tag 'MULTIQC'
    publishDir "$params.outdir" , mode: 'copy',
    saveAs: {filename ->
            if (filename.indexOf("html") > 0)     "multiqc/report/$filename"
        else null            
    }

    input:
    path(qc_zip)
    path(aligned)

    output:
    path("multiqc_report.html"), emit: report
    script:
    """
    multiqc . -f -n multiqc_report.html
    """


}

workflow {

    // Uncomment these lines to see how is structured the channels
    // ch_fasta.view()
    // ch_fastq.view()

    // Fast quality control on reads
    fastqc(ch_fastq)

    // Alignment on reads with bwa
    align(ch_fastq,ch_fasta.collect())

    // Generate a single report with multiqc
    ch_multiqc_input = fastqc.out.qc.collect{ it[1][[1]] },align.out.aligned.collect{ it[1] }
    report(ch_multiqc_input)

}