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
    tuple val(sample_id), path("${read}_fastqc.{zip,html}"), emit: qc

    script:
    """
    ln -s ${reads} ${sample_id}.fastq.gz
    fastqc ${sample_id}.fastq.gz
    """
}

process align {
    echo true
    label 'alignment'
    tag 'BWA'
    publishDir "$params.outdir" , mode: 'copy',
    saveAs: {filename ->
            if (filename.indexOf("sam") > 0)     "align/mapped/$filename"
        else null            
    }

    input:
    tuple val(sample_id), path(read)
    path(fasta)

    output:
    tuple val(sample_id), path("${sample_id}_aligned_reads.sam"), emit: aligned

    script:
    """
    bwa mem ${fasta} ${read} -o ${sample_id}_aligned_reads.sam
    """
}

workflow {

    // ch_fasta.view()
    // ch_fastq.view()

    fastqc(ch_fastq)
    align(ch_fastq,ch_fasta)


}