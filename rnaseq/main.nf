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
    path(read)

    output:
    path("${read}_fastqc.{zip,html}"), emit: qc

    script:
    """
    fastqc ${read}
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
    path(fasta)
    path(reads)

    output:
    path('*.sam'), emit: mapped

    script:
    """
    bwa mem -M ${fasta} ${reads[0]} ${reads[1]} -o ${sample_id}_${rep}_aligned_reads.sam
    """
}

workflow {

    // ch_fasta.view()
    // ch_fastq.view()

    fastqc(ch_fastq)


}