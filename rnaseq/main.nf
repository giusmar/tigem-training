ch_fasta = Channel.fromPath( "/workspace/tigem-training/rnaseq/data/genome/sacCer3.fasta" )
ch_gtf = Channel.fromPath( "/workspace/tigem-training/rnaseq/data/genome/sacCer3.gtf" )
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

process trimming {
    label 'trimming'
    tag 'TRIMGALORE'
    publishDir "$params.outdir" , mode: 'copy',
    saveAs: {filename ->
            if (filename.endsWith("zip"))      "trimming/zip/$filename"    
        else if (filename.endsWith("html"))     "fastqc/html/$filename"
        else null            
    }

    input:
    tuple val(sample_id), path(read)

    output:
    tuple val(sample_id), path("*trimmed.fq.gz"), emit: trimmed

    script:
    """
    trim_galore ${read}
    """
}

process genomegenerate {
    label 'GGenerate'
    tag 'STAR'
    publishDir "$params.outdir" , mode: 'copy',
    saveAs: {filename ->
            if (filename.indexOf("STARindex") > 0)     "STAR/genomegenerate/$filename"
        else null            
    }

    input:
    path(fasta)
    path(gtf)

    output:
    path("STARindex"), emit: starindex

    script:
    """
    mkdir STARindex

    STAR --runMode genomeGenerate \\
    --genomeDir STARindex \\
    --genomeFastaFiles ${fasta} \\
    --sjdbGTFfile ${gtf} \\
    --sjdbOverhang 49
    """
}

process align {
    label 'align'
    tag 'STAR'
    publishDir "$params.outdir" , mode: 'copy',
    saveAs: {filename ->
            if (filename.indexOf("bam") > 0)     "STAR/bams/$filename"
        else null            
    }

    input:
    path(starindex)
    tuple val(sample_id), path(trimmed)

    output:
    tuple val(sample_id), path("*.bam"), emit: sortedBam

    script:
    """
    STAR --genomeDir ${starindex} \\
    --readFilesIn ${trimmed} \\
    --outFileNamePrefix ${sample_id} \\
    --readFilesCommand zcat \\
    --outFilterMultimapNmax 1 \\
    --outReadsUnmapped Fastx \\
    --outSAMtype BAM SortedByCoordinate
    """
}

process featurecounts {
    label 'featurecounts'
    tag 'SUBREAD'
    publishDir "$params.outdir" , mode: 'copy',
    saveAs: {filename ->
            if (filename.indexOf("sam") > 0)     "bwa/align/$filename"
        else null            
    }

    input:
    path(ch_gtf)


    output:
    path("featureCounts_results.txt"), emit: txt
    path("featureCounts_results.log"), emit: txt

    script:
    """
    featureCounts -a ${ch_gtf} -o featureCounts_results.txt . 2> featureCounts_results.log
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

    // Fast quality control on reads
    fastqc(ch_fastq)

    // Trimming on reads
    trimming(ch_fastq)

    // Alignment on reads with STAR
    genomegenerate(ch_fasta,ch_gtf)
    align(genomegenerate.out.starindex.collect(),trimming.out.trimmed)

    // Generate count matrix
    align.out.sortedBam.view()
    //featurecounts()

    // Generate a single report with multiqc
    // report(fastqc.out.qc.collect{ it[1][[1]] },align.out.aligned.collect{ it[1] } )

}
