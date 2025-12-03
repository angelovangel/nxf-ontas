/*
* PARAMS
*/


params.pod5 = "${projectDir}/data/pod5"
params.asfile = "${projectDir}/data/as_decisions.csv"
params.model = "fast"

params.reference = "${projectDir}/data/ref.fasta"
params.bedfile = "${projectDir}/data/regions.bed"
params.outdir = "results"

params.reads = null // Parameter for reads when basecalling is skipped

// barcoded case
params.kit = null
params.samplesheet = null

if (params.kit && !params.samplesheet) {
    error "If --kit is specified, --samplesheet must also be provided."
}

// BASECALL

process DORADO_BASECALL {

    //container 'docker.io/nanoporetech/dorado:latest'

    publishDir "${params.outdir}/00-basecall", mode: 'copy'

    input:
        path decisionfile
        path pod5

    output:
        path "reads.bam"

    script:
    """
    awk -F',' '\$2 == "sequence"' '$decisionfile' | cut -f1 -d, > accepted_reads.txt
    nreads_accept=\$(wc -l < accepted_reads.txt)
    nreads_total=\$(wc -l < '$decisionfile')
    echo -e "Found \$nreads_accept out of \$nreads_total reads to basecall " 

    dorado basecaller -l accepted_reads.txt ${params.model} ${pod5} > reads.bam

    # check if reads.bam has reads and exit if no
    nreads=\$(samtools view -c reads.bam)

    if [ "\$nreads" -eq 0 ]; then
        echo "No reads found in reads.bam, exiting." >&2
        exit 1
    fi

    """
}

// ALIGN
process DORADO_ALIGN {

    container 'docker.io/nanoporetech/dorado:latest'

    publishDir "${params.outdir}/01-align", mode: 'copy', pattern: '*{bam,bai}'
    tag "${reads.simpleName}"

    input:
        tuple path(ref), path(reads)

    output:
        path("*.{bam,bai}")

    script:
    """
    dorado aligner ${ref} ${reads} | samtools sort -o ${reads.simpleName}.align.bam
    samtools index ${reads.simpleName}.align.bam
    """
}

// STATS
process SAMTOOLS_BEDCOV {
    container 'docker.io/aangeloo/nxf-tgs:latest'

    publishDir "${params.outdir}/02-coverage", mode: 'copy', pattern: '*tsv'

    input:
        tuple path(bam), path(bai), path(bed)

    output:
        path "*coverage.tsv"

    script:
    """
    echo -e "chr\tstart\tend\tlabel\tbases\tregion_len\tcoverage" > ${bam.simpleName}.coverage.tsv
    samtools bedcov ${bed} ${bam} | awk '{
        region_length = \$3 - \$2
        coverage = (\$5 / region_length)
        print \$0 "\t" region_length "\t" coverage
    }' >> ${bam.simpleName}.coverage.tsv
    """
}

// REPORT

process REPORT {
    container 'docker.io/aangeloo/nxf-tgs:latest'

    publishDir "${params.outdir}/03-report", mode: 'copy', pattern: '*html'
    
    input:
        path coverage_tsv
    output:
        path "*.html"

    script:
    """
    generate-report.R
    """
}

ch_ref = Channel.fromPath(params.reference, checkIfExists: true)
ch_samplesheet = params.samplesheet ? Channel.fromPath(params.samplesheet, checkIfExists: true) : null

workflow basecall {
    ch_pod5 = Channel.fromPath(params.pod5, checkIfExists: true)
    ch_decisionfile = Channel.fromPath(params.asfile, checkIfExists: true)
    
    main:
    DORADO_BASECALL(ch_decisionfile, ch_pod5)

    emit: 
    ch_bc = DORADO_BASECALL.out
}

workflow {
    
    // If 'reads' parameter is provided create a 
    // channel from that path.
    // also possible to pass a folder with reads, every read is one sample
    if (params.reads) {
        
        if ( file(params.reads).isDirectory() ) {
            pattern = "*.{bam,fastq,fastq.gz,fq,fq.gz}"
            ch_reads = Channel.fromPath(params.reads + "/" + pattern, type: 'file', checkIfExists: true)
            
        } else {
            ch_reads = Channel.fromPath(params.reads, checkIfExists: true)        
        }
    } else {
        // Otherwise, source the channel from the 'basecall' workflow's output.
        ch_reads = basecall().ch_bc
    }
    
    ch_ref \
    .combine( ch_reads ) \
    | DORADO_ALIGN \
    | combine( Channel.fromPath(params.bedfile, checkIfExists: true) ) \
    | SAMTOOLS_BEDCOV \
    
    SAMTOOLS_BEDCOV.out
    .collect()
    | REPORT
}

