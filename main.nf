/*
* PARAMS
*/

params.pod5 = "${projectDir}/data/pod5"
params.asfile = "${projectDir}/data/as_decisions.csv"
params.model = "fast"

params.reference = "${projectDir}/data/ref.fasta"
params.bedfile = "${projectDir}/data/regions.bed"
params.outdir = "${projectDir}/results"

// BASECALL

process DORADO_BASECALL {

    container 'docker.io/nanoporetech/dorado:latest'

    publishDir "${params.outdir}/basecall-${params.model}", mode: 'copy'

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
    """
}

// ALIGN
process DORADO_ALIGN {

    container 'docker.io/nanoporetech/dorado:latest'

    publishDir "${params.outdir}/align-${params.model}", mode: 'copy', pattern: '*{bam,bai}'

    input:
        path ref 
        path bamfile

    output:
        tuple path('align.bam'), path('align.bam.bai')

    script:
    """
    dorado aligner ${ref} ${bamfile} | samtools sort -o align.bam
    samtools index align.bam
    """
}

// STATS
process SAMTOOLS_COV {
    container 'docker.io/aangeloo/nxf-tgs:latest'

    publishDir "${params.outdir}/coverage-${params.model}", mode: 'copy', pattern: '*tsv'

    input:
        path bed
        tuple path(bam), path(bai)

    output:
        path "coverage.tsv"

    script:
    """
    echo -e "chr\tstart\tend\tlabel\tbases\tcoverage" > coverage.tsv
    samtools bedcov ${bed} ${bam} | awk '{
        region_length = \$3 - \$2
        coverage = (\$5 / region_length)
        print \$0 "\t" coverage
    }' >> coverage.tsv
    """
}

workflow {
    pod5_ch = Channel.fromPath(params.pod5, checkIfExists: true)
    decisionfile_ch = Channel.fromPath(params.asfile, checkIfExists: true)
    ref_ch = Channel.fromPath(params.reference, checkIfExists: true)
    
    DORADO_BASECALL(decisionfile_ch, pod5_ch)
    DORADO_ALIGN(ref_ch, DORADO_BASECALL.out)

    SAMTOOLS_COV(params.bedfile, DORADO_ALIGN.out)
}
