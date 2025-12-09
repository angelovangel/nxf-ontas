include {DORADO_BASECALL} from './modules/basecalling.nf'
include {DORADO_BASECALL_BARCODING} from './modules/basecalling.nf'


if (params.help) {
        showHelp()
        exit 0
}

// 
def showHelp() {
        log.info """
=============================================
NXF-ALIGNMENT
basecal, align, and analyze ONT data
=============================================

Usage:
    nextflow run anngelovangel/nxf-alignment [--pod5 <path>] [--reads <path>] --ref <ref.fasta> [options]

Required/important options:
    --ref <path>     Reference FASTA (required unless -entry basecall is used)

Input options:
    --pod5 <dir>           Directory with POD5 files (use when basecalling)
    --reads <file|dir>     BAM/FASTQ file or directory of reads (skips basecalling)
    --asfile <file>        Adaptive sampling CSV (filters reads to basecall)

Processing options:
    --model <name>         Dorado basecalling model (default: fast)
    --kit <name>           Barcoding kit name (required with --samplesheet)
    --samplesheet <file>   CSV with columns: sample,barcode (required with --kit)
    --bed <file>           BED file with regions (auto-generated from reference if omitted)

Output & config:
    --outdir <name>        Output directory name (default: results)
    -profile <name>        Nextflow profile (test)
    -entry <name>          Workflow entry point (basecall, use for basecalling only)

""".stripIndent()
}


if (params.kit && !params.samplesheet) {
    error "If --kit is specified, --samplesheet must also be provided."
}

if (!params.kit && params.samplesheet) {
    error "If --samplesheet is provided, --kit must also be specified."
}


// get run info from bam header
process RUN_INFO {
    
    container 'docker.io/aangeloo/nxf-tgs:latest' 
    
    input:
        // Expects one BAM file as input
        path(bam)

    output:
        // Outputs a CSV file with the extracted info
        path("run_info.csv"), emit: ch_runinfo

    script:
    """
    
    RG_LINE=\$(samtools view -H ${bam} | grep '^@RG' | head -n 1) 
    
    FLOWCELL_ID=\$(echo \$RG_LINE | sed -n 's/.*PU:\\([^[:space:]]*\\).*/\\1/p')
    BASECALL_MODEL=\$(echo \$RG_LINE | sed -n 's/.*basecall_model=\\([^[:space:]]*\\).*/\\1/p')
    RUN_DATE=\$(echo "\$RG_LINE" | sed -n 's/.*DT:\\([^T]*\\).*/\\1/p')
    RUN_ID=\$(echo "\$RG_LINE" | sed -n 's/.*DS:runid=\\([^[:space:]]*\\).*/\\1/p')

    echo "flowcell_id,basecall_model,run_date,run_id" > run_info.csv
    echo "\$FLOWCELL_ID,\$BASECALL_MODEL,\$RUN_DATE,\$RUN_ID" >> run_info.csv
    """
}

process MERGE_READS {
    container 'docker.io/aangeloo/nxf-tgs:latest'
    errorStrategy 'ignore' //because some barcodes defined in the samplesheet might be missing in the data
    
    publishDir "$params.outdir/00-basecall/processed", mode: 'copy', pattern: '*{fastq.gz,fastq,bam}'

    input:
    tuple val(samplename), val(barcode), path(bam_pass)
    
    output: 
    path('*{fastq.gz,fastq,bam}')
    
    script:
    """
    samtools cat ${bam_pass}/${barcode}/*.bam > ${samplename}.bam
    """
}

process READ_STATS {
    container 'docker.io/aangeloo/nxf-tgs:latest'
    publishDir "${params.outdir}/00-basecall", mode: 'copy', pattern: '*readstats.tsv'
    tag "${reads.simpleName} - ${reads.extension} file"

    input:
        path(reads)

    output:
        path("*readstats.tsv")

    script:
    
    """
    echo "file\treads\tbases\tn_bases\tmin_len\tmax_len\tn50\tGC_percent\tQ20_percent" > ${reads.simpleName}.readstats.tsv
    
    if [[ ${reads.extension} == bam ]]; then
        samtools fastq ${reads} | faster2 -ts - >> ${reads.simpleName}.readstats.tsv
    else 
    faster2 -ts ${reads} >> ${reads.simpleName}.readstats.tsv
    fi
    """
}

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

process MAKE_BEDFILE {
    container 'docker.io/aangeloo/nxf-tgs:latest'
    
    input:
        path ref

    output:
        path "fallback.bed"

    script:
    """
    samtools faidx ${ref}
    awk -v OFS='\t' '{print \$1, 0, \$2, \$1}' ${ref}.fai > fallback.bed 
    """
}

process BEDTOOLS_COV {
    container 'docker.io/biocontainers/bedtools:v2.27.1dfsg-4-deb_cv1'

    publishDir "${params.outdir}/02-coverage", mode: 'copy', pattern: '*hist.tsv'
    tag "${bam.simpleName}"

    input:
        tuple path(bam), path(bai), path(bed)

    output:
        path "*hist.tsv"

    script:
    """
    echo -e "chr\tstart\tend\tlabel\tdepth\tbases_at_depth\tsize\tpercent_at_depth" > ${bam.simpleName}.hist.tsv
    bedtools coverage -a ${bed} -b ${bam} -hist >> ${bam.simpleName}.hist.tsv
    """
}


process REPORT {
    container 'docker.io/aangeloo/nxf-tgs:latest'

    publishDir "${params.outdir}", mode: 'copy', pattern: '*html'
    
    input:
       tuple path(hist), path(readstats), path(runinfo)
    output:
        path "*.html"

    script:
    """
    make-report.py --hist $hist --readstats $readstats --runinfo $runinfo -o nxf-alignment-report.html
    """
}


workflow basecall {
    ch_pod5 = Channel.fromPath(params.pod5, checkIfExists: true)
    ch_samplesheet = params.samplesheet ? Channel.fromPath(params.samplesheet, checkIfExists: true) : null

    // if no asfile, use dummy placeholder to still do dorado basecalling without as filtering
    ch_decisionfile = params.asfile ? Channel.fromPath(params.asfile, checkIfExists: true) : Channel.fromPath('EMPTY', type: 'file')
    
    if (params.kit) {
        DORADO_BASECALL_BARCODING(ch_decisionfile, ch_pod5)  
        
        ch_samplesheet
        .splitCsv(header:true)
        .filter{ it -> it.barcode =~ /^barcode*/ }
        .map { row -> tuple( row.sample, row.barcode ) }
        .combine( DORADO_BASECALL_BARCODING.out.ch_bam_pass )
        | MERGE_READS 
    } else {
        DORADO_BASECALL(ch_decisionfile, ch_pod5)
    }
    
    emit: 
    ch_bc = params.kit ? MERGE_READS.out : DORADO_BASECALL.out
}

workflow {
    ch_ref = Channel.fromPath(params.ref)

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
        // Otherwise, source the channel from the 'basecall' (or basecall + merge_reads) workflow's output.
        ch_reads = basecall().ch_bc
    }

    // if no bedfile provided, just use the ref to generate one with the fasta entries
    if ( !params.bed ) {
        // generating bedfile from reference
        MAKE_BEDFILE(Channel.fromPath(params.ref, checkIfExists: true))
        ch_bedfile = MAKE_BEDFILE.out
    } else {
        ch_bedfile = Channel.fromPath(params.bed, checkIfExists: true)
    }

    ch_reads
    .first() \
    | RUN_INFO

    ch_reads
    | READ_STATS

    ch_ref \
    .combine( ch_reads ) \
    //.view()
    | DORADO_ALIGN \
    | combine( ch_bedfile ) \
    | BEDTOOLS_COV \
    
    
    BEDTOOLS_COV.out
    .collect()
    .toList()
    .combine( READ_STATS.out.collect().toList() )
    .combine( RUN_INFO.out )
    //.view()
    | REPORT
    
}
