#!/usr/bin/env nextflow

process REPORT {
    container 'docker.io/aangeloo/nxf-tgs:latest'

    publishDir "${params.outdir}", mode: 'copy', pattern: '*html'
    
    // some inputs are optional, so better as separated inputs
    input:
        path(runinfo) 
        path(wf_props)
        path(readstats)
        path(ref_stats)
        path(hist) 
        path(bedcov)
        path(bedcov_complement) 
        path(flagstat)
        
    output:
        path "*.html"

    script:
    // Check if the file name contains 'empty'
    // If it's empty, we send an empty string; otherwise, we send the full argument
    def hist_arg = hist.name.contains('empty_hist') ? '' : "--hist ${hist}"
    def bedcov_arg = bedcov.name.contains('empty_bedcov') ? '' : "--bedcov ${bedcov}"
    def bedcov_compl_arg = bedcov_complement.name.contains('empty_bedcov_compl') ? '' : "--bedcov-compl ${bedcov_complement}"
    def flagstat_arg = flagstat.name.contains('empty_flagstat') ? '' : "--flagstat ${flagstat}"
    """
    make-report.py \
        --runinfo $runinfo \
        --wfinfo $wf_props \
        --readstats $readstats \
        --refstats $ref_stats \
        $hist_arg \
        $bedcov_arg \
        $bedcov_compl_arg \
        $flagstat_arg \
        -o nxf-alignment-report.html
    """
}