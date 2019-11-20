#!/usr/bin/env nextflow

// DSL2 BRANCH
nextflow.preview.dsl=2

// PRINT HELP AND EXIT
if(params.help){
    println """\

         ===============================================
          E P I D I V E R S E - S N P   P I P E L I N E
         ===============================================
         ~ version ${workflow.manifest.version}

         Usage: 
              nextflow run epidiverse/snp [OPTIONS]...

         Options: GENERAL
              --input [path/to/input/dir]     [REQUIRED] Specify the path to the directory containing each sample output
                                          from the wgbs pipeline to be taken forward for analysis. All the subdirectories must
                                          correspond to sample names in the provided samples file, and contain within them a
                                          bedGraph directories with files in '*.bedGraph' format.

              --reference [path/to/ref.fa]    Path to the input reference genome file in fasta format. Required for the 
                                          variant calling aspect of the pipeline, along with a valid fasta index *.fai file.

              --output [STR]                  A string that will be used as the name for the output results directory, which
                                          will be generated in the working directory. This directory will contain
                                          sub-directories for each set of reads analysed during the pipeline.
                                          [default: dmrs]


         Options: MODIFIERS
              --clusters                      Specify a string that corresponds to a group name in the provided "samples.tsv",
                                          and the pipeline will run DMR comparisons for each group relative to this group.
                                          Otherwise, the pipeline will run all possible pairwise comparisons if no control
                                          group is specified. [default: off]

              --variants                      Specify a string that corresponds to a group name in the provided "samples.tsv",
                                          and the pipeline will run DMR comparisons for each group relative to this group.
                                          Otherwise, the pipeline will run all possible pairwise comparisons if no control
                                          group is specified. [default: off]


         Options: ADDITIONAL
              --help                          Display this help information and exit
              --version                       Display the current pipeline version and exit
              --debug                         Run the pipeline in debug mode    


         Example: 
              nextflow run epidiverse/snp \
              --input /path/to/wgbs/dir \
              --reference /path/to/ref.fa
              --output snps

    """
    ["bash", "${baseDir}/bin/clean.sh", "${workflow.sessionId}"].execute()
    exit 0
}

// PRINT VERSION AND EXIT
if(params.version){
    println """\
         ===============================================
          E P I D I V E R S E - S N P   P I P E L I N E
         ===============================================
         ~ version ${workflow.manifest.version}
    """
    ["bash", "${baseDir}/bin/clean.sh", "${workflow.sessionId}"].execute()
    exit 0
}


// DEFINE PATHS
bam_path = "${params.input}/*/*.bam"

// conditionals for setting --clusters and --variants
if( !params.clusters ){

    variants = true
    clusters = false
} else if( params.variants ){

    variants = true
    clusters = true
} else {

    variants = false
    clusters = true
}

// check reference
if( variants ){

    fasta = file("${params.reference}", checkIfExists: true, glob: false)
    fai = file("${params.reference}.fai", checkIfExists: true, glob: false)
} else {

    fasta = false
    fai = false
}

// PRINT STANDARD LOGGING INFO
log.info ""
log.info "         ================================================"
log.info "          E P I D I V E R S E - S N P    P I P E L I N E"
if(params.debug){
log.info "         (debug mode enabled)"
log.info "         ================================================" }
else {
log.info "         ================================================" }
log.info "         ~ version ${workflow.manifest.version}"
log.info ""
log.info "         input dir     : ${params.input}"
log.info "         reference     : ${params.reference ? "${params.reference}" : "-"}"
log.info "         output dir    : ${params.output}"
log.info "         variant calls : ${variants ? "enabled" : "disabled"}"
log.info "         clustering    : ${clusters ? "enabled" : "disabled"}"
log.info ""
log.info "         ================================================"
log.info "         RUN NAME: ${workflow.runName}"
log.info ""



/////////////////////
// COMMON CHANNELS //
/////////////////////

// STAGE BEDGRAPH CHANNELS FROM TEST PROFILE
if ( workflow.profile.tokenize(",").contains("test") ){

        include check_test_data from './libs/functions.nf' params(BAMPaths: params.BAMPaths)
        BAM = check_test_data(params.BAMPaths)

} else {

    // STAGE BAM CHANNELS
    BAM = Channel.fromFilePairs(bam_path, size: 1)
        .ifEmpty{ exit 1, "ERROR: cannot find valid *.bam files in dir: ${params.input}\n"}
        .map{it.flatten()}

}


////////////////////
// BEGIN PIPELINE //
////////////////////

// INCLUDES
include './libs/snp.nf' params(params)

// WORKFLOWS

// WGBS workflow - primary pipeline
workflow 'SNPS' {

    get:
        BAM
        fasta
        fai
 
    main:
        // samtools sort + index
        preprocessing(BAM, fasta)

        // fork preprocessing.out for clustering and variant calling workflows
        cluster_channel = clusters ? preprocessing.out.map{tuple("clustering", *it)} : Channel.empty()
        variant_channel = variants ? preprocessing.out.map{tuple("variants", *it)} : Channel.empty()

        // position masking
        masking(cluster_channel.mix(variant_channel))

        // clustering workflow
        extracting(masking.out.filter{ it[0] == "clustering" })
        khmer(extracting.out)
        kwip(khmer.out.collect{ it[0] }.mix(khmer.out.collect{ it[1] }).toList())
        clustering(kwip.out)

        // variant calling workflow
        sorting(masking.out.filter{ it[0] == "variants" })
        freebayes_parallel(sorting.out, fasta, fai)
        bcftools(freebayes_parallel.out)

        //freebayes(sorting.out)
        //bcftools(freebayes.out)

        // haplotyping workflow
        //extractHAIRS(sorting.out.combine(freebayes.out, by: 0))
        //HAPCUT2(extractHAIRS.out.combine(freebayes.out, by: 0))
        //bamsplit(sorting.out.combine(HAPCUT2.out, by: 0))
        //MethylDackel(bamsplit.out.transpose())

    emit:
        khmer_publish = khmer.out
        kwip_publish = kwip.out
        clustering_publish = clustering.out
        vcf_unphased = freebayes_parallel.out
        //vcf_unphased = freebayes.out
        vcf_filtered = bcftools.out
        //vcf_phased = HAPCUT2.out
        //bam_haplotypes = bamsplit.out
        //bedGraphs = MethylDackel.out

}


// MAIN workflow
workflow {

    main:
        SNPS(BAM, fasta, fai)

    publish:
        SNPS.out.khmer_publish to: "${params.output}/hashes", mode: 'copy'
        SNPS.out.kwip_publish to: "${params.output}", mode: 'copy'
        SNPS.out.clustering_publish to: "${params.output}", mode: 'move'
        SNPS.out.vcf_unphased to: "${params.output}/vcf", mode: 'copy'

}


//////////////////
// END PIPELINE //
//////////////////



// WORKFLOW TRACING
workflow.onError {
    log.info "Oops... Pipeline execution stopped with the following message: ${workflow.errorMessage}"
}

workflow.onComplete {

    log.info ""
    log.info "         Pipeline execution summary"
    log.info "         ---------------------------"
    log.info "         Name         : ${workflow.runName}${workflow.resume ? " (resumed)" : ""}"
    log.info "         Profile      : ${workflow.profile}"
    log.info "         Launch dir   : ${workflow.launchDir}"    
    log.info "         Work dir     : ${workflow.workDir} ${params.debug ? "" : "(cleared)" }"
    log.info "         Status       : ${workflow.success ? "success" : "failed"}"
    log.info "         Error report : ${workflow.errorReport ?: "-"}"
    log.info ""

    if (params.debug == false && workflow.success) {
        ["bash", "${baseDir}/bin/clean.sh", "${workflow.sessionId}"].execute() }
}