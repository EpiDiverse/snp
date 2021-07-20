#!/usr/bin/env nextflow

// DSL2 BRANCH
nextflow.enable.dsl=2

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
                                          correspond to sample names, and contain within them files in *.bam format.

              --reference [path/to/ref.fa]    Path to the input reference genome file in fasta format. REQUIRED for the 
                                          variant calling aspect of the pipeline, along with a valid fasta index *.fai file.

              --output [STR]                  A string that will be used as the name for the output results directory, which
                                          will be generated in the working directory [default: snps]


         Options: MODIFIERS
              --variants                      Specify to enable variant calling from bisulfite sequencing data. If no MODIFIER
                                          is given, all will run by default. [default: off]

              --phase                         Specify to enable phasing of variants from bisulfite sequencing data. This option 
                                          will also enable --variants. If no MODIFIER is given, all will run by default.
                                          [default: off]

              --clusters                      Specify to enable sample clustering from bisulfite sequencing data. If no MODIFIER
                                          is given, all will run by default. [default: off]


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

// VALIDATE PARAMETERS
ParameterChecks.checkParams(params)

// DEFINE PATHS
bam_path = "${params.input}/*.bam"

// conditionals for setting --variants, --phase and --clusters
if( !params.clusters && !params.phase && !params.variants ){
    variants = true
    clusters = true
    phasings = true
} else {
    variants = (params.variants || params.phase)
    phasings = params.phase
    clusters = params.clusters
}

// check reference
fasta = file("${params.reference}", checkIfExists: true, glob: false)
fai = file("${params.reference}.fai", checkIfExists: true, glob: false)

// determine context for ASM
if ((params.noCpG == true) && (params.noCHH == true) && (params.noCHG == true)) {error "ERROR: please specify methylation context for analysis"}
else if ((params.noCpG == true) && (params.noCHH == true) && (params.noCHG == false)) {context = "--noCpG --CHG "}
else if ((params.noCpG == true) && (params.noCHH == false) && (params.noCHG == true)) {context = "--noCpG --CHH "}
else if ((params.noCpG == true) && (params.noCHH == false) && (params.noCHG == false)) {context = "--noCpG --CHH --CHG "}
else if ((params.noCpG == false) && (params.noCHH == true) && (params.noCHG == true)) {context = " "}
else if ((params.noCpG == false) && (params.noCHH == true) && (params.noCHG == false)) {context = "--CHG "}
else if ((params.noCpG == false) && (params.noCHH == false) && (params.noCHG == true)) {context = "--CHH "}
else {context = "--CHH --CHG "}


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
log.info "         input dir      : ${params.input}"
log.info "         reference      : ${params.reference ? "${params.reference}" : "-"}"
log.info "         output dir     : ${params.output}"
log.info "         clustering     : ${clusters ? "enabled" : "disabled"}"
log.info "         variant calls  : ${variants ? "enabled" : "disabled"}"
log.info "         vcf phasing    : ${phasings ? "enabled" : "disabled"}"
if (phasings) {
log.info "         asm context(s) : ${params.noCpG ? "" : "CpG " }${params.noCHH ? "" : "CHH " }${params.noCHG ? "" : "CHG" }"
}
log.info ""
log.info "         ================================================"
log.info "         RUN NAME: ${workflow.runName}"
log.info ""



/////////////////////
// COMMON CHANNELS //
/////////////////////

// STAGE BEDGRAPH CHANNELS FROM TEST PROFILE
if ( workflow.profile.tokenize(",").contains("test") ){

        include {check_test_data} from './lib/functions.nf' params(BAMPaths: params.BAMPaths)
        BAM = check_test_data(params.BAMPaths)

} else {

    // STAGE BAM CHANNELS
    BAM = Channel.fromPath(bam_path)
        .ifEmpty{ exit 1, "ERROR: cannot find valid *.bam files in dir: ${params.input}\n"}
        .map{ tuple(it.baseName, it) }
        .take(params.take.toInteger())

}

// Error handling when clustering with fewer than three samples
if( clusters ){

    BAM
        .count()
        .subscribe{int c ->
            if( c <= 2 ){
                error "ERROR: clustering is only possible with a minimum of three samples"
                exit 1
            }
        }
}

////////////////////
// BEGIN PIPELINE //
////////////////////

// INCLUDES
include {preprocessing;masking;extracting;khmer;kwip;clustering;sorting;freebayes;bcftools;plot_vcfstats} from './lib/snp.nf' params(params)

// WORKFLOWS

// WGBS workflow - primary pipeline
workflow 'SNPS' {

    take:
        BAM
        fasta
        fai
        context
 
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
        khmer(extracting.out[0])
        kwip(khmer.out.collect())
        clustering(kwip.out)

        // variant calling workflow
        sorting(masking.out.filter{ it[0] == "variants" })
        freebayes(sorting.out[0], fasta, fai)
        if (phasings) { bcftools(WhatsHap_phase.out.groupTuple()) }
        else { bcftools(freebayes.out) }
        plot_vcfstats(bcftools.out)

        // haplotyping workflow (read-based)
        split_scaffolds(freebayes.out)
        WhatsHap_phase(masking.out.filter{ it[0] == "variants" }.map{ it.tail() }.combine(split_scaffolds.out.transpose(), by:0))
        WhatsHap_haplotag(masking.out.filter{ it[0] == "variants" }.map{ it.tail() }.join(bcftools.out))
        WhatsHap_split(preprocessing.out.join(WhatsHap_haplotag.out))

        // allele-specific methylation calling
        MethylDackel(WhatsHap_split.transpose(), fasta, fai, context)


    /*
    emit:
        bam_variants = sorting.out[1]
        //bam_clusters = extracting.out[1]
        khmer_publish = khmer.out
        kwip_publish = kwip.out
        clustering_publish = clustering.out
        vcf_unphased = freebayes.out
        vcf_filtered = bcftools.out
        vcf_vcfstats = plot_vcfstats.out
        
    */

}


// MAIN workflow
workflow {

    main:
        SNPS(BAM, fasta, fai, context)

    /*
    publish:
        //SNPS.out.bam_clusters to: "${params.output}/bam/clusters", mode: 'move'
        SNPS.out.bam_variants to: "${params.output}/bam/variants", mode: 'copy'
        SNPS.out.khmer_publish to: "${params.output}/hashes", mode: 'copy'
        SNPS.out.kwip_publish to: "${params.output}", mode: 'copy'
        SNPS.out.clustering_publish to: "${params.output}", mode: 'move'
        SNPS.out.vcf_unphased to: "${params.output}/vcf", mode: 'copy'
        SNPS.out.vcf_vcfstats to: "${params.output}/stats", mode: 'move'
    */

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
    log.info "         Work dir     : ${workflow.workDir} ${workflow.success && !params.debug ? "(cleared)" : ""}"
    log.info "         Status       : ${workflow.success ? "success" : "failed"}"
    log.info "         Error report : ${workflow.errorReport ?: "-"}"
    log.info ""

    if (workflow.success && !params.debug) {
        ["bash", "${baseDir}/bin/clean.sh", "${workflow.sessionId}"].execute() }
}