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
              --clusters                      Specify to enable sample clustering from bisulfite sequencing data. [default: off]

              --variants                      Specify to enable variant calling from bisulfite sequencing data. [default: off]

              --phase                         Specify to enable phasing of variants from bisulfite sequencing data. This option 
                                          will also enable --variants. [default: off]


         Options: ADDITIONAL
              --help                          Display this help information and exit
              --version                       Display the current pipeline version and exit
              --debug                         Run the pipeline in debug mode    
              --take [INT]                    Limit analysis of input samples to INT. [default: 10]

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
variants = (params.variants || params.phase)
phasings = params.phase
clusters = params.clusters


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

include {
    preprocessing;
    masking;
    extracting;
    khmer;
    kwip;
    clustering;
    sorting;
    freebayes;
    bcftools;
    plot_vcfstats;
    split_scaffolds;
    WhatsHap_phase;
    WhatsHap_haplotag;
    WhatsHap_split;
    MethylDackel;
} from './lib/snp.nf' params(params)

// WORKFLOWS

// SNPS workflow - primary pipeline
workflow 'SNPS' {

    take:
        BAM
        fasta
        fai
 
    main:
        // samtools sort + index
        preprocessing(BAM, fasta)

        // fork preprocessing.out for clustering and variant calling workflows
        cluster_channel = clusters ? preprocessing.out.map{tuple("clustering", *it)} : Channel.empty()
        variant_channel = variants || (!params.variants && !params.clusters && !params.phase) ? preprocessing.out.map{tuple("variants", *it)} : Channel.empty()

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
        
        if(!phasings){ 
        bcftools(freebayes.out)
        plot_vcfstats(bcftools.out) }

    emit:
        preprocessing_out = preprocessing.out
        sorting_out = sorting.out[0]
        freebayes_out = freebayes.out

}


// ASM workflow - secondary pipeline
workflow 'ASM' {

    take:
        preprocessing_out
        sorting_out
        freebayes_out
        fasta
        fai
        context
 
    main:
        // haplotyping workflow (read-based)
        split_scaffolds(freebayes_out)
        WhatsHap_phase(sorting_out.combine(split_scaffolds.out.transpose(), by:0), fasta, fai)
        
        bcftools(WhatsHap_phase.out.groupTuple())
        plot_vcfstats(bcftools.out)
        
        WhatsHap_haplotag(sorting_out.join(bcftools.out))
        WhatsHap_split(preprocessing_out.join(WhatsHap_haplotag.out), fasta, fai)

        // allele-specific methylation calling
        MethylDackel(WhatsHap_split.out.transpose(), fasta, fai, context)

}


// MAIN workflow
workflow {

    main:
        SNPS(BAM, fasta, fai)
        if(phasings){ ASM(SNPS.out.preprocessing_out, SNPS.out.sorting_out, SNPS.out.freebayes_out, fasta, fai, context) }

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