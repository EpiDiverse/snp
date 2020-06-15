#!/usr/bin/env nextflow


// taking input bam files for sorting and indexing
process "preprocessing" {

    label "${params.high ? "high" : "low"}"
    label "finish"
    tag "$sample"

    maxForks "${params.fork}".toInteger()

    input:
    tuple sample, path(bam)
    // eg. [sample, /path/to/sample.bam]
    path fasta

    output:
    tuple sample, path("calmd.bam"), path("calmd.bam.bai")
    // eg. [sample, [/path/to/calmd.bam, /path/to/calmd.bam.bai]]

    script:
    """
    samtools sort -T deleteme -m ${((task.memory.getBytes() / task.cpus) * 0.9).round(0)} -@ ${task.cpus} \\
    -o sorted.bam ${bam} || exit \$?
    samtools calmd -b sorted.bam ${fasta} 1> calmd.bam 2> /dev/null && rm sorted.bam
    samtools index calmd.bam
    """
}


// mask genomic or epigenomic variation according to type
process "masking" {

    label "${params.high ? "high" : "low"}"
    label "finish"
    tag "$sample - $type"

    maxForks "${params.fork}".toInteger()

    input:
    tuple type, sample, path(bam), path(bai)
    // eg. [clustering, sample, /path/to/sample.bam, /path/to/sample.bam.bai]

    output:
    tuple type, sample, path("${type}.bam")
    // eg. [clustering, sample, /path/to/clustering.bam]

    script:
    """
    change_sam_queries.py -Q -T ${task.cpus} -t . ${type == "clustering" ? "-G " : ""}${bam} ${type}.bam || exit \$?
    find -mindepth 1 -maxdepth 1 -type d -exec rm -r {} \\;
    """
}


// extract fastq reads from genomic-masked samples
process "extracting" {

    label "${params.high ? "high" : "low"}"
    label "finish"
    tag "$sample"

    maxForks "${params.fork}".toInteger()

    input:
    tuple type, sample, path(bam)
    // eg. [clustering, sample, /path/to/sample.bam]

    output:
    tuple sample, path("${sample}.fastq.gz")
    // eg. [sample, /path/to/sample.fastq.gz]
    //path "${sample}.bam"

    when:
    params.clusters || (!params.variants && !params.clusters)

    script:
    """
    samtools sort -T deleteme -m ${((task.memory.getBytes() / task.cpus) * 0.9).round(0)} -@ ${task.cpus} \\
    -no ${sample}.bam ${bam} || exit \$?
    samtools fastq ${sample}.bam | gzip -c > ${sample}.fastq.gz && rm ${sample}.bam
    """
}


// run khmer on extracted fastq reads from genomic-masked samples
process "khmer" {

    label "${params.high ? "high" : "low"}"
    label "finish"
    tag "$sample"

    maxForks "${params.fork}".toInteger()

    input:
    tuple sample, path(fastq)
    // eg. [sample, /path/to/sample.fastq.gz]

    output:
    path "${sample}.ct.gz"
    // eg. [/path/to/sample.ct.gz]

    when:
    params.clusters || (!params.variants && !params.clusters)

    script:
    """
    load-into-counting.py -T ${task.cpus} -N 1 -x 1e9 -k 20 -b -f -s tsv ${sample}.ct.gz ${fastq}
    """
}


// run kwip on collected khmer hash tables to get distance matrix
process "kwip" {

    label "${params.high ? "high" : "low"}"
    label "finish"

    input:
    path(hashes)
    // eg. [/path/to/sample_1.ct.gz, ... /path/to/sample_n.ct.gz]

    output:
    tuple path("kern.txt"), path("dist.txt")
    // eg. [/path/to/kern.txt, /path/to/dist.txt]

    when:
    params.clusters || (!params.variants && !params.clusters)

    script:
    """
    kwip -t ${task.cpus} -k kern.txt -d dist.txt ${hashes}
    """
}


// run Rscript on kwip results to generate clustering plots
process "clustering" {

    label "low"
    label "finish"

    input:
    tuple path(kern), path(dist)
    // eg. [/path/to/kern.txt, /path/to/dist.txt]

    output:
    path "*.pdf"

    when:
    params.clusters || (!params.variants && !params.clusters)

    script:
    """
    Rscript ${baseDir}/bin/img.R dist.txt kern.txt clustering 
    """
}



// sorting bam files from bisulfite-masked samples
process "sorting" {

    label "${params.high ? "high" : "low"}"
    label "finish"
    tag "$sample"

    maxForks "${params.fork}".toInteger()

    input:
    tuple type, sample, path(bam)
    // eg. [variant, sample, /path/to/sample.bam]

    output:
    tuple sample, path("${sample}.bam"), path("${sample}.bam.bai")
    // eg. [sample, /path/to/sample.bam, /path/to/sample.bam.bai]
    path "${sample}.bam"

    when:
    params.variants || (!params.variants && !params.clusters)

    script:
    """
    samtools sort -T deleteme -m ${((task.memory.getBytes() / task.cpus) * 0.9).round(0)} -@ ${task.cpus} \\
    -o ${sample}.bam ${bam}
    samtools index ${sample}.bam
    """
}


// variant calling on bisulfite-masked samples
process "freebayes" {

    label "${params.high ? "high" : "low"}"
    label "finish"
    tag "$sample"

    maxForks "${params.fork}".toInteger()

    input:
    tuple sample, path(bam), path(bai)
    // eg. [sample, /path/to/sample.bam, /path/to/sample.bam.bai]
    path fasta
    path fai

    output:
    tuple sample, path("${sample}.vcf")
    // eg. [sample, /path/to/sample.vcf]

    when:
    params.variants || (!params.variants && !params.clusters)

    script:
    """
    fasta_generate_regions.py ${fai} ${params.regions} > regions.txt
    freebayes-parallel regions.txt ${task.cpus} -f ${fasta} ${bam} \\
    --no-partial-observations --report-genotype-likelihood-max --genotype-qualities --min-repeat-entropy 1 > ${sample}.vcf
    """
}



// filtering of variants
process "bcftools" {

    label "${params.high ? "high" : "low"}"
    label "finish"
    tag "$sample"

    input:
    tuple sample, path("raw.vcf")
    // eg. [sample, /path/to/raw.vcf]

    output:
    tuple sample, path("${sample}.bcf")
    // eg. [sample, /path/to/sample.bcf]

    when:
    params.variants || (!params.variants && !params.clusters)

    script:
    """
    #bcftools view -Ob${params.ploidy ? " --max-alleles ${params.ploidy}" : ""} raw.vcf > ${sample}.bcf
    bcftools view -Ob raw.vcf > ${sample}.bcf
    """
}


// generate plots with plot-bamstats
process "plot_vcfstats" {

    label "low"
    label "ignore"
    tag "$sample"

    input:
    tuple sample, path(bcf)
    // eg. [sample, /path/to/sample.bcf]

    output:
    path "${sample}/*"

    when:
    params.variants || (!params.variants && !params.clusters)

    script:
    """
    mkdir ${sample}
    bcftools stats ${bcf} > ${sample}/${bcf}.stats || exit \$?
    plot-vcfstats -P -p ${sample} ${sample}/${bcf}.stats
    """
}


// preprocessing for haplotype phasing
process "extractHAIRS" {

    label "low"
    tag "$sample"

    input:
    tuple sample, path(bam), path(bai), path(vcf)
    // eg. [sample, /path/to/sorted.bam, /path/to/sorted.bam.bai, /path/to/sample.vcf]

    output:
    tuple sample, path("${sample}.txt")
    // eg. [sample, /path/to/sample.txt]

    when:
    params.variants || (!params.variants && !params.clusters)

    script:
    """
    extractHAIRS --bam ${bam} --VCF ${vcf} --out ${sample}.txt
    """
}


// diploid haplotype phasing
process "HAPCUT2" {

    label "low"
    tag "$sample"

    input:
    tuple sample, path(txt), path(vcf)
    // eg. [sample, /path/to/sample.txt, /path/to/sample.vcf]

    output:
    tuple sample, path("phased.${sample}.vcf")
    // eg. [sample, /path/to/sample.vcf]

    when:
    params.variants || (!params.variants && !params.clusters)

    script:
    """
    HAPCUT2 --fragments ${txt} --VCF ${vcf} --output phased.${sample}.vcf
    """
}


// split alignments based on haplotype phasing
process "bamsplit" {

    label "low"
    tag "$sample"

    input:
    tuple sample, path(bam), path(bai), path(vcf)
    // eg. [sample, /path/to/sorted.bam, /path/to/sorted.bam.bai, /path/to/phased.vcf]
    path fasta

    output:
    tuple sample, path("bam/*.bam")
    // eg. [sample, [/path/to/*.bam, ...]]

    when:
    params.variants || (!params.variants && !params.clusters)

    script:
    """
    mkdir bam
    bamsplit.py --ref ${fasta} --reads ${bam} --variants ${vcf} --out_dir bam --ploidy ${params.ploidy}${params.region ? "" : " --region ${params.region}"}
    """
}


// allele-specific Methylation calling
process "MethylDackel" {

    label "low"
    tag "$sample"

    input:
    tuple sample, path(bam)
    // eg. [sample, /path/to/haplotype.bam]
    path fasta
    val context

    output:
    tuple sample, path("${sample}/${bam.baseName}/*.bedGraph")
    // eg. [sample, [/path/to/*_CpG.bedGraph, ...]]
    tuple sample, path("${sample}/${bam.baseName}/logs/*")
    // eg. [sample, [/path/to/logs/*, ...]]

    when:
    params.variants || (!params.variants && !params.clusters)

    script:
    """
    mkdir ${sample} ${sample}/${bam.baseName} ${sample}/logs
    samtools index ${bam}

    STR=\$(echo \$(MethylDackel mbias ${fasta} ${bam} ${sample}/logs/${bam.baseName} ${context} 2>&1 | cut -d ":" -f2))
    MethylDackel extract ${fasta} ${bam} ${context} -o ${sample}/${bam.baseName}/ \$STR \\
    > ${sample}/logs/${bam.baseName}.err 2>&1
    """
}