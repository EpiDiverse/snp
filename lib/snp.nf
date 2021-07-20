#!/usr/bin/env nextflow


// taking input bam files for sorting and indexing
process "preprocessing" {

    label "${params.high ? "high" : "low"}"
    label "finish"
    tag "$sample"

    maxForks "${params.fork}".toInteger()

    input:
    tuple val(sample), path(bam)
    // eg. [sample, /path/to/sample.bam]
    path fasta

    output:
    tuple val(sample), path("calmd.bam"), path("calmd.bam.bai")
    // eg. [sample, /path/to/calmd.bam, /path/to/calmd.bam.bai]

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
    tuple val(type), val(sample), path(bam), path(bai)
    // eg. [clustering, sample, /path/to/sample.bam, /path/to/sample.bam.bai]

    output:
    tuple val(type), val(sample), path("${type}.bam")
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

    publishDir "${params.output}/bam/clusters", mode: 'copy', enabled: true

    input:
    tuple val(type), val(sample), path(bam)
    // eg. [clustering, sample, /path/to/sample.bam]

    output:
    tuple val(sample), path("${sample}.fastq.gz")
    // eg. [sample, /path/to/sample.fastq.gz]
    //path "${sample}.bam"

    when:
    params.clusters || (!params.variants && !params.phase && !params.clusters)

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

    publishDir "${params.output}/hashes", pattern: "${sample}.ct.gz", mode: 'copy', enabled: true

    input:
    tuple val(sample), path(fastq)
    // eg. [sample, /path/to/sample.fastq.gz]

    output:
    path "${sample}.ct.gz"
    // eg. [/path/to/sample.ct.gz]

    when:
    params.clusters || (!params.variants && !params.phase && !params.clusters)

    script:
    """
    load-into-counting.py -T ${task.cpus} -N 1 -x 1e9 -k 20 -b -f -s tsv ${sample}.ct.gz ${fastq}
    """
}


// run kwip on collected khmer hash tables to get distance matrix
process "kwip" {

    label "${params.high ? "high" : "low"}"
    label "finish"

    publishDir "${params.output}", pattern: "*.txt", mode: 'copy', enabled: true

    input:
    path(hashes)
    // eg. [/path/to/sample_1.ct.gz, ... /path/to/sample_n.ct.gz]

    output:
    tuple path("kern.txt"), path("dist.txt")
    // eg. [/path/to/kern.txt, /path/to/dist.txt]

    when:
    params.clusters || (!params.variants && !params.phase && !params.clusters)

    script:
    """
    kwip -t ${task.cpus} -k kern.txt -d dist.txt ${hashes}
    """
}


// run Rscript on kwip results to generate clustering plots
process "clustering" {

    label "low"
    label "finish"

    publishDir "${params.output}", pattern: "*.pdf", mode: 'move', enabled: true

    input:
    tuple path(kern), path(dist)
    // eg. [/path/to/kern.txt, /path/to/dist.txt]

    output:
    path "*.pdf"

    when:
    params.clusters || (!params.variants && !params.phase && !params.clusters)

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

    publishDir "${params.output}/bam/variants", pattern: "${sample}.bam", mode: 'copy', enabled: true

    input:
    tuple val(type), val(sample), path(bam)
    // eg. [variant, sample, /path/to/sample.bam]

    output:
    tuple val(sample), path("${sample}.bam"), path("${sample}.bam.bai")
    // eg. [sample, /path/to/sample.bam, /path/to/sample.bam.bai]
    path "${sample}.bam"

    when:
    params.variants || params.phase || (!params.variants && !params.phase && !params.clusters)

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

    publishDir "${params.output}/vcf", pattern: "${sample}.vcf", mode: 'copy', enabled: true

    input:
    tuple val(sample), path(bam), path(bai)
    // eg. [sample, /path/to/sample.bam, /path/to/sample.bam.bai]
    path fasta
    path fai

    output:
    tuple val(sample), path("${sample}.vcf")
    // eg. [sample, /path/to/sample.vcf]

    when:
    params.variants || params.phase || (!params.variants && !params.phase && !params.clusters)

    script:
    """
    fasta_generate_regions.py ${fai} ${params.regions} > regions.txt
    freebayes-parallel regions.txt ${task.cpus} -f ${fasta} ${bam} \\
    --strict-vcf --no-partial-observations --report-genotype-likelihood-max --genotype-qualities --min-repeat-entropy ${params.entropy} --min-coverage ${params.coverage} > ${sample}.vcf
    """
}



// filtering of variants
process "bcftools" {

    label "${params.high ? "high" : "low"}"
    label "finish"
    tag "$sample"

    input:
    tuple val(sample), path(vcf)
    // eg. [sample, /path/to/raw.vcf]
    // eg. [sample, [/path/to/scaffold1.vcf, /path/to/scaffold1.vcf, ...]]


    output:
    tuple val(sample), path("vcf/${sample}.vcf.gz"), path("vcf/${sample}.vcf.gz.tbi")
    // eg. [sample, /path/to/sample.vcf.gz, /path/to/sample.vcf.gz.tbi]

    when:
    params.variants || params.phase || (!params.variants && !params.phase && !params.clusters)

    script:
    if (params.phase || (!params.variants && !params.phase && !params.clusters))
        """
        mkdir vcf
        bcftools concat -Oz ${vcf} > vcf/${sample}.vcf.gz
        tabix vcf/${sample}.vcf.gz
        """
    else
        """
        mkdir vcf
        bcftools view -Oz ${vcf} > vcf/${sample}.vcf.gz
        tabix vcf/${sample}.vcf.gz
        """

}


// generate plots with plot-bamstats
process "plot_vcfstats" {

    label "low"
    label "ignore"
    tag "$sample"

    publishDir "${params.output}/stats", pattern: "${sample}/*", mode: 'move', enabled: true

    input:
    tuple val(sample), path(vcf), path(tabix)
    // eg. [sample, /path/to/sample.vcf.gz, /path/to/sample.vcf.gz.tbi]

    output:
    path "${sample}/*"

    when:
    params.variants || params.phase || (!params.variants && !params.phase && !params.clusters)

    script:
    """
    mkdir ${sample}
    bcftools stats ${vcf} > ${sample}/${vcf}.stats || exit \$?
    plot-vcfstats -P -p ${sample} ${sample}/${vcf}.stats
    """
}


// split scaffolds for parallel execution of haplotype phasing
process "split_scaffolds" {

    label "low"
    tag "$sample"

    input:
    tuple val(sample), path(vcf)
    // eg. [sample, /path/to/sample.vcf]

    output:
    tuple val(sample), path("scaffolds/*.vcf")
    // eg. [sample, [/path/to/scaffold1.vcf, /path/to/scaffold2.vcf, ...]]

    when:
    params.phase || (!params.variants && !params.phase && !params.clusters)

    script:
    """
    mkdir scaffolds
    bcftools view -H ${vcf} | cut -f1 | sort | uniq | while read line;
    do bcftools view -h ${vcf} > scaffolds/\${line}.vcf
    bcftools view -H ${vcf} | awk -F "\\t" -v line="\$line" 'BEGIN{OFS="\\t"} \$1==line' >> scaffolds/\${line}.vcf;
    done 
    """
}



// read-based haplotype phasing
process "WhatsHap_phase" {

    label "low"
    tag "$sample"

    input:
    tuple val(sample), path(bam), path(bai), path(vcf)
    // eg. [sample, /path/to/masked.bam, /path/to/masked.bam.bai, /path/to/sample.vcf]
    path fasta
    path fai

    output:
    tuple val(sample), path("phased/${sample}.vcf")
    // eg. [sample, /path/to/sample.vcf]

    when:
    params.phase || (!params.variants && !params.phase && !params.clusters)

    script:
    """
    mkdir phased
    whatshap phase -o phased/${sample}.vcf --reference=${fasta} \\
    --tag=PS ${vcf} ${bam} 2> phased/phased.log
    """
}


// phased read assignment
process "WhatsHap_haplotag" {

    label "low"
    tag "$sample"

    input:
    tuple val(sample), path(bam), path(bai), path(vcf), path(tabix)
    // eg. [sample, /path/to/masked.bam, /path/to/masked.bam.bai, /path/to/sample.vcf.gz, /path/to/sample.vcf.gz.tbi]

    output:
    tuple val(sample), path("phased/phased.list")
    // eg. [sample, /path/to/phased.list]

    when:
    params.phase || (!params.variants && !params.phase && !params.clusters)

    script:
    """
    mkdir phased
    whatshap haplotag --ignore-read-groups --output-haplotag-list phased/phased.list \\
    -o phased/phased.bam ${vcf} ${bam}
    """
}


// split read alignments based on phasing
process "WhatsHap_split" {

    label "low"
    tag "$sample"

    input:
    tuple val(sample), path(bam), path(bai), path(phased)
    // eg. [sample, /path/to/sorted.bam, /path/to/sorted.bam.bai, /path/to/phased.list]
    path fasta

    output:
    tuple val(sample), path("bam/*.bam")
    // eg. [sample, [/path/to/phased/*.bam, ...]]

    when:
    params.phase || (!params.variants && !params.phase && !params.clusters)

    script:
    """
    mkdir bam
    whatshap split --output-h1 bam/HP1.sam --output-h2 bam/HP2.sam ${bam} ${phased}
    echo -e "HP1\\nHP2" | xargs -n1 -P2 -i sh -c "samtools view -Sb bam/'{}'.sam > bam/'{}'.bam"
    """
}


// allele-specific Methylation calling
process "MethylDackel" {

    label "low"
    tag "$sample"

    input:
    tuple val(sample), path(bam)
    // eg. [sample, /path/to/haplotype.bam]
    path fasta
    path fai
    val context

    output:
    tuple val(sample), path("${sample}/${bam.baseName}/*.bedGraph")
    // eg. [sample, [/path/to/*_CpG.bedGraph, ...]]
    tuple val(sample), path("${sample}/${bam.baseName}/logs/*")
    // eg. [sample, [/path/to/logs/*, ...]]

    when:
    params.phase || (!params.variants && !params.phase && !params.clusters)

    script:
    """
    mkdir ${sample} ${sample}/${bam.baseName} ${sample}/logs
    samtools index ${bam}

    STR=\$(echo \$(MethylDackel mbias ${fasta} ${bam} ${sample}/logs/${bam.baseName} ${context} 2>&1 | cut -d ":" -f2))
    MethylDackel extract ${fasta} ${bam} ${context} -o ${sample}/${bam.baseName}/ \$STR \\
    > ${sample}/logs/${bam.baseName}.err 2>&1
    """
}