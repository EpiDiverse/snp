#!/usr/bin/env nextflow

// FUNCTION TO LOAD DATASETS IN TEST PROFILE
def check_test_data(BAMPaths) {

    // Set BAM testdata
    BAM = Channel.from(BAMPaths)
        .map { row -> [ row[0], file(row[1]) ] }
        .ifEmpty { exit 1, "test profile BAMPaths was empty - no input files supplied" }

    // Return BAM channel
    return BAM
}