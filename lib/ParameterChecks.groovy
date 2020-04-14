class ParameterChecks {
	static void checkParams(params) {
		assert params.input, "Please specify path to input with --input parameter"
        assert params.reference, "Please specify path to reference fasta with --reference parameter"
        assert params.regions instanceof Integer && params.regions > 0, "--regions parameter must be a positive integer!"
        assert params.ploidy instanceof Integer && params.ploidy > 0, "--ploidy parameter must be a positive integer!"
        assert params.take instanceof Integer && params.take >= 0, "--take parameter must be a non-negative integer!"
        assert params.fork instanceof Integer && params.fork >= 0, "--fork parameter must be a non-negative integer!"
	}
}