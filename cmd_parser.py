from argparse import ArgumentParser
def getCMDArgumentParser():

	"""Parse the command-line arguments for this program."""

	parser = ArgumentParser(
		description='Run the genomics pipeline')

	show_default = ' (default %(default)s)'

	parser.add_argument(
		'-l', '--sample-list',
		action='store',
		dest = 'sample_list',
		help='the sample list')

	parser.add_argument(
		'-d', "--debug",
		dest = 'debug',
		action = 'store_true',
		help='debug the pipeline using default settings' + show_default)

	parser.add_argument(
		"-i", "--ignore-caller-status",
		dest = "ignore_caller_status",
		action = 'store_false',
		help = "Ignore the caller status file.")

	parser.add_argument(
		"-a", "--all-callers",
		dest = 'use_all_callers',
		action = 'store_true',
		help = "Tells the pipeline to use all available callers")

	parser.add_argument(
		"--force-overwrite",
		dest = 'force_overwrite',
		action = 'store_true',
		help = "Whether to delete any caller files and re-do the analysis.")

	return parser

def parseCommandLineOptions(parser):
	dna_callers = []
	rna_callers = []
	copynumber_callers = []
	debug_pipeline = False

	result = {
		'dnaCallers': dna_callers,
		'rnaCallers': rna_callers,
		'copynumberCallers': copynumber_callers,
		'debug': debug_pipeline
	}
	return result
