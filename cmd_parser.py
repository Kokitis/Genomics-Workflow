from argparse import ArgumentParser

def getCmdParser():

	"""Parse the command-line arguments for this program."""

	parser = ArgumentParser(
		description='Run the genomics pipeline')

	show_default = ' (default %(default)s)'

	parser.add_argument(
		'-l', '--sample-list',
		action='store',
		dest = 'sample_list',
		help='the sample list'
	)

	parser.add_argument(
		'-d', "--debug",
		dest = 'debug',
		action = 'store_true',
		help='debug the pipeline using default settings' + show_default
	)

	parser.add_argument(
		"-a", "--all-callers",
		dest = 'use_all_callers',
		action = 'store_true',
		help = "Tells the pipeline to use all available callers"
	)

	parser.add_argument(
		"--dna-callers",
		dest = "dna_callers",
		choices = ['haplotypecaller', 'muse', 'mutect', 'somaticsniper', 'strelka', 'varscan'],
		help = "The callers to use when generating DNA-seq somatic variants"
	)

	parser.add_argument(
		"--rna-callers",
		dest = "rna_callers",
		choices = ['haplotypecaller'],
		help = "The callers to use when generating RNA-seq somatic variants"
	)

	parser.add_argument(
		"--copynumber-callers",
		dest = "copynumber_callers",
		choices = ['cnvkit', 'freec', 'varscan'],
		help = "The callers to use when generating copynumber variants"
	)

	parser.add_argument(
		"--somaticseq-options",
		dest = "somaticseq_options",
		help = "The options to use when running somaticseq. Should be in the form [mode] [truthset] [classifier]"
	)

	parser.add_argument(
		"--force-overwrite",
		dest = 'force_overwrite',
		action = 'store_true',
		help = "Whether to delete any caller files and re-do the analysis."
	)

	return parser

