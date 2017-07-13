from pprint import pprint

from pipelines import MainPipeline
import debug_values
import cmd_parser

# ----------------------------------------------------------------------------------------------------
# --------------------------------------------- Main -------------------------------------------------
# ----------------------------------------------------------------------------------------------------


if __name__ == "__main__" or True:
	
	parser = cmd_parser.getCmdParser()
	parser.parse_args()
	print(parser)
	DEBUG = True


	if DEBUG:
		default_sample_filename = debug_values.default_sample_filename
		default_config_filename = debug_values.default_config_filename
		pipeline_dna_callers 	= debug_values.pipeline_dna_callers
		pipeline_rna_callers 	= debug_values.pipeline_rna_callers
		pipeline_copynumber_callers = debug_values.pipeline_copynumber_callers
		somaticseq_callers 		= debug_values.somaticseq_callers
	else:
		default_sample_filename = ""
		default_config_filename = ""
		pipeline_dna_callers = []
		pipeline_rna_callers = []
		pipeline_copynumber_callers = []
		somaticseq_callers = ""


	pipeline = MainPipeline(
		default_sample_filename,
		default_config_filename,
		dna_callers = pipeline_dna_callers,
		rna_callers = pipeline_rna_callers,
		copynumber_callers = pipeline_copynumber_callers,
		somaticseq_callers = somaticseq_callers,
		debug = True
	)


