from pprint import pprint

from pipelines import MainPipeline
from truthset import TruthsetPipeline
import debug_values
import cmd_parser
import os
# ----------------------------------------------------------------------------------------------------
# --------------------------------------------- Main -------------------------------------------------
# ----------------------------------------------------------------------------------------------------
def runVariantDiscoveryPipeline(cmd_options = None, debug_pipeline = True):
	if debug_pipeline:
		default_sample_filename = debug_values.default_sample_filename
		default_config_filename = debug_values.default_config_filename
		pipeline_dna_callers 	= debug_values.pipeline_dna_callers
		pipeline_rna_callers 	= debug_values.pipeline_rna_callers
		pipeline_copynumber_callers = debug_values.pipeline_copynumber_callers
		somaticseq_callers 		= debug_values.somaticseq_callers
	else:
		PIPELINE_FOLDER = os.path.join(os.getenv("HOME"), "Documents", "Variant_Discovery_Pipeline")
		default_sample_filename = os.path.join(PIPELINE_FOLDER, "sample_list.tsv")
		default_config_filename = os.path.join(PIPELINE_FOLDER, "0_config_files", "pipeline_project_options.txt")
		pipeline_dna_callers = ['all']
		pipeline_rna_callers = ['all']
		pipeline_copynumber_callers = ['all']
		somaticseq_callers = ['table']

	pipeline = MainPipeline(
		default_sample_filename,
		default_config_filename,
		dna_callers = pipeline_dna_callers,
		rna_callers = pipeline_rna_callers,
		copynumber_callers = pipeline_copynumber_callers,
		somaticseq_callers = somaticseq_callers,
		debug = True
	)


def runTruthsetPipeline(debug_pipeline = True):
	
	if debug_pipeline:
		samples = ['TCGA-2H-A9GF']
		default_config_filename = debug_values.default_config_filename
		truthset_type = 'rna'
	else:
		samples = []
		default_config_filename = ""
		truthset_type = 'intersection'

	pipeline = TruthsetPipeline(
		samples = samples,
		options_filename = default_config_filename,
		truthset_type = truthset_type
	)

if __name__ == "__main__" or True:
	
	parser = cmd_parser.getCmdParser()
	parser = parser.parse_args()

	debug_pipeline = parser.debug


	run_variant_discovery = not parser.generate_truthset

	if run_variant_discovery:
		runVariantDiscoveryPipeline(debug_pipeline)
	else:
		runTruthsetPipeline(debug_pipeline)







