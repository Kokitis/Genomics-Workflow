
import datetime
import os
import gdc_api
from settings import Settings
from pprint import pprint

from pipelines import MainPipeline
from github import tabletools

# --------------------------------- Global Variables ----------------------------
now = datetime.datetime.now

# ----------------------------------------------------------------------------------------------------
# --------------------------------------------- Main -------------------------------------------------
# ----------------------------------------------------------------------------------------------------


if __name__ == "__main__" or True:
	
	parser = cmd_parser.getCmdParser()
	parser.parse_args()


	if DEBUG:
		default_sample_filename = debug.default_sample_filename
		default_config_filename = debug.default_config_filename
		pipeline_dna_callers 	= debug.pipeline_dna_callers
		pipeline_rna_callers 	= debug.pipeline_rna_callers
		pipeline_copynumber_callers = debug.pipeline_copynumber_callers
		somaticseq_callers 		= debug.somaticseq_callers
	else:
		default_sample_filename = ""
		default_config_filename = ""
		pipeline_dna_callers = ""
		pipeline_copynumber_callers = ""
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


