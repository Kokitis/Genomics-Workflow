
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
DEBUG = True

if __name__ == "__main__" or True:
	
	if DEBUG:
		from debug import *
	else:
		pass
	pipeline = MainPipeline(
		default_sample_filename,
		default_config_filename,
		dna_callers = pipeline_dna_callers,
		rna_callers = pipeline_rna_callers,
		copynumber_callers = pipeline_copynumber_callers,
		somaticseq_callers = somaticseq_callers,
		debug = True
	)
# /home/upmc/Documents/TCGA-ESCA/RNA-seq/579bce59-438b-4ee2-b199-a91de73bca0e/b33ec9ae-7692-465d-ab40-b9a140df9c2e_gdc_realn_rehead.bam

