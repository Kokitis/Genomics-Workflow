
import datetime
import os
import gdc_api
from settings import Settings
from pprint import pprint

from pipelines import *
from github import tabletools

# --------------------------------- Global Variables ----------------------------
now = datetime.datetime.now

# ----------------------------------------------------------------------------------------------------
# --------------------------------------------- Main -------------------------------------------------
# ----------------------------------------------------------------------------------------------------
DEBUG = True
class GenomicsPipeline:
	def __init__(self, sample_filename, options_filename, dna_callers = None, copynumber_callers = None, rna_callers = None, somaticseq_callers = None, **kwargs):
		if dna_callers is None: dna_callers = []
		if rna_callers is None: rna_callers = []
		if copynumber_callers is None: copynumber_callers = []
		if somaticseq_callers is None: somaticseq_callers = []
		dna_callers = [i.lower() for i in dna_callers]
		rna_callers = [i.lower() for i in rna_callers]
		copynumber_callers = [i.lower() for i in copynumber_callers]

		sample_list = tabletools.Table(sample_filename)

		pipeline_options = Settings(options_filename, **kwargs)

		for index, current_sample in sample_list:
			self._runSample(
				current_sample,
				pipeline_options,
				dna_callers,
				copynumber_callers,
				rna_callers,
				somaticseq_callers
			)

	@staticmethod
	def _runSample(sample, sample_options, dna_callers, copynumber_callers, rna_callers, somaticseq_callers):
		sample = sample.to_dict()
		_use_value = sample.get('Use', False)
		_use_this_sample = _use_value not in {False, 'false', 0, '0', 'no', 'No'}

		if _use_this_sample:
			sample_status = list()
			if dna_callers:
				dna_pipeline_status = DNAWorkflow(sample, sample_options, dna_callers)
				sample_status.append(dna_pipeline_status)
			if rna_callers:
				rna_pipeline_status = RNAWorkflow(sample, sample_options, rna_callers)
				sample_status.append(rna_pipeline_status)
			if copynumber_callers:
				cn_pipeline_status = CopynumberWorkflow(sample, sample_options, copynumber_callers)
				sample_status.append(cn_pipeline_status)
			if somaticseq_callers:
				pass
				#somaticseq_status = somaticseqPipeline()

API = gdc_api.GDCAPI()

if __name__ == "__main__" or True:

	#command_line_options = cmd_parser.getCMDArgumentParser().parse_args()
	PIPELINE_DIRECTORY = "/home/upmc/Documents/Variant_Discovery_Pipeline"
	default_config_filename = os.path.join(PIPELINE_DIRECTORY, "0_config_files", "pipeline_project_options.txt")
	
	default_sample_filename = os.path.join(PIPELINE_DIRECTORY, "rna_sample_list.tsv")
	pipeline_dna_callers = []
	pipeline_rna_callers = ['haplotypecaller']
	pipeline_copynumber_callers = []
	

	pipeline = GenomicsPipeline(
		default_sample_filename,
		default_config_filename,
		dna_callers = pipeline_dna_callers,
		rna_callers = pipeline_rna_callers,
		copynumber_callers = pipeline_copynumber_callers,
		debug = True
	)
# /home/upmc/Documents/TCGA-ESCA/RNA-seq/579bce59-438b-4ee2-b199-a91de73bca0e/b33ec9ae-7692-465d-ab40-b9a140df9c2e_gdc_realn_rehead.bam

