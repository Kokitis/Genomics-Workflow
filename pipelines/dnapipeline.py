from .basepipeline import BasePipeline
from callers import *

class DNAWorkflow(BasePipeline):
	def runWorkflow(self, sample, workflow_options, workflow_callers):
		if len(workflow_callers) == 1 and 'all' in workflow_callers:
			workflow_callers = ['muse', 'mutect', 'somaticsniper', 'strelka', 'varscan']

		if 'haplotypecaller' in workflow_callers:
			haplotypecaller_result = HaplotypeCaller(sample, workflow_options)
		if 'muse' in workflow_callers:
			muse_result = MuSE(sample, workflow_options)
		if 'mutect' in workflow_callers:
			mutect_result = MuTect2(sample, workflow_options)
		if 'somaticsniper' in workflow_callers:
			somaticsniper_result = SomaticSniper(sample, workflow_options)
		if 'strelka' in workflow_callers:
			strelka_result = Strelka(sample, workflow_options)
		if 'varscan' in workflow_callers:
			varscan_result = Varscan(sample, workflow_options)