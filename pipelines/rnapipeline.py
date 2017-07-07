from .basepipeline import BasePipeline
from callers import *
class RNAWorkflow(BasePipeline):
	def runWorkflow(self, sample, options, workflow_callers):

		bqsr = BaseQualityScoreRecalibration(sample, options)
		processed_bam = bqsr.final_output

		if 'haplotypecaller' in workflow_callers:
			haplotypecaller_status = HaplotypeCaller(sample, options, processed_bam)