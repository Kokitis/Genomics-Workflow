from .basepipeline import BasePipeline
from callers import *
class RNAWorkflow(BasePipeline):
	
	def runWorkflow(self, sample, pipeline_options, pipeline_callers):
		if len(pipeline_callers) == 1 and 'all' in pipeline_callers:
			pipeline_callers = ['haplotypecaller']
		bqsr = BaseQualityScoreRecalibration(sample, pipeline_options)
		processed_bam = bqsr.final_output

		if 'haplotypecaller' in pipeline_callers:
			haplotypecaller_status = HaplotypeCaller(sample, pipeline_options, processed_bam)