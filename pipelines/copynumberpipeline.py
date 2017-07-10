from .basepipeline import BasePipeline
from callers import *

class CopynumberWorkflow(BasePipeline):
	def _verifyLocalDependancies(self, options):
		self._checkIfPathExists('bam-readcount', 	options['Programs']['bam-readcount'])
		self._checkIfPathExists('CNVKit', 			options['Programs']['cnvkit'])

		# Verify that other files exist
		self._checkIfPathExists(
			'reference genome index',
			options['Reference Files']['reference genome'] + '.fai'
		) #for FREEC
	def runWorkflow(self, sample, pipeline_options, pipeline_callers):
		if len(pipeline_callers) == 1 and 'all' in pipeline_callers:
			pipeline_callers = ['freec', 'varscan']
		pipeline_status = list()
		if 'cnvkit' in pipeline_callers:
			cnvkit_result = CNVkit(sample, pipeline_options)
			pipeline_status.append(cnvkit_result)
		if 'freec' in pipeline_callers:
			freec_result = FREEC(sample, pipeline_options)
			pipeline_status.append(freec_result)
		if 'varscan' in pipeline_callers:
			varscan_result = VarscanCopynumber(sample, pipeline_options)
			pipeline_status.append(varscan_result)