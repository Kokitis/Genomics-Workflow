from .basepipeline import BasePipeline
from ..callers import *

class CopynumberWorkflow(BasePipeline):
	def runWorkflow(self, sample, options, workflow_callers):
		if len(workflow_callers) == 1 and 'all' in workflow_callers:
			workflow_callers = ['freec', 'varscan']

		if 'cnvkit' in workflow_callers:
			cnvkit_result = CNVkit(sample, options)
		if 'freec' in workflow_callers:
			freec_result = FREEC(sample, options)
		if 'varscan' in workflow_callers:
			varscan_result = VarscanCopynumber(sample, options)