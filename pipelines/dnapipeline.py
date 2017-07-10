from .basepipeline import BasePipeline
from callers import *

class DNAWorkflow(BasePipeline):
	def _verifyLocalDependancies(self, options):
		self._checkIfPathExists('MuSE', 			options['Programs']['muse'])
		self._checkIfPathExists('MuTect2', 			options['Programs']['mutect2'])
		self._checkIfPathExists('SomaticSniper', 	options['Programs']['somaticsniper'])
		self._checkIfPathExists('Strelka', 			options['Programs']['strelka'])
		self._checkIfPathExists('Varscan2', 		options['Programs']['varscan'])
	def runWorkflow(self, sample, pipeline_options, pipeline_callers):
		if len(pipeline_callers) == 1 and 'all' in pipeline_callers:
			pipeline_callers = ['muse', 'mutect', 'somaticsniper', 'strelka', 'varscan']
		pipeline_status = list()
		if 'haplotypecaller' in pipeline_callers:
			haplotypecaller_result = HaplotypeCaller(sample, pipeline_options)
			pipeline_status.append(haplotypecaller_result)
		if 'muse' in pipeline_callers:
			muse_result = MuSE(sample, pipeline_options)
			pipeline_status.append(muse_result)
		if 'mutect' in pipeline_callers:
			mutect_result = MuTect2(sample, pipeline_options)
			pipeline_status.append(mutect_result)
		if 'somaticsniper' in pipeline_callers:
			somaticsniper_result = SomaticSniper(sample, pipeline_options)
			pipeline_status.append(somaticsniper_result)
		if 'strelka' in pipeline_callers:
			strelka_result = Strelka(sample, pipeline_options)
			pipeline_status.append(strelka_result)
		if 'varscan' in pipeline_callers:
			varscan_result = Varscan(sample, pipeline_options)
			pipeline_status.append(varscan_result)