import os


class BasePipeline:

	def __init__(self, sample, pipeline_options, callers):
		if True:
			print("Callers to use: ", callers)
		
		self._verifyPipelineFiles(pipeline_options)

		self.runWorkflow(sample, pipeline_options, callers)

	def _verifyPipelineFiles(self, options):
		""" Verifies that the files required to run the pipeline exist """
		self._verifyCommonDependancies(options)
		self._verifyLocalDependancies(options)

		return options
		
	def _verifyCommonDependancies(self, options):
		# Verify That the Reference Files Exist
		self._checkIfPathExists('GATK', 			options['Programs']['GATK'])
		self._checkIfPathExists('reference genome', options['Reference Files']['reference genome'])
		self._checkIfPathExists('dbSNP', 			options['Reference Files']['dbSNP'])
		self._checkIfPathExists('COSMIC', 			options['Reference Files']['cosmic'])
		self._checkIfPathExists('samtools', 		options['Programs']['samtools'])
		self._checkIfPathExists('samtools (0.1.6)', options['Programs']['samtools-0.1.6'])


	def _verifyLocalDependancies(self, options):
		pass

	@staticmethod
	def _checkIfPathExists(label, path):
		if not os.path.exists(path):
			message = "Missing file for {}: {}".format(label, path)
			raise FileNotFoundError(message)

	def runWorkflow(self, sample, pipeline_options, pipeline_callers):
		pass