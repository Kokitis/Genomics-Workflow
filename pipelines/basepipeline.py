from github import filetools
import os
import gdc_api
API = gdc_api.GDCAPI()

class BasePipeline:

	def __init__(self, sample, pipeline_options, callers):
		if True:
			print("Callers to use: ", callers)
		
		self._verifyPipelineFiles(pipeline_options)
		self._verifySampleFiles(sample, pipeline_options)

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

	def _verifySampleFiles(self, sample, options):
		""" Verifies that all required sample files exist and are not corrupted. """
		#Verify BAM Files
		self._checkIfPathExists('NormalBAM', sample['NormalBAM'])
		self._checkIfPathExists('TumorBAM',  sample['TumorBAM'])

		md5_normal = filetools.generateFileMd5(sample['NormalBAM'])
		expected_md5sum_normal = API(sample['NormalUUID'], 'files')
		if md5_normal != expected_md5sum_normal and not options['globals']['debug']:
			message = "The Normal BAM (ID {}) does not have the correct md5sum!".format(sample['NormalID'])
			raise ValueError(message)

		md5_sample = filetools.generateFileMd5(sample['TumorBAM'])
		expected_md5sum_sample = API(sample['SampleUUID'], 'files')
		if md5_sample != expected_md5sum_sample and not options['globals']['debug']:
			message = "The Tumor BAM (ID {}) does not have the correct md5sum!".format(sample['SampleID'])
			raise ValueError(message)

		#Verify exome targets File
		self._checkIfPathExists("the exome targets", sample['ExomeTargets'])

	@staticmethod
	def _checkIfPathExists(label, path):
		if not os.path.exists(path):
			message = "Missing file for {}: {}".format(label, path)
			raise FileNotFoundError(message)

	def runWorkflow(self, sample, pipeline_options, pipeline_callers):
		pass