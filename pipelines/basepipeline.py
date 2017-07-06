from ..settings import Settings

class BasePipeline:

	def __init__(self, sample, options_filename, callers):
		print("Using callers: ", ",".join(callers))
		pipeline_options = self._verifyPipelineFiles(options_filename)
		self._verifySampleFiles(sample)

		self.runWorkflow(sample, pipeline_options, callers)

	def _verifyPipelineFiles(self, options_filename):
		""" Verifies that the files required to run the pipeline exist """

		# verify that the options file exists and load it.
		self._checkIfPathExists('options file', options_filename)
		options = Settings(options_filename)

		# Verify that the required programs exist
		self._checkIfPathExists('GATK', 			options['Programs']['GATK'])
		self._checkIfPathExists('MuSE', 			options['Programs']['muse'])
		self._checkIfPathExists('MuTect2', 			options['Programs']['mutect2'])
		self._checkIfPathExists('SomaticSniper', 	options['Programs']['somaticsniper'])
		self._checkIfPathExists('Strelka', 			options['Programs']['strelka'])
		self._checkIfPathExists('Varscan2', 		options['Programs']['varscan'])
		self._checkIfPathExists('bam-readcount', 	options['Programs']['bam-readcount'])
		self._checkIfPathExists('CNVKit', 			options['Programs']['cnvkit'])
		self._checkIfPathExists('samtools', 		options['Programs']['samtools'])
		self._checkIfPathExists('samtools (0.1.6)', options['Programs']['samtools-0.1.6'])

		# Verify That the Reference Files Exist
		
		self._checkIfPathExists('reference genome', options['Reference Files']['reference genome'])
		self._checkIfPathExists('dbSNP', 			options['Reference Files']['dbSNP'])
		self._checkIfPathExists('COSMIC', 			options['Reference Files']['cosmic'])

		# Verify that other files exist
		self._checkIfPathExists(
			'reference genome index',
			options['Reference Files']['reference genome'] + '.fai'
		) #for FREEC

		return options

	def _verifySampleFiles(self, sample):
		""" Verifies that all required sample files exist and are not corrupted. """
		#Verify BAM Files
		self._checkIfPathExists('NormalBAM', sample['NormalBAM'])
		self._checkIfPathExists('TumorBAM',  sample['TumorBAM'])

		md5_normal = filetools.generateFileMd5(sample['NormalBAM'])
		expected_md5sum_normal = API(sample['NormalUUID'], 'files')
		if md5_normal != expected_md5sum_normal and not DEBUG:
			message = "The Normal BAM (ID {}) does not have the correct md5sum!".format(sample['NormalID'])
			raise ValueError(message)

		md5_sample = filetools.generateFileMd5(sample['TumorBAM'])
		expected_md5sum_sample = API(sample['SampleUUID'], 'files')
		if md5_sample != expected_md5sum_sample and not DEBUG:
			message = "The Tumor BAM (ID {}) does not have the correct md5sum!".format(sample['SampleID'])
			raise ValueError(message)

		#Verify exome targets File
		self._checkIfPathExists("the exome targets", sample['ExomeTargets'])

	@staticmethod
	def _checkIfPathExists(label, path):
		if not os.path.exists(path):
			message = "Missing file for {}: {}".format(label, path)
			raise FileNotFoundError(message)

	def runWorkflow(self, sample, options, callers):
		pass