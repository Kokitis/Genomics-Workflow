import csv
import datetime
import os
import shutil
import isodate
import gdc_api
import settings
from pprint import pprint
from callers import *

# --------------------------------- Global Variables ----------------------------
now = datetime.datetime.now

# ----------------------------------------------------------------------------------------------------
# ------------------------------------- Variant Caller Pipeline --------------------------------------
# ----------------------------------------------------------------------------------------------------



# ----------------------------------------------------------------------------------------------------
# ------------------------------------------ Pipelines -----------------------------------------------
# ----------------------------------------------------------------------------------------------------

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
		options = settings.read(options_filename)

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


class RNAWorkflow(BasePipeline):
	def runWorkflow(self, sample, options, workflow_callers):

		bqsr = BaseQualityScoreRecalibration(sample, options)
		processed_bam = bqsr.final_output

		if 'haplotypecaller' in workflow_callers:
			haplotypecaller_status = HaplotypeCaller(sample, options, processed_bam)

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

# ----------------------------------------------------------------------------------------------------
# --------------------------------------------- Main -------------------------------------------------
# ----------------------------------------------------------------------------------------------------
DEBUG = True
class GenomicsPipeline:
	def __init__(self, sample_filename, options_filename, dna_callers = [], copynumber_callers = [], rna_callers = []):
		dna_callers = [i.lower() for i in dna_callers]
		rna_callers = [i.lower() for i in rna_callers]
		copynumber_callers = [i.lower() for i in copynumber_callers]

		sample_list = tabletools.Table(sample_filename)

		for index, current_sample in sample_list:
			self._runSample(current_sample, options_filename, dna_callers, copynumber_callers, rna_callers)

	@staticmethod
	def _runSample(sample, options_filename, dna_callers, copynumber_callers, rna_callers):
		sample = sample.to_dict()
		_use_value = sample.get('Use', False)
		_use_this_sample = _use_value not in {False, 'false', 0, '0', 'no', 'No'}

		if _use_this_sample:
			if dna_callers:
				dna_pipeline_status = DNAWorkflow(sample, options_filename, dna_callers)
			if rna_callers:
				rna_pipeline_status = RNAWorkflow(sample, options_filename, rna_callers)
			if copynumber_callers:
				cn_pipeline_status = CopynumberWorkflow(sample, options_filename, copynumber_callers)


API = gdc_api.GDCAPI()


if __name__ == "__main__" or True:

	#command_line_options = cmd_parser.getCMDArgumentParser().parse_args()
	PIPELINE_DIRECTORY = "/home/upmc/Documents/Variant_Discovery_Pipeline"
	default_config_filename = os.path.join(PIPELINE_DIRECTORY, "0_config_files", "pipeline_project_options.txt")
	
	default_sample_filename = os.path.join(PIPELINE_DIRECTORY, "debug_sample_list.tsv")
	pipeline_dna_callers = ['all']
	pipeline_rna_callers = []
	pipeline_copynumber_callers = []
	

	pipeline = GenomicsPipeline(
		default_sample_filename,
		default_config_filename,
		dna_callers = pipeline_dna_callers,
		rna_callers = pipeline_rna_callers,
		copynumber_callers = pipeline_copynumber_callers
	)
# /home/upmc/Documents/TCGA-ESCA/RNA-seq/579bce59-438b-4ee2-b199-a91de73bca0e/b33ec9ae-7692-465d-ab40-b9a140df9c2e_gdc_realn_rehead.bam

