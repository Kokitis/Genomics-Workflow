from .copynumberpipeline import CopynumberWorkflow
from .dnapipeline import DNAWorkflow
from .rnapipeline import RNAWorkflow
from .somaticseqpipeline import SomaticSeqPipeline
from settings import Settings
import os
from github import tabletools
from github import filetools
from github import gdc_api

def _toLower(array):
	return [i.lower() for i in array]

class MainPipeline:
	def __init__(self, sample_filename, options_filename, dna_callers = None, copynumber_callers = None, rna_callers = None, somaticseq_callers = None, **kwargs):

		if True:
			print("sample filename: ", sample_filename)
			print("options filename: ", options_filename)
			print("DNA-seq callers: ", dna_callers)
			print("RNA-seq callers: ", rna_callers)
			print("copynumber callers: ", copynumber_callers)
			print("somaticseq callers: ", somaticseq_callers)
			print("Keyword Arguments:")
			for key, value in kwargs.items():
				print("\t{:<20}\t{}".format(key, value))

		dna_callers = [] if dna_callers is None else _toLower(dna_callers)
		rna_callers = [] if rna_callers is None else _toLower(rna_callers)
		copynumber_callers = [] if copynumber_callers is None else _toLower(copynumber_callers)
		somaticseq_callers = [] if somaticseq_callers is None else _toLower(somaticseq_callers)

		sample_list = tabletools.Table(sample_filename)
		pipeline_options = Settings(options_filename, **kwargs)

		for index, current_sample in sample_list:
			self._runSample(
				current_sample,
				pipeline_options,
				dna_callers,
				copynumber_callers,
				rna_callers,
				somaticseq_callers
			)
	@staticmethod
	def _checkIfPathExists(label, path):
		if not os.path.exists(path):
			message = "Missing file for {}: {}".format(label, path)
			raise FileNotFoundError(message)

	def _runSample(self, sample, sample_options, dna_callers, copynumber_callers, rna_callers, somaticseq_callers):
		sample = sample.to_dict()
		_use_value = sample.get('Use', False)
		_use_this_sample = _use_value not in {False, 'false', 0, '0', 'no', 'No'}
		self._verifySampleFiles(sample, sample_options)

		if _use_this_sample:
			sample_status = list()
			if dna_callers:
				dna_pipeline_status = DNAWorkflow(sample, sample_options, dna_callers)
				sample_status.append(dna_pipeline_status)
			if rna_callers:
				rna_pipeline_status = RNAWorkflow(sample, sample_options, rna_callers)
				sample_status.append(rna_pipeline_status)
			if copynumber_callers:
				cn_pipeline_status = CopynumberWorkflow(sample, sample_options, copynumber_callers)
				sample_status.append(cn_pipeline_status)
			if somaticseq_callers:
				somaticseq_status = SomaticSeqPipeline(sample, sample_options, somaticseq_callers)

	def _verifySampleFiles(self, sample, options):
		""" Verifies that all required sample files exist and are not corrupted. """
		#Verify BAM Files
		self._checkIfPathExists('NormalBAM', sample['NormalBAM'])
		self._checkIfPathExists('TumorBAM',  sample['TumorBAM'])

		md5_normal = filetools.generateFileMd5(sample['NormalBAM'])
		expected_md5sum_normal = gdc_api.request(sample['NormalUUID'], 'files')
		if md5_normal != expected_md5sum_normal and not options['globals']['debug']:
			message = "The Normal BAM (ID {}) does not have the correct md5sum!".format(sample['NormalID'])
			raise ValueError(message)

		md5_sample = filetools.generateFileMd5(sample['TumorBAM'])
		expected_md5sum_sample = gdc_api.request(sample['SampleUUID'], 'files')
		if md5_sample != expected_md5sum_sample and not options['globals']['debug']:
			message = "The Tumor BAM (ID {}) does not have the correct md5sum!".format(sample['SampleID'])
			raise ValueError(message)

		#Verify exome targets File
		self._checkIfPathExists("the exome targets", sample['ExomeTargets'])