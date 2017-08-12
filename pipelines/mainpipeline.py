from .copynumberpipeline import CopynumberWorkflow
from .dnapipeline import DNAWorkflow
from .rnapipeline import RNAWorkflow
from .somaticseqpipeline import SomaticSeqPipeline
from settings import Settings
import os
from github import tabletools
from github import filetools
from github import gdc_api
from github import callertools

caller_classifier = callertools.CallerClassifier()

def _toLower(array):
	return [i.lower() for i in array]

def _getMissingCallers(callset_folder):
	existing_callset = caller_classifier(callset_folder)
	from pprint import pprint

	print(callset_folder)
	pprint(existing_callset)
	missing = list()
	if 'muse' not in existing_callset:
			missing.append('muse')
	if 'mutect' not in existing_callset:
		missing.append('mutect')
	if 'somaticsniper' not in existing_callset:
		missing.append('somaticseq')
	if 'strelka-snp' not in existing_callset or 'strelka-indel' not in existing_callset:
		missing.append('strelka')
	if 'varscan-snp' not in existing_callset:
		missing.append('varscan')
	return missing


class MainPipeline:
	def __init__(self, sample_filename, options_filename, dna_callers = None, copynumber_callers = None, rna_callers = None, somaticseq_callers = None, **kwargs):

		if False:
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
		self._current_index = 0
		for index, current_sample in enumerate(sample_list):
			self._current_index = "{} of {}".format(index + 1, len(sample_list))
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
		print("{} ({})".format(sample['PatientID'], self._current_index))
		print("\tDNA-seq callers: ", dna_callers)
		print("\tRNA-seq callers: ", rna_callers)
		print("\tCopynumber callers: ", copynumber_callers)
		print("\tSomaticSeq callers: ", somaticseq_callers)
		_use_value = sample.get('Use', False)
		_use_this_sample = _use_value not in {False, 'false', 0, '0', 'no', 'No'}
		if 'Histology' in sample:
			_use_this_sample = _use_this_sample and sample['Histology'] == 'Esophagus Adenocarcinoma, NOS'
		if _use_this_sample:
			
			self._verifySampleFiles(sample, sample_options)

		if 'missing' in dna_callers:
			dna_callset_folder = sample_options.getPipelineFolder('callset', sample['PatientID'])
			dna_callers = _getMissingCallers(dna_callset_folder)


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
		else:
			print("SKIPPED SAMPLE")
		print()

	def _verifySampleFiles(self, sample, options, print_errors = True):
		""" Verifies that all required sample files exist and are not corrupted. """
		#Verify BAM Files
		self._checkIfPathExists('NormalBAM', sample['NormalBAM'])
		self._checkIfPathExists('TumorBAM',  sample['TumorBAM'])

		_error_message = message = "\tThe Tumor BAM (ID {}) does not have the correct md5sum!\n\tPath\t{}\n\tExpected\t{}\n\tReceived\t{}"

		md5_normal = filetools.generateFileMd5(sample['NormalBAM'])
		expected_md5sum_normal = gdc_api.request(sample['NormalUUID'], 'files')['md5sum']
		if md5_normal != expected_md5sum_normal and not options['globals']['debug']:
			#message = "\tThe Normal BAM (ID {}) does not have the correct md5sum!".format(sample['NormalID'])
			message = _error_message.format(sample['NormalID'], sample['NormalBAM'], expected_md5sum_normal, md5_normal)
			if not print_errors:
				raise ValueError(message)
			else:
				print(message)

		md5_sample = filetools.generateFileMd5(sample['TumorBAM'])
		expected_md5sum_sample = gdc_api.request(sample['SampleUUID'], 'files')['md5sum']
		if md5_sample != expected_md5sum_sample and not options['globals']['debug']:
			message = _error_message.format(sample['SampleID'], sample['TumorBAM'], expected_md5sum_sample, md5_sample)
			if not print_errors:
				raise ValueError(message)
			else:
				print(message)

		#Verify exome targets File
		self._checkIfPathExists("the exome targets", sample['ExomeTargets'])