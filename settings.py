import configparser
import os
from github import filetools

class Settings:
	"""
		Global Parameters
		-----------------
		* 'debug': 
		* 'overwrite'
	"""
	def __init__(self, filename, **kwargs):
		self.options = configparser.ConfigParser()
		self.options.read(filename)
		kwargs = self._parseKeywordArguments(kwargs)
		self.keyword_arguments = kwargs
		self.base_pipeline_folder = self.options['Pipeline Options']['pipeline folder']

	def __call__(self, *args):
		return self.__getitem__(*args)

	def __getitem__(self, index):

		if isinstance(index, tuple):
			return self.getPipelineFolder(*index)
		elif index == 'globals':
			return self.keyword_arguments
		else:
			return self.options[index]
	@staticmethod
	def _parseKeywordArguments(kwargs):
		kwargs['debug'] = kwargs.get('debug', False)
		kwargs['overwrite'] = kwargs.get('overwrite', False)
		kwargs['verbose'] = kwargs.get('verbose', 'labels')
		return kwargs

	def getPipelineFolder(self, step, patientId = None, caller_name = None):
		if step == 'callset':
			subfolders = ["3_called_variants", patientId]
		elif step == 'variants-somatic':
			subfolders = ["3_called_variants", patientId, caller_name]
		elif step == 'variants-copynumber':
			subfolders = ["4_called_cnvs", patientId, caller_name]
		elif step == 'temporary':
			subfolders = ['5_temporary_files', patientId]
		elif step == 'bam-files':
			subfolders = []
		elif step == 'reference':
			return "/home/upmc/Documents/Reference/"
		elif step == 'variants-rna':
			subfolders = ['7_rna_variants', patientId]

		elif step == 'truthset':
			subfolders = ['truthset']
		elif step == 'somaticseq':
			subfolders = ['somaticseq']
		elif step == 'somaticseq-callset':
			subfolders = ['somaticseq', 'callsets', patientId]
		elif step == 'somaticseq-training':
			subfolders = ['somaticseq', 'training']
		elif step == 'somaticseq-prediction':
			subfolders = ['somaticseq', 'prediction']
		elif step == 'somaticseq-table':
			subfolders = ['somaticseq', 'tables']
		else:
			message = "'{}' is not a valid step in the pipeline!".format(step)
			raise ValueError(message)

		pipeline_folder = os.path.join(self.base_pipeline_folder, *subfolders)
		filetools.checkDir(pipeline_folder, True)
		
		return pipeline_folder
	def getPipelineFile(self, name):
		if name == 'sample log':
			fn = os.path.join(self.base_pipeline_folder, "0_config_files", "benchmark_log.tsv")
		else:
			message = "'{}' is not defined as a pipeline file!".format(name)
			raise FileNotFoundError(message)
		return fn

def read(filename):
	return Settings(filename)

if __name__ == "__main__":
	filename = "C:\\Users\\Deitrickc\\Google Drive\\Genomics\\api_files\\pipeline_project_options.txt"
	options = Settings(filename)

	print("Pipeline Folder: ", options['Pipeline Options']['pipeline folder'])
	print("variants-somatic", options['pipelineFolder', 'variants-somatic', 'A', 'B'])
	print(options['Programs']['GATK'])
	print(options['Globals'])


