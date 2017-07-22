import csv

import os
import shutil
import isodate

from github import systemtools
from github import filetools
from github import gdc_api

import datetime


class Workflow:
	def __init__(self, sample, options):

		##### Define commonly-used variables
		self.caller_name = self.__class__.__name__
		self.sample_log = options.getPipelineFile('sample log')
		self.targets    = sample['ExomeTargets']
		self.reference  = options['Reference Files']['reference genome']
		self.dbSNP      = options['Reference Files']['dbSNP']
		self.cosmic     = options['Reference Files']['COSMIC']

		self.verbose_level = options['globals']['verbose']

		self.program             = options['Programs'].get(self.caller_name.lower())
		self.gatk_program        = options['Programs']['GATK']
		self.max_cpu_threads     = options['Parameters']['MAX_CORES']
		self.max_memory_usage    = options['Parameters']['JAVA_MAX_MEMORY_USAGE']
		self.min_base_quality    = options['Parameters']['MIN_NUCLEOTIDE_QUALITY']
		self.min_mapping_quality = options['Parameters']['MIN_MAPPING_QUALITY']
		self.min_somatic_quality = options['Parameters']['SOMATIC_QUALITY']
		self.min_coverage        = options['Parameters']['MIN_COVERAGE']

		##### Define the paths and common partial filenames
		self.output_folder = options['variants-somatic', sample['PatientID'], self.caller_name]
		#self.output_folder = options['Pipeline Options']['somatic pipeline folder']

		self.base_prefix = "{normal}_vs_{tumor}.{prefix}".format(
			tumor   = sample['SampleID'], 
			normal  = sample['NormalID'],
			prefix = self.caller_name.lower()
		)
		self.abs_prefix = os.path.join(self.output_folder, self.base_prefix)
	
		self.temp_folder = options['temporary', sample['PatientID']]
		#self.temp_folder = options['Pipeline Options']['temporary folder']
		self.temp_files = list()
		self.full_output = []

		self.setCustomEnvironment(sample, options)
		filetools.checkDir(self.output_folder, True)
		filetools.checkDir(self.temp_folder, True)
		
		self.console_file = os.path.join(
			self.output_folder,
			self.caller_name + ".console_log.txt"
		)
		self.readme_filename = os.path.join(
			self.output_folder,
			self.caller_name + ".readme.txt"
		)


		if options['globals']['overwrite'] and os.path.exists(self.output_folder):
			shutil.rmtree(self.output_folder)
		print("\tRunning ", self.caller_name)
		if False:
			
			print("\tprogram location: ", self.program)
			print("\toutput folder: ", self.output_folder)
			print("\ttemp folder: ", self.temp_folder)
			print("\tcaller output prefix: ", self.base_prefix)

		self.createReadmeFile(sample)

		self.runCallerWorkflow(sample)

		self.renameOutputFiles()

		self.verifyOutputFiles(self.full_output)


		#self.updateSampleLog(sample, program_start, program_stop)

	def runCallerCommand(self, command, label, expected_output = None, **kwargs):
		"""
			Keyword Arguments:
				command:
				label:
				expected_output:
				output_filename:
		"""
		#expected_output = kwargs.get('expected_output')
		#output_filename = kwargs.get('output_filename')
		#label = kwargs.get('label', self.caller_name + ".runCallerCommand")

		if expected_output is None:
			expected_output = []
		elif isinstance(expected_output, str):
			expected_output = [expected_output]
	
		if 'verbose' in kwargs: kwargs.pop('verbose')
		self.addToReadme(command, label, expected_output)
		caller_session = systemtools.Terminal(
			command,
			label = label,
			expected_output = expected_output,
			command_filename = self.console_file,
			#output_filename = output_filename,
			verbose = ['all'],
			**kwargs
		)

		#self._verifySessionStatus(caller_session)
		self.verifyOutputFiles(expected_output)

		command_status = caller_session.getStatus()
		return command_status

	def generatePileup(self, bam_file, bam_name):
		""" Generates pileup files. If single is True, only
			one file will be generated.
		"""
		output_file = os.path.join(
			self.temp_folder,
			'{}.mpileup'.format(
				bam_name
			)
		)
		self.temp_files.append(output_file)
		#pileup_command = "samtools mpileup -q 1 -B -f {reference} {sample} > {output}".format(
		pileup_command = "samtools mpileup -q 1 -B -f {reference} {sample}".format(
			reference = self.reference,
			sample = bam_file,
			output = output_file,
			filename = output_file
		)

		label = "Generate MPileup File"
		self.runCallerCommand(pileup_command, label, output_file, output_filename = output_file)
		return output_file

	def renameOutputFiles(self):
		pass

	def createReadmeFile(self, sample):

		with open(self.readme_filename, 'a') as readme_file:
			readme_file.write(self.caller_name + '\n')
			readme_file.write("Started the caller at {0}\n".format(datetime.datetime.now().isoformat()))

			for key, value in sample.items():
				readme_file.write("{0:<15}{1}\n".format(key + ':', value))

			return self.readme_filename

	def addToReadme(self, command, label, expected_output):
		if not os.path.exists(self.readme_filename): return None
		current_datetime = datetime.datetime.now()
		line_len = 100
		num = line_len - (len(label) + len(current_datetime.isoformat()))
		linebreak = '#' * int(num / 2) + "{0} {1}".format(current_datetime.isoformat(), label) + "#" * int(num/2) + '\n'

		with open(self.readme_filename, 'a') as readme:
			readme.write(linebreak)
			readme.write(command + '\n')
			if isinstance(expected_output, list):
				expected_output = '|'.join(expected_output)
			readme.write("Expected Output: " + expected_output + '\n')
			readme.write('\n')

	def updateSampleLog(self, sample, program_start, program_stop):
		status = self.getCallerStatus()
		duration = program_stop - program_start
		iso_duration = isodate.duration_isoformat(duration)
		caller_log = {
			'patientId':  sample['PatientID'],
			'caller':    self.caller_name,
			'startTime': program_start.isoformat(),
			'stopTime':  program_stop.isoformat(),
			'duration':   duration,
			'isoDuration': iso_duration,
			'intermediateFiles': ', '.join(self.temp_files),
			'outputFiles':    ', '.join(self.full_output),
			'notes':      "",
			'status':     status
		}
		writeheaders = not os.path.exists(self.sample_log) or os.path.getsize(self.sample_log) == 0

		if not writeheaders:
			with open(self.sample_log, 'r') as file1:
				reader = csv.DictReader(file1, delimiter = '\t')
				fieldnames = reader.fieldnames
		else:
			fieldnames = sorted(caller_log.keys())

		with open(self.sample_log, 'a', newline = '') as file1:
			writer = csv.DictWriter(file1, fieldnames = fieldnames, delimiter = '\t')
			if writeheaders:
				writer.writeheader()
			writer.writerow(caller_log)

	def setCustomEnvironment(self, sample, options):
		pass
	
	def runCallerWorkflow(self, sample):
		pass

	def getCallerStatus(self):
		caller_failed = any(not os.path.exists(fn) for fn in self.full_output)
		status = not caller_failed

		return status

	def _verifySessionStatus(self, session):
		""" Verifies that the command completed properly 
			Parameters
			----------
				session: systemtools.Terminal
		"""

		session_status = session.status
		if not session_status:
			print(self.caller_name)
			print(session)
			raise FileNotFoundError("{}: The expected output files are not present!".format(self.caller_name))

	def verifyOutputFiles(self, expected_output):
		if isinstance(expected_output, str): expected_output = [expected_output]

		output_status = [(os.path.exists(fn), fn) for fn in expected_output]
		if any(not i[0] for i in output_status):
			print("The output files from {} were not generated correctly.".format(self.caller_name))
			print("Some output files were not created: ")
			for s, f in output_status:
				print("\t{}\t{}".format(s, f))
			raise FileNotFoundError()

	def _downloadGdcFile(self, sample, output_filename):
		case_id = sample['CaseID']
		case_info = gdc_api.request(case_id, 'cases')

		caller_name = self.caller_name.lower()
		if caller_name.endswith('2'): caller_name = caller_name[:-1]

		gdc_file_info = gdc_api.methods.extract.extractFile(case_info, 'raw', 'caller', caller_name)

		file_id = gdc_file_info['fileId']
		file_name = gdc_file_info['fileName'] #will be .gz

		# Download The Files
		download_output = os.path.join(self.output_folder, file_id, file_name)
		download_command = gdc_api.generateCommand(file_id, folder = self.output_folder)
		download_status = self.runCallerCommand(
			download_command, 'Downloading GDC File', download_output, verbose = [])

		# Extract the Files
		extract_command = "gunzip {fn}".format(fn = download_output)
		extract_output = os.path.splitext(download_output)[0]
		extract_status = self.runCallerCommand(
			extract_command, "Extracting GDC Files", extract_output)

		# Rename Files

		os.rename(extract_output, output_filename)

		self.verifyOutputFiles(output_filename)

		return extract_status

	def _overwriteExistingFiles(self):
		""" Deletes any existing files """

		shutil.rmtree(self.output_folder)
		filetools.checkDir(self.output_folder)
