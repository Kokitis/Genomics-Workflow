import csv
import datetime
import os
import shutil
import isodate
import gdc_api
import cmd_parser
import configparser
from pprint import pprint
GITHUB_FOLDER = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
import sys
sys.path.append(GITHUB_FOLDER)
print(GITHUB_FOLDER)
import pytools.systemtools as systemtools
import pytools.filetools as filetools
import pytools.tabletools as tabletools
import pytools.timetools as timetools
# --------------------------------- Global Variables ----------------------------
now = datetime.datetime.now

# The parent folder the pipeline will be run in.
PIPELINE_DIRECTORY = "/home/upmc/Documents/Variant_Discovery_Pipeline"
# File to save the console output to. Only used when the console output is supressed.
CONSOLE_LOG_FILE = ""
# File containing a test sample.
SAMPLE_LOG_FILE = os.path.join(PIPELINE_DIRECTORY, "0_config_files", "sample_logV2.tsv")
README_FILE     = os.path.join(PIPELINE_DIRECTORY, "0_readme_files", "readme.{0}.txt".format(now().isoformat()))
# Whether to use backwards-compatible filenames
BACKWARDS_COMPATIBLE = False
# Whether to overwrite any existing files
FORCE_OVERWRITE = False

DEBUG = True
VERBOSE_LEVEL = 4
# ----------------------------------------------------------------------------------------------------
# ------------------------------------ Set Up Global Functions ---------------------------------------
# ----------------------------------------------------------------------------------------------------

def generateTimestamp():
	timestamp = now().isoformat()
	timestamp = timestamp.split('.')[0]
	return timestamp


def getsize(path, total = True):
	sizes = []
	if os.path.isfile(path):
		sizes += [os.path.getsize(path)]
	elif os.path.isdir(path):
		dirsize = list()
		for fn in os.listdir(path):
			new_path = os.path.join(path, fn)
			dirsize += getsize(new_path, total = False)
		sizes += dirsize
	else:
		pass
	if total: sizes = sum(sizes)
	return sizes

def getPipelineFolder(step, patientId = None, caller_name = None):

	if step == 'variants-somatic':
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
	else:
		message = "'{}' is not a valid step in the pipeline!".format(step)
		raise ValueError(message)

	pipeline_folder = os.path.join(PIPELINE_DIRECTORY, *subfolders)
	filetools.checkDir(pipeline_folder, True)
	
	return pipeline_folder

# ----------------------------------------------------------------------------------------------------
# ------------------------------------- Variant Caller Pipeline --------------------------------------
# ----------------------------------------------------------------------------------------------------


class Caller:
	def __init__(self, sample, options):

		##### Define commonly-used variables
		program_start = now()
		self.caller_name = self.__class__.__name__
		self.reference  = options['Reference Files']['reference genome']
		self.dbSNP      = options['Reference Files']['dbSNP']
		self.cosmic     = options['Reference Files']['COSMIC']

		self.program             = options['Programs'].get(self.caller_name.lower())
		self.gatk_program        = options['Programs']['GATK']
		self.max_cpu_threads     = options['Parameters']['MAX_CORES']
		self.max_memory_usage    = options['Parameters']['JAVA_MAX_MEMORY_USAGE']
		self.min_base_quality    = options['Parameters']['MIN_NUCLEOTIDE_QUALITY']
		self.min_mapping_quality = options['Parameters']['MIN_MAPPING_QUALITY']
		self.min_somatic_quality = options['Parameters']['SOMATIC_QUALITY']
		self.min_coverage        = options['Parameters']['MIN_COVERAGE']

		##### Define the paths and common partial filenames
		self.output_folder = getPipelineFolder('variants-somatic', sample['PatientID'], self.caller_name)
		#self.output_folder = options['Pipeline Options']['somatic pipeline folder']
		
		self.base_prefix = "{normal}_vs_{tumor}.{prefix}".format(
			tumor   = sample['SampleID'], 
			normal  = sample['NormalID'],
			prefix = self.caller_name.lower()
		)
		self.abs_prefix = os.path.join(self.output_folder, self.base_prefix)
	
		self.temp_folder = getPipelineFolder('temporary', sample['PatientID'])
		#self.temp_folder = options['Pipeline Options']['temporary folder']
		self.temp_files = list()

		filetools.checkDir(self.output_folder, True)
		filetools.checkDir(self.temp_folder, True)
		
		self.setCustomEnvironment(sample, options)
		self.console_file = os.path.join(
			self.output_folder,
			self.caller_name + ".console_log.txt"
		)
		self.readme_filename = os.path.join(
			self.output_folder,
			self.caller_name + ".readme.txt"
		)
		
		if FORCE_OVERWRITE and os.path.exists(self.output_folder):
			shutil.rmtree(self.output_folder)

		if DEBUG:
			print("Running ", self.caller_name)
			print("\tprogram location: ", self.program)
			print("\toutput folder: ", self.output_folder)
			print("\ttemp folder: ", self.temp_folder)
			print("\tcaller output prefix: ", self.base_prefix)

		self.createReadmeFile(sample)

		self.runCallerWorkflow(sample)

		self.renameOutputFiles()

		self.verifyOutputFiles()
		
		##### Update the caller log
		program_stop = now()

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
		output_filename = kwargs.get('output_filename')
		#label = kwargs.get('label', self.caller_name + ".runCallerCommand")

		if expected_output is None:
			expected_output = []
		elif isinstance(expected_output, str):
			expected_output = [expected_output]
	

		self.addToReadme(command, label, expected_output)
		caller_session = systemtools.Terminal(
			command,
			label = label,
			expected_output = expected_output,
			command_filename = self.console_file,
			output_filename = output_filename,
			verbose = VERBOSE_LEVEL,
		)

		self._verifySessionStatus(caller_session)

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
			readme_file.write("Started the caller at {0}\n".format(now().isoformat()))

			for key, value in sample.items():
				readme_file.write("{0:<15}{1}\n".format(key + ':', value))

			return self.readme_filename

	def addToReadme(self, command, label, expected_output):
		if not os.path.exists(self.readme_filename): return None
		current_datetime = now()
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
			'status':     status,
			'commands':   '|'.join(self.command_list)
		}
		writeheaders = not os.path.exists(SAMPLE_LOG_FILE) or os.path.getsize(SAMPLE_LOG_FILE) == 0

		if not writeheaders:
			with open(SAMPLE_LOG_FILE, 'r') as file1:
				reader = csv.DictReader(file1, delimiter = '\t')
				fieldnames = reader.fieldnames
		else:
			fieldnames = sorted(caller_log.keys())

		with open(SAMPLE_LOG_FILE, 'a', newline = '') as file1:
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

	def verifyOutputFiles(self):

		output_status = [(os.path.exists(fn), fn) for fn in self.full_output]
		if any(not i[0] for i in output_status):
			print("Some output files were not created: ")
			for s, f in output_status:
				print("\t{}\t{}".format(s, f))
			raise FileNotFoundError()


class DepthOfCoverage(Caller):

	def setCustomEnvironment(self, sample, options):
		self.caller_name = "DepthOfCoverage"
		self.gene_list = os.path.join(
			getPipelineFolder('reference'),
			"GRCh38-hg38-NCBIRefSeq-UCSCRefSeq-allFieldsFromTable-WholeGene.txt.sorted.tsv"
		)
		self.final_output = self.abs_prefix + ".sample_gene_summary"
		full_output_suffixes = [
			"", 
			'.sample_cumulative_coverage_counts', 
			".sample_cumulative_coverage_proportions",
			".sample_gene_summary",
			".sample_interval_statistics", 
			".sample_interval_summary",
			".sample_statistics", 
			".sample_summary"
		]
		self.full_output = [self.abs_prefix + i for i in full_output_suffixes]

	def runCallerWorkflow(self, sample):
		self.status = self.determineDepthOfCoverage(sample, output_filename = self.final_output)

	def determineDepthOfCoverage(self, sample, output_filename):

		command = """java -jar {GATK} \
			--analysis_type DepthOfCoverage \
			--reference_sequence {reference} \
			--out {prefix} \
			--input_file {normal} \
			--input_file {tumor} \
			--calculateCoverageOverGenes {genes} \
			--intervals {targets} \
			--out {outfile}""".format(
				GATK =      self.gatk_program,
				reference = self.reference,
				prefix =    self.abs_prefix,
				normal =    sample['NormalBAM'],
				tumor =     sample['TumorBAM'],
				targets =   sample['ExomeTargets'],
				genes =     self.gene_list,
				outfile =	output_filename
			)
		label = "DepthofCoverage"
		result = self.runCallerCommand(command, label, output_filename)

		return result


class MuSE(Caller):

	def setCustomEnvironment(self, sample, options):
		self.caller_name = "MuSE"
		#Define relevant filenames
		self.call_output = '.'.join(self.abs_prefix.split('.')[:-1]) + '.MuSE.txt' #remove caller name form prefix
		self.sump_output = os.path.splitext(self.call_output)[0] + ".vcf"

		self.full_output = [self.call_output, self.sump_output]
		self.final_output= self.sump_output

	def runCallerWorkflow(self, sample):
		#Run the caller commands
		self.call_output_status = self.runMuseCall(sample)
		self.sump_output_status = self.runMuseSump()

	def runMuseCall(self, sample):

		call_command = "{program} call -O {prefix} -f {reference} {tumor} {normal}".format(
			program 	= self.program,
			reference 	= self.reference,
			prefix 		= self.call_output.replace('.MuSE.txt', ''),
			tumor 		= sample['TumorBAM'],
			normal 		= sample['NormalBAM'])
		label = "MuSE Call"
		output_result = self.runCallerCommand(call_command, label, self.call_output)

		return output_result

	def runMuseSump(self):
		
		sump_command = "{program} sump -I {call_output} -E -D {dbSNP} -O {output}".format(
			program 	= self.program,
			prefix 		= self.abs_prefix,
			dbSNP 		= self.dbSNP,
			call_output = self.call_output,
			output 		= self.sump_output)
		label = "MuSE Sump"
		output_result = self.runCallerCommand(sump_command, label, self.sump_output)

		return output_result


class MuTect(Caller):

	def setCustomEnvironment(self, sample, options):
		self.caller_name = "MuTect"
		self.java_program = "/home/upmc/Downloads/jre1.7.0_80/bin/java"
		self.variant_file = self.abs_prefix + '.mutect117.vcf'
		self.coverage_file= self.abs_prefix + '.mutect117.coverage.wiggle.txt'
		self.final_output = self.coverage_file
		self.full_output = [self.variant_file, self.coverage_file]

	def runCallerWorkflow(self, sample):
		caller_status = self.runMutect(sample)

	def runMutect(self, sample):
		# MuTect 1.1.7 requires java 7, rather than the currently installed java 8
		
		command = """{java} -jar {program} \
			--analysis_type MuTect \
			--reference_sequence {reference} \
			--dbsnp {dbsnp} \
			--cosmic {cosmic} \
			--intervals {targets} \
			--input_file:normal {normal} \
			--input_file:tumor {tumor} \
			--out {output} \
			--coverage_file {coverage}""".format(
				java =      self.java_program,
				program =   self.program,
				reference = self.reference,
				dbsnp =     self.dbSNP,
				cosmic =    self.cosmic,
				targets =   sample['ExomeTargets'],
				output =    self.variant_file,
				coverage =  self.coverage_file,
				normal =    sample['NormalBAM'],
				tumor =     sample['TumorBAM'])
		label = "Running the original mutect..."
		output_result = self.runCallerCommand(command, label, self.variant_file)
		return output_result


class MuTect2(Caller):

	def setCustomEnvironment(self, sample, options):
		self.caller_name = "MuTect2"
		self.final_output = self.abs_prefix + '.vcf'
		self.full_output = [self.final_output]

	def runCallerWorkflow(self, sample):
		self.status = self.runMutect2(sample)
		

	def runMutect2(self, sample):
		mutect2_command = """java {memory} -jar {GATK} \
			-T MuTect2 \
			-R {reference} \
			-L {targets} \
			-I:normal {normal} \
			-I:tumor {tumor} \
			--dbsnp {dbSNP} \
			--cosmic {cosmic} \
			--out {output}""".format(
			GATK        = self.gatk_program, 
			memory      = self.max_memory_usage,
			dbSNP       = self.dbSNP,
			cosmic      = self.cosmic,
			reference   = self.reference, 
			normal      = sample['NormalBAM'], 
			tumor       = sample['TumorBAM'], 
			targets     = sample['ExomeTargets'],
			output      = self.final_output)
		label = "MuTect2 Variant Discovery"
		output_result = self.runCallerCommand(mutect2_command, label, self.final_output)

		return output_result


class SomaticSniper(Caller):

	def setCustomEnvironment(self, sample, options):
		# Define the relevant filenames
		self.caller_name = "SomaticSniper"
		somaticsniper_folder= self.program # for readability
		script_folder 		= os.path.join(somaticsniper_folder, 'src', 	'scripts')
		self.program 		= os.path.join(somaticsniper_folder, 'build', 	'bin', 	'bam-somaticsniper')
		self.readcount_program 	= options['Programs']['bam-readcount']
		self.samtools_program 	= options['Programs']['samtools-0.1.6']


		self.snpfilter_script= os.path.join(script_folder, 'snpfilter.pl')
		self.readcount_script= os.path.join(script_folder, 'prepare_for_readcount.pl')
		self.hc_script 		=  os.path.join(script_folder, 'highconfidence.pl')
		self.fpfilter  		=  os.path.join(script_folder, 'fpfilter.pl')

		self.raw_variants 	= self.abs_prefix + '.vcf'
		self.hq_variants 	= self.abs_prefix + '.hq.vcf'
		self.lq_variants 	= self.abs_prefix + '.lq.vcf'		

	def runCallerWorkflow(self, sample):
		# Will print "Couldn't find single-end mapping quality. Check to see if the SM tag is in BAM."
		# This doesn't invalidate results, but try not to use single-end mapping quality in output

		self.variant_discovery_status = self.runVariantDiscovery(sample)
		# Generate pileup files
		normal_pileup_file = self.generatePileup(sample['NormalBAM'], 	'normal') #returns pileup filename
		tumor_pileup_file  = self.generatePileup(sample['TumorBAM'], 	'tumor')

		# Filter LOH
		_intermediate_loh_filtered_output = self.abs_prefix + ".SNPfilter.intermediate"
		loh_filtered_output = self.abs_prefix + ".SNPfilter.final"
		_intermediate_file  = self.removeLOH(self.raw_variants, normal_pileup_file, _intermediate_loh_filtered_output)
		loh_filtered_output = self.removeLOH(_intermediate_loh_filtered_output, tumor_pileup_file, loh_filtered_output)
		# loh_filtered_output = _intermediate_file
		
		readcount_status 			= self.readcounts(loh_filtered_output, sample['TumorBAM'])
		readcounts 					= readcount_status['outputFiles']
		false_positive_status 		= self.removeFalsePositives(loh_filtered_output, readcounts)
		false_positive_output 		= false_positive_status['outputFiles']
		high_confidence_status	 	= self.calculateConfidence(false_positive_output)

		self.full_output = [self.raw_variants, self.hq_variants, self.lq_variants]
		self.final_output = self.hq_variants
	
	def runVariantDiscovery(self, sample):
		# -------------------------------- Variant Discovery Command -----------------------------
		somaticsniper_command = """{program} \
			-q 1 \
			-Q 15 \
			-s 0.01 \
			-T 0.85 \
			-N 2 \
			-r 0.001 \
			-G \
			-L \
			-n NORMAL \
			-t TUMOR \
			-F vcf \
			-f {reference} {tumor} {normal} {outputfile}""".format(
				program     = self.program,
				reference   = self.reference,
				tumor       = sample['TumorBAM'],
				normal      = sample['NormalBAM'],
				outputfile  = self.raw_variants
			)
		label = "SomaticSniper Variant Discovery"
		output_result = self.runCallerCommand(somaticsniper_command, label, self.raw_variants)
		return output_result

	def generatePileup(self, bam_file, bam_name):

		raw_pileup_file = os.path.join(self.temp_folder, "{0}.somaticsniper.pileup".format(bam_name))
		filtered_output_file = raw_pileup_file + '.pileup'
		
		#samtools_command = "{samtools} pileup -cvi -f {reference} {bam} > {pileup}".format(
		samtools_command = "{samtools} pileup -cvi -f {reference} {bam}".format(
			samtools = self.samtools_program,
			reference = self.reference,
			bam = bam_file,
			pileup = raw_pileup_file)

		#filter_command = """samtools.pl varFilter -Q {basequality} {inputpileup} > {outputpileup}""".format(
		filter_command = """samtools.pl varFilter -Q {basequality} {inputpileup}""".format(
			basequality = self.min_base_quality,
			inputpileup = raw_pileup_file,
			outputpileup = filtered_output_file)
		# samtools.pl varFilter raw.pileup | awk '$6>=20' > final.pileup
		label = "SomaticSniper: Generate Pileup File" 
		self.runCallerCommand(samtools_command, label, raw_pileup_file, output_filename = raw_pileup_file)

		label = "SomaticSniper: Filter Pileup File"
		self.runCallerCommand(filter_command, label, filtered_output_file, output_filename = filtered_output_file)

		return filtered_output_file

	def removeLOH(self, vcf_file, pileup_file, output_file):
		# -------------------------- Filter and remove LOH --------------------------------
		command = """perl {snpfilter} \
			--snp-file {vcf} \
			--indel-file {pileup} \
			--out-file {output}""".format(
				snpfilter = self.snpfilter_script,
				vcf 	= vcf_file,
				pileup 	= pileup_file,
				output 	= output_file
			)
		label = "SomaticSniper: removeLOH"
		self.runCallerCommand(command, label, output_file)
		return output_file

	def readcounts(self, loh_file, tumor_bam):
		""" Expected Output: 
		"""
		prepare_readcount_output 	= loh_file + '.pos'
		readcount_output 			= loh_file + '.readcounts.rc'
		# Prepare readcounts
		pr_command = """perl {script} \
			--snp-file {inputfile} \
			--out-file {output}""".format(
				script 		= self.readcount_script,
				inputfile 	= loh_file,
				output 		= prepare_readcount_output
			)
		label = "SomaticSniper: Prepare Readcounts"
		preparation_result = self.runCallerCommand(pr_command, label, prepare_readcount_output)
		# Readcounts

		readcount_command = """{program} \
			-b {mbq} \
			-q 1 \
			-f {reference} \
			-l {proutput} \
			{tumor} > {output}""".format(
				program 	= self.readcount_program, # different from readcount script
				mbq 		= self.min_base_quality,
				reference 	= self.reference,
				proutput 	= prepare_readcount_output,
				tumor 		= tumor_bam,
				output 		= readcount_output
			)
		label = "SomaticSniper: Generate Readcounts"
		# Note: The output for this command will be saved in the console log file.
		readcount_result = self.runCallerCommand(readcount_command, label, expected_output = readcount_output, output_filename = readcount_output)
		#self._captureReadcountOutput(readcount_output)

		return readcount_result

	def _captureReadcountOutput(self, output_file):
		""" Since the terminal output is saved to the console file now,
			the readcount output must be scraped from that file.
		"""
		if os.path.exists(output_file): return None
		readcounts = list()
		console_file = self._searchForConsoleFile()
		with open(console_file, 'r') as file1:
			for line in file1.read().splitlines():
				if line[:3] == 'chr': readcounts.append(line)
		with open(output_file, 'w') as file2:
			for line in readcounts:
				file2.write(line + '\n')

	def _searchForConsoleFile(self):
		if not os.path.exists(self.console_file):
			candidates = sorted(os.listdir(self.output_folder))
			if len(candidates) == 0:
				candidate = self.console_file
			else:
				candidate = candidates[-1]
		else:
			candidate = self.console_file

		return candidate

	def removeFalsePositives(self, loh_file, readcounts):
		false_positive_output = loh_file + '.fp_pass'
		fp_command = """perl {script} \
			--snp-file {snpfile} \
			--readcount-file {readcounts}""".format(
			script 		= self.fpfilter,
			snpfile 	= loh_file,
			readcounts 	= readcounts
		)

		label = "SomaticSniper: Remove False Positives"
		output_status = self.runCallerCommand(fp_command, label, false_positive_output)
		return output_status

	def calculateConfidence(self, vcf_file):
		command = """perl {script} \
			--snp-file {snpfile} \
			--min-mapping-quality {mmq} \
			--min-somatic-score {ss} \
			--lq-output {lq} \
			--out-file {hq}""".format(
				script 	= self.hc_script,
				mmq 	= self.min_mapping_quality,
				ss 		= self.min_somatic_quality,
				snpfile	= vcf_file,
				hq 		= self.hq_variants,
				lq 		= self.lq_variants
		)

		label = "SomaticSniper: Filter lq Variants"
		output_result = self.runCallerCommand(command, label, [self.hq_variants, self.lq_variants])
		return output_result


class Strelka(Caller):

	def setCustomEnvironment(self, sample, options):
		self.caller_name = "Strelka"
		patient_folder = os.path.dirname(self.output_folder)

		strelka_folder 						= os.path.join(patient_folder, 'Strelka')
		self.results_folder 				= os.path.join(patient_folder, 'Strelka', 	'results')
		self.strelka_project_config_file 	= os.path.join(patient_folder, 'strelka_configuration.ini')
		self.config_script 					= os.path.join(self.program,   'bin', 		'configureStrelkaWorkflow.pl')
		
		full_output = [
			'all.somatic.indels.vcf',
			'all.somatic.snvs.vcf',
			'passed.somatic.indels.vcf',
			'passed.somatic.snvs.vcf'
		]
		self.full_output = [os.path.join(self.results_folder, fn) for fn in full_output]
		self.final_output = self.full_output[0]
	
	def runCallerWorkflow(self, sample):
		# Strelka won't run if the output folder already exists.
		if not os.path.exists(self.final_output) and os.path.exists(self.output_folder):
			shutil.rmtree(self.output_folder)

		# --------------------------------- Configure Strelka --------------------------------
		self.generateStrelkaConfigFile(self.strelka_project_config_file)
		
		self.configureStrelka(sample, self.strelka_project_config_file)
		
		self.runStrelka()


	def generateStrelkaConfigFile(self, configuration_file):
		strelka_configuration_options = self.getStrelkaConfiguration()
		print("Generating Strelka Config File...")
		with open(configuration_file, 'w') as file1:
			# strelka_configuration_options = [i.split('\t')[0].strip() for i in strelka_configuration_options.s]
			file1.write('[user]\n')
			for row in strelka_configuration_options.split('\n'):
				#row = [i for i in row.split('\t') if '=' in i]
				#if len(row) > 0:
				file1.write(row.strip() + '\n')
		return configuration_file

	def configureStrelka(self, sample, config_ini):
		strelka_config_command = """perl {script} \
			--tumor={tumor} \
			--normal={normal} \
			--ref={reference} \
			--config={config} \
			--output-dir={output_dir}""".format(
				script 		= self.config_script,
				tumor 		= sample['TumorBAM'],
				normal 		= sample['NormalBAM'],
				reference 	= self.reference,
				config 		= config_ini,
				output_dir 	= self.output_folder
		)
		#print(strelka_config_command)
		output_result = self.runCallerCommand(strelka_config_command, "Configure Strelka", [self.output_folder])
		return output_result

	def runStrelka(self):
		#print(os.path.exists(self.output_folder), '\t', self.output_folder)
		output_suffixes = [
			'all.somatic.indels.vcf',
			'all.somatic.snvs.vcf',
			'passed.somatic.indels.vcf',
			'passed.somatic.snvs.vcf'
		]
		output_files = [os.path.join(self.results_folder, i) for i in output_suffixes]
		#print(self.output_folder)
		strelka_run_command = "make -j {threads} -C {outputdir}".format(
			outputdir = self.output_folder, 
			threads = self.max_cpu_threads
		)

		label = "Strelka Run"
		output_result = self.runCallerCommand(strelka_run_command, label, output_files)

	@staticmethod
	def getStrelkaConfiguration():
		strelka_configuration_options = """
			isSkipDepthFilters = 1						isSkipDepthFilters should be set to 1 to skip depth filtration for whole exome or other targeted sequencing data
			maxInputDepth = 10000						strelka will not accept input reads above this depth. Set this value <= 0 to disable this feature.
			depthFilterMultiple = 3.0					If the depth filter is not skipped, all variants which occur at a depth greater than depthFilterMultiple*chromosome mean depth will be filtered out.
			snvMaxFilteredBasecallFrac = 0.4            Somatic SNV calls are filtered at sites where greater than this fraction of basecalls have been removed by the mismatch density filter in either sample.
			snvMaxSpanningDeletionFrac = 0.75           Somatic SNV calls are filtered at sites where greater than this fraction of overlapping reads contain deletions which span the SNV call site.
			indelMaxRefRepeat = 8                       Somatic indel calls are filtered if they represent an expansion or contraction of a repeated pattern with a repeat count greater than indelMaxRefRepeat in the reference
			indelMaxWindowFilteredBasecallFrac = 0.3    Somatic indel calls are filtered if greater than this fraction of basecalls in a window extending 50 bases to each side of an indel's call position have been removed by the mismatch density filter.
			indelMaxIntHpolLength = 14                  Somatic indels are filtered if they overlap ’interrupted homopolymers’ greater than this length. 
			ssnvPrior = 0.000001                        prior probability of a somatic snv or indel
			sindelPrior = 0.000001
			ssnvNoise = 0.0000005                       probability of an snv or indel noise allele 
			sindelNoise = 0.000001
			ssnvNoiseStrandBiasFrac = 0.5               Fraction of snv noise attributed to strand-bias. It is not recommended to change this setting.
			minTier1Mapq = 20                           minimum MAPQ score for PE reads at tier1
			minTier2Mapq = 0                            minimum MAPQ score for PE and SE reads at tier2
			ssnvQuality_LowerBound = 15                 Somatic quality score (QSS_NT, NT=ref) below which somatic SNVs are marked as filtered
			sindelQuality_LowerBound = 30               Somatic quality score (QSI_NT, NT=ref) below which somatic indels are marked as filtered
			isWriteRealignedBam = 0                     Optionally write out read alignments which were altered during the realignment step. Location: ${ANALYSIS_DIR}/realigned/{normal,tumor}.realigned.bam
			binSize = 25000000                          Jobs are parallelized over segments of the reference genome no larger than this size"""

		strelka_configuration_options = """
			isSkipDepthFilters = 1
			maxInputDepth = 10000
			depthFilterMultiple = 3.0
			snvMaxFilteredBasecallFrac = 0.4
			snvMaxSpanningDeletionFrac = 0.75
			indelMaxRefRepeat = 8
			indelMaxWindowFilteredBasecallFrac = 0.3
			indelMaxIntHpolLength = 14
			ssnvPrior = 0.000001
			sindelPrior = 0.000001
			ssnvNoise = 0.0000005
			sindelNoise = 0.000001
			ssnvNoiseStrandBiasFrac = 0.5
			minTier1Mapq = 20
			minTier2Mapq = 0
			ssnvQuality_LowerBound = 15
			sindelQuality_LowerBound = 30
			isWriteRealignedBam = 0
			binSize = 25000000"""
		return strelka_configuration_options

	def renameOutputFiles(self):
		for source in self.full_output:
			folder, basename = os.path.split(source)
			if "all.somatic.indels" in basename:
				basename = "all.somatic.indels"
			elif "all.somatic.snvs" in basename:
				basename = "all.somatic.snps"
			elif "passed.somatic.indels" in basename:
				basename = "passed.somatic.indels"
			else:
				basename = "passed.somatic.snps"

			destination = self.abs_prefix + basename + ".strelka.vcf"

			shutil.copyfile(source, destination)


class Varscan(Caller):

	def setCustomEnvironment(self, sample, options):
		self.caller_name = "Varscan"
		self.raw_snps   = self.abs_prefix + '.snp.vcf'
		self.raw_indels = self.abs_prefix + '.indel.vcf'

		self.somatic_hc = self.abs_prefix + '.snp.Somatic.hc.vcf'
		self.somatic_lc = self.abs_prefix + '.snp.Somatic.lc.vcf'
		self.germline   = self.abs_prefix + '.snp.Germline.vcf'
		self.germline_hc= self.abs_prefix + '.snp.Germline.hc.vcf'
		self.loh        = self.abs_prefix + '.snp.LOH.vcf'
		self.loh_hc     = self.abs_prefix + '.snp.LOH.hc.vcf'

		self.full_output = [
			self.raw_snps, 
			self.raw_indels,
			self.somatic_hc,
			#self.somatic_lc,
			self.germline,
			self.germline_hc,
			self.loh,
			self.loh_hc]

		self.final_output = self.somatic_hc
	
	def runCallerWorkflow(self, sample):

		pileup_status = self.generateSinglePileup(sample)
		pileup_file = pileup_status['outputFiles']
		variant_discovery_status = self.runSingleVariantDiscovery(pileup_file)
		#self.runDoubleVariantDiscovery(sample)
		processing_status = self.postProcessing()

	def generateSinglePileup(self, sample):
		pileup = os.path.join(
			self.temp_folder,
			"{}_vs_{}.intermediate.mpileup".format(
				sample['NormalID'],
				sample['SampleID']
			)
		)
		# command = """samtools mpileup -f {reference} -q 1 -B {normal} {tumor} > {pileup}""".format(
		command = """samtools mpileup \
					-f {reference} \
					-q 1 \
					-B {normal} \
					{tumor}""".format(
			reference = self.reference,
			normal 	= sample['NormalBAM'],
			tumor 	= sample['TumorBAM'],
			pileup 	= pileup
		)
		#output_result = self.runCallerCommand(command, "GenerateSinglePileup", pileup, output_filename = pileup, show_output = True)
		#Need to use os.system due to how samtools generates mpileups and the size of the output.
		if not os.path.exists:
			os.system(command + " > {}".format(pileup))
		output_result = {'outputFiles': pileup}
		pprint(output_result)
		return output_result

	def runSingleVariantDiscovery(self, pileup):
		expected_output = [self.raw_snps, self.raw_indels]
		command = """java -jar {program} somatic {pileup} {output} \
			--mpileup 1 \
			--min-coverage 8 \
			--min-coverage-normal 8 \
			--min-coverage-tumor 6 \
			--min-var-freq 0.10 \
			--min-freq-for-hom 0.75 \
			--normal-purity 1.0 \
			--tumor-purity 1.00 \
			--p-value 0.99 \
			--somatic-p-value 0.05 \
			--strand-filter 0 \
			--output-vcf""".format(
				program = self.program,
				pileup  = pileup,
				output  = self.abs_prefix)
		print(command)
		output_result = self.runCallerCommand(command, "RunSingleVariantDiscovery", expected_output)

		return output_result
	
	def runDoubleVariantDiscovery(self, sample):

		normal_pileup = self.generatePileup(sample['NormalBAM'], sample['NormalID'])
		tumor_pileup = self.generatePileup(sample['TumorBAM'], sample['SampleID'])
		command = """java {memory} -jar {varscan} somatic {normal} {tumor} \
			--output-snp {snp} \
			--output-indel {indel} \
			--output-vcf 1""".format(
				varscan = self.program,
				memory  = self.max_memory_usage,
				normal  = normal_pileup,
				tumor   = tumor_pileup,
				snp     = self.raw_snps,
				indel   = self.raw_indels,
				mc      = self.min_coverage)
		output_files = [self.raw_snps, self.raw_indels]
		label = "Variant Discovery"
		status = self.runCallerCommand(command, label, output_files)

		return status

	def postProcessing(self):
		expected_output = [self.somatic_hc, self.germline, self.germline_hc, self.loh, self.loh_hc]

		process_command = """java -jar {varscan} processSomatic {vcf} \
			--min-tumor-freq 0.10 \
			--max-normal-freq 0.05 \
			--p-value 0.07""".format(
				memory  = self.max_memory_usage,
				varscan = self.program,
				vcf     = self.raw_snps
			)
		label = "Varscan Postprocessing"
		output_result = self.runCallerCommand(process_command, label, expected_output)
		return output_result

	def renameOutputFiles(self):
		return None
		for source in os.listdir(self.output_folder):
			#destination = os.path.join(self.output_folder, source.replace('varscan', 'raw'))
			destination = os.path.splitext(source)[0] + ".varscan.vcf"
			abs_source = os.path.join(self.output_folder, source)
			shutil.copyfile(abs_source, destination)


class HaplotypeCaller(Caller):
	def __init__(self, sample, options, input_bam = None):
		self.input_bam = input_bam
		self.options = options
		super().__init__(sample, options)

	def setCustomEnvironment(self, sample, options):
		self.caller_name = "HaplotypeCaller"

	def runCallerWorkflow(self, sample):

		if self.input_bam is None:
			call_status = self.dnaWorkflow(sample, self.options)
		else:
			call_status = self.rnaWorkflow(sample, self.options)

		self.final_output = call_status['outputFiles']
		self.full_output = [self.final_output]

	def dnaWorkflow(self, sample, options):
		raw_variant_file = self.abs_prefix + ".RNA.raw_variants.vcf"
		call_status = self.dnaVariantDiscovery(sample, raw_variant_file)
		return call_status

	def rnaWorkflow(self, sample, options):
		raw_variant_file 	= self.abs_prefix + ".RNA.raw_variants.vcf"
		variant_file 		= self.abs_prefix + ".RNA.filtered_variants.vcf"
		call_status = self.rnaVariantDiscovery(self.input_bam, raw_variant_file)
		filter_status = self.filterVariants(call_status['outputFiles'], variant_file)
		return filter_status

	def rnaVariantDiscovery(self, bam_file, output_filename):
		command = """java -jar {GATK} \
			--analysis_type HaplotypeCaller \
			--reference_sequence {reference} \
			--input_file {sample} \
			--dbsnp {dbSNP} \
			--dontUseSoftClippedBases \
			--out {output}""".format(
				GATK        = self.gatk_program,
				reference   = self.reference,
				sample      = bam_file,
				dbSNP       = self.dbSNP,
				output      = output_filename)
		label = 'RNA Variant Discovery'
		output_result = self.runCallerCommand(command, label, output_filename)
		return output_result
	
	def filterVariants(self, variant_file, output_filename):
		command = """java -jar {GATK} \
			--analysis_type VariantFiltration \
			--reference_sequence {reference} \
			--variant {inputfile} \
			--clusterSize 3 \
			--clusterWindowSize 35 \
			--filterName FS --filterExpression \"FS > 30.0\" \
			--filterName QD --filterExpression \"QD < 2.0\" \
			--out {output}""".format(
				GATK        = self.gatk_program,
				reference   = self.reference,
				inputfile   = variant_file,
				output      = output_filename)
		label = 'Variant Filtering'
		output_result = self.runCallerCommand(command, label, output_filename)
		return output_result
	
	def dnaVariantDiscovery(self, sample, output_filename):
		command = """java -jar {GATK} \
			--analysis_type HaplotypeCaller \
			--reference_sequence {reference} \
			--input_file {normal} \
			--input_file {tumor} \
			--intervals {targets} \
			--num_cpu_threads_per_data_thread {threads} \
			--dbsnp {dbSNP} \
			--out {output}""".format(
				GATK        = self.gatk_program,
				reference   = self.reference,
				normal      = sample['NormalBAM'],
				tumor       = sample['TumorBAM'],
				targets     = sample['ExomeTargets'],
				dbSNP       = self.dbSNP,
				output      = output_filename,
				threads     = self.max_cpu_threads
			)
		label = 'DNA Variant Discovery'
		output_result = self.runCallerCommand(command, label, output_filename)
		return output_result


class BaseQualityScoreRecalibration(Caller):

	def setCustomEnvironment(self, sample, options):
		self.caller_name = "BaseQualityScoreRecalibration"
		
		self.output_folder = os.path.join(
			getPipelineFolder('variants-rna', sample['PatientID'])
		)

		self.input_bam = sample['RNABAM']

		self.recalibration_table    = os.path.join(
			self.output_folder, sample['SampleID'] + ".RNA.recalibration_data.table"
		)
		
		self.covariate_table        = os.path.join(
			self.output_folder, sample['SampleID'] + ".RNA.covariate_data.table"
		)
		
		self.recalibration_plots    = os.path.join(
			self.output_folder, sample['SampleID'] + ".RNA.recalibration_plots.pdf"
		)
		
		self.cigar_bam              = os.path.join(
			self.output_folder, sample['SampleID'] + ".RNA.cigar.bam"
		)
		
		self.realigned_bam          = os.path.join(
			self.output_folder, sample['SampleID'] + ".RNA.recalibrated.bam"
		)

		self.final_output   = self.realigned_bam
		self.full_output    = [
			self.recalibration_table,
			self.covariate_table,
			self.recalibration_plots,
			self.cigar_bam,
			self.realigned_bam
		]

	def runCallerWorkflow(self, sample):

		
		raw_bam 		     = self._getRNABAM(sample)
		cigar__status        = self.splitCigarReads()
		recalibration_status = self.generateRecalibrationTable()
		realign_status       = self.recalibrateBAM()
		covariate_status     = self.generateCovariateTable()
		recalibration_status = self.generateRecalibrationPlots()

	@staticmethod
	def _getRNABAM(sample):
		return sample['RNABAM']
	
	def splitCigarReads(self):
		command = """java -jar {GATK} \
			--analysis_type SplitNCigarReads \
			--reference_sequence {reference} \
			--input_file {inputbam} \
			--out {outputbam} \
			--read_filter ReassignOneMappingQuality -RMQF 255 -RMQT 60 \
			--unsafe ALLOW_N_CIGAR_READS""".format(
				GATK        = self.gatk_program,
				reference   = self.reference,
				inputbam    = self.input_bam,
				outputbam   = self.cigar_bam
			)
		label = "Split Cigar Reads"
		output_result = self.runCallerCommand(command, label, self.cigar_bam)
		return output_result

	def generateRecalibrationTable(self):
		command = """java -jar {GATK} \
			--analysis_type BaseRecalibrator \
			--reference_sequence {reference} \
			--input_file {bam} \
			--knownSites {dbSNP} \
			--out {output}""".format(
				GATK        = self.gatk_program,
				reference   = self.reference,
				dbSNP       = self.dbSNP,
				bam         = self.cigar_bam,
				output      = self.recalibration_table)
		label = "Generate Recalibration Table"
		output_result = self.runCallerCommand(command, label, self.recalibration_table)
		return output_result

	def recalibrateBAM(self):
		command = """java -jar {GATK} \
			--analysis_type PrintReads \
			--reference_sequence {reference} \
			--input_file {inputfile} \
			--BQSR {table} \
			--out {output}""".format(
				GATK 		= self.gatk_program,
				reference 	= self.reference,
				inputfile 	= self.cigar_bam,
				output 		= self.realigned_bam,
				table 		= self.recalibration_table)
		label = "Recalibrate BAM"
		output_result = self.runCallerCommand(command, label, self.realigned_bam)
		return output_result
	
	def generateCovariateTable(self):
		command = """java -jar {GATK} \
			--analysis_type BaseRecalibrator \
			--reference_sequence {reference} \
			--input_file {realigned_bam} \
			--BQSR {table} \
			--knownSites {dbSNP} \
			--out {output}""".format(
				GATK        = self.gatk_program,
				reference   = self.reference,
				dbSNP       = self.dbSNP,
				table       = self.recalibration_table,
				realigned_bam = self.realigned_bam,
				output      = self.covariate_table)
		label = "Generate Covariate Table"
		output_result = self.runCallerCommand(command, label, self.covariate_table)
		return output_result
	
	def generateRecalibrationPlots(self):
		command = """ java -jar {GATK} \
			--analysis_type AnalyzeCovariates \
			--reference_sequence {reference} \
			--beforeReportFile {before} \
			--afterReportFile {after} \
			--plotsReportFile {plots}""".format(
				GATK        = self.gatk_program,
				reference   = self.reference,
				before      = self.recalibration_table,
				after       = self.covariate_table,
				plots       = self.recalibration_plots)
		label = "Generate Recalibration Plots"
		output_result = self.runCallerCommand(command, label, self.recalibration_plots)
		return output_result


# ----------------------------------------------------------------------------------------------------
# ------------------------------------ Copynumber Tool Functions -------------------------------------
# ----------------------------------------------------------------------------------------------------


class VarscanCopynumber(Caller):

	def setCustomEnvironment(self, sample, options):
		self.caller_name = "Varscan"

		self.output_folder = getPipelineFolder('variants-copynumber', sample['PatientID'], self.caller_name)

		self.base_prefix = "{normal}_vs_{tumor}.{prefix}".format(
			tumor=sample['SampleID'],
			normal=sample['NormalID'],
			prefix=self.caller_name.lower()
		)

		self.abs_prefix = os.path.join(
			self.output_folder,
			self.base_prefix
		)

		self.rscript_filename   = os.path.join(self.output_folder, "{0}.varscan_CBS.r".format(sample['PatientID']))
		self.copynumber_output  = self.abs_prefix + '.copynumber'           # [prefix].copynumber
		self.called_copynumbers = self.copynumber_output + '.called'    # [prefix].copynumber.called
		self.called_homdels     = self.called_copynumbers + '.homdel'   # [prefix].copynumber.called.homdel
		self.copynumber_segments =self.called_copynumbers + '.segments' # [prefix].copynumber.called.segments
		self.full_output = [
			self.copynumber_output,
			self.called_copynumbers,
			self.called_homdels,
			self.copynumber_segments
		]
		self.final_output = self.copynumber_segments

	def runCallerWorkflow(self, sample):
		normal_pileup   = self.generatePileup(sample['NormalBAM'], sample['NormalID'])
		tumor_pileup    = self.generatePileup(sample['TumorBAM'], sample['SampleID'])

		copynumber_ratio_status = self.callCopynumberRatios(normal_pileup, tumor_pileup)
		copynumber_caller_status  = self.copyCaller()

		segmentation_status = self.circularBinarySegmentation()

	def circularBinarySegmentation(self):
		# NEED to strip first row of table before r script
		"""
				circular_binary_segmentation(
			varscan_prefix + '.copynumber.called', 
			varscan_prefix + '.copynumber.called.segments',
			rscript_file)
		"""
		script_contents = """library(DNAcopy)
			cn <- read.table("{input_file}",header=TRUE)
			CNA.object <-CNA( genomdat = cn[,6], chrom = cn[,1], maploc = cn[,2], data.type = 'logratio')
			CNA.smoothed <- smooth.CNA(CNA.object)
			segs <- segment(CNA.smoothed, verbose=0, min.width=2)
			segs2 = segs$output
			write.table(segs2[,2:6], file="{output_file}", row.names=F, col.names=F, quote=F, sep="\\t")""".format(
			input_file = self.called_copynumbers,
			output_file = self.copynumber_segments)

		with open(self.rscript_filename, 'w') as file1:
			for line in script_contents.splitlines():
				file1.write(line + '\n')

		command = """Rscript {script}""".format(script = self.rscript_filename)
		label = "Circulr Binary Segmentation"
		output_result = self.runCallerCommand(command, label, self.copynumber_segments)
		return output_result

	def callCopynumberRatios(self, normal_pileup, tumor_pileup):
		# -------------------------- Varscan (copynumber - caller) -----------------------------
		copynumber_command = """java {memory} -jar {varscan} copynumber {normal} {tumor} {prefix} \
			--min-base_qual {mbq} \
			--min-map-qual {mmq}""".format(
				mmq 	= self.min_mapping_quality,
				mbq 	= self.min_base_quality,
				normal 	= normal_pileup,
				tumor 	= tumor_pileup,
				varscan = self.program,
				memory 	= self.max_memory_usage,
				prefix 	= self.abs_prefix
		)
		label = "Varscan Copynumber Ratios"
		output_result = self.runCallerCommand(copynumber_command, label, self.copynumber_output)

		return output_result

	def copyCaller(self):
		command = """java {memory} -jar {varscan} copyCaller {copynumber} \
		--output-file {called} \
		--output-homdel-file {homdels} \ 
		--min-coverage {mq}""".format(
			memory 		= self.max_memory_usage,
			varscan 	= self.program,
			copynumber 	= self.copynumber_output,
			called 		= self.called_copynumbers,
			homdels 	= self.called_homdels,
			mq 			= self.min_coverage)
		label = "Varscan CopyCaller"
		output_result = self.runCallerCommand(command, label, self.called_copynumbers)

		return output_result

class CNVkit(Caller):

	def setCustomEnvironment(self, sample, options):
		self.caller_name = "CNVKit"
		self.output_folder = getPipelineFolder('variants-copynumber', sample['PatientID'], self.caller_name)
		self.final_output = os.path.join(self.output_folder, 'reference_cnv.cnn')
		self.full_output = [
			sample['SampleID'] + '.cns',
			sample['SampleID'] + '.cnr',
			sample['SampleID'] + '.targetcoverage.cnn'
		]
		self.full_output = [os.path.join(self.output_folder, fn) for fn in self.full_output]
	
	def runCallerWorkflow(self, sample):
		call_status = self.runBatchCommand(sample)


	def runBatchCommand(self, sample):
		reference_cnn = os.path.join(self.output_folder, "reference_cnv.cnn")

		batch_command = """{cnvkit} batch {tumor} \
			--normal {normal} \
			--targets {targets} \
			--fasta {reference} \
			--output-reference {refcnn} \
			--output-dir {results} \
			--diagram \
			--scatter""".format(
				cnvkit = self.program,
				tumor = sample['TumorBAM'],
				normal = sample['NormalBAM'],
				targets = sample['ExomeTargets'],
				reference = self.reference,
				refcnn = reference_cnn,
				results = self.output_folder
			)

		label = "Batch Command"
		output_result = self.runCallerCommand(batch_command, label, reference_cnn)

		return output_result


class FREEC(Caller):
	def setCustomEnvironment(self, sample, options):
		self.caller_name = "FREEC"
		self.output_folder = getPipelineFolder('variants-copynumber', sample['PatientID'], self.caller_name)
		
		self.base_prefix = "{normal}_vs_{tumor}.{prefix}".format(
			tumor  = sample['SampleID'],
			normal = sample['NormalID'],
			prefix = self.caller_name.lower()
		)
		
		self.abs_prefix = os.path.join(
			self.output_folder,
			self.base_prefix
		)

		self.samtools_program = options['Programs']['samtools']

		self.program_folder = self.program # for readability
		self.program        = os.path.join(self.program_folder, 'src', 'freec')
		self.script_folder  = os.path.join(self.program_folder, 'scripts')

		self.assess_significance_script = os.path.join(self.script_folder, "assess_significance.R")
		self.plot_script                = os.path.join(self.script_folder, "makeGraph.R")
		self.config_file                = self.abs_prefix + "_config.txt"
		self.chrlenfile                 = os.path.join(self.output_folder, "chromosome_lengths.txt")

		pileup_normal_output_filename = os.path.join(
			self.temp_folder,
			'{}.{}.mpileup'.format(
				sample['NormalID'],
				self.caller_name.lower()
			)
		)

		pileup_sample_output_filename = os.path.join(
			self.temp_folder,
			'{}.{}.mpileup'.format(
				sample['NormalID'],
				self.caller_name.lower()
			)
		)
		nbasename = os.path.basename(pileup_normal_output_filename)
		tbasename = os.path.basename(pileup_sample_output_filename)
		# files for the normal sample
		self.copynumber_normal 	= os.path.join(self.output_folder, tbasename + "_normal_CNVs")
		self.ratio_normal  		= os.path.join(self.output_folder, tbasename + "_normal_ratio.txt")
		self.baf_normal 		= os.path.join(self.output_folder, nbasename + "_BAF.txt")
		self.plot_normal 		= os.path.join(self.output_folder, self.ratio_normal + ".png")
		# files for the tumor sample
		self.copynumber_sample  = os.path.join(self.output_folder, tbasename + "_CNVs")
		self.ratio_sample      	= os.path.join(self.output_folder, tbasename + "_ratio.txt")
		self.baf_sample    		= os.path.join(self.output_folder, tbasename + "_BAF.txt")
		self.plot_sample 		= os.path.join(self.output_folder, self.ratio_sample + ".png")
		
		self.full_output = [
			self.copynumber_normal, self.ratio_normal,  self.baf_normal,
			self.copynumber_sample,	self.ratio_sample,  self.baf_sample
		]

	def runCallerWorkflow(self, sample):

		# --------------------- Generate Pileup Files -----------------------------
		# Generate Pileup Files
		pileup_normal   = self.generatePileup(sample['NormalBAM'], sample['NormalID'])
		pileup_sample   = self.generatePileup(sample['TumorBAM'],  sample['SampleID'])
		
		# --------------------- Generate ChrLen File  -----------------------------
		# Detect Chromosome Lengths
		self.createChrLenFile(sample['ExomeTargets'])

		# --------------------- Generate Config File ------------------------------
		# Configure FREEC

		self.configureFREEC(pileup_normal, pileup_sample, sample['ExomeTargets'])

		# ------------------------------------- Run FREEC --------------------------------
		# Main Analysis
		# Filenames are automatically generated by the program based on the pileup filename.

		self.runFREEC()

		# ------------------------------------ Add Log2Ratios ----------------------------
		# --------------------------------- Calculate Significance -----------------------
		
		# Calculate Significance
		# files for the normal sample

		significance_normal = self.calculateSignificance(self.copynumber_normal, self.ratio_normal)
		significance_tumor  = self.calculateSignificance(self.copynumber_sample, self.ratio_sample)
	
		# Generate Plots

		plot_normal_status = self.generatePlots(self.ratio_normal, self.baf_normal)
		plot_sample_status = self.generatePlots(self.ratio_sample, self.baf_sample)
	
	def createChrLenFile(self, targets_file):
		if os.path.exists(self.chrlenfile): return None
		genome_index_file = self.reference + '.fai'
		exclude_chroms = [ #manually selected
			'chr1_KI270706v1_random',
			'chr4_GL000008v2_random',
			'chr14_GL000009v2_random',
			'chrUn_KI270742v1'
		]
		# generated as the intersection of the genome index and exome targets files
		# existence of targets file is verified at the start of the pipeline
		with open(targets_file, 'r') as bedfile:
			bedchrs = set([i[0] for i in csv.reader(bedfile, delimiter = '\t')])

		with open(genome_index_file, 'r') as file1:
			chrlens = ['\t'.join(i[:3]) + '\n'
					   for i in csv.reader(file1, delimiter = '\t')
					   if (i[0] in bedchrs and i[0] not in exclude_chroms)]
		
		with open(self.chrlenfile, 'w') as file1:
			[file1.write(i) for i in chrlens]

	def calculateSignificance(self, cnvs, ratios):
		expected_output = cnvs + '.p.value.txt'

		#Command used in the documentation
		command = """cat {script} | R --slave --args {cnvs} {ratios}""".format(
			script 	= self.assess_significance_script,
			cnvs 	= cnvs,
			ratios 	= ratios
		)

		#Reformmatted command
		command = """Rscript {script} --slave --args {cnvs} {ratios}""".format(
			script = self.assess_significance_script,
			cnvs = cnvs,
			ratios = ratios
		)

		label = "Calculate Significance"
		output_result = self.runCallerCommand(command, label, expected_output)
		return output_result

	def generatePlots(self, ratios, baf_file):
		# --------------------------------- Generate Plots --------------------------------
		expected_output = ratios + '.png'

		#Command used in the documentation
		command = """cat {script} | R --slave --args {ploidy} {ratios} {baf}""".format(
			script = self.plot_script,
			ploidy = "2",
			ratios = ratios,
			baf = baf_file)

		#reformatted command
		command = """Rscript {script} --slave --args {ploidy} {ratios} {baf}""".format(
			script = self.plot_script,
			ploidy = "2",
			ratios = ratios,
			baf = baf_file
		)

		label = "Generate Plots"
		output_result = self.runCallerCommand(command, label, expected_output)
		return output_result

	def configureFREEC(self, normal_pileup, tumor_pileup, targets):
		general_options = {
			'bedtools': "/usr/bin/bedtools",
			# Default: 0.8 use something like 0.6 to get more segments (and thus more predicted CNVs)
			'breakPointThreshold': 0.8,
			'chrFiles': os.path.join(os.path.dirname(self.reference), 'chromosomes'),
			# a list of chromosomes and chromosome lengths. Basically the reference dict.
			'chrLenFile': self.chrlenfile,
			'maxThreads': '6',
			'noisyData': 'TRUE',
			'outputDir': self.output_folder,
			'ploidy': '2',
			# set FALSE to avoid printing "-1" to the _ratio.txt files Useful for exome-seq or targeted sequencing data
			'printNA': 'FALSE',
			# Default: 10, recommended value >= 50 for for exome data
			'readCountThreshold': '50',
			'samtools': self.samtools_program,
			# 'sex': "",
			# for whole exome sequencing: "window=0"
			'window': "0"
		}
		general_options = ["{0} = {1}".format(k, v) for k, v in general_options.items()]
		
		# [sample]
		sample_options = {
			'mateFile': tumor_pileup,
			'inputFormat': 'pileup',
			'mateOrientation': '0'
		}
		sample_options = ["{0} = {1}".format(k, v) for k, v in sample_options.items()]

		# [control]
		control_options = {
			'mateFile': normal_pileup,
			'inputFormat': 'pileup',
			'mateOrientation': '0'
		}
		
		control_options = ["{0} = {1}".format(k, v) for k, v in control_options.items()]

		# [BAF]
		baf_options = {
			'SNPfile': self.dbSNP,
			"minimalCoveragePerPosition": '5'
		}
		baf_options = ["{0} = {1}".format(k, v) for k, v in baf_options.items()]

		# [target]
		target_options = {
			"captureRegions": targets,
		}
		target_options = ["{0} = {1}".format(k, v) for k, v in target_options.items()]

		all_options  = ["[general]"]+ general_options
		all_options += ["[sample]"] + sample_options
		all_options += ["[control"] + control_options
		all_options += ["[BAF]"]    + baf_options
		all_options += ["[target]"] + target_options

		with open(self.config_file, 'w') as file1:
			for line in all_options:
				file1.write(line + '\n')
		return os.path.exists(self.config_file)

	def runFREEC(self):
		# expected_output = [self.normal_CNVs, self.sample_CNVs]
		freec_command = "{freec} -conf {config}".format(
			freec = self.program,
			config = self.config_file
		)

		label = "Call Variants"
		output_result = self.runCallerCommand(freec_command, label, self.full_output)
		return output_result

# ----------------------------------------------------------------------------------------------------
# ------------------------------------------ Pipelines -----------------------------------------------
# ----------------------------------------------------------------------------------------------------

class BasePipeline:

	def __init__(self, sample, options_filename, callers):

		pipeline_options = self._verifyPipelineFiles(options_filename)
		self._verifySampleFiles(sample)

		self.runWorkflow(sample, pipeline_options, callers)

	def _verifyPipelineFiles(self, options_filename):
		""" Verifies that the files required to run the pipeline exist """

		# verify that the options file exists and load it.
		self._checkIfPathExists('options file', options_filename)
		options = configparser.ConfigParser()
		options.read(options_filename)

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
	def _runWorkflow(self, sample, options, workflow_callers):

		processed_bam = BaseQualityScoreRecalibration(sample, options)

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


if __name__ == "__main__":

	command_line_options = cmd_parser.getCMDArgumentParser().parse_args()

	default_config_filename = os.path.join(PIPELINE_DIRECTORY, "0_config_files", "pipeline_project_options.txt")
	
	default_sample_filename = os.path.join(PIPELINE_DIRECTORY, "debug_sample_list.tsv")
	pipeline_dna_callers = ['all']
	pipeline_copynumber_callers = []
	pipeline_rna_callers = []


	pipeline = GenomicsPipeline(
		default_sample_filename,
		default_config_filename,
		dna_callers = pipeline_dna_callers,
		rna_callers = pipeline_rna_callers,
		copynumber_callers = pipeline_copynumber_callers
	)

