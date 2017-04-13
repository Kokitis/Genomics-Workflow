import csv
import logging
import datetime
import os
import hashlib
import shutil
import isodate
#import gdc_api
import gdc_api
import subprocess
import math
import configparser
import hashlib
import shlex
from pprint import pprint
from argparse import ArgumentParser

#--------------------------------- Global Variables ----------------------------
now = datetime.datetime.now
GLOBAL_START = now() #Used to log when a series of samples were run together

#The parent folder the pipeline will be run in.
PIPELINE_DIRECTORY = "/home/upmc/Documents/Variant_Discovery_Pipeline"
#File to save the console output to. Only used when the console output is supressed.
CONSOLE_LOG_FILE = ""
#File containing a test sample.
SAMPLE_LOG_FILE = os.path.join(PIPELINE_DIRECTORY, "0_config_files", "sample_logV2.tsv")
README_FILE = os.path.join(PIPELINE_DIRECTORY, "0_readme_files", "readme.{0}.txt".format(now().isoformat()))
#Whether to use backwards-compatible filenames
BACKWARDS_COMPATIBLE = True
#initial_working_directory = PIPELINE_DIRECTORY

"""Set up the LOGGER"""
def configurePipelineLogger():

	logger_filename = os.path.join(PIPELINE_DIRECTORY, '0_config_files', 'pipeline_log.log')
	LOGGER = logging.getLogger('genome_pipeline')
	hdlr = logging.FileHandler(logger_filename)
	formatter = logging.Formatter('%(asctime)s %(levelname)s %(message)s')
	hdlr.setFormatter(formatter)
	LOGGER.addHandler(hdlr) 
	LOGGER.setLevel(logging.INFO)

	LOGGER.info("#" * 120)

	LOGGER.info('#'*30 + 'Starting the genomics pipeline at ' + GLOBAL_START.isoformat() + '#'*30)

	LOGGER.info("#" * 120)
	return LOGGER

#----------------------------------------------------------------------------------------------------
#------------------------------------ Set Up Global Functions ---------------------------------------
#----------------------------------------------------------------------------------------------------
#configuration = configparser.ConfigParser()
#config.read()
def generate_file_md5(filename, blocksize=2**20):
    m = hashlib.md5()
    with open( filename , "rb" ) as f:
        while True:
            buf = f.read(blocksize)
            if not buf:
                break
            m.update( buf )
    return m.hexdigest()

def checkdir(path, full = False):
	if not os.path.exists(path): 
		if full:
			os.makedirs(path)
		else:
			os.mkdir(path)

def format_options(options):
	""" Transforms a dictionary of options into a string.
		Use '|' to separate values with the same option key
	"""
	options = options.copy()
	option_string = list()
	for param, value in sorted(options.items()):
		if param in {'output', 'program', 'notes'}: continue
		if value in {None}: continue
		if value == "":
			string = "-{0}".format(param)
		elif isinstance(value, str) and '|' in value:
			value = value.split('|')
			string = ["-{0} {1}".format(param, i) for i in value]
			string = ' '.join(string)
		else:
			string = "-{0} {1}".format(param,value)
		option_string.append(string)
	option_string = " ".join(option_string)
	return option_string

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

def Terminal(command, label = None, show_output = True):
	""" Calls the system shell """
	terminal_log = os.path.join(PIPELINE_DIRECTORY, "0_config_files", "terminal_log.log")
	label = "{0} {1}".format(label if label is not None else "", now().isoformat())
	terminal_label  = '--'*15 + label + '--'*15 + '\n'
	terminal_label += '..'*40 + '\n'
	terminal_label += command + '\n'
	terminal_label += '..'*40 + '\n'
	terminal_label += '--'*40 + '\n'

	#Try using exceptions to catch timeout errors
	logging.info("System Command: " + str(command))
	if show_output:
		print(terminal_label)
		process = os.system(command)
		output = ""
	else:
		command = shlex.split(command)
		process = subprocess.Popen(command, stdout = subprocess.PIPE, stderr = subprocess.STDOUT)
		output = str(process.stdout.read(),'utf-8')
		with open(CONSOLE_LOG_FILE, 'a') as console_file:
			console_file.write(now().isoformat() + '\n')
			console_file.write(command + '\n')
			console_file.write(output + '\n\n')

	return output

def readTSV(filename, headers = False):
	with open(filename, 'r') as file1:
		reader = csv.DictReader(file1, delimiter = '\t')
		fieldnames = reader.fieldnames
		reader = list(reader)

	if headers: return reader, fieldnames
	else:       return reader

#----------------------------------------------------------------------------------------------------
#------------------------------------- Variant Caller Pipeline --------------------------------------
#----------------------------------------------------------------------------------------------------
class Caller:
	def __init__(self, sample, options, DEBUG = False):
		#self.caller_name = caller_name
		program_start = now()
		print("{0}(DEBUG = {1})".format(self.__name__, DEBUG))
		self.setCallerEnvironment(sample, options)

		if DEBUG:
			print("Running ", self.caller_name)
			print("\tprogram location: ", self.program)
			print("\toutput folder: ", self.output_folder)
			print("\ttemp folder: ", self.temp_folder)
			print("\tcaller output prefix: ", self.prefix)
		self.createReadMe(sample)
		self.runCallerWorkflow(sample, options)
		if BACKWARDS_COMPATIBLE:
			#Rename the output files to make them backwards compatible 
			#with previous versions of the pipeline.
			self.runBackwardsCompatibility()
		
		program_stop = now()
		self.updateSampleLog(sample, program_start, program_stop)
	
	def setCallerEnvironment(self, sample, options):
		self.caller_name = self.__name__
		self.reference 	= options['Reference Files']['reference genome']
		self.dbSNP 		= options['Reference Files']['dbSNP']
		self.cosmic 	= options['Reference Files']['COSMIC']


		self.program 			= options['Programs'].get(self.caller_name)
		self.gatk_program 		= options['Programs']['GATK']
		self.max_cpu_threads 	= options['Parameters']['MAX_CORES']
		self.max_memory_usage 	= options['Parameters']['JAVA_MAX_MEMORY_USAGE']
		self.min_base_quality 	= options['Parameters']['MIN_NUCLEOTIDE_QUALITY']
		self.min_mapping_quality = options['Parameters']['MIN_MAPPING_QUALITY']
		self.min_somatic_quality = options['Parameters']['SOMATIC_QUALITY']
		self.min_coverage 		 = options['Parameters']['MIN_COVERAGE']

		self.output_folder = os.path.join(
			options['Pipeline Options']['somatic pipeline folder'],
			sample['PatientID'], self.caller_name)
		
		self.prefix = os.path.join(self.output_folder, 
			"{normal}_vs_{tumor}.{prefix}".format(
			tumor   = sample['SampleID'], 
			normal  = sample['NormalID'],
			prefix = self.caller_name.lower()))
		
		self.temp_folder = os.path.join(
			options['Pipeline Options']['temporary folder'], 
			sample['PatientID'])

		self.command_list = list()
		self.temp_files = list()
		
		checkdir(self.output_folder)
		checkdir(self.temp_folder, True)

	def runCallerCommand(self, command, label, expected_output):
		self.command_list.append(command)
		if isinstance(expected_output, str): expected_output = [expected_output]

		files_missing = any([not os.path.exists(fn) for fn in expected_output])

		if len(expected_output) == 0 or files_missing:
			self.addToReadMe(command, label, expected_output)
			Terminal(command)
		else:
			basenames = [os.path.basename(fn) for fn in expected_output]
			print("The following output files already exist:")
			for basename in basenames:
				print("\t", basename)
		command_status = len(expected_output) == 0 or not any([not os.path.exists(fn) for fn in expected_output])
		return command_status

	def generatePileup(self, bam_file, bam_name):
		""" Generates pileup files. If single is True, only
			one file will be generated.
		"""
		output_file = os.path.join(self.temp_folder, bam_name + '.mpileup')
		pileup_command = "samtools mpileup -q 1 -B -f {reference} {sample} > {output}".format(
			reference = self.reference,
			sample = bam_file,
			output = output_file)

		self.runCallerCommand(pileup_command, output_file)
		return output_file

	def runBackwardsCompatibility(self):
		""" Renames the output files to be compatible with previous versions of the pipeline. """
		pass

	def createReadMe(self, sample):
		current_datetime =  now().isoformat()

		readme_filename = "{0}.{1}.readme.txt".format(self.caller_name, current_datetime)
		readme_filename = os.path.join(self.output_folder, readme_filename)

		with open(readme_filename, 'w') as readme_file:
			readme_file.write(self.caller_name + '\n')
			readme_file.write("Started the caller at {0}\n".format(current_datetime))
			for key, value in sample.items():
				readme_file.write("{0}:\t{1}\n".format(key, value))

	def addToReadMe(self, command, label, expected_output):
		current_datetime = now()
		line_len = 100
		num = line_len - (len(label) + len(current_datetime.isoformat()))
		linebreak = '#' * int(num / 2) + "{0} {1}".format(current_datetime.isoformat(), label) + "#" * int(num/2) + '\n'

		with open(self.readme_file, 'a') as readme:
			readme.write(linebreak)
			readme.write(command + '\n')
			readme.write("Expected Output: ", expected_output)

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
			'outputFiles':    ', '.join(self.full_output),#'|'.join(outputs),
			'notes': 	  "",
			'status':	  status,
			'commands':   '|'.join(self.command_list)#'|'.join(commands)
		}
		writeheaders = not os.path.exists(SAMPLE_LOG_FILE) or os.path.getsize(SAMPLE_LOG_FILE) == 0
		#line 	= sorted(caller_log.items())
		#line 	= list(reversed(line))
		#headers = '\t'.join([i[0] for i in line]) + '\n'
		#line 	= '\t'.join([str(i[1]) for i in line]) + '\n'

		#with open(SAMPLE_LOG_FILE, 'a') as file1:
		#	if writeheaders: file1.write(headers)
		#	file1.write(line)
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

	def getCallerStatus(self):
		caller_failed = any(not os.path.exists(fn) for fn in self.full_output)
		status = not caller_failed

		return status
	@staticmethod
	def _fileNotFound(filename):
		print("Could not locate ", filename)


class MuSE(Caller):
	__name__ = "MuSE"
	def runCallerWorkflow(self, sample, options):
		#self.caller_name = 'MuSE'
		#super().__init__(self, sample, options, caller_name)

		self.call_output = self.runMuseCall(sample)
		self.sump_output = os.path.splitext(self.call_output)[0] + '.vcf'
		sump_output = self.runMuseSump(self.call_output)

		self.full_output = [self.call_output, self.sump_output]
		self.final_output= self.sump_output

	def runMuseCall(self, sample):
		prefix = '.'.join(self.prefix.split('.')[:-1])
		muse_call_output = prefix + '.MuSE.txt'
		call_command = "{program} call -O {prefix} -f {reference} {tumor} {normal}".format(
			program = self.program,
			reference = self.reference,
			prefix = prefix,
			tumor = sample['TumorBAM'],
			normal = sample['NormalBAM'])
		
		self.runCallerCommand(call_command, muse_call_output)

		return muse_call_output

	def runMuseSump(self, call_output):
		
		#muse_sump_output = output_prefix + ".MuSe.vcf"
		sump_command = "{program} sump -I {call_output} -E -D {dbSNP} -O {output}".format(
			program = self.program,
			prefix = self.prefix,
			dbSNP = self.dbSNP,
			call_output = call_output,
			output = self.sump_output)

		status = self.runCallerCommand(sump_command, self.sump_output)

		return status

class MuTect2(Caller):
	__name__ = "MuTect2"
	def runCallerWorkflow(self, sample, options):
		#super().__init__(self, sample, options, 'MuTect2')
		self.final_output = self.prefix + '.vcf'
		mutect2_output = self.runMutect2(sample)
		self.full_output = [self.final_output]

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
			cosmic 		= self.cosmic,
			reference   = self.reference, 
			normal      = sample['NormalBAM'], 
			tumor       = sample['TumorBAM'], 
			targets     = sample['ExomeTargets'],
			output      = self.final_output)

		status = self.runCallerCommand(mutect2_command, self.final_output)

		return status

class SomaticSniper(Caller):
	__name__ = "SomaticSniper"
	def runCallerWorkflow(self, sample, options):
		#Will print "Couldn't find single-end mapping quality. Check to see if the SM tag is in BAM."
		#This doesn't invalidate results, but try not to use single-end mapping quality in output
		#super().__init__(self, sample, options, 'SomaticSniper')
		
		self.readcount_program = options['Programs']['bam-readcount']
		self.samtools_program = options['Programs']['samtools-0.1.6']
		self.somaticsniper_folder = self.program
		self.program = os.path.join(self.somaticsniper_folder, 'build', 'bin', 'bam-somaticsniper')
		self.script_folder = os.path.join(self.somaticsniper_folder, 'src', 'scripts')

		self.snpfilter_script= os.path.join(self.script_folder, 'snpfilter.pl')
		self.readcount_script =       os.path.join(self.script_folder, 'prepare_for_readcount.pl')
		self.hc_script =       os.path.join(self.script_folder, 'highconfidence.pl')
		self.fpfilter  =       os.path.join(self.script_folder, 'fpfilter.pl')

		self.raw_variants = self.prefix + '.vcf'
		self.hq_variants = self.prefix + '.hq.vcf'
		self.lq_variants = self.prefix + '.lq.vcf'


		raw_output = self.runVariantDiscovery(sample)
		#Generate pileup files
		normal_pileup_file = self.generatePileupFile(sample['NormalBAM'], 'normal')
		tumor_pileup_file  = self.generatePileupFile(sample['TumorBAM'], 'tumor')

		#Filter LOH
		_intermediate_loh_filtered_output = self.prefix + ".SNPfilter.intermediate"
		loh_filtered_output  = self.prefix + ".SNPfilter.final"
		_intermediate_file  = self.removeLOH(self.raw_variants, normal_pileup_file, _intermediate_loh_filtered_output)
		loh_filtered_output = self.removeLOH(_intermediate_loh_filtered_output, tumor_pileup_file, loh_filtered_output)

		
		readcounts = self.readcounts(loh_filtered_output, sample['TumorBAM'])
		false_positive_output = self.removeFalsePositives(loh_filtered_output, readcounts)
		high_confidence_variants = self.calculateConfidence(false_positive_output)

		self.full_output = [self.raw_variants, self.hq_variants, self.lq_variants]
		self.final_output = self.hq_variants
	
	def _setDefultOptions(self):
		default_options = {
			'q': self.min_mapping_quality, #filtering reads with mapping quality less than INT [0]
			'Q': self.min_somatic_quality, #filtering somatic snv output with somatic quality less than INT [15]
			'L': None,# FLAG do not report LOH variants as determined by genotypes

			'G': None,	#do not report Gain of Referene variants as determined by genotypes

			'p': None, 	#disable priors in the somatic calculation. Increases sensitivity for solid tumors.

			'J': None, 	#Use prior probabilities accounting for the somatic mutation rate

			's': 0.01, 	#FLOAT prior probability of a somatic mutation (implies -J) [0.01]

			'T': 0.850000, #FLOAT theta in maq consensus calling model (for -c/-g) [0.850000]

			'N': 2, 		#INT number of haplotypes in the sample (for -c/-g) [2]

			'r': 0.001000, #FLOAT prior of a diﬀerence between two haplotypes (for -c/-g) [0.001000]
			'F': 'vcf' 	#STRING select output format (vcf or classic) [classic]
		}
		
		default_options = format_options(default_options)
		return default_options
	def runVariantDiscovery(self, sample):
		#-------------------------------- Variant Discovery Command -----------------------------
		#output_file = self.prefix + '.vcf'
		#default_options = self._setDefultOptions()
		somaticsniper_command = "{program} -q {mmq} -Q {mss} -F vcf -f {reference} {tumor} {normal} {outputfile}".format(
			program 	= self.program,
			reference 	= self.reference,
			mmq 		= self.min_mapping_quality,
			mss 		= self.min_somatic_quality,
			tumor 		= sample['TumorBAM'],
			normal 		= sample['NormalBAM'],
			outputfile 	= self.raw_variants)
		label = "SomaticSniper Variant Discovery"
		status = self.runCallerCommand(somaticsniper_command, label, self.raw_variants)
		return status

	def generatePileupFile(self, bam_file, bam_name):
		pileup_file = os.path.join(self.temp_folder, "somaticsniper.{0}_indel".format(bam_name))
		output_file = pileup_file + '.pileup'
		samtools_command = "{samtools} pileup -cvi -f {reference} {bam} > {pileup}".format(
			samtools = self.samtools_program,
			reference = self.reference,
			bam = bam_file,
			pileup = pileup_file)
		#base_filter_command = """samtools.pl varFilter {inputpileup} | awk '$6>={basequality}' | grep -P "\\t\\*\\t" > {outputpileup}.pileup"""
		filter_command = """samtools.pl varFilter -Q {basequality} {inputpileup} > {outputpileup}""".format(
			basequality = self.min_base_quality,
			inputpileup = pileup_file,
			outputpileup = output_file)
		#samtools.pl varFilter raw.pileup | awk '$6>=20' > final.pileup
		label = "SomaticSniper: Generate Pileup File" 
		self.runCallerCommand(samtools_command, label, pileup_file)
		label = "SomaticSniper: Filter Pileup File"
		self.runCallerCommand(filter_command, label, output_file)

		return output_file

	def removeLOH(self, vcf_file, pileup_file, output_file):
		#-------------------------- Filter and remove LOH --------------------------------
		command = "perl {snpfilter} --snp-file {vcf} --indel-file {pileup} --out-file {output}".format(
			snpfilter = self.snpfilter_script,
			vcf = vcf_file,
			pileup = pileup_file,
			output = output_file)
		label = "SomaticSniper: removeLOH"
		self.runCallerCommand(command, label, output_file)
		return output_file

	def readcounts(self, loh_file, tumor_bam):
		""" Expected Output: 
		"""
		prepare_readcount_output = loh_file + '.pos'
		readcount_output = loh_file + '.readcounts.rc'
		#Prepare readcounts
		pr_command = "perl {script} --snp-file {inputfile} --out-file {output}".format(
			script = self.readcount_script,
			inputfile = loh_file,
			output = prepare_readcount_output)
		label = "SomaticSniper: Prepare Readcounts"
		self.runCallerCommand(pr_command, label, prepare_readcount_output)
		#Readcounts

		readcount_command = "{program} -b {mbq} -q 1 -f {reference} -l {proutput} {tumor} > {output}".format(
			program = self.readcount_program,
			mbq = self.min_base_quality,
			reference = self.reference,
			proutput = prepare_readcount_output,
			tumor = tumor_bam,
			output = readcount_output)
		label = "SomaticSniper: Generate Readcounts"
		self.runCallerCommand(readcount_command, label, readcount_output)

		return readcount_output

	def removeFalsePositives(self, loh_file, readcounts):
		false_positive_output = loh_file + '.fp_pass'
		fp_command = "perl {fpfilter} --snp-file {snpfilter} -readcount-file {readcounts}".format(
			fpfilter = self.fpfilter,
			snpfilter = loh_file,
			readcounts = readcounts)
		label = "SomaticSniper: Remove False Positives"
		self.runCallerCommand(fp_command, label, false_positive_output)
		return false_positive_output

	def calculateConfidence(self, vcf_file):
		command = "perl {script} --snp-file {vcf} --min-mapping-quality {mmq} --min-somatic-score {ss} --lq-output {lq} --out-file {hq}".format(
			script = self.hc_script,
			mmq = self.min_mapping_quality,
			ss = self.min_somatic_quality,
			vcf = vcf_file,
			hq = self.hq_variants,
			lq = self.lq_variants)
		label = "SomaticSniper: Filter lq Variants"
		status = self.runCallerCommand(command, label, [self.hq_variants, self.lq_variants])
		return status

	def runBackwardsCompatibility(self):

		for fn in os.listdir(folder):
			abs_fn = os.path.join(folder, fn)
			if "" in fn:
				pass
			elif "" in fn:
				pass
	def writeReadMe(self):
		pass

class Strelka(Caller):
	__name__ = "Strelka"
	def runCallerWorkflow(self, sample, options):
		#super().__init__(self, sample, options, 'Strelka')
		
		patient_folder = os.path.dirname(self.output_folder)

		self.results_folder = os.path.join(patient_folder, 'Strelka', 'results')
		self.config_script = os.path.join(self.program, 'bin', 'configureStrelkaWorkflow.pl')
		self.strelka_project_config_file = os.path.join(patient_folder, "strelka_configuration.ini")
		self.final_output = os.path.join(self.results_folder, 'passed.somatic.snvs.vcf')
		#Strelka won't run if the output folder already exists.
		if not os.path.exists(self.final_output) and os.path.exists(self.output_folder):
			shutil.rmtree(self.output_folder)


		strelka_folder = os.path.join(patient_folder, 'Strelka')
		
		#--------------------------------- Configure Strelka --------------------------------
		#change working directory
		config_script = self.generateStrelkaConfigFile(self.strelka_project_config_file)
		self.configureStrelka(sample, self.strelka_project_config_file)
		output_files = self.runStrelka()

		self.full_output = output_files
		
	def generateStrelkaConfigFile(self, configuration_file):
		strelka_configuration_options = self.getStrelkaConfiguration()
		
		with open(configuration_file, 'w') as file1:
			#strelka_configuration_options = [i.split('\t')[0].strip() for i in strelka_configuration_options.s]
			file1.write('[user]\n')
			for row in strelka_configuration_options.split('\n'):
				row = [i for i in row.split('\t') if '=' in i]
				if len(row) > 0:
					file1.write(row[0] + '\n')
		return configuration_file


	def configureStrelka(self, sample, config_ini):
		strelka_config_command = "{script} --tumor={tumor} --normal={normal} --ref={reference} --config={config} --output-dir={output_dir}".format(
			script = self.config_script,
			tumor = sample['TumorBAM'],
			normal= sample['NormalBAM'],
			reference = self.reference,
			config = config_ini,
			output_dir = self.output_folder)
		self.runCallerCommand(strelka_config_command, self.results_folder)

	def runStrelka(self):
		output_suffixes = ['all.somatic.indels.vcf', 'all.somatic.snvs.vcf', 'passed.somatic.indels.vcf', 'passed.somatic.snvs.vcf']
		output_files = [os.path.join(self.results_folder, i) for i in output_suffixes]
		strelka_run_command = "make -j {threads} -C {outputdir}".format(
			outputdir = self.output_folder, threads = self.max_cpu_threads)
		self.runCallerCommand(strelka_run_command, output_files)
		return output_files

	def getStrelkaConfiguration(self):
		strelka_configuration_options = """
			isSkipDepthFilters = 1 						isSkipDepthFilters should be set to 1 to skip depth filtration for whole exome or other targeted sequencing data
			maxInputDepth = 10000						strelka will not accept input reads above this depth. Set this value <= 0 to disable this feature.
			depthFilterMultiple = 3.0 					If the depth filter is not skipped, all variants which occur at a depth greater than depthFilterMultiple*chromosome mean depth will be filtered out.
			snvMaxFilteredBasecallFrac = 0.4 			Somatic SNV calls are filtered at sites where greater than this fraction of basecalls have been removed by the mismatch density filter in either sample.
			snvMaxSpanningDeletionFrac = 0.75 			Somatic SNV calls are filtered at sites where greater than this fraction of overlapping reads contain deletions which span the SNV call site.
			indelMaxRefRepeat = 8 						Somatic indel calls are filtered if they represent an expansion or contraction of a repeated pattern with a repeat count greater than indelMaxRefRepeat in the reference
			indelMaxWindowFilteredBasecallFrac = 0.3	Somatic indel calls are filtered if greater than this fraction of basecalls in a window extending 50 bases to each side of an indel's call position have been removed by the mismatch density filter.
			indelMaxIntHpolLength = 14 					Somatic indels are filtered if they overlap ’interrupted homopolymers’ greater than this length. 
			ssnvPrior = 0.000001 						prior probability of a somatic snv or indel
			sindelPrior = 0.000001
			ssnvNoise = 0.0000005 						probability of an snv or indel noise allele 
			sindelNoise = 0.000001
			ssnvNoiseStrandBiasFrac = 0.5 				Fraction of snv noise attributed to strand-bias. It is not recommended to change this setting.
			minTier1Mapq = 20 							minimum MAPQ score for PE reads at tier1
			minTier2Mapq = 0 							minimum MAPQ score for PE and SE reads at tier2
			ssnvQuality_LowerBound = 15 				Somatic quality score (QSS_NT, NT=ref) below which somatic SNVs are marked as filtered
			sindelQuality_LowerBound = 30				Somatic quality score (QSI_NT, NT=ref) below which somatic indels are marked as filtered
			isWriteRealignedBam = 0 					Optionally write out read alignments which were altered during the realignment step. Location: ${ANALYSIS_DIR}/realigned/{normal,tumor}.realigned.bam
			binSize = 25000000 							Jobs are parallelized over segments of the reference genome no larger than this size"""
		return strelka_configuration_options

class Varscan(Caller):
	__name__ = "Varscan"
	def runCallerWorkflow(self, sample, options):
		#super().__init__(self, sample, options, 'Varscan')

		self.raw_snps = self.prefix + '.raw.snp.vcf'
		self.raw_indels=self.prefix + '.raw.indel.vcf'

		self.somatic_hc = self.prefix + '.snp.Somatic.hc'
		self.somatic_lc = self.prefix + '.snp.Somatic.lc'
		self.germline 	= self.prefix + '.snp.Germline'
		self.germline_hc= self.prefix + '.snp.Germline.hc'
		self.loh 		= self.prefix + '.snp.LOH'
		self.loh_hc 	= self.prefix + '.snp.LOH.hc'

		normal_pileup 	= self.generatePileup(sample['NormalBAM'], sample['NormalID'])
		tumor_pileup 	= self.generatePileup(sample['TumorBAM'], sample['SampleID'])

		self.runVarscan(normal_pileup, tumor_pileup)

		processed_variants = self.postProcessing()
		
		self.full_output = [
			self.raw_snps, 
			self.raw_indels,
			self.somatic_hc,
			self.somatic_lc,
			self.germline,
			self.germline_hc,
			self.loh,
			self.loh_hc]
		self.final_output = self.somatic_hc

	def runVarscan(self, normal_pileup, tumor_pileup):
		command = "java {memory} -jar {varscan} somatic {normal} {tumor} --output-snp {snp} --output-indel {indel} --output-vcf 1".format(
			varscan = self.program,
			memory 	= self.max_memory_usage,
			normal 	= normal_pileup,
			tumor 	= tumor_pileup,
			snp 	= self.raw_snps,
			indel 	= self.raw_indels,
			mc 		= self.min_coverage)
		output_files = [self.raw_snps, self.raw_indels]
		status = self.runCallerCommand(command, output_files)

		return status

	def postProcessing(self):
		expected_output = [self.somatic_hc, self.somatic_lc, self.germline, self.germline_hc, self.loh, self.loh_hc]

		process_command = "java {memory} -jar {varscan} processSomatic {vcf}"
		process_command = process_command.format(
			memory  = self.max_memory_usage,
			varscan = self.program,
			vcf  	= self.raw_snps)
		status = self.runCallerCommand(process_command, expected_output)
		return status

class HaplotypeCaller(Caller):
	__name__ = "HaplotypeCaller"
	def runCallerWorkflow(self, sample, options, workflow = 'DNA-seq'):
		#super().__init__(self, sample, options, "HaplotypeCaller")

		#self.output_folder = os.path.join(self.output_folder, workflow)
		#self.prefix = os.path.join(self.output_folder, "{0}_vs_{1}.haplotypecaller".format(sample['NormalID'], sample['SampleID']))
		checkdir(self.output_folder, True)

		self.dna_output = self.prefix + ".DNA.raw_snps_indels.vcf"
		self.rna_output = self.prefix + ".RNA.raw_snps_indels.vcf"
		self.filtered_variants = os.path.splitext(self.rna_output)[0] + ".filtered.vcf"

		if workflow == 'DNA-seq':
			output = self.callDNAVariants(sample, options)
			self.full_output = [self.dna_output]
			self.final_output = self.dna_output
		else:
			bqsr = BaseQualityScoreRecalibration(sample, options)
			raw_variants = self.callRNAVariants(bqsr.bam)
			filtered_variants = self.filterVariants(self.rna_output)
			self.full_output = [self.rna_output, self.filtered_variants]
			self.final_output = self.filtered_variants

	def callRNAVariants(self, bam_file):
		command = """java -jar {GATK} \
			-T HaplotypeCaller \
			-R {reference} \
			-I {sample} \
			--dbsnp {dbSNP} \
			--dontUseSoftClippedBases \
			-o {output}""".format(
				GATK 		= self.gatk_program,
				reference 	= self.reference,
				sample 		= bam_file,
				dbSNP 		= self.dbSNP,
				output 		= self.rna_output)
		status = self.runCallerCommand(command, self.rna_output)
		return status
	
	def filterVariants(self, vcf_file):
		command = """java -jar {GATK} \
			-T VariantFiltration \
			-R {reference} \
			-V {inputfile} \
			-window 35 \
			-cluster 3 \
			--filterName FS --filterExpression \"FS > 30.0\" \
			--filterName QD --filterExpression \"QD < 2.0\" \
			-o {output}""".format(
				GATK 		= self.gatk_program,
				reference 	= self.reference,
				inputfile 	= vcf_file,
				output 		= self.filtered_variants)
		status = self.runCallerCommand(command, self.filtered_variants)
		return status
	
	def callDNAVariants(self, sample, options):
		command = """java -jar {GATK} \
			-T HaplotypeCaller \
			-R {reference} \
			-I {normal} \
			-I {tumor} \
			-L {targets} \
			-nct {threads} \
			--dbsnp {dbSNP} \
			-o {output}""".format(
				GATK 		= self.gatk_program,
				reference 	= self.reference,
				normal 		= sample['NormalBAM'],
				tumor 		= sample['TumorBAM'],
				targets 	= sample['ExomeTargets'],
				dbSNP 		= self.dbSNP,
				output 		= self.dna_output,
				threads 	= self.max_cpu_threads)

		status = self.runCallerCommand(command, self.dna_output)
		return status


class BaseQualityScoreRecalibration(Caller):
	__name__ = "BaseQualityScoreRecalibration"
	def runCallerWorkflow(self, sample, options):
		#super().__init__(self, sample, options, 'BaseQualityScoreRecalibration')
		self.output_folder = os.path.join("/media/upmc/WD_Partition_2/RNA-seq", 'recalibrated_genomes', sample['PatientID'])

		self.recalibration_table 	= os.path.join(self.output_folder, sample['SampleID'] + ".RNA.recalibration_data.table")
		self.covariate_table 		= os.path.join(self.output_folder, sample['SampleID'] + ".RNA.covariate_data.table")
		self.recalibration_plots 	= os.path.join(self.output_folder, sample['SampleID'] + ".RNA.recalibration_plots.pdf")
		self.cigar_bam 				= os.path.join(self.output_folder, sample['SampleID'] + ".RNA.cigar.bam")
		self.realigned_bam 			= os.path.join(self.output_folder, sample['SampleID'] + ".RNA.recalibrated.bam")
		self.bam 					= self.realigned_bam
		
		raw_rna_bam_file 	= self._getRNABAM(sample)
		cigar_bam 			= self.splitCigarReads(raw_rna_bam_file)
		recalibration_table = self.generateRecalibrationTable(cigar_bam)
		realigned_bam 		= self.recalibrateBAM(cigar_bam)
		covariate_table 	= self.generateCovariateTable(realigned_bam)
		recalibration_plots = self.generateRecalibrationPlots()

		self.final_output 	= realigned_bam
		self.full_output 	= [self.recalibration_table, self.covariate_table, self.recalibration_plots, self.cigar_bam, self.realigned_bam]

	def _getRNABAM(self, sample):
		return sample['RNABAM']
	
	def splitCigarReads(self, bam_file):
		command = """java -jar {GATK} \
			-T SplitNCigarReads \
			-R {reference} \
			-I {inputbam} \
			-o {outputbam} \
			-rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 \
			-U ALLOW_N_CIGAR_READS""".format(
				GATK 		= self.gatk_program,
				reference 	= self.reference,
				inputbam 	= bam_file,
				outputbam 	= self.cigar_bam)
		status = self.runCallerCommand(command, self.cigar_bam)
		return status

	def generateRecalibrationTable(self, bam_file):
		command = """java -jar {GATK} \
			-T BaseRecalibrator \
			-R {reference} \
			-I {bam} \
			-knownSites {dbSNP} \
			-o {output}""".format(
				GATK 		= self.gatk_program,
				reference 	= self.reference,
				dbSNP 		= self.dbSNP,
				bam 		= bam_file,
				output 		= self.recalibration_table)
		status = self.runCallerCommand(command, self.recalibration_table)
		return status

	def recalibrateBAM(self, bam_file):
		command  = """java -jar {GATK} \
			-T PrintReads \
			-R {reference} \
			-I {bam} \
			-BQSR {table} \
			-o {output}""".format(
				GATK = self.gatk_program,
				reference = self.reference,
				bam = bam_file,
				output = self.realigned_bam,
				table = self.recalibration_table)
		status = self.runCallerCommand(command, self.realigned_bam)
		return status
	
	def generateCovariateTable(self, bam_file):
		command = """java -jar {GATK} \
			-T BaseRecalibrator \
			-R {reference} \
			-I {realigned_bam} \
			-BQSR {table} \
			-knownSites {dbSNP} \
			-o {output}""".format(
				GATK 		= self.gatk_program,
				reference 	= self.reference,
				dbSNP 		= self.dbSNP,
				table 		= self.recalibration_table,
				realigned_bam = bam_file,
				output 		= self.covariate_table)
		status = self.runCallerCommand(command, self.covariate_table)
		return status
	
	def generateRecalibrationPlots(self):
		command = """ java -jar {GATK} \
			-T AnalyzeCovariates \
			-R {reference} \
			-before {before} \
			-after {after} \
			-plots {plots}""".format(
				GATK 		= self.gatk_program,
				reference 	= self.reference,
				before 		= self.recalibration_table,
				after 		= self.covariate_table,
				plots 		= self.recalibration_plots)
		status = self.runCallerCommand(command, self.recalibration_plots)
		return status

class UnifiedGenotyper(Caller):
	__name__ = "UnifiedGenotyper"
	def runCallerWorkflow(self, sample, options, workflow = 'DNA-seq'):
		#super().__init__(self, sample, options, 'UnifiedGenotyper')
		#self.output_folder = os.path.join(self.output_folder, workflow)
		#self.prefix = os.path.join(self.output_folder, "{0}_vs_{1}.unifiedgenotyper".format(sample['NormalID'], sample['SampleID']))

		self.dna_variants = self.prefix + '.DNA.raw_snps_indels.vcf'
		self.rna_variants = self.prefix + '.RNA.raw_snps_indels.vcf'
		self.filtered_variants = os.path.splitext(self.rna_variants)[0] + '.filtered.vcf'
		
		if workflow == 'DNA-seq':
			dna_variants = self.callDNAVariants(sample)
			self.full_output = [self.dna_variants]
			self.final_output = self.dna_variants
		else:
			raw_rna_variants = self.callRNAVariants(sample)
			filtered_rna_variants = self.filterVariants(sample)
			self.full_output += [self.rna_variants, self.filtered_variants]
			self.final_output = self.filtered_variants

	def callDNAVariants(self, sample):
		command = """java -jar {GATK} \
			-T UnifiedGenotyper \
			-R {reference} \
			-I {normal} \
			-I {tumor} \
			-L {targets} \
			-nct 6 \
			--dbsnp {dbSNP} \
			-o {output}""".format(
				GATK = self.gatk_program,
				reference = self.reference,
				normal = sample['NormalBAM'],
				tumor = sample['TumorBAM'],
				targets = sample['ExomeTargets'],
				dbSNP = self.dbSNP,
				output = self.dna_variants)
		status = self.runCallerCommand(command, self.dna_variants)
		return status

	def callRNAVariants(self, bam_file):
		command = """java -jar {GATK} \
			-T UnifiedGenotyper \
			-R {reference} \
			-I {sample} \
			--dbsnp {dbSNP} \
			-nct 6 \
			--maxRuntime 30 \
			--filter_reads_with_N_cigar \
			-o {output}""".format(
				GATK = self.gatk_program,
				reference = self.reference,
				sample = bam_file,
				dbSNP = self.dbSNP,
				output = self.rna_output)
		status = self.runCallerCommand(command, self.rna_output)
		return status

	def filterVariants(self, vcf_file):
		command = """java ‐jar {GATK} \
			‐T VariantFiltration 
			‐R {reference} \
			‐V {inputfile} \
			‐window 35 \
			‐cluster 3 \
			‐filterName FS ‐filter "FS > 30.0" \
			‐filterName QD ‐filter "QD < 2.0" \
			‐o {output}""".format(
				GATK = self.gatk_program,
				reference = self.reference,
				inputfile = vcf_file,
				output = self.filtered_variants)
		status = self.runCallerCommand(command, self.filtered_variants)
		return status

#----------------------------------------------------------------------------------------------------
#------------------------------------ Copynumber Tool Functions -------------------------------------
#----------------------------------------------------------------------------------------------------
class VarscanCopynumber(Caller):
	__name__ = "Varscan"
	def runCallerWorkflow(self, sample, options):
		#super().__init__(self, sample, options, 'Varscan')
		self.output_folder = os.path.join(
			options['Pipeline Options']['copynumber pipeline folder'],
			sample['PatientID'], self.caller_name)
		
		
		self.prefix = os.path.join(self.output_folder, 
			"{normal}_vs_{tumor}.{prefix}".format(
			tumor   = sample['SampleID'], 
			normal  = sample['NormalID'],
			prefix = self.caller_name.lower()))
		
		self.rscript_filename 	= os.path.join(self.output_folder, "{0}.varscan_CBS.r".format(sample['PatientID']))
		self.copynumber_output 	= self.prefix + '.copynumber'			#[prefix].copynumber
		self.called_copynumbers = self.copynumber_output + '.called'	#[prefix].copynumber.called
		self.called_homdels 		= self.called_copynumbers + '.homdel'	#[prefix].copynumber.called.homdel
		self.copynumber_segments =self.called_copynumbers + '.segments' #[prefix].copynumber.called.segments

		normal_pileup 	= self.generatePileup(sample['NormalBAM'], sample['NormalID'])
		tumor_pileup 	= self.generatePileup(sample['TumorBAM'], sample['SampleID'])

		copynumber_ratios = self.callCopynumberRatios(normal_pileup, tumor_pileup)
		copynumber_calls  = self.copyCaller()

		copynumber_segments = self.circularBinarySegmentation()


		self.full_output = [
			self.copynumber_output,
			self.called_copynumbers,
			self.called_homdels,
			self.copynumber_segments]
		self.final_output = self.copynumber_segments

	def circularBinarySegmentation(self):
		#NEED to strip first row of table before r script
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

		status = self.runCallerCommand(command, self.copynumber_segments)
		return status

	def callCopynumberRatios(self, normal_pileup, tumor_pileup):
		#-------------------------- Varscan (copynumber - caller) -----------------------------
		copynumber_command = "java {memory} -jar {varscan} copynumber {normal} {tumor} {prefix} --min-base_qual {mbq} --min-map-qual {mmq}".format(
			mmq = self.min_mapping_quality,
			mbq = self.min_base_quality,
			normal = normal_pileup,
			tumor = tumor_pileup,
			varscan = self.program,
			memory = self.max_memory_usage,
			prefix = self.prefix)
		self.runCallerCommand(copynumber_command, self.copynumber_output)

		return self.copynumber_output

	def copyCaller(self):
		command = """java {memory} -jar {varscan} copyCaller {copynumber} \
		--output-file {called} \
		--output-homdel-file {homdels} \ 
		--min-coverage {mq}""".format(
			memory = self.max_memory_usage,
			varscan = self.program,
			copynumber = self.copynumber_output,
			called = self.called_copynumbers,
			homdels = self.called_homdels,
			mq = self.min_coverage)

		self.runCallerCommand(command, self.called_copynumbers)

		return self.called_copynumbers

def varscan_copynumber(sample, options):
	"""
		Output
		----------
			Varscan
				{sample}.copynumber
				{sample}.copynumber.called
				{sample}.copynumber.called.gc
				{sample}.copynumber.called.homdel
				{sample}.copynumber.called.segments
		Options
		----------
			--output-file	Output file to contain the calls
			--min-coverage	Minimum read depth at a position to make a call [8]
			--amp-threshold	Lower bound for log ratio to call amplification [0.25]
			--del-threshold	Upper bound for log ratio to call deletion (provide as positive number) [0.25]
			--min-region-size	Minimum size (in bases) for a region to be counted [10]
			--recenter-up	Recenter data around an adjusted baseline > 0 [0]
			--recenter-down	Recenter data around an adjusted baseline < 0 [0]
	"""
	#samtools mpileup -q 1 -f ref.fa normal.bam tumor.bam | 
	#java -jar VarScan.jar copynumber varScan --mpileup 1
	#
	#samtools mpileup -f reference.fasta -q 1 -B normal.bam tumor.bam >normal-tumor.mpileup 
	#java -jar VarScan.jar somatic normal-tumor.mpileup normal-tumor.varScan.output --mpileup 1
	print("\tVarscan (copynumber)")
	patientID = sample['PatientID']
	logging.info("{0}: {1} Running Varscan (copynumber) {1}".format(patientID, '*'*30))
	default_options = {
			'--min-coverage': options['Parameters']['MIN_COVERAGE']	#Minimum read depth at a position to make a call [8]
			#'--amp-threshold':			#Lower bound for log ratio to call amplification [0.25]
			#'--del-threshold':	0.25	#Upper bound for log ratio to call deletion (provide as positive number) [0.25]
			#'--min-region-size': 10 	#	Minimum size (in bases) for a region to be counted [10]
			#'--recenter-up':			#Recenter data around an adjusted baseline > 0 [0]
			#'--recenter-down':			#Recenter data around an adjusted baseline < 0 [0]
	}



	varscan = options['Programs']['Varscan']
	samtools= options['Programs']['samtools']

	reference = options['Reference Files']['reference genome']

	output_folder = os.path.join(
		options['Pipeline Options']['copynumber pipeline folder'],
		sample['PatientID'], 'Varscan')

	varscan_prefix = os.path.join(output_folder, "{normal}_vs_{tumor}.varscan".format(
		normal = sample['NormalID'],
		tumor = sample['SampleID']))

	if not os.path.isdir(output_folder):
		os.makedirs(output_folder)
	varscan_basename = os.path.basename(varscan_prefix)
	varscan_output = varscan_prefix + '.copynumber.called'

	#-------------------------- Varscan (copynumber - Pileup) ----------------------------
	#Generates {varscan_prefix}.copynumber
	print("\t\t(1/4) Generating Pileup Files...")
	logging.info("{0} Varscan 1 of 4 (Pileup): Generating pileup files...".format(patientID))
	sample, pileup_log = Pileup(sample, options)
	logging.info("{0} Generated the Pileup files.".format(patientID))

	#-------------------------- Varscan (copynumber - caller) -----------------------------
	copynumber_command = "java {memory} -jar {varscan} copynumber {normal} {tumor} {prefix} --min-base_qual {mbq} --min-map-qual {mmq}"
	copynumber_command = copynumber_command.format(
		mmq = options['Parameters']['MIN_MAPPING_QUALITY'],
		mbq = options['Parameters']['MIN_NUCLEOTIDE_QUALITY'],
		normal = sample['NormalPileup'],
		tumor = sample['TumorPileup'],
		varscan = varscan,
		memory = options['Parameters']['JAVA_MAX_MEMORY_USAGE'],
		prefix = varscan_basename)

	print("\t\t(2/4) Calculating Ratios...")
	if not os.path.isfile(varscan_basename + '.copynumber'):
		logging.info("{0} Varscan 2 of 4 (copynumber command): {1}".format(patientID, copynumber_command))
		Terminal(copynumber_command, 
			label = 'Varscan (copynumber)', 
			show_output = False)
	else: 
		logging.info("{0} Varscan 2 of 4 (copynumber command): THe output already exists as {1}".format(patientID, varscan_basename + '.copynumber'))
		print('\t\tThe copynumber file already exists!')
		print("\t\t", varscan_basename)
	#varscan_prefix += '.copynumber'

	#--------------------------- Varscan Copycaller ----------------------------------
	copycaller_command = """java {memory} -jar {varscan} copyCaller {prefix}.copynumber \
		--output-file {prefix}.copynumber.called \
		--output-homdel-file {prefix}.copynumber.called.homdel \ 
		--min-coverage {mq}""".format(
			memory = options['Parameters']['JAVA_MAX_MEMORY_USAGE'],
			varscan = varscan,
			prefix = varscan_basename,
			mq = options['Parameters']['MIN_COVERAGE'])
	
	print("\t\t(3/4) Processing Ratios...")
	if not os.path.isfile(varscan_basename + '.copynumber.called'):
		logging.info("{0} Varscan 3 of 4 (copycaller): {1}".format(patientID, copycaller_command))
		Terminal(copycaller_command, 
			label = "Varscan (copycaller)", 
			show_output = False)
	else: 
		logging.info("{0} Varscan 3 of 4 (copycaller): The output already exists as {1}".format(patientID, varscan_basename + '.copynumber.called'))
		print('\t\tThe copycaller file already exists!')

	#--- Move the generated files from  the lcoal directory to the varscan directory----
	varscan_outputs = list()
	for suffix in ['.copynumber', '.copynumber.called', '.copynumber.called.gc', '.copynumber.called.homdel']:
		shutil.move(os.path.join(os.getcwd(), varscan_basename + suffix), varscan_prefix + suffix)
		varscan_outputs.append(varscan_prefix + suffix)
		
	
	#--------------------------- Circular Binary Segmentation --------------------------
	try:
		print("\t\t(4/4) Generating Segments...")
		rscript_file = os.path.join(
			options['Pipeline Options']['temporary folder'],
			sample['PatientID'],
			"{patient}.varscan_CBS.r".format(patient = sample['PatientID']))

		circular_binary_segmentation(
			varscan_prefix + '.copynumber.called', 
			varscan_prefix + '.copynumber.called.segments',
			rscript_file)
		logging.info("{0} Varscan 4 of 4 (Generate Segements): Successfully generated the copynumber segments.".format(patientID))
	except Exception as exception:
		logging.error("{0} Varscan 4 or 4 (Generate Segments): The segments could not be generated ({1})".format(patientID, exception))

	vc_log = {
		'Inputs': [sample['NormalPileup'], sample['TumorPileup']],
		'Outputs': [output_folder],
		'Commands': [copynumber_command, copycaller_command]
	}

	return sample, vc_log

class CNVkit(Caller):
	__name__ = "CNVkit"
	def runCallerWorkflow(self, sample, options):
		#super().__init__(self, sample, options, 'CNVkit')
		self.output_folder = os.path.join(
			options['Pipeline Options']['copynumber pipeline folder'],
			sample['PatientID'], self.caller_name)
		
		self.final_output = os.path.join(self.output_folder, 'reference_cnv.cnn')
		self.full_output = [
			sample['SampleID'] + '.cns',
			sample['SampleID'] + '.cnr',
			sample['SampleID'] + '.targetcoverage.cnn'
		]
		self.full_output = [os.path.join(self.output_folder, fn) for fn in self.full_output]
		reference_cnn = self.runBatchCommand(sample)

	def runBatchCommand(self, sample):
		reference_cnn = os.path.join(self.output_folder, "reference_cnv.cnn")
		expected_output = []
		batch_command = """{cnvkit} batch {tumor} --normal {normal} \
			--targets {targets} \
			--fasta {reference} \
			--output-reference {refcnn} \
			--output-dir {results} \
			--diagram --scatter""".format(
				cnvkit = self.program,
				tumor = sample['TumorBAM'],
				normal = sample['NormalBAM'],
				targets = sample['ExomeTargets'],
				reference = self.reference,
				refcnn = reference_cnn,
				results = self.output_folder)
		self.runCallerCommand(batch_command, self.final_output)

		expected_output = [
		]
		return expected_output

class FREEC(Caller):
	__name__ = "FREEC"
	def runCallerWorkflow(self, sample, options):
		self.output_folder = os.path.join(
			options['Pipeline Options']['copynumber pipeline folder'],
			sample['PatientID'], self.caller_name)
		
		
		self.prefix = os.path.join(self.output_folder, 
			"{normal}_vs_{tumor}.{prefix}".format(
			tumor   = sample['SampleID'], 
			normal  = sample['NormalID'],
			prefix = self.caller_name.lower()))
		
		self.samtools_program = options['Programs']['samtools']

		self.program_folder = self.program
		self.program 		= os.path.join(self.program_folder, 'src', 'freec')
		self.script_folder 	= os.path.join(self.program_folder, 'scripts')

		self.assess_significance_script = os.path.join(self.script_folder, "assess_significance.R")
		self.plot_script 				= os.path.join(self.script_folder, "makeGraph.R")
		self.config_file 				= self.prefix + "_config.txt"
		self.chrlenfile 				= os.path.join(self.output_folder, "chromosome_lengths.txt")
		#--------------------- Generate Pileup Files -----------------------------
		print("\t\t(1/6) Generate Pileup Files")
		normal_pileup 	= self.generatePileup(sample['NormalBAM'], sample['NormalID'])
		tumor_pileup 	= self.generatePileup(sample['TumorBAM'], sample['SampleID'])
		
		#--------------------- Generate ChrLen File  -----------------------------
		print("\t\t(2/6) Detect Chromosome Lengths")
		self.createChrLenFile(sample['ExomeTargets'])

		#--------------------- Generate Config File ------------------------------
		print("\t\t(3/6) Configure FREEC")
		#[general]
		#_configure_freec(sample, options, reference, chrlenfile, output_folder):
		all_options = self.configureFREEC(sample, normal_pileup, tumor_pileup, sample['ExomeTargets'])

		#------------------------------------- Run FREEC --------------------------------
		print("\t\t(4/6) Main Analysis")
		self.normal_CNVs 	= os.path.join(self.output_folder, "{0}.pileup_normal_CNVs".format(sample['SampleID']))
		self.normal_ratios	= os.path.join(self.output_folder, "{0}.pileup_normal_ratio.txt".format(sample['SampleID']))
		self.normal_baf_file = os.path.join(self.output_folder, "{0}.pileup_BAF.txt".format(sample['SampleID']))
		
		#files for the tumor sample
		self.sample_CNVs 		= os.path.join(self.output_folder, "{0}.pileup_CNVs".format(sample['SampleID']))
		self.sample_ratios 		= os.path.join(self.output_folder, "{0}.pileup_ratio.txt".format(sample['SampleID']))
		self.sample_baf_file 	= os.path.join(self.output_folder, "{0}.pileup_BAF.txt".format(sample['SampleID']))
		self.full_output = [self.normal_CNVs, self.normal_ratios, self.normal_baf_file,
							self.sample_CNVs, self.sample_ratios, self.sample_baf_file]
		self.runFREEC()

		#------------------------------------ Add Log2Ratios ----------------------------
		#--------------------------------- Calculate Significance -----------------------
		print("\t\t(5/6) Calculate Significance")
		#files for the normal sample

		self.normal_significance = self.calculateSignificance(self.normal_CNVs, self.normal_ratios)
		self.sample_significance = self.calculateSignificance(self.sample_CNVs, self.sample_ratios)
	
		print("\t\t(6/6) Generate Plots")

		self.normal_plot = os.path.join(self.output_folder, self.normal_ratios + ".png")
		self.sample_plot = os.path.join(self.output_folder, self.sample_ratios + ".png")
		self.generatePlots(self.normal_ratios, self.normal_baf_file)
		self.generatePlots(self.sample_ratios, self.sample_baf_file)
	
	def createChrLenFile(self, targets_file):
		
		genome_index_file = self.reference + '.fai'
		exclude_chroms = ['chr1_KI270706v1_random', 'chr4_GL000008v2_random', 'chr14_GL000009v2_random', 'chrUn_KI270742v1']
		#generated as the intersection of the genome index and exome targets files
		if os.path.exists(targets_file):
			with open(targets_file, 'r') as bedfile:
				bedchrs = set([i[0] for i in csv.reader(bedfile, delimiter = '\t')])
		else:
			self._fileNotFound(targets_file)
		
		if os.path.exists(genome_index_file):
			with open(genome_index_file, 'r') as file1:
				chrlens = ['\t'.join(i[:3]) + '\n' for i in csv.reader(file1, delimiter = '\t') if (i[0] in bedchrs and i[0] not in exclude_chroms)]
		else:
			self._fileNotFound(genome_index_file)
		
		with open(self.chrlenfile, 'w') as file1:
			[file1.write(i) for i in chrlens]

	
	def calculateSignificance(self, cnvs, ratios):
		expected_output = cnvs + '.p.value.txt'
		command = """cat {script} | R --slave --args {CNVs} {ratios}""".format(
			script = self.assess_significance_script,
			CNVs = cnvs,
			ratios = ratios)
		
		self.runCallerCommand(command, expected_output)
		return expected_output

	def generatePlots(self, ratios, baf_file):
		#--------------------------------- Generate Plots --------------------------------
		expected_output = ratios + '.png'
		command = """cat {script} | R --slave --args {ploidy} {ratios} {baf}""".format(
			script = self.plot_script,
			ploidy = "2",
			ratios = ratios,
			baf = baf_file)
		self.runCallerCommand(command, expected_output)

	def configureFREEC(self, sample, normal_pileup, tumor_pileup, targets):
		general_options = {
			'bedtools': "/usr/bin/bedtools",
			'breakPointThreshold': 0.8, #Default: 0.8 use something like 0.6 to get more segments (and thus more predicted CNVs)
			'chrFiles': os.path.join(os.path.dirname(self.reference), 'chromosomes'),
			'chrLenFile': self.chrlenfile, #a list of chromosomes and chromosome lengths. Basically the reference dict.
			'maxThreads': '6',
			'noisyData': 'TRUE',
			'outputDir': self.output_folder,
			'ploidy': '2',
			'printNA': 'FALSE', #set FALSE to avoid printing "-1" to the _ratio.txt files Useful for exome-seq or targeted sequencing data
			'readCountThreshold': '50', #Default: 10, recommended value >=50 for for exome data
			'samtools': self.samtools_program,
			#'sex': "",
			'window': "0" #for whole exome sequencing: "window=0"
		}
		general_options = ["{0} = {1}".format(k, v) for k, v in general_options.items()]
		
		#[sample]
		sample_options = {
			'mateFile': tumor_pileup,
			'inputFormat': 'pileup',
			'mateOrientation': '0'
		}
		sample_options = ["{0} = {1}".format(k, v) for k, v in sample_options.items()]

		#[control]
		control_options = {
			'mateFile': normal_pileup,
			'inputFormat': 'pileup',
			'mateOrientation': '0'
		}
		
		control_options = ["{0} = {1}".format(k, v) for k, v in control_options.items()]

		#[BAF]
		baf_options = {
			'SNPfile': self.dbSNP,
			"minimalCoveragePerPosition": '5'
		}
		baf_options = ["{0} = {1}".format(k, v) for k, v in baf_options.items()]

		#[target]
		target_options = {
			"captureRegions": targets,
		}
		target_options = ["{0} = {1}".format(k, v) for k, v in target_options.items()]

		all_options = ["[general]"] + general_options
		all_options +=["[sample]"] 	+ sample_options
		all_options +=["[control"] 	+ control_options
		all_options +=["[BAF]"]		+ baf_options
		all_options +=["[target]"] 	+ target_options

		with open(self.config_file, 'w') as file1:
			for line in all_options:
				file1.write(line + '\n')
		return os.path.exists(self.config_file)

	def runFREEC(self):
		#expected_output = [self.normal_CNVs, self.sample_CNVs]
		freec_command = "{freec} -conf {config}".format(
			freec = self.program,
			config = self.config_file)
		self.runCallerCommand(freec_command, self.full_output)



#----------------------------------------------------------------------------------------------------
#------------------------------------------ Pipelines -----------------------------------------------
#----------------------------------------------------------------------------------------------------
class PipelineOBS:
	def __init__(self, sample, options, callers):
		self.pipeline_name = type(self).__name__
		print("   Running ", self.pipeline_name)
		#LOGGER.info("{0}: Running {1}".format(sample['PatientID'], self.pipeline_name))
		self.start = now()
		self.caller_logs = list()
		#self.pipeline_folder = self._get_pipeline_folder(options)

		
		callers = self._get_callers(callers) #Dictionary of caller methods

		self._run_pipeline(sample, options, callers)
		
		self.stop = now()

		self._update_log(sample, options, callers, self.start, self.stop)

	
	@staticmethod
	def _get_callers():
		return []

	def _run_caller(self, caller, sample, options):
		start = now()
		sample, caller_log = run_program(caller, sample, options)
		self.caller_logs.append(caller_log)
		
		stop = now()
		duration = stop - start


		

		csv_log([caller_log])
		return sample

	def _run_pipeline(self, sample, options, callers):
		
		for caller_name, caller in callers.items():
			pstart = now()
			sample = self._run_caller(caller, sample, options)
			pstop = now()
			duration = pstop - pstart
			LOGGER.info("{0}: Ran {1} in {2}".format(sample['PatientID'], caller_name, duration))
		
	def _update_log(self, sample, options, callers, start, stop):
		intermediate_files = list()
		input_files = list()
		output_files = list()
		for caller_log in self.caller_logs:
			intermediate_files += caller_log.get('Intermediate Files', [])
			input_files += caller_log.get('Input Files', [])
			output_files += caller_log.get('Output Files', [])

		intermediate_files = sorted(set(intermediate_files))
		input_files = sorted(set(input_files))
		output_files = sorted(set(output_files))
		try:
			input_size = sum([getsize(i) for i in input_files])
			output_size = sum([getsize(i) for i in output_files])
		except:
			input_size = output_size = 0
		status = True
		pipeline_log = {
			'PatientID':  sample['PatientID'],
			'Program':    self.pipeline_name,
			'Start Time': start.isoformat(),
			'Stop Time':  stop.isoformat(),
			'Duration':   stop - start,
			'Inputs':     os.path.dirname(sample['TumorBAM']),
			'Input Size': input_size,
			'Intermediate Files': intermediate_files,
			'Outputs':    options['Pipeline Options']['somatic pipeline folder'],
			'Output Size':output_size,
			'Notes': 	  "",
			'Status':	  status,
			'Commands':   ",".join(sorted(callers.keys()))
		}
		csv_log([pipeline_log])



class Pipeline:
	def __init__(self, sample, options, sample_callers, DEBUG = False):
		pipeline_callers = self._getCallers(sample_callers)

		for caller_name, callerClass in sorted(pipeline_callers):
			print(caller_name)
			callset = callerClass(sample, options, DEBUG)

class SomaticPipeline(Pipeline):
	@staticmethod
	def _getCallers(callers):
		available_callers = {
			#'pon': Mutect_pon_detection,
			#'gdc': GDC_somatic,
			'muse': MuSE,
			'mutect2': MuTect2,
			'somaticsniper': SomaticSniper,
			'strelka': Strelka,
			'varscan': Varscan,
			'haplotypecaller': HaplotypeCaller,
			'unifiedgenotyper': UnifiedGenotyper
		}
		callers = [(i,j) for i, j in available_callers.items() if i in callers]
		return callers
	@staticmethod
	def _get_pipeline_folder(options):
		folder = os.path.join(options['working directory'], "3_called_variants")
		return folder

class CopynumberPipeline(Pipeline):
	@staticmethod
	def _getCallers(callers):
		available_callers = {
			'varscan': VarscanCopynumber,
			'cnvkit': CNVkit,
			'freec': FREEC
		}
		callers = [(i,j) for i, j in available_callers.items() if i in callers]
		return callers
	@staticmethod
	def _get_pipeline_folder(options):
		return os.path.join(options['working directory'], "4_copynumber_pipeline")


#----------------------------------------------------------------------------------------------------
#--------------------------------------------- Main -------------------------------------------------
#----------------------------------------------------------------------------------------------------
class SampleBAMFiles:
	""" A class for locating and validating the BAM files for each sample.
	"""
	def __init__(self, sample, config, verify_md5sum = False, DEBUG = False):
		"""
		"""
		self.DEBUG = DEBUG
		if self.DEBUG:
			print("SampleBAMFiles({0}, verify_md5sum = {1}, DEBUG = {2})".format(sample['PatientID'], verify_md5sum, DEBUG))
		sample['NormalBAM'] = self._fetchFilename(sample.get('NormalBAM'), sample['NormalUUID'])
		sample['TumorBAM'] = self._fetchFilename(sample.get('TumorBAM'), sample['SampleUUID'])

		if verify_md5sum:
			normal_md5sum = API(sample['NormalUUID'], 'files')
			tumor_md5sum = API(sample['SampleUUID'], 'files')
		else:
			normal_md5sum = None
			tumor_md5sum = None

		normal_file_status = self._verifyFile(sample['NormalBAM'], normal_md5sum)
		tumor_file_status  = self._verifyFile(sample['TumorBAM'], tumor_md5sum)

		if self.DEBUG:
			self.status = True
		else:
			self.status = normal_file_status and tumor_file_status

	def _fetchFilename(self, filename, file_id):
		"""
			Parameters
			----------
				filename: string [PATH]
					Path to a BAM file. if the filename is None or "", the program will
					attempt to locate the file on the current computer.
				file_id: string
					UUID of the file.
				options
		"""
		if filename is None or not os.path.exists(filename):
			filename = gdc_api.getFileLocation(file_id)['file_location']

		return filename

	@staticmethod
	def _verifyFile(filename, expected_md5sum):
		""" Verifies a single file
			If expected_md5sum is None, the program will skip checking it.
		"""
		file_exists = os.path.isfile(filename)
		if file_exists and expected_md5sum is not None:
			file_md5sum = generate_file_md5(filename)
			file_md5sum_status = file_md5sum == expected_md5sum
		else: file_md5sum_status = True

		file_is_valid = file_exists and file_md5sum_status

		return file_is_valid

class GenomicsPipeline:
	def __init__(self, sample_filename, **kwargs):
		"""
			Required Arguments
			------------------
				sample_list
			Optional Arguments
			------------------
				somatic
				copynumber
				caller_status_filename
				config_filename
		"""
		if kwargs['parser'].debug:
			print("GenomicsPipeline")
			for k, v in kwargs.items(): print('\t', k, '\t', v)


		self.parser = kwargs['parser']
		#Parse the provided arguments
		somatic_callers = [i.lower() for i in kwargs.get("somatic", [])]
		copynumber_callers = [i.lower() for i in kwargs.get("copynumber", [])]

		config_filename = kwargs.get('config_filename', os.path.join(PIPELINE_DIRECTORY, "0_config_files", "pipeline_project_options.txt"))
		caller_status_filename = kwargs.get('caller_status_filename', os.path.join(PIPELINE_DIRECTORY, "0_config_files", "caller_status.tsv"))

		sample_list, config = self._loadPipelineConfiguration(sample_filename, config_filename)		
		
		LOGGER.info("Somatic Callers: " + ', '.join(somatic_callers))
		LOGGER.info("Copynumber Callers: " + ', '.join(copynumber_callers))
		LOGGER.info("Running through the genomics pipeline with {0} samples.".format(len(sample_list)))
		#sample_list = []
		for index, sample in enumerate(sample_list):
			print("({0}/{1}) {2}\t{3}".format(index+1, len(sample_list), sample['PatientID'], now().isoformat()), flush = True)

			sample_status = self.runSample(sample, config, somatic_callers, copynumber_callers)



	@staticmethod
	def _loadPipelineConfiguration(sample_filename, config_filename):
		if not os.path.isfile(sample_filename):
			message = "The sample list does not exists at " + sample_filename
		elif not os.path.isfile(config_filename):
			message = "The config file does not exist at " + config_filename
		else:
			message = None
		if message is not None:
			raise FileNotFoundError(message)

		sample_list = readTSV(sample_filename)

		#----------------------------------------Process Config ------------------------------------------
		config = configparser.ConfigParser()
		config.read(config_filename)
		
		return sample_list, config
	
	@staticmethod
	def _makePatientFolders(patientID, options):

		snv_folder = options['Pipeline Options']['somatic pipeline folder']
		cnv_folder = options['Pipeline Options']['copynumber pipeline folder']
		temp_folder= options['Pipeline Options']['temporary folder']

		snv_folder = os.path.join(snv_folder,patientID)
		cnv_folder = os.path.join(cnv_folder, patientID)
		temp_folder = os.path.join(temp_folder, patientID)

		checkdir(snv_folder)
		checkdir(cnv_folder)
		checkdir(temp_folder)

	@staticmethod
	def _useSample(sample):
		use_value = sample.get('Use', True)

		if isinstance(use_value, str): use_value = use_value.lower()
		
		manually_skipped = use_value in {False, 'false', 0, '0', 'no'}

		already_completed = False

		invalid_targets = not os.path.isfile(sample['ExomeTargets'])



		use_status = {
			'status': not (manually_skipped or already_completed or invalid_targets),
			'manual': manually_skipped,
			'completed': already_completed,
			'targets': invalid_targets
		}

		return use_status

	def generate_readme(self, options):
		readmefile = "readme.txt"

		general_notes = [
			'Time Started                = ' + "" 
			'Time Ended                  = ' + now().isoformat(),
			'Pipeline Directory          = ' + PIPELINE_DIRECTORY]

		reference_files = [
			'\n------------------- Reference Files --------------------',
			]
		reference_files += sorted(["{0:<20} = {1}".format(k, v) for k, v, in options['Reference Files'].items()])

		global_parameter_values = [
			'\n------------------ Common Parameters -------------------',
			"MAX_CORES                   = {0}".format(options['Parameters']['MAX_CORES']),
			"MIN_MAPPING_QUALITY         = {0}".format(options['Parameters']['MIN_MAPPING_QUALITY']),
			"MIN_BASE_POSITION_FREQUENCY = {0}".format(options['Parameters']['MIN_BASE_POSITION_FREQUENCY']),
			"MIN_COVERAGE                = {0}".format(options['Parameters']['MIN_COVERAGE']),
			"MIN_NUCLEOTIDE_QUALITY      = {0}".format(options['Parameters']['MIN_NUCLEOTIDE_QUALITY']),
			"\n-------------- GATK Specific parameters -----------------",
			"GATK_NUM_THREADS            = {0}".format(options['Parameters']['GATK_NUM_THREADS']),

			"\n------------- Bambino Specific parameter ----------------",
			"MIN_FLANKING_QUALITY        = {0}".format(options['Parameters']['MIN_FLANKING_QUALITY']),
			"MIN_ALT_ALLELE_COUNT        = {0}".format(options['Parameters']['MIN_ALT_ALLELE_COUNT']),
			"MIN_MINOR_FREQUENCY         = {0}".format(options['Parameters']['MIN_MINOR_FREQUENCY']),
			"MMF_MAX_HQ_MISMATCHES       = {0}".format(options['Parameters']['MMF_MAX_HQ_MISMATCHES']),
			"MMF_MIN_HQ_THRESHOLD        = {0}".format(options['Parameters']['MMF_MIN_HQ_THRESHOLD']),
			"MMF_MAX_LQ_MISMATCHES       = {0}".format(options['Parameters']['MMF_MAX_LQ_MISMATCHES']),
			"UNIQUE_FILTER_COVERAGE      = {0}".format(options['Parameters']['UNIQUE_FILTER_COVERAGE']),

			"\n-------------- Mpileup specific parameter ---------------",
			"BWA_DOWNGRADE_COFF          = {0}".format(options['Parameters']['MIN_FLANKING_QUALITY']),
			"NO_OF_READS_TO_CONSIDER_REALIGNMENT = {0}".format(options['Parameters']['MIN_FLANKING_QUALITY']),
			"FREQ_OF_READS               = {0}".format(options['Parameters']['MIN_FLANKING_QUALITY']),
			"MPILEUP_QUALITY_THRESHOLD   = {0}".format(options['Parameters']['MIN_FLANKING_QUALITY']),

			"\n----------- Somatic Sniper specific parameter ----------",
			"SOMATIC_QUALITY             = {0}".format(options['Parameters']['SOMATIC_QUALITY']),

			"\n--------------- Java specific parameters ---------------"
			"JAVA_MAX_MEMORY_USAGE       = {0}".format(options['Parameters']['JAVA_MAX_MEMORY_USAGE']),
			"JAVA_GARBAGE_COLLECTION     = {0}".format(options['Parameters']['JAVA_GARBAGE_COLLECTION'])]

		with open(readmefile, 'w') as file1:
			for line in general_notes + reference_files + global_parameter_values:
				file1.write(line + '\n')
	
	def _getSampleCompletedCallers(self, sample):

		pass
	def runSample(self, sample, config, somatic_callers, copynumber_callers):
		print("#"*180)
		print("#"*90 + sample['PatientID'] + '#'*90)
		print("#"*180)

		if self.parser.debug:
			print("PatientID: ", sample['PatientID'])
			print("\tSomatic Callers: ", somatic_callers)
			print("\tCopynumber Callers: ", copynumber_callers)

		sample_start = now()
		
		use_this_sample = self._useSample(sample)['status']
		sample_caller_status = self._getSampleCompletedCallers(sample)

		self._makePatientFolders(sample['PatientID'], config)

		#normal_file_api, tumor_file_api = self._get_file_info(sample)
			
		#--------------------------- Download and verify the BAM files ---------------------------------
		try:
			prepared_files = SampleBAMFiles(sample, config, DEBUG = self.parser.debug)
			file_status = prepared_files.status
			if file_status:
				LOGGER.info("{0}: The BAM files exist and are valid. NormalBAM={1}. TumorBAM={2}".format(sample['PatientID'],sample['NormalBAM'], sample['TumorBAM']))
			else:
				LOGGER.error("{0}: The BAM files are invalid! NormalBAM={1}. TumorBAM={2}".format(sample['PatientID'],sample['NormalBAM'], sample['TumorBAM']))
		except Exception as exception:
			file_status = False
			message = "{0}: GenomicsPipeline.run_sample: The BAM files could not be loaded ({1})".format(sample['PatientID'], str(exception))
			print(message)
			LOGGER.critical(message)

		if self.parser.debug:
			print("\tFile Status: ", file_status)



		if file_status:
			if not self.parser.ignore_caller_status and not self.parser.debug:
				somatic_callers = [i for i in somatic_callers if i not in sample_completed_callers]
			somatic_pipeline = SomaticPipeline(sample, config, somatic_callers, DEBUG = self.parser.debug)
			if self.parser.debug:
				print("\t\tSomatic Callers: ", somatic_callers)

			if not self.parser.ignore_caller_status and not self.parser.debug:
				copynumber_callers = [i for i in somatic_callers if i not in sample_completed_callers]
			copynumber_pipeline = CopynumberPipeline(sample, config, copynumber_callers)
			if self.parser.debug:
				print("Copynumber Callers: ", copynumber_callers)
		else:
			logging.critical("{0}: The BAM files are invalid!".format(sample['PatientID']))
			print("\tThe BAM files did not download correctly!")
		sample_stop = now()


		sample_log = {
			'PatientID':  sample['PatientID'],
			'Program':    "Genomics Pipeline",
			'Start Time': sample_start.isoformat(),
			'Stop Time':  sample_stop.isoformat(),
			'Duration':   sample_stop - sample_start,
			'Inputs':     os.path.dirname(sample['TumorBAM']),
			'Input Size': None,
			'Intermediate Files': None,
			'Outputs':    PIPELINE_DIRECTORY,
			'Output Size':None,
			'Notes': 	  "",
			'Status':	  file_status,
			'Commands':   "Genomics Pipeline({0})".format(sample['PatientID'])
		}

		return file_status

def getCMDArgumentParser():

	"""Parse the command-line arguments for this program."""

	parser = ArgumentParser(
		description='Run the genomics pipeline')

	show_default = ' (default %(default)s)'

	parser.add_argument('-l', '--sample-list', 
		action='store',
		dest = 'sample_list',
		help='the sample list')

	parser.add_argument('-d', "--debug",
		dest = 'debug',
		action = 'store_true',
		help='debug the pipeline using default settings' + show_default)

	parser.add_argument("-i", "--ignore-caller-status",
		dest = "ignore_caller_status",
		action = 'store_false',
		help = "Ignore the caller status file.")

	return parser

LOGGER = configurePipelineLogger()
API = gdc_api.GDCAPI()

if __name__ == "__main__":

	CMD_PARSER = getCMDArgumentParser().parse_args()

	config_filename 		= os.path.join(PIPELINE_DIRECTORY, "0_config_files", "pipeline_project_options.txt")
	caller_status_filename 	= os.path.join(PIPELINE_DIRECTORY, "0_config_files", "caller_status.tsv")
	
	if CMD_PARSER.debug:
		sample_filename = os.path.join(PIPELINE_DIRECTORY, "sample_list.tsv")
		somatic_callers = ['MuSE', 'Varscan', 'Strelka', 'SomaticSniper', 'Mutect2', "HaplotypeCaller", "UnifiedGenotyper"]
		copynumber_callers = ['Varscan', 'CNVkit', 'FREEC']
	else:
		sample_filename = CMD_PARSER.sample_list
		somatic_callers = ['MuSE', 'Varscan', 'Strelka', 'Somaticsniper', 'Mutect']
		copynumber_callers = ['varscan', 'cnvkit', 'freec']	

	pipeline = GenomicsPipeline(sample_filename, config_filename = config_filename, somatic = somatic_callers, copynumber = copynumber_callers, parser = CMD_PARSER)
else:
	pass
#3872