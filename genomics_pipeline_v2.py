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

PIPELINE_DIRECTORY = "/home/upmc/Documents/Variant_Discovery_Pipeline"
CONSOLE_LOG_FILE = ""
SAMPLE_LOG_FILE = os.path.join(PIPELINE_DIRECTORY, "0_config_files", "sample_logV2.tsv")
README_FILE = os.path.join(PIPELINE_DIRECTORY, "0_readme_files", "readme.{0}.txt".format(now().isoformat()))
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

def Terminal(command, label = None, show_output = False, timeout = None):
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
	else:
		command = shlex.split(command)
		process = subprocess.Popen(command, stdout = subprocess.PIPE, stderr = subprocess.STDOUT)
		output = str(process.stdout.read(),'utf-8')
		with open(CONSOLE_LOG_FILE, 'a') as console_file:
			console_file.write(now().isoformat() + '\n')
			console_file.write(output + '\n\n')

	return process
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
	def __init__(self, sample, options):
		#self.caller_name = caller_name
		program_start = now()
		self.caller_name = self.__name__
		self.setCallerEnvironment()
		self.runCallerWorkflow()
		program_stop = now()
		self.updateSampleLog(sample, program_start, program_stop)
	
	def setCallerEnvironment(self):
		self.caller_name = caller_name
		self.reference 	= options['Reference Files']['reference genome']
		self.dbSNP 		= options['Reference Files']['dbSNP']
		self.cosmic 	= options['Reference Files']['COSMIC']


		self.program 			= options['Programs'].get(caller_name)
		self.gatk_program 		= options['Programs']['GATK']
		self.max_cpu_threads 	= options['Parameters']['MAX_CORES']
		self.max_memory_usage 	= options['Parameters']['JAVA_MAX_MEMORY_USAGE']
		self.min_base_quality 	= options['MIN_NUCLEOTIDE_QUALITY']
		self.min_mapping_quality = options['Parameters']['MIN_MAPPING_QUALITY'],
		self.min_somatic_quality = options['Parameters']['SOMATIC_QUALITY']
		
		self.output_folder = os.path.join(
			options['Pipeline Options']['somatic pipeline folder'],
			sample['PatientID'], caller_name)
		self.temp_folder = os.path.join(
			options['Pipeline Options']['temporary folder'], 
			sample['PatientID'], caller_name)

		self.command_list = list()
		self.temp_files = list()
		
		checkdir(self.output_folder)
		checkdir(self.temp_folder, True)

		self.prefix = os.path.join(output_folder, 
			"{normal}_vs_{tumor}.{prefix}".format(
			tumor   = sample['SampleID'], 
			normal  = sample['NormalID']),
			prefix = caller_name.lower())

	def runCallerCommand(self, command, expected_output):
		self.command_list.append(command)
		if isinstance(expected_output, str): expected_output = [expected_output]
		if expected_output is None or any([not os.path.exists(fn) for fn in expected_output]):
			Terminal(command)
		else:
			print(os.path.basename(expected_output), 'already exists!')
	def runCallerWorkflow(self):
		pass

	def generatePileup(self, bam_file, bam_name):
		""" Generates pileup files. If single is True, only
			one file will be generated.
		"""
		output_file = os.path.join(self.temp_folder, bam_name + '.mpileup')
		pileup_command = "samtools mpileup -q 1 -B -f {reference} {sample} > {output}".format(
			reference = self.reference,
			sample = bam_file,
			output = outpuit_file)

		self.runCallerCommand(pileup_command, output_file)
		return output_file

	def updateSampleLog(self, sample, program_start, program_stop):
		status = self.getCallerStatus()
		duration = program_stop - program_start
		caller_log = {
			'PatientID':  sample['PatientID'],
			'Program':    self.caller_name,
			'Start Time': program_start.isoformat(),
			'Stop Time':  program_stop.isoformat(),
			'Duration':   duration,
			'Intermediate Files': ', '.join(self.temp_files),
			'Outputs':    ', '.join(self.full_output),#'|'.join(outputs),
			'Notes': 	  "",
			'Status':	  status,
			'Commands':   '|'.join(self.command_list)#'|'.join(commands)
		}
		writeheaders = os.getsize(SAMPLE_LOG_FILE) == 0
		line 	= sorted(caller_log.items())
		line 	= reversed(line)
		headers = '\t'.join([i[0] for i in line]) + '\n'
		line 	= '\t'.join([i[1] for i in line]) + '\n'

		with open(SAMPLE_LOG_FILE, 'a') as file1:
			if writeheaders: file1.write(headers)
			file1.write(line)

	def getCallerStatus(self):
		caller_failed = any(not os.path.exists(fn) for fn in self.full_output)
		status = not Failed

		return status

class cCaller(Caller):
	""" Base class for the copynumber callers """
	def __init__(self, sample, options, caller_name):
		super().__init__(self, sample, options, caller_name)

		self.output_folder = os.path.join(
			options['Pipeline Options']['copynumber pipeline folder'],
			sample['PatientID'], caller_name.lower())

class MuSE(Caller):
	def runCallerWorkflow(self, sample, options):
		#self.caller_name = 'MuSE'
		#super().__init__(self, sample, options, caller_name)
		call_output = self.runMuseCall(sample)
		sump_output = self.runMuseSump(call_output)

		self.full_output = [call_output, sump_output]
		self.final_output= sump_output

	def runMuseCall(self, sample):
		muse_call_output = self.prefix + '.txt'
		call_command = "{program} call -O {prefix}.MuSE.txt -f {reference} {tumor} {normal}".format(
			program = self.program,
			reference = self.reference,
			prefix = '.'.join(self.prefix.split('.')[:-1]),
			tumor = sample['TumorBAM'],
			normal = sample['NormalBAM'])
		
		self.runCallerCommand(call_command, muse_call_output)

		return muse_call_output

	def runMuseSump(self, call_output):
		
		muse_sump_output = output_prefix + ".MuSe.vcf"
		sump_command = "{program} sump -I {call_output} -E -D {dbSNP} -O {output}".format(
			program = self.program,
			prefix = self.prefix,
			dbSNP = self.dbSNP,
			call_output = call_output,
			output = muse_sump_output)

		self.runCallerCommand(sump_command, muse_sump_output)

		return muse_sump_output

class MuTect2(Caller):
	def runCallerWorkflow(self, sample, options):
		#super().__init__(self, sample, options, 'MuTect2')

		mutect2_output = self.runMutect2(sample)
		self.full_output = [mutect2_output]
		self.final_output = mutect2_output

	def runMutect2(self, sample):
		mutect2_output = self.prefix + '.vcf'

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
			memory      = self.memory_usage,
			dbSNP       = self.dbSNP,
			cosmic 		= self.cosmic,
			reference   = self.reference, 
			normal      = sample['NormalBAM'], 
			tumor       = sample['TumorBAM'], 
			targets     = sample['ExomeTargets'],
			output      = mutect2_output)

		self.runCallerCommand(mutect2_command, mutect2_output)

		return mutect2_output

class SomaticSniper(Caller):
	def runCallerWorkflow(self, sample, options):
		#Will print "Couldn't find single-end mapping quality. Check to see if the SM tag is in BAM."
		#This doesn't invalidate results, but try not to use single-end mapping quality in output
		#super().__init__(self, sample, options, 'SomaticSniper')
		
		self.readcount_program = options['Programs']['bam-readcount']
		self.samtools_program = options['Programs']['samtools-0.1.6']
		self.somaticsniper_folder = self.program
		self.program = os.path.join(self.somaticsniper_folder, 'build', 'bin', 'bam-somaticsniper')
		self.scripts_folder = os.path.join(self.somaticsniper_folder, 'src', 'scripts')

		self.snpfilter_script= os.path.join(self.script_folder, 'snpfilter.pl')
		self.readcount =       os.path.join(self.script_folder, 'prepare_for_readcount.pl')
		self.hc_script =       os.path.join(self.script_folder, 'highconfidence.pl')
		self.fpfilter  =       os.path.join(self.script_folder, 'fpfilter.pl')

		raw_output = self.runVariantDiscovery(sample)
		#Generate pileup files
		normal_pileup_file = self.generatePileupFile(sample['NormalBAM'], 'normal')
		tumor_pileup_file  = self.generatePileupFile(sample['TumorBAM'], 'tumor')

		#Filter LOH
		_intermediate_file = self.prefix + ".SNPfilter.intermediate"
		loh_filtered_output  = self.prefix + ".SNPfilter.final"
		_intermediate_file  = self.removeLOH(raw_output, normal_pileup_file, _intermediate_file)
		loh_filtered_output = self.removeLOH(_intermediate_file, tumor_pileup_file, loh_filtered_output)

		
		readcounts = self.readcounts(loh_filtered_output, sample['TumorBAM'])
		false_positive_output = self.removeFalsePositives(loh_filtered_output, readcounts)
		high_confidence_variants = self.calculateConfidence(self)

		suffixes = ['.vcf', ".SNPfilter.intermediate", ".SNPfilter.final", ".SNPfilter.final.pos", ".SNPfilter.final.readcounts.rc",
			".SNPfilter.final.fp_pass", ".lq.vcf", ".hq.vcf"]
		self.full_output = [self.prefix + s for s in suffixes]
		self.final_output = high_confidence_variants
	@staticmethod
	def _setDefultOptions():
		default_options = {
			'q': options['Parameters']['MIN_MAPPING_QUALITY'], #filtering reads with mapping quality less than INT [0]
			'Q': options['Parameters']['SOMATIC_QUALITY'], #filtering somatic snv output with somatic quality less than INT [15]
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
		output_file = self.prefix + '.vcf'
		default_options = self._setDefultOptions()
		somaticsniper_command = "{program} {options} -f {reference} {tumor} {normal} {outputfile}".format(
			program 	= self.program,
			options 	= default_options,
			reference 	= self.reference,
			tumor 		= sample['TumorBAM'],
			normal 		= sample['NormalBAM'],
			outputfile 	= output_file)
		self.runCallerCommand(somaticsniper_command, output_file)
		return output_file

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
		self.runCallerCommand(samtools_command, pileup_file)
		self.runCallerCommand(filter_command, output_file)

		return output_file

	def removeLOH(self, vcf_file, pileup_file, output_file):
		#-------------------------- Filter and remove LOH --------------------------------
		loh_filter_command = "perl {snpfilter} --snp-file {vcf} --indel-file {pileup} --out-file {output}".format(
			snpfilter = self.snpfilter_script,
			vcf = vcf_file,
			pileup = pileup_file,
			output = output_file)

		self.runCallerCommand(loh_filter_command, output_file)
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
		self.runCallerCommand(pr_command, prepare_readcount_output)
		#Readcounts

		readcount_command = "{program} -b {mbq}  -q 1 -f {reference} -l {proutput} {tumor} > {output}".format(
			program = self.readcount_program,
			mbq = self.min_base_quality,
			reference = self.reference,
			proutput = prepare_readcount_output,
			tumor = sample['TumorBAM'],
			output = readcount_output)
		self.runCallerCommand(readcount_command, readcount_output)

		return readcount_output

	def removeFalsePositives(self, loh_file, readcounts):
		false_positive_output = loh_file + '.fp_pass'
		fp_command = "perl {fpfilter} --snp-file {snpfilter} -readcount-file {readcounts}".format(
			fpfilter = self.fpfilter,
			snpfilter = loh_file,
			readcounts = readcounts)
		self.runCallerCommand(fp_command, false_positive_output)
		return false_positive_output

	def calculateConfidence(self, vcf_file):
		output_file = self.prefix + '.hq.vcf'
		command = "perl {script} --snp-file {vcf} --min-mapping-quality {mmq} --min-somatic-score {ss} --lq-output {prefix}.lq.vcf --out-file {prefix}.hq.vcf".format(
			script = self.hc_script,
			mmq = self.min_mapping_quality,
			ss = self.min_somatic_quality,
			vcf = vcf_file,
			prefix = self.prefix)
		self.runCallerCommand(command, output_file)
		return output_file

class Strelka(Caller):
	def runCallerWorkflow(self, sample, options):
		#super().__init__(self, sample, options, 'Strelka')
		
		patient_folder = os.path.dirname(self.output_folder)
		self.results_folder = os.path.join(patient_folder, 'Strelka', 'results')
		self.config_script = os.path.join(self.program, 'bin', 'configureStrelkaWorkflow.pl')


		strelka_folder = os.path.join(patient_folder, 'Strelka')
		
		#--------------------------------- Configure Strelka --------------------------------
		#change working directory
		config_script = self.generateConfigScript(self.config_script)
		self.configureStrelka(sample, config_script)
		output_files = self.runStrelka()

		self.full_output = output_files
		self.full_output = self.prefix + '.passed.somatic.snvs.vcf'

	
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
		
		strelka_run_command = "make -j {threads} -C {outputdir}".format(
			outputdir = self.output_folder, threads = self.max_cpu_threads)
		self.runCallerCommand(strelka_run_command, self.results_folder)

		output_suffixes = ['.all.somatic.indels.vcf', '.all.somatic.snvs.vcf', '.passed.somatic.indels.vcf', '.passed.somatic.snvs.vcf']
		output_files = [self.prefix + i for i in output_suffixes]

		return output_files

	def getStrelkaConfiguration():
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
	def runCallerWorkflow(self, sample, options):
		#super().__init__(self, sample, options, 'Varscan')

		normal_pileup 	= self.generatePileup(sample['NormalBAM'], sample['NormalID'])
		tumor_pileup 	= self.generatePileup(sample['TumorBAM'], sample['SampleID'])

		raw_snp, raw_indel = self.runVarscan(normal_pileup, tumor_pileup)

		processed_variants = self.postProcessing(raw_snp)
		self.full_output = [raw_snp, raw_indel] + processed_variants
		self.final_output = processed_variants[0]

	def runVarscan(self, normal_pileup, tumor_pileup):

		output_files = [self.prefix +'.raw.snp.vcf', self.prefix + '.raw.indel.vcf']
		command = "java {memory} -jar {varscan} somatic {normal} {tumor} --output-snp {prefix}.raw.snp.vcf --output-indel {prefix}.raw.indel.vcf --output-vcf 1".format(
			varscan = self.prgoram,
			memory = self.max_memory_usage,
			normal = normal_pileup,
			tumor = tumor_pileup,
			prefix =self.prefix,
			mc = self.min_coverage)
		self.runCallerCommand(command, output_files)

		return output_files

	def postProcessing(self, snp_vcf):
		output_file = self.prefix + '.snp.Somatic.hc'

		expected_output = [self.prefix + i for i in ['.snp.Somatic.hc', '.snp.Somatic.lc', '.snp.Germline', '.snp.LOH', '.snp.Germline.hc', '.snp.LOH.hc']]

		process_command = "java {memory} -jar {varscan} processSomatic {snp_vcf}"
		process_command = process_command.format(
			memory  = self.max_memory_usage,
			varscan = self.program,
			vcf  = snp_vcf)
		self.runCallerCommand(process_command, output_file)
		return output_files



class HaplotypeCaller(Caller):
	def runCallerWorkflow(self, sample, options, workflow = 'RNA-seq'):
		#super().__init__(self, sample, options, "HaplotypeCaller")

		self.output_folder = os.path.join(self.output_folder, workflow)
		self.prefix = os.path.join(self.output_folder, "{0}_vs_{1}.haplotypecaller".format(sample['NormalID'], sample['SampleID']))
		checkdir(self.output_folder, True)

		self.dna_output = self.prefix + ".DNA.raw_snps_indels.vcf"
		self.rna_output = self.prefix + ".RNA.raw_snps_indels.vcf"
		self.rna_filtered_output = self.prefix + ".RNA.final_snps_indels.vcf"

		if workflow == 'DNA-seq':
			output = self.callDNAVariants(sample)
			self.full_output = [self.dna_output]
			self.final_output = self.dna_output
		else:
			bqsr = BaseQualityScoreRecalibration(sample, options)
			raw_variants = self.callRNAVariants(bqsr.bam)
			filtered_variants = self.filterVariants(raw_variants)
			self.full_output = [self.rna_output, self.rna_filtered_output]
			self.final_output = self.rna_filtered_output

	def callRNAVariants(self, bam_file):
		rna_output = self.rna_output

		command = """java -jar {GATK} \
			-T HaplotypeCaller
			-R {reference} \
			-I {sample} \
			--dbsnp {dbSNP} \
			--dontUseSoftClippedBases \
			-o {output}""".format(
				GATK 		= self.gatk_program,
				reference 	= self.reference,
				sample 		= bam_file,
				dbSNP 		= self.dbSNP,
				output 		= rna_output)
		self.runCallerCommand(command, rna_output)
		return rna_output
	
	def filterVariants(self, vcf_file):
		rna_filter_command = """java -jar {GATK} \
			-T VariantFiltration \
			-R {reference} \
			-V {inputfile} \
			-window 35 \
			-cluster 3 \
			--filterName FS --filterExpression \"FS > 30.0\" \
			--filterName QD --filterExpression \"QD < 2.0\" \
			-o {output}""".format(
				GATK = self.gatk_program,
				reference = self.reference,
				inputfile = vcf_file,
				output = self.rna_filtered_output)
		self.runCallerCommand(rna_filter_command, filtered_output)
		return filtered_output
	
	def callDNAVariants(self, sample, options):
		output_file = self.dna_output
		command = """java -jar {GATK} \
			-T HaplotypeCaller
			-R {reference} \
			-I {normal} \
			-I {tumor} \
			-L {targets} \
			-nct {threads} \
			--dbsnp {dbSNP} \
			-o {output}""".format(
				GATK = self.gatk_program,
				reference = self.reference,
				normal = sample['NormalBAM'],
				tumor = sample['TumorBAM'],
				targets = sample['ExomeTargets'],
				dbSNP = self.dbSNP,
				output = output_file,
				threads = self.max_cpu_threads)

		self.runCallerCommand(command, output_file)
		return output_file


class BaseQualityScoreRecalibration(Caller):
	def runCallerWorkflow(self, sample, options):
		#super().__init__(self, sample, options, 'BaseQualityScoreRecalibration')
		self.output_folder = os.path.join("/media/upmc/WD_Partition_2/RNA-seq", 'recalibrated_genomes', sample['PatientID'])

		self.recalibration_table 	= os.path.join(self.output_folder, sample['SampleID'] + ".RNA.recalibration_data.table")
		self.covariate_table 		= os.path.join(self.output_folder, sample['SampleID'] + ".RNA.covariate_data.table")
		self.recalibration_plots 	= os.path.join(self.output_folder, sample['SampleID'] + ".RNA.recalibration_plots.pdf")
		self.cigar_bam 				= os.path.join(self.output_folder, sample['SampleID'] + ".RNA.cigar.bam")
		self.realigned_bam 			= os.path.join(self.output_folder, sample['SampleID'] + ".RNA.recalibrated.bam")
		self.bam = self.realigned_bam
		raw_rna_bam_file = self._getRNABAM(sample)
		cigar_bam = self.splitCigarReads(raw_rna_bam_file)
		recalibration_table = self.generateRecalibrationTable(cigar_bam)
		realigned_bam = self.recalibrateBAM(cigar_bam, recalibration_table)
		covariate_table = self.generateCovariateTable(realigned_bam, recalibration_table)
		recalibration_plots = self.generateRecalibrationPlots(recalibration_table, covariate_table)

		self.final_output = realigned_bam
		self.full_output = [self.recalibration_table, self.covariate_table, self.recalibration_plots, self.cigar_bam, self.realigned_bam]

	def _getRNABAM(self, sample):
		return sample['RNABAM']
	def splitCigarReads(self, rna_bam):
		output_file = self.cigar_bam
		command = """java -jar {GATK} \
			-T SplitNCigarReads \
			-R {reference} \
			-I {inputbam} \
			-o {outputbam} \
			-rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 \
			-U ALLOW_N_CIGAR_READS""".format(
				GATK 		= self.gatk_program,
				reference 	= self.reference,
				inputbam 	= rna_bam,
				outputbam 	= cigar_bam)
		self.runCallerCommand(command, output_file)
		return output_file

	def generateRecalibrationTable(self, bam_file):
		recalibration_table = self.recalibration_table
		command = """java -jar {GATK} \
			-T BaseRecalibrator \
			-R {reference} \
			-I {bam} \
			-knownSites {dbSNP} \
			-o {output}""".format(
				GATK = self.gatk_program,
				reference = self.reference,
				dbSNP = self.dbSNP,
				bam = bam_file,
				output = recalibration_table)
		self.runCallerCommand(command, recalibration_table)
		return recalibration_table

	def recalibrateBAM(self, bam_file, recalibration_table):
		output_file = self.realigned_bam
		command  = """java -jar {GATK} \
			-T PrintReads \
			-R {reference} \
			-I {bam} \
			-BQSR {table} \
			-o {output}""".format(
				GATK = self.gatk_program,
				reference = self.reference,
				bam = bam_file,
				output = realigned_bam,
				table = recalibration_table)
		self.runCallerCommand(command, output_file)
		return output_file
	
	def generateCovariateTable(self, bam_file, recalibration_table):
		covariate_table = self.covariate_table
		command = """java -jar {GATK} \
			-T BaseRecalibrator \
			-R {reference} \
			-I {realigned_bam} \
			-BQSR {table} \
			-knownSites {dbSNP} \
			-o {output}""".format(
				GATK = self.gatk_program,
				reference =self.reference,
				dbSNP = self.dbSNP,
				table = recalibration_table,
				realigned_bam = bam_file,
				output = covariate_table)
		self.runCallerCommand(command, covariate_table)
		return covariate_table
	
	def generateRecalibrationPlots(self, recalibration_table, covariate_table):
		recalibration_plots = self.recalibration_plots
		command = """ java -jar {GATK} \
			-T AnalyzeCovariates \
			-R {reference} \
			-before {before} \
			-after {after} \
			-plots {plots}""".format(
				GATK = self.gatk_program,
				reference = self.reference,
				before = recalibration_table,
				after = covariate_table,
				plots = recalibration_plots)
		self.runCallerCommand(command, recalibration_plots)
		return recalibration_plots

class UnifiedGenotyper(Caller):
	def runCallerWorkflow(self, sample, options, workflow = 'RNA-seq'):
		#super().__init__(self, sample, options, 'UnifiedGenotyper')
		self.output_folder = os.path.join(self.output_folder, workflow)
		self.prefix = os.path.join(self.output_folder, "{0}_vs_{1}.unifiedgenotyper".format(sample['NormalID'], sample['SampleID']))

		if workflow == 'DNA-seq':
			dna_variants = self.cannDNAVariants(sample)
			self.full_output = [dna_variants]
			self.final_output = dna_variants
		else:
			raw_rna_variants = self.callRNAVariants(sample)
			filtered_rna_variants = self.filterVariants(sample)
			self.full_output += [raw_rna_variants, filtered_rna_variants]
			self.final_output = filtered_rna_variants

	def callDNAVariants(self, sample):
		dna_output = self.prefix + '.DNA.raw_snps_indels.vcf'
		dna_command = """java -jar {GATK} \
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
				output = dna_output)
		self.runCallerCommand(dna_command, dna_output)
		return dna_output

	def callRNAVariants(self, bam_file):
		rna_output = self.prefix + '.RNA.raw_snps_indels.vcf'
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
				output = rna_output)
		self.runCallerCommand(command, rna_output)
		return rna_output

	def filterVariants(self, vcf_file):
		filtered_output = os.path.splitext(vcf_file)[0] + '.filtered.vcf'
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
				output = filtered_output)
		self.runCallerCommand(command, filtered_output)
		return filtered_output

#----------------------------------------------------------------------------------------------------
#------------------------------------ Copynumber Tool Functions -------------------------------------
#----------------------------------------------------------------------------------------------------
class VarscanCopynumber(cCaller):
	def runCallerWorkflow(self, sample, options):
		#super().__init__(self, sample, options, 'Varscan')

		normal_pileup 	= self.generatePileup(sample['NormalBAM'], sample['NormalID'])
		tumor_pileup 	= self.generatePileup(sample['TumorBAM'], sample['SampleID'])

		copynumber_ratios = self.callCopynumberRatios(normal_pileup, tumor_pileup)
		copynumber_calls  = self.copyCaller()

		copynumber_segments = self.circularBinarySegmentation()

	def circularBinarySegmentation(input_file, output_file, rscript):
		#NEED to strip first row of table before r script
		script = """library(DNAcopy)
			cn <- read.table("{input_file}",header=TRUE)
			CNA.object <-CNA( genomdat = cn[,6], chrom = cn[,1], maploc = cn[,2], data.type = 'logratio')
			CNA.smoothed <- smooth.CNA(CNA.object)
			segs <- segment(CNA.smoothed, verbose=0, min.width=2)
			segs2 = segs$output
			write.table(segs2[,2:6], file="{output_file}", row.names=F, col.names=F, quote=F, sep="\\t")""".format(
			input_file = input_file,
			output_file = output_file)
		with open(rscript, 'w') as file1:
			for line in script.splitlines():
				file1.write(line + '\n')

		command = "Rscript {script}".format(script = rscript)

		self.runCallerCommand(command, output_file)
		return output_file

	def callCopynumberRatios(self, normal_pileup, tumor_pileup):
		copynumber_output = varscan_prefix + '.copynumber.called'
		#-------------------------- Varscan (copynumber - caller) -----------------------------
		copynumber_command = "java {memory} -jar {varscan} copynumber {normal} {tumor} {prefix} --min-base_qual {mbq} --min-map-qual {mmq}".format(
			mmq = self.min_mapping_quality,
			mbq = self.min_base_quality,
			normal = normal_pileup,
			tumor = tumor_pileup,
			varscan = self.program,
			memory = self.max_memory_usage,
			prefix = self.prefix)
		self.runCallerCommand(copynumber_command, copynumber_output)

		return copynumber_output

	def copyCaller(self, copynumber_ratios):
		copycaller_command = """java {memory} -jar {varscan} copyCaller {prefix}.copynumber \
		--output-file {prefix}.copynumber.called \
		--output-homdel-file {prefix}.copynumber.called.homdel \ 
		--min-coverage {mq}""".format(
			memory = self.max_memory_usage,
			varscan = self.program,
			prefix = self.prefix,
			mq = self.min_coverage)

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

class CNVkit(cCaller):
	def runCallerWorkflow(self, sample, options):
		#super().__init__(self, sample, options, 'CNVkit')
		reference_cnn = self.runBatchCommand(sample)

	def runBatchCommand(self, sample):
		reference_cnn = os.path.join(self.output_folder, "reference_cnv.cnn")

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
		self.runCallerCommand(batch_command, reference_cnn)

		expected_output = [
			os.path.join()
		]
		return reference_cnn

class FREEC(cCaller):
	def runCallerWorkflow(sself, sample, options):
		pass
def _configure_freec(sample, options, reference, chrlenfile, output_folder):
	general_options = {
		'bedtools': "/usr/bin/bedtools",
		'breakPointThreshold': 0.8, #Default: 0.8 use something like 0.6 to get more segments (and thus more predicted CNVs)
		'chrFiles': os.path.join(os.path.dirname(reference), 'chromosomes'),
		'chrLenFile': chrlenfile, #a list of chromosomes and chromosome lengths. Basically the reference dict.
		'maxThreads': '6',
		'noisyData': 'TRUE',
		'outputDir': output_folder,
		'ploidy': '2',
		'printNA': 'FALSE', #set FALSE to avoid printing "-1" to the _ratio.txt files Useful for exome-seq or targeted sequencing data
		'readCountThreshold': '50', #Default: 10, recommended value >=50 for for exome data
		'samtools': options['Programs']['samtools'],
		#'sex': "",
		'window': "0" #for whole exome sequencing: "window=0"
	}
	general_options = ["{0} = {1}".format(k, v) for k, v in general_options.items()]
	
	#[sample]
	sample_options = {
		'mateFile': sample['TumorPileup'],
		'inputFormat': 'pileup',
		'mateOrientation': '0'
	}
	sample_options = ["{0} = {1}".format(k, v) for k, v in sample_options.items()]

	#[control]
	control_options = {
		'mateFile': sample['NormalPileup'],
		'inputFormat': 'pileup',
		'mateOrientation': '0'
	}
	
	control_options = ["{0} = {1}".format(k, v) for k, v in control_options.items()]

	#[BAF]
	baf_options = {
		'SNPfile': options['Reference Files']['dbSNP'],
		"minimalCoveragePerPosition": '5'
	}
	baf_options = ["{0} = {1}".format(k, v) for k, v in baf_options.items()]

	#[target]
	target_options = {
		"captureRegions": sample['ExomeTargets'],
	}
	target_options = ["{0} = {1}".format(k, v) for k, v in target_options.items()]

	all_options = ["[general]"] + general_options
	all_options +=["[sample]"] 	+ sample_options
	all_options +=["[control"] 	+ control_options
	all_options +=["[BAF]"]		+ baf_options
	all_options +=["[target]"] 	+ target_options

	return all_options

def FREEC(sample, options):
	#http://boevalab.com/FREEC/tutorial.html
	print("\tControl-FREEC")
	freec_location = options['Programs']['freec']
	reference = options['Reference Files']['reference genome']
	
	output_folder = os.path.join(options['Pipeline Options']['copynumber pipeline folder'], sample['PatientID'], 'FREEC')
	
	if not os.path.isdir(output_folder):
		os.makedirs(output_folder)
	
	freec_config_file = os.path.join(output_folder, "{0}_vs_{1}_config.txt".format(sample['NormalID'], sample['SampleID']))
	chrlenfile = os.path.join(output_folder, "chromosome_lengths.txt")
	genome_index_file = options['Reference Files']['reference genome'] + '.fai'
	
	script_folder = os.path.dirname(freec_location)
	script_folder = os.path.join(os.path.dirname(script_folder), 'scripts')

	assess_significance_script = os.path.join(script_folder, "assess_significance.R")
	plot_script = os.path.join(script_folder, "makeGraph.R")
	


	print("Output Folder: ", output_folder)
	print("FREEC Config File: ", freec_config_file)
	print("Genome Index File: ", genome_index_file)
	print("Script Folder: ", script_folder)
	print("Significance Script: ", assess_significance_script)
	print("Plot Script: ", plot_script)
	

	#--------------------- Generate Pileup Files -----------------------------
	print("\t\t(1/6) Generate Pileup Files")
	sample, pileup_log = Pileup(sample, options)
	pprint(sample)
	#--------------------- Generate ChrLen File  -----------------------------
	print("\t\t(2/6) Detect Chromosome Lengths")
	#generated as the intersection of the genome index and exome targets files
	exclude_chroms = ['chr1_KI270706v1_random', 'chr4_GL000008v2_random', 'chr14_GL000009v2_random', 'chrUn_KI270742v1']
	
	
	with open(sample['ExomeTargets'], 'r') as bedfile:
		bedchrs = set([i[0] for i in csv.reader(bedfile, delimiter = '\t')])
	
	with open(genome_index_file, 'r') as file1:
		chrlens = ['\t'.join(i[:3]) + '\n' for i in csv.reader(file1, delimiter = '\t') if (i[0] in bedchrs and i[0] not in exclude_chroms)]
	
	with open(chrlenfile, 'w') as file1:
		[file1.write(i) for i in chrlens]

	#--------------------- Generate Config File ------------------------------
	print("\t\t(3/6) Configure FREEC")
	#[general]
	#_configure_freec(sample, options, reference, chrlenfile, output_folder):
	all_options = _configure_freec(sample, options, reference, chrlenfile, output_folder)
	with open(freec_config_file, 'w') as file1:
		for line in all_options:
			file1.write(line + '\n')
	#------------------------------------- Run FREEC --------------------------------
	print("\t\t(4/6) Main Analysis")
	freec_command = "{freec} -conf {config}".format(
		freec = freec_location,
		config = freec_config_file)
	Terminal(freec_command, show_output = True)

	#------------------------------------ Add Log2Ratios ----------------------------
	#--------------------------------- Calculate Significance -----------------------
	print("\t\t(5/6) Calculate Significance")
	
	#files for the normal sample
	normal_CNVs = os.path.join(output_folder, "{0}.pileup_normal_CNVs".format(sample['SampleID']))
	normal_ratios= os.path.join(output_folder, "{0}.pileup_normal_ratio.txt".format(sample['SampleID']))
	normal_baf_file = os.path.join(output_folder, "{0}.pileup_BAF.txt".format(sample['SampleID']))
	
	#files for the tumor sample
	CNVs = os.path.join(output_folder, "{0}.pileup_CNVs".format(sample['SampleID']))
	ratios = os.path.join(output_folder, "{0}.pileup_ratio.txt".format(sample['SampleID']))
	baf_file = os.path.join(output_folder, "{0}.pileup_BAF.txt".format(sample['SampleID']))

	normal_sig_command = "cat {script} | R --slave --args {CNVs} {ratios}".format(
		script = assess_significance_script,
		CNVs = normal_CNVs,
		ratios = normal_ratios)
	other_sig_command = "cat {script} | R --slave --args {CNVs} {ratios}".format(
		script = assess_significance_script,
		CNVs = CNVs,
		ratios = ratios)
	
	Terminal(normal_sig_command, show_output = True)
	Terminal(other_sig_command, show_output = True)
	#--------------------------------- Generate Plots --------------------------------
	print("\t\t(6/6) Generate Plots")
	

	print("Plot Script: ", plot_script)
	
	plot_command = "cat {script} | R --slave --args {ploidy} {ratios} {baf}"
	normal_plot_command = plot_command.format(
		script = plot_script,
		ploidy = "2",
		ratios = normal_ratios,
		baf = normal_baf_file)
	plot_command = plot_command.format(
		script = plot_script,
		ploidy = "2",
		ratios = ratios,
		baf = baf_file)
	
	Terminal(normal_plot_command, show_output = True)
	Terminal(plot_command, show_output = True)

	freec_log = {
		'Inputs': [sample['NormalBAM'], sample['TumorBAM'], sample['ExomeTargets'], chrlenfile],
		'Outputs': [output_folder],
		'Intermediate Files': [],
		'Commands': [freec_command]
	}

	return sample, freec_log

#----------------------------------------------------------------------------------------------------
#------------------------------------------ Pipelines -----------------------------------------------
#----------------------------------------------------------------------------------------------------
class Pipeline:
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
	def __init__(self, sample, options, sample_callers):
		pipeline_callers = self._getCallers(sample_callers)

		for caller_name, caller in pipeline_callers:
			print(caller_name)
			callset = caller(sample, options)

class SomaticPipeline(Pipeline):
	@staticmethod
	def _getCallers(callers):
		available_callers = {
			#'pon': Mutect_pon_detection,
			#'gdc': GDC_somatic,
			'muse': MuSE,
			'mutect': MuTect2,
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
		sample['NormalBAM'] = self._fetchFilename(sample.get('NormalBAM'), sample['NormalUUID'])
		sample['TumorBAM'] = self._fetchFilename(sample.get('TumorBAM'), sample['SampleUUID'])

		if verify_md5sum:
			normal_md5sum = API(sample['NormalUUID'], 'files')
			tumor_md5sum = API(sample['SampleUUID'], 'files')
		else:
			normal_md5sum = None
			tumor_md5sum = None

		normal_file_status = self.verifyFile(sample['NormalBAM'], normal_md5sum)
		tumor_file_status  = self.verifyFile(sample['TumorBAM'], tumor_md5sum)

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
		if self.DEBUG:
			print("\t\tExpected File: ", filename)
			print("\t\tFile Exists: ", file_exists)
			print("\t\tExpected md5sum: ", expected_md5sum)
			print("\t\tFile md5sum: ", file_md5sum)
			print("\t\tFile is Valid: ", file_is_valid)

		return file_is_valid

class GenomicsPipeline:
	def __init__(self, **kwargs):
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
		if parser.debug:
			print("GenomicsPipeline")
			for k, v in kwargs.items(): print('\t', k, '\t', v)


		self.parser = parser
		#Parse the provided arguments
		somatic_callers = [i.lower() for i in kwargs.get("copynumber", [])]
		copynumber_callers = [i.lower() for i in kwargs.get("somatic", [])]
		config_filename = kwargs.get('config_filename', os.path.join(PIPELINE_DIRECTORY, "0_config_files", "pipeline_project_options.txt"))
		caller_status_filename = kwargs.get('caller_status_filename', os.path.join(PIPELINE_DIRECTORY, "0_config_files", "caller_status.tsv"))

		sample_list, config = self._loadPipelineConfiguration(sample_filename, config_filename)		
		
		LOGGER.info("Somatic Callers: " + ', '.join(somatic_callers))
		LOGGER.info("Copynumber Callers: " + ', '.join(copynumber_callers))
		LOGGER.info("Running through the genomics pipeline with {0} samples.".format(len(sample_list)))
		#sample_list = []
		for index, sample in enumerate(sample_list):
			print("({0}/{1}) {2}\t{3}".format(index+1, len(sample_list), sample['PatientID'], now().isoformat()), flush = True)

			sample_status = self.runSample(sample, config, somatic_callers, copynumber_callers, caller_status)



	@staticmethod
	def _loadPipelineConfiguration(sample_filename, config_filename, caller_status_filename):
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

		snv_folder = os.path.join(snv_folder,PatientID)
		cnv_folder = os.path.join(cnv_folder, PatientID)
		temp_folder = os.path.join(temp_folder, PatientID)

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

		if os.path.exists(self.caller_status_filename):
			pass
		else:
			pass
		return []
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

		self._makePatientFolders(sample, config)

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
			if not self.parser.ignore_caller_status or self.parser.debug:
				somatic_callers = [i for i in somatic_callers if i not in sample_completed_callers]
			somatic_pipeline = SomaticPipeline(sample, config, somatic_callers)
			if self.parser.debug:
				print("Somatic Pipeline: ", type(somatic_pipeline))
				print("\t\tSomatic Callers: ", somatic_callers)

			if not self.parser.ignore_caller_status or self.parser.debug:
				copynumber_callers = [i for i in somatic_callers if i not in sample_completed_callers]
			copynumber_pipeline = CopynumberPipeline(sample, config, copynumber_callers)
			if self.parser.debug:
				print("Copynumber Pipeline: ", type(copynumber_pipeline))
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
			'Notes': 	  "Normal Only" if normal_only else "",
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

	parser.add_argument('-d' ,"--debug",
		dest = 'debug',
		action = 'store_true',
		help='debug the pipeline using default settings' + show_default)

	parser.add_argument('-i', "--ignore-caller-status",
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
		somatic_callers = ['MuSE', 'Varscan', 'Strelka', 'SomaticSniper', 'Mutect2']
		copynumber_callers = []
	else:
		sample_filename = CMD_PARSER.sample_list
		somatic_callers = ['MuSE', 'Varscan', 'Strelka', 'Somaticsniper', 'Mutect']
		copynumber_callers = ['varscan', 'cnvkit', 'freec']	

	pipeline = GenomicsPipeline(sample_filename, config_filename, caller_status_filename, somatic_callers, copynumber_callers, parser = CMD_PARSER)
else:
	pass
#3872