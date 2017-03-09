import csv
import logging
import datetime
import os
import shutil
import isodate
import gdc_api
import subprocess
import math
import configparser
import hashlib
import shlex
from pprint import pprint

now = datetime.datetime.now
global_start = now() #Used to log when a series of samples were run together

PIPELINE_DIRECTORY = "/home/upmc/Documents/Variant_Discovery_Pipeline"
#initial_working_directory = PIPELINE_DIRECTORY
file_exists_string = "SKIPPING {program}: THE FILE {f} ALREADY EXISTS!"
console_log_file_filename = os.path.join(PIPELINE_DIRECTORY, "console_log_file.txt")
"""Set up the logger"""

logger = logging.getLogger('genome_pipeline')
hdlr = logging.FileHandler('pipeline_log_file.log')
formatter = logging.Formatter('%(asctime)s %(levelname)s %(message)s')
hdlr.setFormatter(formatter)
logger.addHandler(hdlr) 
logger.setLevel(logging.INFO)

logger.info("#" * 120)

logger.info('#'*30 + 'Starting the genomics pipeline at ' + global_start.isoformat() + '#'*30)

logger.info("#" * 120)
#----------------------------------------------------------------------------------------------------
#-------------------------------------Set Up Global Parameters---------------------------------------
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
"""
class DefaultConfig:
	def __init__(self, filename = None):

		self.working_directory = PIPELINE_DIRECTORY
		if filename:
			print("Loading config file from ", filename)
			self.config = self.load_file(filename)
		else:
			print("Loading default options...")
			self.config = self.generate_default_config(self.working_directory, "/home/upmc/Programs/")
	
		self.validate_options(self.config)
	@staticmethod
	def load_file(filename):
		_, ext = os.path.splitext(filename)
		if ext == '.txt' or ext == '.ini': #python config file
			configuration = configparser.ConfigParser()
			configuration.read(filename)
		elif ext == '.json':
			pass
		else:
			raise FileError("The configuration file is not valid.")
		return configuration

	@staticmethod
	def set_program_outputs(path):
		""" Sets the default output templates for each program.
			Pipelines:
				1: Aligned Reads
					1_aligned_reads
				2: Pre-process the files
					2_processing_pipeline
				3: Call Somatic Variants
					3_variant_calling_pipeline
				4: Call Copynumber Alterations
					4_copynumber_pipeline
				5: Intermediate Files
					5_intermediary_files
			Setup:
				Pipeline_folder
					patient_folder
						variant_caller_folder
						temporary_files_folder
						final_output_folder

		"""
		outputs = {
			#--------------------------------- Preprocessing Programs -----------------------------------
			'RealignerTargetCreator': 	"2_processing_pipeline/{patient}/{sample}.{ref}.targets.intervals",
			'IndelRealigner': 			"2_processing_pipeline/{patient}/{sample}.{ref}.realigned.bam",
			'BaseRecalibrator': 		"2_processing_pipeline/{patient}/{sample}.{ref}.recalibrated_data_table.grp",
			'PrintReads': 				"2_processing_pipeline/{patient}/{sample}.{ref}.recalibrated.bam",
			#----------------------------------- Variant callers ---------------------------------------
			'Bambino':					"3_called_variants/{patient}/Bambino/", 	#[CUSTOM FORMAT]
			'GDC (somatic)':			"3_called_variants/{patient}/GDC/",
			'mpileup':					"3_called_variants/{patient}/mpileup/",
			'HaplotypeCaller':			"3_called_variants/{patient}/HaplotypeCaller/", 	#[VCF]
			'MuSE':						"3_called_variants/{patient}/MuSE/", 			#[VCF]
			'MuTect2':					"3_called_variants/{patient}/MuTect2/",		#[VCF]
			'Seurat':					"3_called_variants/{patient}/Seurat/{normal}_vs_{tumor}.vcf",		#[VCF]
			'SomaticSniper':			"3_called_variants/{patient}/SomaticSniper/",#[VCF]
			'Strelka':					"3_called_variants/{patient}/Strelka/", 							#[VCF]
			'Varscan (somatic)':		"3_called_variants/{patient}/Varscan/", 		#[CUSTOM FORMAT|VCF?]
			'pileup':					"5_temporary_files/{patient}/pileup/{sample}.{ref}.pileup", #manually defined in the varscan (somatic funciton)
			#----------------------------------- Variant Evaluation -------------------------------------
			'VariantEval': 				"F_Evaluated_Variants/{sample}.{ref}.varianteval.vcf",
			#---------------------------------Copynumber Pipeline Tools----------------------------------
			
			'fixfile':					"4_called_cnvs/{patient}/{sample}.pileup.fixed", #Custom function to fix errors caused by mpileup
			'CoNIFER':					"4_called_cnvs/{patient}/CoNIFER/",
			'FREEC': 					"4_called_cnvs/{patient}/FREEC/",
			'GDC (copynumber)':			"4_called_cnvs/{patient}/GDC/",
			'Varscan (copynumber)':		"4_called_cnvs/{patient}/Varscan/",
			'ADTEx':					"4_called_cnvs/{patient}/ADTEx/",
			#-------------------------------------------Other--------------------------------------------
			'Intermediate Folder': 		"5_temporary_files/{patient}/",
			'Panel of Normals (Mutect)':"D_panel_of_normals/{sample}.{ref}.artifact.vcf",
			'Panel of Normals (Combine)':"D_panel_of_normals/panel_of_normals.{ref}.{length}.vcf"
		}

		for key, value in outputs.items():
			outputs[key] = os.path.join(path, value)
		return outputs
	@staticmethod
	def set_program_paths(program_folder):
		#programs = "/home/upmc/Programs/"
		gatk_path = os.path.join(program_folder, "GATK", "GATK3.6.jar")
		program_paths = {
			#------------------------------------- Preprocessing ----------------------------------------
			'RealignerTargetCreator': 	gatk_path, #'GATK' indicates it is run through GATK
			'IndelRealigner': 			gatk_path,
			'BaseRecalibrator': 		gatk_path,
			'PrintReads': 				gatk_path,
			#----------------------------------- Variant Discovery --------------------------------------
			'Bambino':					os.path.join(program_folder, 'bambino_bundle_1.06.jar'),
			'MuSE':						os.path.join(program_folder, "MuSE", "MuSEv1.0rc_submission_c039ffa"),
			'MuTect2':					gatk_path,
			'Seurat':					os.path.join(program_folder, "Seurat-2.5.jar"),
			'SomaticSniper':			os.path.join(program_folder, "somatic-sniper", "build", "bin", "bam-somaticsniper"),
			'Strelka':					os.path.join(program_folder, "Strelka"),
			#----------------------------------- Variant Evaluation -------------------------------------
			'VariantEval': 				gatk_path,
			#-------------------------------- Copynumber Pipeline Tools ---------------------------------
			'ADTEx':					os.path.join(program_folder, 'ADTEx.v.2.0/ADTEx.py'),
			'cnvkit':					os.path.join(program_folder, 'cnvkit', 'cnvkit.py'),
			'CoNIFER':					os.path.join(program_folder, 'conifer_v0.2.2', "conifer.py"),
			'FREEC': 					os.path.join(program_folder, 'FREEC-9.6', 'src', 'freec'),
			'XHMM': 					os.path.join(program_folder, 'xhmm', 'xhmm'),
			#---------------------------------- Other Pipeline Tools ------------------------------------
			'GATK':						os.path.join(program_folder, "GATK", "GATK3.6.jar"),
			'samtools':					os.path.join(program_folder, "samtools-0.1.18", "samtools"),
			'bcftools':					os.path.join(program_folder, "samtools-0.1.18", "bcftools", "bcftools"),
			'fixfile':					"", #Custom function to fix errors caused by mpileup
			'plinkseq':					os.path.join(program_folder, "plinkseq-0.10", "pseq"),
			'Varscan':					os.path.join(program_folder, 'VarScan.v2.3.9.jar'),
			#--------------------------------------- Evaluation -----------------------------------------
			'VariantEffectPredictor':   os.path.join(program_folder,'ensembl-tools-release-86','scripts', 'variant_effect_predictor', 'variant_effect_predictor.pl'),
			#-------------------------------------------Other--------------------------------------------
			'samtools 0.1.6':			os.path.join(program_folder, "samtools-0.1.6", "samtools"),
			'bam-readcount':			os.path.join(program_folder, 'bam-readcount', 'bin', 'bam-readcount'),
			'GDC Download Tool':		os.path.join(program_folder, 'gdc_data_transfer_tool', 'gdc-client')
		}
		return program_paths
	@staticmethod
	def set_pipeline_parameters():
		somatic_quality = 20
		parameters = {
			'MAX_ALT_ALLELES':		   	 8,
			'MAX_CORES': 			   	 6,
			'MIN_MAPPING_QUALITY':    	15, 
			'MIN_BASE_POSITION_FREQUENCY': 0.6,
			'MIN_COVERAGE':            	 8, #originally 4, changed to 8 to match GDC
			'MIN_COVERAGE_NORMAL':		 8, #From GDC, used in Varscan
			'MIN_COVERAGE_TUMOR':		 6,
			'MIN_NUCLEOTIDE_QUALITY': 	10,
			'MIN_CALL_CONF_THRESHOLD':	30, #The minimum phred-scaled confidence threshold at which variants should be called
			'MIN_EMIT_CONF_THRESHOLD':	30,	#The minimum phred-scaled confidence threshold at which variants should be emitted 
			'INITIAL_NORMAL_LOD':		 0.5,#Initial LOD threshold for calling normal variant
			'INITIAL_TUMOR_LOD':		 4.0, #Initial LOD threshold for calling normal variant
			## -- GATK Specific parameters
			'GATK_NUM_THREADS':        	 6,

			## -- Bambino Specific parameter --
			'MIN_FLANKING_QUALITY':   	15,
			'MIN_ALT_ALLELE_COUNT':    	 2,
			'MIN_MINOR_FREQUENCY':     	 0,
			'MMF_MAX_HQ_MISMATCHES':   	 5,
			'MMF_MIN_HQ_THRESHOLD':   	15,
			'MMF_MAX_LQ_MISMATCHES':   	 6,
			'UNIQUE_FILTER_COVERAGE':  	 2,

			## -- Mpileup specific parameter --
			'BWA_DOWNGRADE_COFF':     	50,
			'NO_OF_READS_TO_CONSIDER_REALIGNMENT': 3,
			'FREQ_OF_READS':           	0.0002,
			'MPILEUP_QUALITY_THRESHOLD':10,

			## -- Caveman specific parameter
			'NO_OF_BASES_TO_INCREMENT': 250000,


			## Somatic Sniper specific parameter ##
			'SOMATIC_QUALITY':        somatic_quality,

			## Mutect2 specific parameters
			'CONTAMINATION_FRACTION': 0.0,

			## Java specific parameters
			'JAVA_MAX_MEMORY_USAGE': "-Xmx{0}m".format(12288), #12GB
			'JAVA_GARBAGE_COLLECTION':'-XX:+UseSerialGC'
			}
		return parameters

	def validate_options(self, options):
		""" Confirms that the programs and files are exitst where their respective paths point to """
		#-------------------------------- Validate Program Locations ------------------------------------

		print("--"*20 + "Valid Files" + "--"*20)
		pprint(options.keys())
		for name, path in sorted(options['Reference Files'].items()):
			isvalid = os.path.isfile(path)
			print("{name:<25} {valid}".format(name = name, valid = isvalid))


	def generate_default_config(self, working_directory, program_folder):
		""" An updated version of the default config file
		"""
		#cwd = "/home/upmc/Documents/Variant_Discovery_Pipeline/"

		known_indels = [
			"/home/upmc/Documents/Reference/Known_Indels/resources_broad_hg38_v0_Homo_sapiens_assembly38.known_indels.vcf",
			"/home/upmc/Documents/Reference/Known_Indels/resources_broad_hg38_v0_Mills_and_1000G_gold_standard.indels.hg38.vcf"
		]

		parameters = {
			#The directory to save all intermediate files in
			'working directory': 	working_directory,
			'reference directory': 	"/home/upmc/Documents/Reference/",	#The directory with all relevant reference files.
			'program directory': 	"/home/upmc/Programs/", 			#The folder where the relevant programs are saved.
			'log file':				os.path.join(working_directory, "pipeline_log.tsv"),
			'sample list':			"",									#A tsv file of the samples to use.
			#'program options':		set_program_options(),				#optional commands for each program
			'output':				self.set_program_outputs(path = working_directory),	#Template filenames for each tool output
			'programs':				self.set_program_paths(program_folder),
			'reference': {
				'GRCh37': "/home/upmc/Documents/Reference/GRCh37/Homo_sapiens_assembly19.fasta",
				'GRCh38': "/home/upmc/Documents/Reference/GRCh38/GRCh38_full_analysis_set_plus_decoy_hla.fa",
				'GRCh38.d1.vd1': "/home/upmc/Documents/Reference/GRCh38/GRCh38.d1.vd1.fa", 
				'dbSNP':  "/home/upmc/Documents/Reference/dbSNP/dbsnp_144.hg38.vcf.gz",
				#'dbSNP': "/home/upmc/Documents/Reference/GRCh38_Resources/dbsnp_144.hg38.vcf.gz",
				'hapmap':"/home/upmc/Documents/Reference/GRCh38_Resources/hapmap_3.3.hg38.vcf",
				'omni':  "/home/upmc/Documents/Reference/GRCh38_Resources/1000G_omni2.5.hg38.vcf",
				'1000G': "/home/upmc/Documents/Reference/GRCh38_Resources/1000G_phase1.snps.high_confidence.hg38.vcf",
				'cosmic':"/home/upmc/Documents/Reference/Cosmic/CosmicCodingMuts_v78.vcf",
				'known indels': '|'.join(known_indels),
				'user token': "/home/upmc/Programs/gdc_data_transfer_tool/gdc-user-token.2016-10-03T14_49_48-04_00.txt"
				},
			'parameters': self.set_pipeline_parameters(),
			'GATK': {
				'-logging_level': 'ERROR',
				'-log_to_file': None,
				'-num_threads': 2, #Number of data threads to allocate to this analysis
				'-performanceLog': None #Write GATK runtime performance log to this file
			}
		}

		return parameters
	def to_json(self):
		return self.config
	def save(self, ext = '.json'):
		import json
		filename = "pipeline_options.json"
		with open(filename, 'w') as file1:
			file1.write(json.dumps(self.config, indent = 4, sort_keys = True))
"""
def checkdir(path):
	if not os.path.isdir(path): os.makedirs(path)
#----------------------------------------------------------------------------------------------------
#---------------------------------------- Logging functions -----------------------------------------
#----------------------------------------------------------------------------------------------------
def csv_log(rows, filename = None):
	"""row should be a list of dictionaries"""
	if not isinstance(rows, list):
		rows = [rows]
	if len(rows) == 0: return None
	if filename is None:
		filename = "pipeline_log.tsv"

	fieldnames = sorted((rows[0].keys())) + ['Global Start']

	with open(filename, 'a') as file1:
		writer = csv.DictWriter(file1, fieldnames = fieldnames, delimiter = '\t')

		if os.path.getsize(filename) < 1:
			writer.writeheader()
		for row in rows:
			for field in fieldnames:
				value = row.get(field, "")
				if isinstance(value, (list, set)): value = '|'.join(value)
				if isinstance(value, str): value = value.replace(PIPELINE_DIRECTORY, 'cwd/')

			row['Global Start'] = global_start
			writer.writerow(row)

	print("Wrote {0} rows.".format(len(rows)))

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


def print_update(label, level):
	indent = '\t' * level
	print(indent + label, flush = True)

def Terminal(command, label = None, show_output = False, timeout = None):
	""" Calls the system shell """
	terminal_log = "terminal_log.txt"
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
		#from subprocess import STDOUT, check_output
		#output = subprocess.check_output(command, stderr=subprocess.STDOUT, timeout=timeout)
		command = shlex.split(command)
		#command += ['>', 'console_output.txt']
		#print(command)
		try:
			with open(console_log_file_filename, 'a') as console_log_file:
				process = subprocess.call(command, timeout=timeout, shell=False, stdout = console_log_file)
		except Exception as exception:
			print("EXCEPTION: ", exception)
			logging.error(exception)
			process = None

	return process
#----------------------------------------------------------------------------------------------------
#------------------------------------ Other Pipeline Functions --------------------------------------
#----------------------------------------------------------------------------------------------------
def cleanup(sample, options):
	"""	Removes files generated by the pipelines that are no longer needed.
		These files are stored in the 'Remove' field of the sample
		Remove
			Intermediate Files
			PileupFiles
			BAM files
	"""
	print("\tCleanup")
	pprint(sample)

	for filename in sample.get('Remove', []):
		try:
			os.remove(filename)
		except Exception as failure:
			print("\t\tCould not remove {0} ({1})".format(filename, failure))
	

	for filename in sample.get('Intermediate Files', []):
		try:
			os.remove(filename)
		except Exception as failure:
			print("\t\tCould not remove {0} ({1})".format(filename, failure))

	try:
		pass
		#shutil.rmtree(os.path.dirname(sample['NormalBAM']))
	except Exception as failure:
		print("\t\tCould not remove the BAM directory: ", failure)
	try:
		shutil.rmtree(os.path.join(options['Pipeline Options']['temporary_files'], sample['PatientID']))
	except Exception as failure:
		print("\t\tCould not remove the intermediary files: ", failure)

	cleanup_log = {
		'Inputs': sample['Remove'],
		'Outputs': [],
		'Commands': ""
	}

	return sample, cleanup_log

def GDC_DNAseq(sample, options, file_types, output_dir):
	print("\tGDC ({0})".format(file_types[0]))
	patient_id = sample['CaseID']
	destination = os.path.join(output_dir, sample['PatientID'], 'GDC')
	if not os.path.exists(destination):
		os.makedirs(destination)
	#destination = output_dir.format(patient = sample['PatientID'])
	#destination = "3_called_variants/{patient}/GDC/".format(patient = sample['PatientID'])

	case = gdc_api.case_api(patient_id)
	cnvs = [f for f in case['files'] if f['data_category'] in file_types]
	cnvs = [gdc_api.file_api(f['file_id']) for f in cnvs]
	

	gdc_output = list()
	notes = list()
	submitter_ids = dict()
	for index, f in enumerate(cnvs):
		print("\t\t({0}/{1}) Downloading {2}".format(index+1, len(cnvs), f['file_id']))

		file_id = f['file_id']
		destination_file = os.path.join(destination, f['basic_info']['file_name'])
		submitter_ids[destination_file] = f['submitter_id']

		source = os.path.join(PIPELINE_DIRECTORY, file_id, f['file_name'])
		tries = 5
		while tries > 0 and not os.path.isfile(destination_file):
			tries -= 1
			gdc_api.download_file(
				file_id = file_id,
				destination = destination)
			if os.path.isfile(destination_file):
				gdc_output.append(destination)
				break
		else:
			notes.append("The file {0} failed to download.".format(
				f['basic_info']['file_name']))

	#pprint(submitter_ids)
	#----------------- Unpack the compressed files and rename them -------------------------
	for fn in os.listdir(destination):
		abspath = os.path.join(destination, fn)
		#print(abspath)

		dirname, basename = os.path.split(abspath)   #path, filename.vcf.gz
		file_name, ext = os.path.splitext(basename)	 #filename.vcf, .gz
		
		if ext == '.gz':
			_, new_ext = os.path.splitext(file_name)  #filename, .vcf
			file_name, _ = os.path.splitext(basename) #filname.vcf, .gz
			unzip_command = "gunzip {filename}".format(filename = abspath)
			#print(unzip_command)
			Terminal(unzip_command, show_output = True)
		else:
			file_name += ext
			new_ext = ext

		
		if abspath in submitter_ids.keys():
			os.rename(os.path.join(dirname, file_name), 
					  os.path.join(dirname, submitter_ids[abspath] + new_ext))

	
	gdc_log = {
		'Inputs': [],
		'Outputs': [output_dir],
		'Commands': ['gdc_api.download_file(destination = {0})'.format(file_id, destination)],
		'Notes': '|'.join(notes)
	}

	return sample, gdc_log

def Pileup(sample, options, single = False):
	""" Generates pileup files. If single is True, only
		one file will be generated.
	"""
	patientID = sample['PatientID']
	logging.info("*"*20 + 'Generating Pileup Files' + '*'*20)
	
	reference = options['Reference Files']['reference genome']
	samtools = options['Programs']['samtools']
	pileup_command = "{samtools} mpileup -q 1 -B -f {reference} {sample} > {output}"
	output_folder = os.path.join(
		options['Pipeline Options']['temporary folder'],
		sample['PatientID'], 'Pileup')
	
	if not os.path.isdir(output_folder):
		os.makedirs(output_folder)
	#------------------------ Normal Pileup ----------------------------
	"""
	normal_pileup_output = options['output']['pileup'].format(
		patient = sample['PatientID'],
		sample = sample['NormalID'],
		ref = sample['Reference'])"""
	normal_pileup_output = os.path.join(
		output_folder, sample['NormalID'] + '.pileup')
	tumor_pileup_output = os.path.join(
		output_folder, sample['SampleID'] + '.pileup')

	normal_pileup_command = pileup_command.format(
		samtools = samtools,
		reference = reference,
		sample = sample['NormalBAM'],
		output = normal_pileup_output)
	tumor_pileup_command = pileup_command.format(
		samtools = samtools,
		reference = reference,
		sample = sample['TumorBAM'],
		output = tumor_pileup_output)
	

	sample['NormalPileup'] = normal_pileup_output
	sample['TumorPileup'] = tumor_pileup_output
	#sample['Remove'] += [normal_pileup_output, tumor_pileup_output]


	if not os.path.isfile(sample['NormalPileup']):
		print("Generating the Normal Pileup File...")
		logging.info("{0}: (Pileup 1 of 2 - Normal Pileup): {1}".format(patientID, normal_pileup_command))
		Terminal(normal_pileup_command, label = 'Pileup (normal)',
			show_output = True)
		
	else:
		print("\t\tThe Normal Pileup file already exists")
		logging.info("{0}: (Pileup 1 of 2 - Normal Pileup): The file already exists at {1}".format(patientID, normal_pileup_output))

	if not os.path.isfile(sample['TumorPileup']):
		print("Generating the Tumor Pileup File...")
		logging.info("{0}: (Pileup 2 of 2 - Tumor Pileup): {1}".format(patientID, tumor_pileup_command))
		Terminal(tumor_pileup_command, label = 'Pileup (tumor)',
			show_output = True)
		
	else:
		print("The Tumor Pileup already exists!")
		logging.info("{0}: (Pileup 2 of 2 - Normal Pileup): The file already exists at {1}".format(patientID, tumor_pileup_output))

	logging.info("{0} (Normal Pileup Check): {1}".format(patientID, os.path.isfile(normal_pileup_output)))
	logging.info("{0} (Tumor Pileup Check): {1}".format(patientID, os.path.isfile(tumor_pileup_output)))

	pileup_log = {
		'Inputs': [sample['NormalBAM'], sample['TumorBAM']],
		'Outputs':[output_folder],
		'Commands':[normal_pileup_command, tumor_pileup_command]
	}

	return sample, pileup_log

#----------------------------------------------------------------------------------------------------
#----------------------------------------- Common Functions -----------------------------------------
#----------------------------------------------------------------------------------------------------
def run_program(function, sample_var, option_var, DEBUG = True, logs = None):
	""" Runs each caller
		Parameters
		----------
			function: function
				The function that runs the specified caller
			sample_var: dict<>
				A dictionary containing the relevant variables
				for the sample
			option_var: dict<>
				A dictionary holding relevant variables
			DEBUG: bool; default True
				If False, the output in the terminal is
				suppressed and any exceptions thrown by the function 
				are caught.
		Returns
		----------
			sample: dict<>
				A dictionary containing the sample variables. Some
				functions add or modify these.
			function_log: dict<>
				A log of some benchmarking variable for this run.

	"""
	notes = ""
	program_start = now()
	sample = sample_var.copy()
	function_name = function.__name__
	#function_log = dict()
	if DEBUG:
		#print("Debugging", function.__name__)
		sample, function_log = function(sample_var, option_var)
	else:
		try:
			#print("trying to run...")
			sample, function_log = function(sample = sample_var, options = option_var)
			#print('Debugged successfully')
		except Exception as failure:
			print("\t\t\tException: ", failure)
			#print("Have an exception: ", fail)
			sample = sample_var
			function_log = {
				'Inputs': [],
				'Outputs': [],
				'Commands': [],
				'Notes': failure,
				'Status': 'Failed'
			}		

	program_stop = now()
	inputs   = function_log.get('Inputs', [])
	intermediate_files = function_log.get('Intermediate Files', [])
	outputs  = function_log.get('Outputs', [])
	notes    = function_log.get('Notes', "")
	status   = function_log.get('Status', 'Passed')
	commands = function_log.get('Commands', [])

	if not isinstance(inputs, list):     inputs = []
	if not isinstance(outputs, list):   outputs = []
	if not isinstance(commands, list): commands = []
	inputs = [i for i in inputs   if isinstance(i, str)]
	outputs= [i for i in outputs  if isinstance(i, str)]
	command= [i.replace('\t\t', "") for i in commands if isinstance(i, str)]

	try:
		input_size = sum([getsize(fn) for fn in inputs])
	except Exception as Ex:
		print("\t\t\t", Ex)
		input_size = 0

	try:
		output_size = sum([getsize(fn) for fn in outputs])
	except Exception as Ex:
		print("\t\t\t", Ex)
		output_size = 0

	function_log = {
		'PatientID':  sample_var['PatientID'],
		'Program':    function.__name__,
		'Start Time': program_start.isoformat(),
		'Stop Time':  program_stop.isoformat(),
		'Duration':   program_stop - program_start,
		'Inputs':     inputs,#'|'.join(inputs),
		'Input Size': input_size,
		'Intermediate Files': intermediate_files,
		'Outputs':    outputs,#'|'.join(outputs),
		'Output Size':output_size,
		'Notes': 	  notes,
		'Status':	  status,
		'Commands':   commands#'|'.join(commands)
	}

	#csv_log([function_log])

	return sample, function_log

def verify_output(output):
	print("Verifying the existence of:", output)
	exists = [not os.path.exists(f) for f in output]
	status = any(exists)
	return status

def run_commands(commands, start = 0, patientID = "", command_label = ""):
	logging.info("{0}: Running {1} commands for {2}".format(patientID, len(commands), command_label))
	for index, process in enumerate(commands, start = start):
		label, command, output_file = process\
		
		print("\t\t({0}/{1}) {2}".format(index+1, len(commands), label), flush = True)
		print("\t\t\tExpected Output: ", output_file)
		
		if not os.path.isfile(output_file):
			logging.info("{0}: Running {1} of {2} commands ({3}): {4}".format(patientID, index+1, len(commands), label, command))
			Terminal(command, show_output = False)
		else:
			logging.info("{0}: Skipping {1} of {2} commands ({3}) (the output already exists) : {4}".format(patientID, index+1, len(commands), label, output_file))
			print("\t\t\t{0} already exists.".format(output_file))

def get_program_parameters(sample, options, program):
	output_folder = os.path.join(options['Programs'][program],
		sample['PatientID'])
		
	if '(' in program:
		program = program.split('(')[0].strip()
	program = options['Programs'][program]

	reference = options['Reference Files']['reference genome']

	return program, output_folder, reference
#----------------------------------------------------------------------------------------------------
#------------------------------------- Variant Caller Pipeline --------------------------------------
#----------------------------------------------------------------------------------------------------
class Workflow:
	def __init__(self, sample, pipeline_type):
		self._make_folders(sample, "3_called_variants")
	def _default_options(self):
		pass
	def _make_folders(self):
		pipeline_folder = options['Pipeline Options'][pipeline_type]
		folder = os.path.join(pipeline_folder, sample['PatientID'], 'Strelka')

def GDC_somatic(sample, options):
	output_dir = options['Pipeline Options']['somatic pipeline folder']
	print("GDC_somatic output dir: ", output_dir)
	return GDC_DNAseq(sample, options, ['Simple Nucleotide Variation'], output_dir)

def MuSE(sample, options):
	"""
		Output
		---------
			"3_called_variants/{patient}/MuTect2/{sample}.{ref}.vcf" [VCF]
		Optional Inputs
		---------------
			dbSNP
		Notes
		----------
			The dbSNP file must be compressed ('gz')
	"""
	logging.info("{0}: {1} Running MuSE {1}".format(sample['PatientID'], '*'*30))
	
	print("\tMuSE", flush = True)
	output_folder = os.path.join(
		options['Pipeline Options']['somatic pipeline folder'],
		sample['PatientID'], 'MuSE')

	output_prefix = os.path.join(output_folder, 
		"{normal}_vs_{tumor}".format(
			tumor   = sample['SampleID'], 
			normal  = sample['NormalID']))
	if not os.path.exists(output_folder): os.mkdir(output_folder)
	muse_call_output = output_prefix# + ".MuSE.txt"


	muse_call_output = output_prefix + '.MuSE.txt'
	muse_sump_output = output_prefix + ".MuSe.vcf"
	reference = options['Reference Files']['reference genome']
	location = options['Programs']['muse']
	call_command = "{program} call -O {prefix} -f {reference} {tumor} {normal}"
	sump_command = "{program} sump -I {prefix}.MuSE.txt -E -D {dbSNP} -O {prefix}.Muse.vcf"
	output_file = os.path.join(output_folder, output_prefix + '.vcf')


	#MuSE call
	print("\t\t(1/2) Calling Variants...")
	print("\t\t\tExpected Output: ", output_prefix + '.MuSE.txt')
	if not os.path.isfile(output_prefix + '.MuSE.txt'):
		call_command = call_command.format(	
			program 	= location,
			prefix 		= output_prefix,
			reference 	= reference,
			tumor 		= sample['TumorBAM'],
			normal 		= sample['NormalBAM'])

		logger.info("{0} (MuSE 1 of 2 - Call Command): {1}".format(sample['PatientID'], call_command))
		Terminal(call_command,
			label = 'MuSE (call)', show_output = True)
	else:
		logger.info("{0} (MuSE 1 of 2 - Call Command): The file already exists as ".format(sample['PatientID'], output_prefix + '.MuSE.txt'))
		print("\t\tThe call file already exists!")

	#-------------------------- Muse sump ----------------------------
	print("\t\t(2/2) Filtering Variants", flush = True)
	print("\t\t\tExpected Output: ", muse_sump_output)
	if not os.path.isfile(muse_sump_output):
		sump_command = sump_command.format(	
			program = location,
			prefix = output_prefix,
			dbSNP = options['Reference Files']['dbSNP'])
		logger.info("{0} (MuSE 2 of 2 - Sump Command): {1}".format(sample['PatientID'], sump_command))
		Terminal(sump_command,
			label = 'MuSE (sump)', show_output = True)
	else:
		logger.info("{0} (MuSE 2 of 2 - Sump Command): The file already exists as ".format(sample['PatientID'], muse_sump_output))
		print("\t\tThe sump vcf file already exists!")

	muse_log = {
		'Inputs': [sample['NormalBAM'], sample['TumorBAM'], reference, options['Reference Files']['dbSNP']],
		'Intermediate Files': [muse_call_output],
		'Outputs': [output_folder],
		'Commands': [call_command, sump_command],
		'Status': verify_output([output_file])
	}

	return sample, muse_log

def Mutect2(sample, options):
	""" Runs Mutect2 to call snp variants
		Outputs
		----------
			"3_called_variants/{patient}/MuTect2/{normal}_vs_{tumor}.vcf" [VCF]
		Optional Inputs
		---------------
			targets
			dbSNP
			cosmic
		Options
		----------
			-num_threads
			-max_base_quality_score
			-max_alternate_alleles
		GDC Command
		-----------
		MuSE call \
		-f <reference> \
		-r <region> \                                   
		<tumor.bam> \
		<normal.bam> \
		-O <intermediate_muse_call.txt>

		MuSE sump \
		-I <intermediate_muse_call.txt> \                           
		-E \                                
		-D <dbsnp_known_snp_sites.vcf> \
		-O <muse_variants.vcf>  
		Notes
		-----
			If you get a contig error, check to make sure the chromosomes start with 'chr'

	"""
	print("\tMuTect2", flush = True)
	patientID = sample['PatientID']
	logging.info("{0}: {1} Running MuTect2 {1}".format(patientID, '*'*30))

	#--------------------------------------------- Set Mutect2 Options --------------------------------------------
	default_options = {
		'-dbsnp': options['Reference Files']['dbSNP'],
		'-cosmic': options['Reference Files']['cosmic'],
		#'L': "",
		'min_base_quality_score': options['Parameters']['MIN_MAPPING_QUALITY']
	}
	default_options = format_options(default_options)
	

	reference = options['Reference Files']['reference genome']
	gatk_location = options['Programs']['GATK']

	output_folder = os.path.join(
		options['Pipeline Options']['somatic pipeline folder'],
		sample['PatientID'], 'MuTect2')
	mutect_output = os.path.join(
		output_folder,
		"{normal}_vs_{tumor}.mutect2.vcf".format(
			tumor  = sample['SampleID'], 
			normal = sample['NormalID']))

	if not os.path.exists(output_folder): os.makedirs(output_folder)


	#------------------------------------------ Call MuTect2 -------------------------------------------------
	print("\t\t(1/1) Calling Variants", flush = True)
	#Modified_command
	mutect2_command = """java {memory} -jar {GATK} \
		-T MuTect2 \
		-R {reference} \
		-L {targets} \
		-I:normal {normal} \
		-I:tumor {tumor} \
		--dbsnp {dbSNP} \
		--contamination_fraction_to_filter {cff} \
		--standard_min_confidence_threshold_for_calling {smctc} \
		--standard_min_confidence_threshold_for_emitting {smcte} \
		--min_base_quality_score {mbq} \
		--max_alternate_alleles {maa} \
		--initial_normal_lod {inl} \
		--initial_tumor_lod {itl} \
		--output_mode EMIT_VARIANTS_ONLY \
		--out {output}""".format(
		GATK        = gatk_location, 
		memory      = options['Parameters']['JAVA_MAX_MEMORY_USAGE'],
		NUM_THREADS = options['Parameters']['MAX_CORES'],
		dbSNP       = options['Reference Files']['dbSNP'],
		cff 		= options['Parameters']['CONTAMINATION_FRACTION'],
		mbq         = options['Parameters']['MIN_NUCLEOTIDE_QUALITY'],
		maa         = options['Parameters']['MAX_ALT_ALLELES'],
		smctc       = options['Parameters']['MIN_EMIT_CONF_THRESHOLD'],
		smcte       = options['Parameters']['MIN_CALL_CONF_THRESHOLD'],
		inl         = options['Parameters']['INITIAL_NORMAL_LOD'],
		itl         = options['Parameters']['INITIAL_TUMOR_LOD'],
		reference   = reference, 
		normal      = sample['NormalBAM'], 
		tumor       = sample['TumorBAM'], 
		targets     = sample['ExomeTargets'],
		output      = mutect_output, 
		options     = default_options)

	mutect2_command = """java {memory} -jar {GATK} \
		-T MuTect2 \
		-R {reference} \
		-L {targets} \
		-I:normal {normal} \
		-I:tumor {tumor} \
		--dbsnp {dbSNP} \
		--cosmic {cosmic} \
		--out {output}""".format(
		GATK        = gatk_location, 
		memory      = options['Parameters']['JAVA_MAX_MEMORY_USAGE'],
		NUM_THREADS = options['Parameters']['MAX_CORES'],
		dbSNP       = options['Reference Files']['dbSNP'],
		cosmic 		= options['Reference Files']['cosmic'],
		cff 		= options['Parameters']['CONTAMINATION_FRACTION'],
		mbq         = options['Parameters']['MIN_NUCLEOTIDE_QUALITY'],
		maa         = options['Parameters']['MAX_ALT_ALLELES'],
		smctc       = options['Parameters']['MIN_EMIT_CONF_THRESHOLD'],
		smcte       = options['Parameters']['MIN_CALL_CONF_THRESHOLD'],
		inl         = options['Parameters']['INITIAL_NORMAL_LOD'],
		itl         = options['Parameters']['INITIAL_TUMOR_LOD'],
		reference   = reference, 
		normal      = sample['NormalBAM'], 
		tumor       = sample['TumorBAM'], 
		targets     = sample['ExomeTargets'],
		output      = mutect_output, 
		options     = default_options)

	if not os.path.isfile(mutect_output):
		logging.info("{0} (MuTect2 1 of 1 - variant calling): {1}".format(patientID, mutect2_command))
		Terminal(mutect2_command, label = 'Mutect2', timeout = None)
		
	else:
		logging.info("{0} (MuTect2 1 of 1 - variant calling): The output already exists as {1}".format(patientID, mutect_output))
	status = verify_output([mutect_output])
	if not status:
		print("\t\tOUTPUT ERROR: Could not locate one or more of: ")
		print("\t\t", mutect_output)
	mutect2_log = {
		'Inputs': [sample['NormalBAM'], sample['TumorBAM'], reference],
		'Outputs': [output_folder],
		'Commands': [mutect2_command],
		'Status': status
	}

	return sample, mutect2_log

def get_somaticsniper_commands(sample, options):
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

	#--------------------------------- Constants --------------------------------------------
	program = os.path.join(
		options['Programs']['somaticsniper'],
		'build', 'bin', 'bam-somaticsniper')
	scriptdir = os.path.join(
		options['Programs']['somaticsniper'],
		'src', 'scripts')
	samtools = options['Programs']['samtools']
	output_folder = os.path.join(
		options['Pipeline Options']['somatic pipeline folder'],
		sample['PatientID'],
		'SomaticSniper')
	
	intermediate_folder = os.path.join(options['Pipeline Options']['temporary folder'], sample['PatientID'], 'SomaticSniper')
	if not os.path.isdir(intermediate_folder): os.mkdir(intermediate_folder)
	
	reference = options['Reference Files']['reference genome']
	output_file = os.path.join(output_folder, "{normal}_vs_{tumor}.somaticsniper.vcf".format(
		tumor = sample['SampleID'], 
		normal = sample['NormalID']))

	#Scripts
	#----------------------------------- Script Paths ---------------------------------------
	
	#samtools_script = os.path.join(scriptdir, '')
	
	snpfilter_script= os.path.join(scriptdir, 'snpfilter.pl')
	readcount =       os.path.join(scriptdir, 'prepare_for_readcount.pl')
	hc_script =       os.path.join(scriptdir, 'highconfidence.pl')
	fpfilter  =       os.path.join(scriptdir, 'fpfilter.pl')

	#-------------------------------- Variant Discovery Command -----------------------------
	somaticsniper_command = "{program} {options} -f {reference} {tumor} {normal} {outputfile}"
	somaticsniper_command = somaticsniper_command.format(
		program 	= program,
		options 	= default_options,
		reference 	= reference,
		tumor 		= sample['TumorBAM'],
		normal 		= sample['NormalBAM'],
		outputfile 	= output_file)

	#----------------------------------- Generate Pileup Strings -----------------------------------
	
	normal_output = os.path.join(intermediate_folder, "somaticsniper.normal_indel")
	tumor_output  = os.path.join(intermediate_folder, "somaticsniper.tumor_indel")
	
	base_samtools_command = "{samtools} pileup -cvi -f {reference} {tumorbam} > {pileup}"
	base_filter_command = """{samtools}.pl varFilter {inputpileup} | awk  '$6>={basequality}' | grep -P "\\t\\*\\t" > {outputpileup}.pileup"""

	normal_samtools_command = base_samtools_command.format(
		samtools = samtools,
		reference = reference,
		tumorbam = sample['NormalBAM'],
		pileup = normal_output)

	tumor_samtools_command = base_samtools_command.format(
		samtools = samtools,
		reference = reference,
		tumorbam = sample['TumorBAM'],
		pileup = tumor_output)

	normal_filter_command = base_filter_command.format(
		samtools = samtools,
		basequality = options['Parameters']['MIN_NUCLEOTIDE_QUALITY'],
		inputpileup = normal_output,
		outputpileup= normal_output)

	tumor_filter_command = base_filter_command.format(
		samtools = samtools,
		basequality = options['Parameters']['MIN_NUCLEOTIDE_QUALITY'],
		inputpileup = tumor_output,
		outputpileup= tumor_output)

	#-------------------------- Filter and remove LOH --------------------------------
	
	normal_loh_filter_command = "perl {snpfilter} --snp-file {inputvcf} --indel-file {inputpileup} --out-file {inputvcf}.SNPfilter.intermediate".format(
		snpfilter = snpfilter_script,
		inputvcf = output_file,
		inputpileup = normal_output,
		normal = sample['NormalID'])
	tumor_loh_filter_command = "perl {snpfilter} --snp-file {inputvcf}.SNPfilter.intermediate --indel-file {inputpileup} --out-file {output}.SNPfilter.final".format(
		snpfilter = snpfilter_script,
		inputvcf = output_file,
		inputpileup = tumor_output,
		normal = sample['NormalID'],
		output = output_file)
	filter_intermediate_output = output_file + ".SNPfilter.intermediate"
	filter_output = output_file + '.SNPfilter.final'

	#---------------------------- Readcount Filter --------------------------------
	#Will print warning about undefined value due to the first ~17 lines being a single-column of comments.
	prepare_readcount_output = filter_output + '.pos'
	readcount_output = filter_output + 'readcounts.rc'
	pr_command = "perl {script} --snp-file {inputfile} --out-file {output}".format(
		script = readcount,
		inputfile = filter_output,
		output = prepare_readcount_output)
	#Will print "Couldn't find single-end mapping quality. Check to see if the SM tag is in BAM."
	#This doesn't invalidate results, but try not to use single-end mapping quality in output
	readcount_command = "{program} -b {mbq}  -q 1 -f {reference} -l {proutput} {tumor} > {output}".format(
		program = options['Programs']['bam-readcount'],
		mbq = options['Parameters']['MIN_NUCLEOTIDE_QUALITY'],
		reference = reference,
		proutput = prepare_readcount_output,
		tumor = sample['TumorBAM'],
		output = readcount_output)

	#------------------------------ Remove False Positives -------------------------
	#print("\t\t(4/4) Removing False Positives")
	basename, _ = os.path.splitext(os.path.basename(output_file))
	final_output = os.path.join(output_folder, basename)
	false_positive_output = filter_output + '.fp_pass'
	
	fp_command = "perl {fpfilter} --snp-file {snpfilter} -readcount-file {readcounts}".format(
		fpfilter = fpfilter,
		snpfilter = filter_output,
		readcounts = readcount_output)
	hc_command = "perl {script} --snp-file {snpfilter}.fp_pass --min-mapping-quality {mmq} --min-somatic-score {ss} --lq-output {output}.lq.vcf --out-file {output}.hq.vcf".format(
		script = hc_script,
		mmq = options['Parameters']['MIN_MAPPING_QUALITY'],
		ss = options['Parameters']['SOMATIC_QUALITY'],
		snpfilter = filter_output,
		output = final_output)

	commands = [
		('Variant Discovery', somaticsniper_command, output_file),
		#Manage Pileup Files
		('Normal Pileup File', normal_samtools_command, normal_output),
		('Tumor Pileup File', tumor_samtools_command, tumor_output),
		('Normal Pileup Filter', normal_filter_command, normal_output + '.pileup'),
		('Tumor Pileup Filter', tumor_filter_command, tumor_output + '.pileup'),
		#Filter LOH
		('Filter LOH from Normal', normal_loh_filter_command,"{0}_intermediate.SNPfilter".format(output_file)),
		('Filter LOH from Tumor', tumor_loh_filter_command,"{0}_final.SNPfilter".format(output_file)),
		#Readcounts
		('Prepare Readcounts', pr_command, prepare_readcount_output),
		('Calculate Readcounts', readcount_command, readcount_output),
		#Remove False Positives
		('Remove False Positives', fp_command, ""),
		('Generate HC Files', hc_command, final_output + '.hq.vcf')
	]

	return commands

def SomaticSniper(sample, options):
	print("\tSomaticSniper")
	patientID = sample['PatientID']

	logging.info("{0}: {1} Running SomaticSniper {1}".format(patientID, '*'*30))
	output_folder  = os.path.join(
		options['Pipeline Options']['somatic pipeline folder'],
		sample['PatientID'],
		'SomaticSniper')
	
	if not os.path.isdir(output_folder):
		os.makedirs(output_folder)

	output_file = os.path.join(output_folder, 
		"{0}_vs_{1}.somaticsniper.hq.vcf".format(sample['NormalID'], 
		sample['SampleID']))
	
	try:
		commands = get_somaticsniper_commands(sample, options)
		_labels, _commands, _outputs = zip(*commands)
		logging.info("{0}: SomaticSniper - successfully prepared the commands.".format(patientID))

		run_commands(commands, patientID = patientID, command_label = 'SomaticSniper')
	except Exception as exception:
		logging.critical("{0}: Could not run SomaticSniper: {1}".format(patientID, exception))


	somaticsniper_log = {
		'Inputs': [sample['NormalBAM'], sample['TumorBAM']],
		'Intermediate Files': [],
		'Outputs': [output_folder],
		'Commands': _commands,
		'Status': verify_output([output_file])
	}
	return sample, somaticsniper_log


def Strelka(sample, options):
	""" Calls vairants via Strelka
		Docs
		----------
			https://sites.google.com/site/strelkasomaticvariantcaller/home/configuration-and-analysis
		Output
		----------
			"3_called_variants/{patient}/Strelka/"
			Analysis
				chromosomes
					...
				config
					run.config.ini
				results
					all.somatic.indels.vcf   [VCF]
					all.somatic.snvs.vcf     [VCF]
					passed.somatic.indels.vcf[VCF]
					passed.somatic.snvs.vcf  [VCF]
	"""
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
	
	print('\tStrelka', flush = True)
	patientID = sample['PatientID']
	logging.info("{0}: {1} Running Strelka {1}".format(patientID, '*'*30))
	reference = options['Reference Files']['reference genome']
	strelka_location = options['Programs']['strelkafolder']
	patient_folder = os.path.join(
		options['Pipeline Options']['somatic pipeline folder'],
		sample['PatientID'])

	if not os.path.isdir(patient_folder): os.mkdir(patient_folder)

	strelka_folder = os.path.join(patient_folder, 'Strelka')
	strelka_results_folder = os.path.join(strelka_folder, 'results')
	
	#if not os.path.exists(strelka_output_dir):
	#	os.makedirs(strelka_output_dir)
	#--------------------------------- Configure Strelka --------------------------------
	#change working directory
	#base_config_file = os.path.join(strelka, "etc/strelka_config_isaac_default.ini")
	#case_config_file = os.path.join(strelka_output_dir, 'config.ini')
	#strelka_config = shutil.copyfile(base_config_file, case_config_file)
	strelka_configuration_ini = os.path.join(
		patient_folder, 'strelka_configuration.ini')

	try:
		with open(strelka_configuration_ini, 'w') as file1:
			#strelka_configuration_options = [i.split('\t')[0].strip() for i in strelka_configuration_options.s]
			file1.write('[user]\n')
			for row in strelka_configuration_options.split('\n'):
				row = [i for i in row.split('\t') if '=' in i]
				if len(row) > 0:
					file1.write(row[0] + '\n')
		logging.info("{0} (strelka 1 of 3 - configuration ini): Successfully created the ini file.".format(patientID))
	except Exception as exception:
		logging.critical("{0} (strelka 1 of 3 - configuration ini): The strelka configuration file could not be created at {1} ({2})".format(patientID, strelka_configuration_ini, exception))

	strelka_config_script = os.path.join(strelka_location, "bin", "configureStrelkaWorkflow.pl")

	
	strelka_config_command = "{script} --tumor={tumor} --normal={normal} --ref={reference} --config={config} --output-dir={output_dir}"
	strelka_config_command = strelka_config_command.format(
		script = strelka_config_script,
		tumor = sample['TumorBAM'],
		normal= sample['NormalBAM'],
		reference = reference,
		config = strelka_configuration_ini,
		output_dir = strelka_folder)

	print("\t\t(2/3) Calling Variants...", flush = True)
	if not os.path.exists(strelka_results_folder):
		print(strelka_config_command)
		Terminal(strelka_config_command, label = 'Strelka (configure)', show_output = True)

		#--------------------------------- Run Strelka --------------------------------
		print("Running Strelka...")
		
		_num_threads = options['Parameters']['MAX_CORES']
		strelka_run_command = "make -j {threads} -C {outputdir}".format(
			outputdir = strelka_folder, threads = _num_threads)
		Terminal(strelka_run_command, label = "Strelka (run)", show_output = True)
	else:
		logging.info("{0} (strelka 2 of 3 - variant calling): The output already exists.".format(patientID))

	output_prefixes = ['all.somatic.indels.vcf', 'all.somatic.snvs.vcf', 'passed.somatic.indels.vcf', 'passed.somatic.snvs.vcf']
	output_files = [os.path.join(strelka_results_folder, i) for i in output_prefixes]
	renamed_output_files = list()
	for of in output_files:
		output_folder, output_filename = os.path.split(of) #path, basename.vcf
		output_basename, _ = os.path.splitext(output_filename)# basename, vcf
		destination = "{0}_vs_{1}.".format(sample['NormalID'], sample['SampleID']) + output_filename + '.strelka.vcf' #_vs_.basename.strelka
		destination = os.path.join(output_folder, destination)
		try:
			os.rename(of, destination)
			renamed_output_files.append(destination)
		except:
			print("\t\tCould not rename {0} to {1}".format(of, destination))
	logging.info("{0} (strelka 3 of 3 - renamed output)")
	strelka_log = {
		'Inputs': [sample['NormalBAM'], sample['TumorBAM'], reference],
		'Outputs': [strelka_results_folder],
		'Commands': [strelka_config_command],
		'Status': verify_output(renamed_output_files)
	}
	return sample, strelka_log


def get_varscan_somatic_commands(sample, options):
	#varscan_location, output_folder, reference = get_program_parameters(sample, options, 'Varscan (somatic)')
	varscan_location = options['Programs']['varscan']
	output_folder = os.path.join(
		options['Pipeline Options']['somatic pipeline folder'],
		sample['PatientID'], 'Varscan')
	if not os.path.isdir(output_folder):
		os.makedirs(output_folder)
	#----------------------------------- Generate Pileup Files ------------------------------------------
	output_prefix = os.path.join(output_folder, "{0}_vs_{1}".format(sample['NormalID'], sample['SampleID']))
	#output_file = os.path.join(output_folder, output_prefix)
	#varscan_command = "java {memory} -jar {varscan} somatic {normal} {tumor} {output} --min-coverage {mc} --output-vcf 1"
	varscan_command = "java {memory} -jar {varscan} somatic {normal} {tumor} --output-snp {prefix}.raw.snp.vcf --output-indel {prefix}.raw.indel.vcf --output-vcf 1"
	varscan_command = varscan_command.format(
		varscan = varscan_location,
		memory = options['Parameters']['JAVA_MAX_MEMORY_USAGE'],
		normal = sample['NormalPileup'],
		tumor = sample['TumorPileup'],
		prefix = output_prefix,
		mc = options['Parameters']['MIN_COVERAGE'])

	#--------------------------------------- Postprocessing ---------------------------------------------
	pp_output_file = output_prefix + '.snp.Somatic.hc'
	process_command = "java {memory} -jar {varscan} processSomatic {vcf}.raw.snp.vcf"
	process_command = process_command.format(
		memory  = options['Parameters']['JAVA_MAX_MEMORY_USAGE'],
		varscan = varscan_location,
		vcf  = output_prefix)

	commands = [
		('Variant Discovery', varscan_command, output_prefix + '.raw.snp.vcf'),
		('Postprocessing', process_command, pp_output_file)
	]
	return commands

def varscan_somatic(sample, options):
	print("\tVarscan (somatic)")
	patientID = sample['PatientID']
	logging.info("{0}: {1} Varscan (somatic) {1}".format(patientID, '*'*30))
	#varscan_location, output_folder, reference = get_program_parameters(sample, options, 'Varscan (somatic)')
	output_folder = os.path.join(
		options['Pipeline Options']['somatic pipeline folder'],
		sample['PatientID'], 'Varscan')
	if not os.path.isdir(output_folder): os.makedirs(output_folder)

	print("\t\t(1/2) Generating pileup files...")
	sample, pileup_log = Pileup(sample, options)
	print("\t\tFinished generating pileup files...")
	try:
		commands = get_varscan_somatic_commands(sample, options)
		_, _commands, _ = zip(*commands)
		logging.info("{0}: Successfully prepared the varscan commands.".format(patientID))
	except Exception as exception:
		logging.critical("{0}: Could not prepare the Varscan commands.".format(patientID))

	run_commands(commands, 1, patientID = patientID, command_label = 'Varscan Somatic')

	varscan_output = _, _, varscan_output = zip(*commands)

	varscan_log = {
		'Inputs': [sample['NormalPileup'], sample['TumorPileup']],
		'Intermediate Files': [sample['NormalPileup'], sample['TumorPileup']],
		'Outputs': [output_folder],
		'Commands': _commands,
		'Status': verify_output(varscan_output)
	}
	return sample, varscan_log


def Mutect_pon_detection(sample, options):
	print("\tPanel of Normals (Artifact Detection via Mutect2)", flush = True)
	reference = options['Reference Files']['reference genome']
	gatk_program = options['Programs']['GATK']
	output_folder = os.path.join(
		options['Pipeline Options']['panel of normals'], 
		"6_panel_of_normals/artifacts/", 'Artifacts')
	pon_output = os.path.join(output_folder, "{sample}.artifact.vcf".format(sample = sample['NormalID']))
	#pon_output = os.path.join(PIPELINE_DIRECTORY, pon_output)

	if not os.path.isdir(_output_dir): os.makedirs(_output_dir)

	pon_command = "java {memory} -jar {GATK} -T MuTect2 -R {reference} -I:tumor {normal} --dbsnp {dbSNP} --artifact_detection_mode -L {targets} -o {output}"
	pon_command = pon_command.format(
		memory    = options['Parameters']['JAVA_MAX_MEMORY_USAGE'],
		GATK      = gatk_program,
		reference = reference,
		normal    = sample['NormalBAM'],
		dbSNP     = options['Reference Files']['dbSNP'],
		#cosmic    = options['Reference Files']['cosmic'],
		targets   = sample['ExomeTargets'],
		output    = pon_output)

	if not os.path.isfile(pon_output):
		Terminal(pon_command, show_output = False)
	else:
		print("\t\tThe PON file for {0} already exists.".format(sample['NormalID']))

	pon_log = {
		'Inputs': [sample['NormalBAM'], reference, options['Reference Files']['dbSNP']],
		'Intermediate Files': [],
		'Outputs': [pon_output],
		'Commands': [pon_command],
		'Status': verify_output([pon_output])
	}

	return sample, pon_log

def HaplotypeCaller(sample, options):
	sample, bqsr_log = BQSR(sample, options)
	#UnifiedGenotyper(sample, options)
	output_folder = os.path.join("/media/upmc/WD_Partition_2/RNA-seq/output/", sample['PatientID'], 'HaplotypeCaller')
	
	reference = options['Reference Files']['reference genome']
	gatk_location = options['Programs']['GATK']
	dna_output = os.path.join(output_folder, 'DNA-seq', "{0}_vs_{1}.raw_snps_indels.vcf".format(sample['NormalID'], sample['SampleID']))
	rna_output = os.path.join(output_folder, 'RNA-seq', "{0}.RNA.raw_snps_indels.vcf".format(sample['SampleID']))
	rna_filtered_output = os.path.join(output_folder, 'RNA-seq', "{0}.RNA.final_snps_indels.vcf".format(sample['SampleID']))
	
	checkdir(os.path.dirname(rna_output))
	checkdir(os.path.dirname(dna_output))
	
	dna_command = """java -jar {GATK} \
		-T HaplotypeCaller
		-R {reference} \
		-I {normal} \
		-I {tumor} \
		-L {targets} \
		-nct 6 \
		--dbsnp {dbSNP} \
		-o {output}""".format(
			GATK = gatk_location,
			reference = reference,
			normal = sample['NormalBAM'],
			tumor = sample['TumorBAM'],
			targets = sample['ExomeTargets'],
			dbSNP = options['Reference Files']['dbSNP'],
			output = dna_output)

	rna_format_command = """java -jar {GATK} \
		-T SplitNCigarReads \
		-R {reference} \
		-I {sample} \
		-o {outputbam} \
		-rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 \
		-UALLOW_N_CIGAR_READS
		"""
	rna_command = """java -jar {GATK} \
		-T HaplotypeCaller
		-R {reference} \
		-I {sample} \
		--dbsnp {dbSNP} \
		--dontUseSoftClippedBases \
		-o {output}""".format(
			GATK = gatk_location,
			reference = reference,
			sample = sample['RNABAM'],
			dbSNP = options['Reference Files']['dbSNP'],
			output = rna_output)
	#--filter_reads_with_N_cigar \
	#-drf mappingqualityunavailable\
	#-rf ReassignMappingQuality -DMQ 255 \
	rna_filter_command = """java -jar {GATK} \
		-T VariantFiltration \
		-R {reference} \
		-V {inputfile} \
		-window 35 \
		-cluster 3 \
		--filterName FS --filterExpression \"FS > 30.0\" \
		--filterName QD --filterExpression \"QD < 2.0\" \
		-o {output}""".format(
			GATK = gatk_location,
			reference = reference,
			inputfile = rna_output,
			output = rna_filtered_output)

	if not os.path.isfile(rna_output):
		Terminal(rna_command)
		
	else:
		print("RNA output already exists.")

	if not os.path.isfile(rna_filtered_output):
		print("FILTERING RNASEQ DATA")
		Terminal(rna_filter_command, show_output = True)
	else:
		print("RNASEQ DATA ALREADY EXISTS")


	if not os.path.isfile(dna_output) and False:
		Terminal(dna_command)
	else:
		print('DNA Output already exists.')

	hc_log = {
		'Inputs': [],
		'Outputs': [output_folder],
		'Intermediate Files': [],
		'Commands': []
	}

	return sample, hc_log

def UnifiedGenotyper(sample, options):

	output_folder = os.path.join("/media/upmc/WD_Partition_2/RNA-seq/output/", sample['PatientID'], 'UnifiedGenotyper')
	
	reference = options['Reference Files']['reference genome']
	gatk_location = options['Programs']['GATK']
	dna_output = os.path.join(output_folder, 'DNA-seq', "{0}_vs_{1}.raw_snps_indels.unifiedgenotyper.vcf".format(sample['NormalID'], sample['SampleID']))
	rna_output = os.path.join(output_folder, 'RNA-seq', "{0}.RNA.raw_snps_indels.unifiedgenotyper.vcf".format(sample['SampleID']))
	rna_filtered_output = os.path.join(output_folder, 'RNA-seq', "{0}.RNA.final_snps_indels.unifiedgenotyper.vcf".format(sample['SampleID']))
	checkdir(os.path.dirname(rna_output))
	checkdir(os.path.dirname(dna_output))
	dna_command = """java -jar {GATK} \
		-T UnifiedGenotyper \
		-R {reference} \
		-I {normal} \
		-I {tumor} \
		-L {targets} \
		-nct 6 \
		--dbsnp {dbSNP} \
		-o {output}""".format(
			GATK = gatk_location,
			reference = reference,
			normal = sample['NormalBAM'],
			tumor = sample['TumorBAM'],
			targets = sample['ExomeTargets'],
			dbSNP = options['Reference Files']['dbSNP'],
			output = dna_output)

	rna_command = """java -jar {GATK} \
		-T UnifiedGenotyper \
		-R {reference} \
		-I {sample} \
		--dbsnp {dbSNP} \
		-nct 6 \
		--maxRuntime 30 \
		--filter_reads_with_N_cigar \
		-o {output}""".format(
			GATK = gatk_location,
			reference = reference,
			sample = sample['RNABAM'],
			dbSNP = options['Reference Files']['dbSNP'],
			output = rna_output)

	rna_filter_command = """java ‐jar {GATK} \
		‐T VariantFiltration 
		‐R {reference} \
		‐V {inputfile} \
		‐window 35 \
		‐cluster 3 \
		‐filterName FS ‐filter "FS > 30.0" \
		‐filterName QD ‐filter "QD < 2.0" \
		‐o {output}""".format(
			GATK = gatk_location,
			reference = reference,
			inputfile = rna_output,
			output = rna_filtered_output)

	if not os.path.isfile(rna_output):
		Terminal(rna_command)
		#Terminal(rna_filter_command)
	else:
		print("RNA output already exists.")

	if not os.path.isfile(dna_output):
		Terminal(dna_command)
	else:
		print('DNA Output already exists.')

	ug_log = {
		'Inputs': [],
		'Outputs': [output_folder],
		'Intermediate Files': [],
		'Commands': []
	}

	return sample, ug_log

def BQSR(sample, options):
	reference = options['Reference Files']['reference Genome']
	dbSNP = options['Reference Files']['dbSNP']
	gatk_location = options['Programs']['GATK']
	output_folder = os.path.join("/media/upmc/WD_Partition_2/RNA-seq", 'recalibrated_genomes', sample['PatientID'])
	checkdir(output_folder)

	recalibration_table = os.path.join(output_folder, "{0}.RNA.recalibration_data.table".format(sample['SampleID']))
	covariate_table = os.path.join(output_folder, "{0}.RNA.covariate_data.table".format(sample['SampleID']))
	recalibration_plots = os.path.join(output_folder, "{0}.RNA.recalibration_plots.pdf".format(sample['SampleID']))
	cigar_bam = os.path.join(output_folder, "{0}.RNA.cigar.bam".format(sample['SampleID']))
	realigned_bam = os.path.join(output_folder, "{0}.RNA.recalibrated.bam".format(sample['SampleID']))
	_MESSAGE = """-L RESTRICTS ANALYSIS TO CHROM 1"""
	
	cigar_command = """java -jar {GATK} \
		-T SplitNCigarReads \
		-R {reference} \
		-I {inputbam} \
		-o {outputbam} \
		-rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 \
		-U ALLOW_N_CIGAR_READS""".format(
			GATK = gatk_location,
			reference = reference,
			inputbam = sample['RNABAM'],
			outputbam = cigar_bam)

	br_command = """java -jar {GATK} \
		-T BaseRecalibrator \
		-R {reference} \
		-I {bam} \
		-knownSites {dbSNP} \
		-o {output}""".format(
			GATK = gatk_location,
			reference = reference,
			dbSNP = dbSNP,
			bam = cigar_bam,
			output = recalibration_table)

	re_command  = """java -jar {GATK} \
		-T PrintReads \
		-R {reference} \
		-I {bam} \
		-BQSR {table} \
		-o {output}""".format(
			GATK = gatk_location,
			reference = reference,
			bam = cigar_bam,
			output = realigned_bam,
			table = recalibration_table)

	co_command = """java -jar {GATK} \
		-T BaseRecalibrator \
		-R {reference} \
		-I {realigned_bam} \
		-BQSR {table} \
		-knownSites {dbSNP} \
		-o {output}""".format(
			GATK = gatk_location,
			reference = reference,
			dbSNP = dbSNP,
			table = recalibration_table,
			realigned_bam = realigned_bam,
			output = covariate_table)

	plot_command = """ java -jar {GATK} \
		-T AnalyzeCovariates \
		-R {reference} \
		-before {before} \
		-after {after} \
		-plots {plots}""".format(
			GATK = gatk_location,
			reference = reference,
			before = recalibration_table,
			after = covariate_table,
			plots = recalibration_plots)
	#print(br_command)
	#print(re_command)
	#print(co_command)
	#print(plot_command)
	if not os.path.exists(cigar_bam):
		Terminal(cigar_command)
	if not os.path.exists(recalibration_table):
		Terminal(br_command)
	if not os.path.exists(realigned_bam):
		Terminal(re_command)
	if not os.path.exists(covariate_table):
		Terminal(co_command)
	if not os.path.exists(recalibration_plots):
		Terminal(plot_command)
	
	bqsr_log = None
	sample['RNABAM'] = realigned_bam
	return sample, bqsr_log

#----------------------------------------------------------------------------------------------------
#------------------------------------ Evaluate Somatic Variants -------------------------------------
#----------------------------------------------------------------------------------------------------
"""
	1. Annotate variants to genes, variant severity, and transcript (Oncotator, Annovar?, snpEFF?)
		a. Use Variant Effect Predictor, as the output is more standardized to the The Sequence Ontology Project
	1.5 Recalibrate quality scores?
	2. Filter against panel of normals
	3. Combine variants into single file per sample
	4. Use mutsigcv to assess mutation significance (want q-values less than 0.1)
"""

def GetVariantList(sample, options):

	patientID = sample['PatientID']
	normalID  = sample['NormalID']
	tumorID   = sample['SampleID']

	variants   = {
		#os.path.join(options['output']['Bambino'].format(patient = patientID),"{0}_vs_{1}.bambino.vcf".format(normalID, tumorID)),
		#os.path.join(options['output']['Haplotypecaller'].format(patient = patientID), "{0}_vs_{1}.raw.snps.indels.vcf".join(normalID, tumorID)),
		'muse':          os.path.join(options['output']['MuSE'].format(patient = patientID), "{0}_vs_{1}.MuSE.vcf".format(normalID, tumorID)),
		'mutect2':       os.path.join(options['output']['MuTect2'].format(patient = patientID), "{0}_vs_{1}.mutect2.vcf".format(normalID, tumorID)),
		'somaticsniper': os.path.join(options['output']['SomaticSniper'].format(patient = patientID), "{0}_vs_{1}.somaticsniper.hq.vcf".format(normalID, tumorID)),
		'strelka indel': os.path.join(options['output']['Strelka'].format(patient = patientID),"results", "{0}_vs_{1}.passed.somatic.indels.strelka.vcf".format(normalID, tumorID)),
		'strelka snv':   os.path.join(options['output']['Strelka'].format(patient = patientID),"results", "{0}_vs_{1}.passed.somatic.snvs.strelka.vcf".format(normalID, tumorID)),
		'varscan':       os.path.join(options['output']['Varscan (somatic)'].format(patient = patientID), "{0}_vs_{1}.snp.Somatic.hc".format(normalID, tumorID))
	}

	variants = {k:v for k, v in variants.items() if os.path.isfile(v)}

	return variants

	
def VariantEffectPredictor(sample, options, variants, output_folder):
	""" VariantEffectPredictor uses terminology as defined by the Sequence Ontology Project
	"""
	#perl variant_effect_predictor.pl --cache -i input.txt -o output.txt	

	reference = options['Reference Files']['reference genome']
	program = options['Programs']['varianteffectpredictor']

	pprint(variants)
	for index, source in enumerate(variants):
		_, fn = os.path.split(source)
		destination = os.path.join(output_folder, fn + '.VEP_annotated.vcf')

		command = """perl {vep} \
			--input_file {inputfile} \
			--output_file {output} \
			--fasta {reference} \
			--species homo_sapiens \
			--assembly GRCh38 \
			--format vcf \
			--cache \
			--html \
			--symbol \
			--biotype \
			--total_length \
			--numbers \
			--fields Consequence,Codons,Aminoacids,Gene,SYMBOL,Feature,EXON,PolyPhen,SIFT,Protein_position,BIOTYPE \
			--vcf""".format(
				vep = program,
				inputfile = source,
				output = destination,
				reference = reference)
		#print(command)
		Terminal(command, show_output = True)

	vep_log = {

	}
	return vep_log

def vcftomaf(vcf_file = None):
	vcftomaf_script = "/home/upmc/Programs/vcf2maf-1.6.12/vcf2maf.pl"
	output_folder = "/home/upmc/Documents/Variant_Discovery_Pipeline/8_combined_variants"
	output_file = os.path.join(output_folder, "maf_test.maf")

	if vcf_file is None:
		vcf_file = os.path.join(output_folder, "TCGA-2H-A9GF-01A-11D-A37C-09_TCGA-2H-A9GF-11A-11D-A37F-09_mutect_annotated.vcf")

	#perl vcf2maf.pl --input-vcf tests/test.vcf --output-maf tests/test.vep.maf --tumor-id WD1309 --normal-id NB1308

	command = """perl {script} \
				--input-vcf {vcf_file} \
				--output-maf {output} \
				--tumor-id {tumor} \
				--normal-id {normal}""".format(
					script = vcftomaf_script,
					vcf_file = vcf_file,
					output = output_file)

def VariantAnnotator(sample, options):
	"""
		 java -jar GenomeAnalysisTK.jar \
		-R reference.fasta \
		-T VariantAnnotator \
		-I input.bam \
		-o output.vcf \
		-A Coverage \
		-V input.vcf \
		-L input.vcf \
   --dbsnp dbsnp.vcf
	"""

	va_command =  """java {memory} -jar {GATK} \
		-R {reference} \
		-T VariantAnnotator \
		-I {bam} \
		-o output.vcf \
		-A Coverage \
		-V input.vcf \
		-L input.vcf \
		--dbsnp dbsnp.vcf""".format(
		memory = options['Parameters']['JAVA_MAX_MEMORY_USAGE'],
		GATK = gatk_program,
		reference = reference,
		bam = "")

def VariantRecalibrator(sample, options, variant_type):
	"""	Recalibrates variant quality scores.
		Example Usage
		-------------
		java -jar GenomeAnalysisTK.jar \ 
		-T VariantRecalibrator \ 
		-R reference.fa \ 
		-input raw_variants.vcf \ 
		-resource:hapmap,known=false,training=true,truth=true,prior=15.0 hapmap.vcf \ 
		-resource:omni,known=false,training=true,truth=true,prior=12.0 omni.vcf \ 
		-resource:1000G,known=false,training=true,truth=false,prior=10.0 1000G.vcf \ 
		-resource:dbsnp,known=true,training=false,truth=false,prior=2.0 dbsnp.vcf \ 
		-an DP \ 
		-an QD \ 
		-an FS \ 
		-an SOR \ 
		-an MQ \
		-an MQRankSum \ 
		-an ReadPosRankSum \ 
		-an InbreedingCoeff \
		-mode SNP \ 
		-tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 \ 
		-recalFile recalibrate_SNP.recal \ 
		-tranchesFile recalibrate_SNP.tranches \ 
		-rscriptFile recalibrate_SNP_plots.R 

		Reference
		---------
			https://software.broadinstitute.org/gatk/guide/article?id=2805
	"""
	#-------------------------- Generate the sample-specific options ---------------------------
	prefix = "{0}_vs_{1}".format(sample['NormalID'], sample['SampleID'])
	reference = options['Reference Files']['reference genome']

	#------------------------ Build the recalibration model ----------------------

	snp_build_command = """java {memory} -jar {GATK} \
		-T VariantRecalibrator \
		-R {reference} \
		-input {raw_variants} \
		-an QD \
		-an FS \
		-an SOR \
		-an MQ \
		-mode SNP \
		-resource:hapmap,known=false,training=true,truth=true,prior=15.0 {hapmap} \
		-resource:omni,known=false,training=true,truth=true,prior=12.0 {OMNI} \
		-resource:1000G,known=false,training=true,truth=false,prior=10.0 {OTG} \
		-resource:dbsnp,known=true,training=false,truth=false,prior=2.0 {dbSNP} \
		-tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 \
		-recalFile {prefix}.recalibrate_SNP.recal \
		-tranchesFile {prefix}.recalibrate_SNP.tranches \
		-rscriptFile {prefix}.recalibrate_SNP_plots.R""".format(
			memory       = options['Parameters']['JAVA_MAX_MEMORY_USAGE'],
			GATK         = options['Programs']['GATK'],
			reference    = reference,
			prefix       = prefix,
			raw_variants = "",
			mode         = variant_type.upper(),

			hapmap       = options['Reference Files']['hapmap'],
			OMNI         = options['Reference Files']['omni'],
			OTG          = options['Reference Files']['1000G'],
			dbSNP        = options['Reference Files']['dbSNP'])

	#------------------------- Recalibrate the SNP quality scores ------------------------
	#ts_filter_level = 99.9 recommended by GATK best practices
	#Retrieves 99.9% of true sites, includes an amount of false positives
	snp_apply_command = """java {memory} -jar {gatk} \
		-T ApplyRecalibration \
		-R {reference} \
		-input {raw_variants} \
		-mode SNP \
		--ts_filter_level 99.9 \
		-recalFile {prefix}.recalibrate_SNP.recal \
		-tranchesFile {prefix}.recalibrate_SNP_plots.R \
		-o {prefix}.recalibrated_snps_raw_indels.vcf""".format(
			memory = options['Parameters']['JAVA_MAX_MEMORY_USAGE'],
			gatk = gatk_program,
			reference = reference,
			raw_variants = "",
			prefix = prefix)
	#----------------------------- Recalibrate Indels ----------------------------------
	indel_build_command = """java {memory} -jar {gatk} \
		-T VariantRecalibrator \
		-R {reference} \
		-input {prefix}.recalibrated_snps_raw_indels.vcf \
		-an QD \
		-an FS \
		-an SOR \
		-an MQ \
		-mode INDEL \
		--maxGaussians 4 \
		-resource:dbsnp,known=true,training=false,truth=false,prior=2.0 {dbSNP} \
		-resource:mills,known=false,training=true,truth=true,prior=12.0 {mills} \
		-tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 \
		-recalFile {prefix}.recalibrate_INDEL.recal \
		-tranchesFile {prefix}.recalibrate_INDEL.tranches \
		-rscriptFile {prefix}.recalibrate_INDEL_plots.R""".format(
			memory       = options['Parameters']['JAVA_MAX_MEMORY_USAGE'],
			GATK         = options['Programs']['GATK'],
			reference    = reference,
			prefix       = prefix,
			raw_variants = "",

			mills        = "",
			dbSNP        = options['Reference Files']['dbSNP'])
	
	indel_apply_command = """java {memory} -jar {gatk} \
		-T ApplyRecalibration \
		-R {reference} \
		-input {raw_variants} \
		-mode INDEL \
		--ts_filter_level 99.9 \
		-recalFile {prefix}.recalibrate_INDEL.recal \
		-tranchesFile {prefix}.recalibrate_INDEL_plots.R \
		-o {prefix}.recalibrated_variants.vcf""".format(
			memory = options['Parameters']['JAVA_MAX_MEMORY_USAGE'],
			gatk = gatk_program,
			reference = reference,
			raw_variants = "",
			prefix = prefix)

def ProcessPipelineMethod(sample, options):
	#-------------------------- make all of the relevant folders ----------------------------
	pipeline_folder = os.path.join(PIPELINE_DIRECTORY, '8_combined_variants') #Holds folders for each patient
	patient_folder = os.path.join(pipeline_folder, sample['PatientID']) #holds folders for an individual patient
	raw_variant_folder = os.path.join(patient_folder, 'raw_variants')
	annotated_variants_folder = os.path.join(patient_folder, 'annotated_variants')

	for _dir in [pipeline_folder, patient_folder, raw_variant_folder, annotated_variants_folder]:
		if not os.path.isdir(_dir):
			os.mkdir(_dir)
	
	#-------------------------- Get a list of all the variants used in this analysis ---------------------
	sample_variants = GetVariantList(sample, options)
	pprint(sample_variants)
	
	#------------------------------ Annotate Variants ----------------------------------------
	raw_variants = list()

	for key, source in sample_variants.items():
		folder, file_name = os.path.split(source)
		destination = os.path.join(raw_variant_folder, file_name)
		shutil.copy2(source, destination)
		raw_variants.append(destination)


	vep_log = VariantEffectPredictor(sample, options, raw_variants, annotated_variants_folder)

class ProcessPipeline:
	def __init__(self, sample, options):
		#------------- Move Variant Files --------------------
		output_folder = ""
		variants = self.GetVariantList(sample, options)

	def move_files(sample, options, output_folder):
		variants = self.GetVariantList(sample, options)

	def run(self, sample, options):

		sample, log = VariantEffectPredictor(sample, options)


#----------------------------------------------------------------------------------------------------
#------------------------------------ Copynumber Tool Functions -------------------------------------
#----------------------------------------------------------------------------------------------------

def circular_binary_segmentation(input_file, output_file, rscript):
	"""
	library(DNAcopy)
	cn <- read.table("your.cn.file",header=F)
	CNA.object <-CNA( genomdat = cn[,6], chrom = cn[,1], maploc = cn[,2], data.type = 'logratio')
	CNA.smoothed <- smooth.CNA(CNA.object)
	segs <- segment(CNA.smoothed, verbose=0, min.width=2)
	segs2 = segs$output
	write.table(segs2[,2:6], file="out.file", row.names=F, col.names=F, quote=F, sep="\t")
	"""
	#NEED to strip first row of table before r script
	script = """
		library(DNAcopy)
		cn <- read.table("{input_file}",header=TRUE)
		CNA.object <-CNA( genomdat = cn[,6], chrom = cn[,1], maploc = cn[,2], data.type = 'logratio')
		CNA.smoothed <- smooth.CNA(CNA.object)
		segs <- segment(CNA.smoothed, verbose=0, min.width=2)
		segs2 = segs$output
		write.table(segs2[,2:6], file="{output_file}", row.names=F, col.names=F, quote=F, sep="\\t")
	""".format(
		input_file = input_file,
		output_file = output_file)
	with open(rscript, 'w') as file1:
		file1.write(script)
	Terminal("Rscript {fn}".format(fn = rscript),
		show_output = False)
	return output_file

def GDC_copynumber(sample, options):
	output_dir = options['Pipeline Options']['copynumber pipeline folder']
	return GDC_DNAseq(sample, options, ['Copy Number Variation'], output_dir)

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
		shutil.move(os.path.join(PIPELINE_DIRECTORY, varscan_basename + suffix), varscan_prefix + suffix)
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


def CNVkit(sample, options):
	print("\tCNVkit")
	patientID = sample['PatientID']
	logging.info("{0}: {1} Running CNVkit {1}".format(patientID, '*'*30))
	cnvkit_location = options['Programs']['cnvkit']
	output_folder = os.path.join(
		options['Pipeline Options']['copynumber pipeline folder'],
		sample['PatientID'],
		"CNVkit")
	reference = options['Reference Files']['reference genome']
	if not os.path.isdir(output_folder):
		os.makedirs(output_folder)

	reference_cnn = os.path.join(output_folder, "reference_cnv.cnn")

	batch_command = """{cnvkit} batch {tumor} --normal {normal} \
		--targets {targets} \
		--fasta {reference} \
		--output-reference {refcnn} \
		--output-dir {results} \
		--diagram --scatter""".format(
			cnvkit = cnvkit_location,
			tumor = sample['TumorBAM'],
			normal = sample['NormalBAM'],
			targets = sample['ExomeTargets'],
			reference = reference,
			refcnn = reference_cnn,
			results = output_folder)
	if not os.path.isfile(reference_cnn):
		logging.info("{0} CNVkit Command: {1}".format(patientID, batch_command))
		Terminal(batch_command, show_output = False)
	else:
		logging.info("{0}: CNVkit Command: THe output already exists as {1}".format(patientID, output_folder))
		print("The output already exists!")

	#I modified cnvkit/coverage.py, line 213 to add support for files with 6 rows (inc. '+' column)

	#-------------------------- Generate Coverage Files ---------------------------------

	#-------------------------- Genrate normal reference --------------------------------
	#Recommended to use all normal samples from the full cohort if they were sequenced on the same platform.
	#For each normal sample...

	#-------------------------- Generate CNV Plots --------------------------------------
	"""Batch Pipeline:
	cnvkit.py target baits.bed --split [--annotate --short-names] -o my_targets.bed
	cnvkit.py antitarget my_targets.bed [--access] -o my_antitargets.bed

	# For each sample...
	cnvkit.py coverage Sample.bam my_targets.bed -o Sample.targetcoverage.cnn
	cnvkit.py coverage Sample.bam my_antitargets.bed -o Sample.antitargetcoverage.cnn

	# With all normal samples...
	cnvkit.py reference *Normal.bam -t my_targets.bed -a my_antitargets.bed \
		[--fasta hg19.fa --male-reference] -o my_reference.cnn

	# For each tumor sample...
	cnvkit.py fix Sample.targetcoverage.cnn Sample.antitargetcoverage.cnn my_reference.cnn -o Sample.cnr
	cnvkit.py segment Sample.cnr -o Sample.cns

	# Optionally, with --scatter and --diagram
	cnvkit.py scatter Sample.cnr -s Sample.cns -o Sample-scatter.pdf
	cnvkit.py diagram Sample.cnr -s Sample.cns [--male-reference] -o Sample-diagram.pdf
	"""

	prepare_targets_command = "{CNVkit} target {targets} --split [--annotate --short-names] -o {localtargets}"

	#Add Heatmap Option


	cnvkit_log = {
		'Inputs': [],
		'Outputs': [output_folder],
		'Intermediate Files': [],
		'Commands': []
	}

	return sample, cnvkit_log

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
		logger.info("{0}: Running {1}".format(sample['PatientID'], self.pipeline_name))
		self.start = now()
		self.caller_logs = list()
		#self.pipeline_folder = self._get_pipeline_folder(options)

		try:
			callers = self._get_callers(callers) #Dictionary of caller methods
		except Exception as exception:
			logger.critical("{0}: Could not get the callers for pipeline {1}: {2}".format(sample['PatientID'], self.pipeline_name, exception))
		
		self._run_pipeline(sample, options, callers)
		
		self.stop = now()
		logger.info("{0}: {1} Ran from {2} to {3} for a duration of {4}.".format(
			sample['PatientID'],
			self.pipeline_name,
			self.start.isoformat(),
			self.stop.isoformat(),
			self.stop - self.start))

		self._update_log(sample, options, callers, self.start, self.stop)

	
	@staticmethod
	def _get_callers():
		return {}

	def _get_status(self):
		pass
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
			logger.info("{0}: Ran {1} in {2}".format(sample['PatientID'], caller_name, duration))
		
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

class SomaticPipeline(Pipeline):
	@staticmethod
	def _get_callers(callers):
		available_callers = {
			'pon': Mutect_pon_detection,
			'gdc': GDC_somatic,
			'muse': MuSE,
			'mutect2': Mutect2,
			'somaticsniper': SomaticSniper,
			'strelka': Strelka,
			'varscan': varscan_somatic,
			'haplotypecaller': HaplotypeCaller,
			'unifiedgenotyper': UnifiedGenotyper
		}
		callers = [i.lower() for i in callers]
		callers = {caller: available_callers[caller] for caller in callers if caller in available_callers}
		return callers
	@staticmethod
	def _get_pipeline_folder(options):
		folder = os.path.join(options['working directory'], "3_called_variants")
		return folder

class CopynumberPipeline(Pipeline):
	@staticmethod
	def _get_callers(callers):
		available_callers = {
			'gdc': GDC_copynumber,
			'varscan': varscan_copynumber,
			'cnvkit': CNVkit,
			'freec': FREEC
		}
		callers = [i.lower() for i in callers]
		callers = {caller: available_callers[caller] for caller in callers if caller in available_callers}
		return callers
	@staticmethod
	def _get_pipeline_folder(options):
		return os.path.join(options['working directory'], "4_copynumber_pipeline")


#----------------------------------------------------------------------------------------------------
#--------------------------------------------- Main -------------------------------------------------
#----------------------------------------------------------------------------------------------------
class SampleFiles:
	def __init__(self, sample, config, normal_api, tumor_api, normal_only = False):
		import hashlib

		#Appended to the beginning of the filename
		self.sample = sample
		self.normal_only = normal_only
		print("Finding NormalBAM")
		if 'NormalBAM' not in sample.keys():
			self.sample['NormalBAM'] = self._fetch_file(sample, config, normal_api)['path']
		
		if not self.normal_only and 'TumorBAM' not in sample.keys():
			self.sample['TumorBAM'] = self._fetch_file(sample, config, tumor_api)['path']
		

		self.status = self.verify_files(self.sample, normal_api, tumor_api)
		#self.status = True
		#print("SAMPLEFILES STATUS MANUALLY SET FOR TESTING PURPOSES!")
		#self.status = True

	def _fetch_file(self, sample, config, file_data):
		
		bam_folder = config['Pipeline Options']['bam folder']
		file_id = file_data['basic_info']['file_id']
		file_name= file_data['basic_info']['file_name']
		#source = os.path.join(PIPELINE_DIRECTORY, file_id, file_name) #File is automatically downloaded to here
		
		bam_file = os.path.join(bam_folder, file_id, file_name)


		_file_already_exists_and_is_valid = self._verify_file(bam_file, file_data['md5sum'])

		if not _file_already_exists_and_is_valid:
			print("\t\tDownloading ", file_id)
			response = gdc_api.download_file(file_id, output_folder, basename = None)
		else:
			response = {
				'path': bam_file,
				'output folder': os.path.dirname(bam_file),
				'status': True
			}
			print("\t\tThe file already exists and is valid")

		return response
	
	@staticmethod
	def _verify_file(filename, expected_md5sum):
		""" Verifies a single file
		def generate_file_md5(filename, blocksize=2**20):
		"""
		file_exists = os.path.isfile(filename)
		if file_exists:
			file_md5sum = generate_file_md5(filename)
		else: file_md5sum = False

		file_is_valid = file_exists and (file_md5sum == expected_md5sum)
		if False:
			print("\t\tExpected File: ", filename)
			print("\t\tFile Exists: ", file_exists)
			print("\t\tExpected md5sum: ", expected_md5sum)
			print("\t\tFile md5sum: ", file_md5sum)
			print("\t\tFile is Valid: ", file_is_valid)

		if 'chr1' in filename:
			file_is_valid = True

		return file_is_valid

	def verify_files(self, sample, normal_api, tumor_api):
		""" Verifies that both the normal and tumor bams are valid and completly downloaded.
		"""
		normal_file_valid = self._verify_file(sample['NormalBAM'], normal_api['md5sum'])
		print("\t\tThe Normal BAM file {0} verified.".format('was' if normal_file_valid else 'was not'))
		
		tumor_file_valid = self._verify_file(sample['TumorBAM'], tumor_api['md5sum']) or self.normal_only
		if not self.normal_only:
			print("\t\tThe Tumor BAM file {0} verified.".format('was' if tumor_file_valid else 'was not'))
		return normal_file_valid and tumor_file_valid


class GenomicsPipeline:
	def __init__(self, sample_filename, config_filename = None, somatic_callers = None, copynumber_callers = None):
		
		logger.info("Somatic Callers: " + ', '.join(somatic_callers))
		logger.info("Copynumber Callers: " + ', '.join(copynumber_callers))
		
		try:
			sample_list, config, self.completed_samples_list = self._load_files(sample_filename, config_filename)
		except Exception as exception:
			print("The was an error in GenomicsPipeline._load_files: ", exception)
			logger.critical("The was an error in GenomicsPipeline._load_files: " + str(Exception))
		print("Running the pipeline with {0} samples.".format(len(sample_list)), flush = True)
		pprint(self.completed_samples_list)

		logger.info("Running through the genomics pipeline with {0} samples.".format(len(sample_list)))

		for index, sample in enumerate(sample_list):
			print("({0}/{1}) {2}\t{3}".format(index+1, len(sample_list), sample['PatientID'], now().isoformat()), flush = True)

			use_this_sample = self._use_sample(sample, self.completed_samples_list)

			if use_this_sample['status']:
				logger.info("Analyzing {0}".format(sample['PatientID']))
				sample_status = self.run_sample(sample, config, somatic_callers, copynumber_callers)
			else:
				if use_this_sample['manual']: use_this_sample_message = 'Manually Skipped'
				elif use_this_sample['completed']: use_this_sample_message = 'Already Analyzed'
				elif use_this_sample['targets']: use_this_sample_message = 'Invalid Targets File: {0}'.format(sample['ExomeTargets'])
				else:
					use_this_sample_message = 'Unknown Error'
				logger.info("Skipped sample {0} ({1}).".format(sample['PatientID'], use_this_sample_message))


	@staticmethod
	def _load_files(sample_filename, config_filename):
		try:
			with open(sample_filename, 'r') as samplefile:
				sample_list = list(csv.DictReader(samplefile, delimiter = '\t'))
			logger.info("Loaded the sample file list at " + sample_filename)
		except:
			logger.critical("Could not load the sample file located at " + sample_filename)
		#----------------------------------------Process Config ------------------------------------------
		try:
			#config = DefaultConfig(filename = config_filename).to_json()
			config = configparser.ConfigParser()
			config.read(config_filename)
			logger.info("Loaded the config file at " + config_filename)
		except:
			logger.critical("Could not load the config file at " + config_filename)
		
		completed_samples = config['General Options']['completed sample log']
		try:
			with open(completed_samples, 'r') as file1:
				reader = list(csv.DictReader(file1, delimiter = '\t'))
				logger.info("Loaded the completed samples file at " + completed_samples)
		except:
			reader = list()
			logger.critical("Loaded the completed samples file at " + completed_samples)

		spf = config['Pipeline Options']['somatic pipeline folder']

		is_empty = lambda p: len( list( os.listdir( os.path.join(spf, p) ) ) ) == 0

		reader += [{'Barcode': i, 'Duration': "", 'Callers': ""} for i in os.listdir(spf) if not is_empty(i)]
		reader = []
		return sample_list, config, reader
	
	def _get_file_info(self, sample):
		normal_info = gdc_api.file_api(sample['NormalUUID'])
		tumor_info  = gdc_api.file_api(sample['SampleUUID'])
		return normal_info, tumor_info
	@staticmethod
	def _make_patient_folders(sample, options):

		snv_folder = options['Pipeline Options']['somatic pipeline folder']
		cnv_folder = options['Pipeline Options']['copynumber pipeline folder']
		temp_folder= options['Pipeline Options']['temporary folder']

		snv_folder = os.path.join(snv_folder, sample['PatientID'])
		cnv_folder = os.path.join(cnv_folder, sample['PatientID'])
		temp_folder = os.path.join(temp_folder, sample['PatientID'])

		if not os.path.isdir(snv_folder): os.mkdir(snv_folder)
		if not os.path.isdir(cnv_folder): os.mkdir(cnv_folder)
		if not os.path.isdir(temp_folder): os.mkdir(temp_folder)

	@staticmethod
	def _use_sample(sample, completed_samples):
		use_value = sample.get('Use', True)

		if isinstance(use_value, str): use_value = use_value.lower()
		
		manually_skipped = use_value in {False, 'false', 0, '0', 'no'}

		already_completed = sample['PatientID'] in [i['Barcode'] for i in completed_samples]

		invalid_targets = not os.path.isfile(sample['ExomeTargets'])

		use = {
			'status': not (manually_skipped or already_completed or invalid_targets),
			'manual': manually_skipped,
			'completed': already_completed,
			'targets': invalid_targets
		}

		return use

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
	

	def run_sample(self, sample, config, somatic_callers, copynumber_callers):
		print("#"*200)
		print("#"*90 + sample['PatientID'] + '#'*90)
		print("#"*200)
		sample_start = now()
		self._make_patient_folders(sample, config)
		normal_only = all(c == 'pon' for c in (somatic_callers + copynumber_callers))

		normal_file_api, tumor_file_api = self._get_file_info(sample)
			
		pprint(sample)
		#--------------------------- Download and verify the BAM files ---------------------------------
		try:
			prepared_files = SampleFiles(sample, config, normal_file_api, tumor_file_api, normal_only)
			file_status = prepared_files.status
			if file_status:
				logger.info("{0}: The BAM files exist and are valid. NormalBAM={1}. TumorBAM={2}".format(sample['PatientID'],sample['NormalBAM'], sample['TumorBAM']))
			else:
				logger.error("{0}: The BAM files are invalid! NormalBAM={1}. TumorBAM={2}".format(sample['PatientID'],sample['NormalBAM'], sample['TumorBAM']))
		except Exception as exception:
			file_status = False
			logger.critical("{0}: GenomicsPipeline.run_sample: The BAM files could not be loaded ({1})".format(sample['PatientID'], exception))



		if file_status:
			with open(console_log_file_filename, 'a') as file1:
				file1.write("#" * 52 + sample['PatientID'] + "#" * 52)
			if somatic_callers is not None:

				somatic_pipeline = SomaticPipeline(sample, config, somatic_callers)

			if copynumber_callers is not None:

				copynumber_pipeline = CopynumberPipeline(sample, config, copynumber_callers)

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
		csv_log(sample_log)

		#mark the sample as completed.
		
		self.save_finished_sample(
			sample = sample,
			options = config,
			duration = sample_stop - sample_start, 
			callers = somatic_callers + copynumber_callers)

		return file_status

	def save_finished_sample(self, sample, options, duration = "", callers = []):
		try:
			sample_row = {
				'Barcode': sample['PatientID'],
				'Duration': duration,
				'Callers': ','.join(callers)
			}

			with open(options['General Options']['completed sample log'], 'w', newline = "") as file1:
				writer = csv.DictWriter(file1, delimiter = '\t', fieldnames = ['Barcode', 'Duration', 'Callers'])
				writer.writeheader()
				writer.writerows(self.completed_samples_list)
			logging.info("{0}: Successfully logged the results of the completed analysis.".format(sample['PatientID']))
		except Exception as exception:
			logging.error("{0} completed successfully, but the program was unable to log the result.".format(sample['PatientID']))



def debug():
	fn = "C:\\Users\\Deitrickc\\Google Drive\\Genomics\\Cake_1.0.tar\\Cake\\trunk\\scripts\\variant_effect_predictor.pl"
	print(os.path.getsize(fn))

if __name__ == "__main__" and True:
	#run_pipelines(samples = "/home/upmc/Documents/Variant_Discovery_Pipeline/sample_list.tsv")
	config_filename = os.path.join(PIPELINE_DIRECTORY, "api_files", "pipeline_project_options.txt")
	if True:
		sample_filename = os.path.join(PIPELINE_DIRECTORY, "sample_list_WD.tsv")
		somatic_callers = ['muse']
		copynumber_callers = ['']
		#somatic_callers = ['MuSE', 'Varscan', 'Strelka', 'Somaticsniper', 'Mutect2']
		#copynumber_callers = ['varscan', 'cnvkit', 'freec']
	else:
		sample_filename = os.path.join("/media/upmc/09B504C10EA957B6/DNA-seq", "DNA-seq_Sample_List.tsv")
		#sample_filename = os.path.join(PIPELINE_DIRECTORY, "tcga_esca_sample_list.adenocarcinoma.DELL.2017-02-09.tsv")
		somatic_callers = ['muse']#'MuSE', 'Varscan', 'Strelka', 'Somaticsniper', 'Mutect2']
		copynumber_callers = []
	pipeline = GenomicsPipeline(sample_filename, config_filename, somatic_callers, copynumber_callers)
else:
	with open("sample_list.tsv", 'r') as samplefile:
		sample_list = list(csv.DictReader(samplefile, delimiter = '\t'))
	#----------------------------------------Process Config ------------------------------------------
	config = DefaultConfig(filename = None).to_json()

	pprint(sample_list)
	ProcessPipelineMethod(sample_list[0], config)
#3872