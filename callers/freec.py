import os
import csv
from .basicworkflow import Workflow

class FREEC(Workflow):
	def setCustomEnvironment(self, sample, options):
		self.caller_name = "FREEC"
		self.output_folder = options['variants-copynumber', sample['PatientID'], self.caller_name]
		
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

		workflow_status = [
			significance_normal,
			significance_tumor,
			plot_normal_status,
			plot_sample_status
		]

	
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
		#command = """cat {script} | R --slave --args {cnvs} {ratios}""".format(
		#	script 	= self.assess_significance_script,
		#	cnvs 	= cnvs,
		#	ratios 	= ratios
		#)

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
		#command = """cat {script} | R --slave --args {ploidy} {ratios} {baf}""".format(
		#	script = self.plot_script,
		#	ploidy = "2",
		#	ratios = ratios,
		#	baf = baf_file)

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