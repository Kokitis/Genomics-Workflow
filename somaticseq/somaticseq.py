#from truthset import Truthset 
from callers import Workflow
import os
import shutil
from pprint import pprint

from github import callertools
from github import vcftools

class SomaticSeq(Workflow):
	def __init__(self, sample, options, mode, truthset = None, classifier = None):
		"""
			Parameters
			----------
				mode: {'trainer', 'prediction', 'table'}
					* 'trainer': Trains the classifier.
						Requires 'truthset'
					* 'prediction': Calculated the probability that a variant is a true variant.
						Requires 'snp_classifier', the output from the 'trainer' step.
					* 'table': generates the table used for the training and prediction steps.
				truthset: Truthset; default None
					The truthset object to use. Provides a reference to a ground truth file containing
						true variant positions.
				classifier: string; default None
					Path to an .RData file the was generated from the "trainer" mode.

		"""
		self.mode = mode 
		self.ss_classifier = classifier # Should only be provided in training mode.
		self.truthset = truthset # Should be None if in prediction/table mode

		super().__init__(sample, options)

	def setCustomEnvironment(self, sample, options):
		self.somaticseq_folder 		= options['Programs']['somaticseq']
		self.somaticseq_program 	= os.path.join(self.somaticseq_folder, "SomaticSeq.Wrapper.sh")
		self.modify_vjsd_script 	= os.path.join(self.somaticseq_folder, "modify_VJSD.py")
		self.ada_trainer_script 	= os.path.join(self.somaticseq_folder, "r_scripts", "ada_model_builder.R")
		self.ada_prediction_script	= os.path.join(self.somaticseq_folder, "r_scripts", "ada_model_predictor.R")
		self.tsv_to_vcf_script  	= os.path.join(self.somaticseq_folder, "SSeq_tsv2vcf.py")
		self.merged_vcf2tsv_script  = os.path.join(self.somaticseq_folder, "SSeq_merged.vcf2tsv.py")

		self.output_folder 				= options.getPipelineFolder('somaticseq-' + self.mode, sample['PatientID'])
		self.original_callset_folder 	= options.getPipelineFolder('callset', sample['PatientID'])
		self.split_callset_folder 		= options.getPipelineFolder('somaticseq-callset', sample['PatientID'], 'original-split')

	def runCallerWorkflow(self, sample):
		print("Running Workflow...")
		patientId = sample['PatientID']
		callset = self._getRawCallset()
		#pprint(callset)
		#process_output_folder = getPipelineFolder('somaticseq-' + self.mode, patientId)
		
		processed_callset = self._processVJSDFiles(
			callset, 
			self.output_folder,
			patientId
		)
		
		merged_raw_variant_file = self._mergeVariantFiles(
			processed_callset, 
			self.output_folder,
			patientId
		)
		
		filtered_raw_variant_file= self._filterVariantTargets(
			merged_raw_variant_file, 
			self.output_folder,
			patientId
		)
		
		self.trained_snp_table = self._generateCovariateTable(
			sample, 
			callset, 
			filtered_raw_variant_file, 
			self.output_folder
		)

		if self.mode == 'table':
			pass #The trained_snp_table is the final result.
		elif self.mode == 'trainer':
			self.classifier = self.buildTrainer(self.trained_snp_table)
		elif self.mode == 'prediction':
			self.classifier = None
			prediction_table = self.runPredictor(self.trained_snp_table)
			prediction_vcf = self._convertToVcf(sample, prediction_table)
		else:
			message = "'{}' is not a supported mode for SomaticSeq! {'trainer', 'prediction', 'table'}".format(self.mode)
			raise ValueError(message)

		return self.trained_snp_table

	def _getRawCallset(self):
		classifier = callertools.CallerClassifier()

		original_callset_folder = self.original_callset_folder

		split_callset_folder = self.split_callset_folder

		original_callset = classifier(original_callset_folder)
		#pprint(original_callset)
		vcftools.splitCallset(original_callset, split_callset_folder)
		
		callset = classifier(split_callset_folder, type = 'snp')
		#pprint(callset)	

		return callset

	def _processVJSDFiles(self, callset, output_folder, patientId):
		"""
			Parameters
			----------
				callset: dict<sting:string)
					A dictionary with references to the patient's callset.
				output_folder: string
					The folder to place all output files in.
				patientId: string
					The patient's barcode
		"""
		processed_callset = dict()
		for caller, input_file in callset.items():
			print("\tProcessing {}...".format(caller))
			if caller == 'muse':
				method = 'MuSE'
			elif caller == 'somaticsniper':
				method = 'SomaticSniper'
			elif 'varscan' in caller:
				method = 'VarScan2'
			else:
				method = None

			output_filename = os.path.join(
				output_folder, 
				"{}.training.modified.vcf".format(patientId)
			)

			if method:
				command = """python3 {program} \
					--call-method {method} \
					--input-vcf {infile} \
					--output-vcf {outfile}""".format(
					program = self.modify_vjsd_script,
					method = method,
					infile = input_file,
					outfile = output_filename)
				#systemtools.Terminal(command, use_system = True)
				self.runCallerCommand(command, 'VJSD ({})'.format(method), output_filename)
			else:
				shutil.copy2(input_file, output_filename)

			processed_callset[caller] = output_filename


		return processed_callset

	def _mergeVariantFiles(self, callset, output_folder, patientId):
		"""
			Note: The only records that are merged are those that are
			unfiltered in at least one caller.
		"""
		print("Merging files...")
		output_filename = os.path.join(
			output_folder,
			"{}.{}.modified.merged.vcf".format(patientId, self.mode)
		)

		command = """java -jar {gatk} \
			--analysis_type CombineVariants \
			--reference_sequence {reference} \
			--genotypemergeoption UNSORTED \
			--filteredrecordsmergetype KEEP_IF_ANY_UNFILTERED \
			--variant {muse} \
			--variant {mutect2} \
			--variant {somaticsniper} \
			--variant {strelka} \
			--variant {varscan} \
			--out {output}""".format(
			gatk 		= self.gatk_program,
			reference 	= self.reference,
			muse 		= callset['muse-snp'],
			mutect2 	= callset['mutect2-snp'],
			somaticsniper = callset['somaticsniper-snp'],
			strelka 	= callset['strelka-snp'],
			varscan 	= callset['varscan-snp'],
			output 		= output_filename
		)

		self.runCallerCommand(command, 'CombineVariants', output_filename)
		return output_filename
	
	def _filterVariantTargets(self, input_filename, output_folder, patientId):
		print("Excluding non-exome targets...")
		output_filename = os.path.join(
			output_folder,
			"{}.{}.modified.merged.excluded.vcf".format(patientId, self.mode)
		)

		command = """intersectBed -header -a {infile} -b {targets} > {output}""".format(
			infile = input_filename,
			targets = self.targets,
			output = output_filename
		)
		command = """intersectBed -header -a {infile} -b {targets}""".format(
			infile = input_filename,
			targets = self.targets,
			output = output_filename
		)

		#systemtools.Terminal(command, use_system = True)
		self.runCallerCommand(
			command = command, 
			label = "Filtering Targets", 
			expected_output = output_filename, 
			output_filename = output_filename,
			verbose = ['command', 'status']
		)
		#print("\tResult: {}\t{}".format(os.path.exists(output_filename), output_filename))
		return output_filename
	
	def _generateCovariateTable(self, sample, callset, merged_callset, output_folder):
		print("Generating the covariate table...")
		#start_time = time.time()
		output_filename = os.path.join(
			output_folder,
			"{}.{}.modified.merged.excluded.snp.tsv".format(sample['PatientID'], self.mode)
		)
		#print("Merged_callset: {}\t{}".format(os.path.exists(merged_callset), merged_callset))
		command = """python3 {script} \
			--p-scale phred \
			--genome-reference {reference} \
			--normal-bam-file {normal} \
			--tumor-bam-file {tumor} \
			--dbsnp-vcf {dbSNP} \
			--cosmic-vcf {cosmic} \
			--vcf-format {merged} \
			--muse-vcf {muse} \
			--mutect-vcf {mutect2} \
			--somaticsniper-vcf {somaticsniper} \
			--strelka-vcf {strelka} \
			--varscan-vcf {varscan} \
			--output-tsv-file {output}""".format(
				script 			= self.merged_vcf2tsv_script,
				reference 		= self.reference,
				normal 			= sample['NormalBAM'],
				tumor  			= sample['TumorBAM'],
				dbSNP  			= self.dbSNP,
				cosmic 			= self.cosmic,
				truth  			= self.truthset,
				merged 			= merged_callset,
				muse   			= callset['muse-snp'],
				mutect2 		= callset['mutect2-snp'],
				somaticsniper 	= callset['somaticsniper-snp'],
				strelka 		= callset['strelka-snp'],
				varscan 		= callset['varscan-snp'],
				output 			= output_filename
			)

		#print(command)
		if self.mode == 'trainer':
			command += " --ground-truth-vcf " + self.truthset

		#if not os.path.exists(output_filename):
		#print(command)
		self.runCallerCommand(command, "Generating Table", output_filename)
		#os.system(command)
		#stop_time = time.time()

		#duration = timetools.Duration(seconds = stop_time - start_time)
		#print("Converted the table in ", duration.isoformat())
		return output_filename

	def _convertToVcf(self, sample, input_filename):
		print("Converting to a VCF file...")

		output_filename = os.path.splitext(input_filename)[0] + ".vcf"

		command = """python3 {script} \
			--tsv-in {infile} \
			--vcf-out {outfile} \
			--normal-sample-name {normalid} \
			--tumor-sample-name {tumorid} \
			--emit-all \
			--individual-mutation-tools {tools} \
			--phred-scale""".format(
				script = self.tsv_to_vcf_script,
				infile = input_filename,
				outfile = output_filename,
				normalid = sample['NormalID'],
				tumorid = sample['SampleID'],
				tools = "MuSE CGA SomaticSniper Strelka VarScan2"
			)

		self.runCallerCommand(command, "Generating VCF", output_filename)

		

	def buildTrainer(self, input_filename):
		print("Building model...")
		command = "{script} {infile}".format(
			script = self.ada_trainer_script,
			infile = input_filename
		)
		expected_output = input_filename + '.Classifier.RData'
		if not os.path.exists(expected_output):
			#systemtools.Terminal(command, use_system = True)
			self.runCallerCommand("Training the Model", expected_output)
		else:
			print("The Somaticseq classifier already exists.")
		return expected_output

	def runPredictor(self, input_filename):

		"""
		"""
		basename = os.path.splitext(os.path.basename(input_filename))[0]
		output_folder = os.path.dirname(input_filename)
		output_filename = "{}.predicted_scores.tsv".format(basename)
		output_filename = os.path.join(output_folder, output_filename)

		print("Predicting Scores...")
		command = "{script} {classifier} {infile} {outfile}".format(
			script = self.ada_prediction_script,
			classifier = self.ss_classifier,
			infile = input_filename,
			outfile = output_filename
		)
		print(command)
		#systemtools.Terminal(command, use_system = True)
		self.runCallerCommand(command, "Calculating Predictions", output_filename)

		return output_filename
