import os
import shutil
from .basicworkflow import Workflow

class Strelka(Workflow):

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