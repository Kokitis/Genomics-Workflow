import os
from .basicworkflow import Workflow, getPipelineFolder

class BaseQualityScoreRecalibration(Workflow):

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
		#	self.recalibration_plots,
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