import os
from .basicworkflow import Workflow, getPipelineFolder

class MuTect(Workflow):

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


class MuTect2(Workflow):

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