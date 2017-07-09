import os
from .basicworkflow import Workflow

class HaplotypeCaller(Workflow):
	def __init__(self, sample, options, input_bam = None):
		self.input_bam = input_bam
		self.options = options
		print("input_bam: ", input_bam)
		self.full_output = []
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