import os
from .basicworkflow import Workflow, getPipelineFolder

class DepthOfCoverage(Workflow):

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