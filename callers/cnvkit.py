import os
from .basicworkflow import Workflow

class CNVkit(Workflow):

	def setCustomEnvironment(self, sample, options):
		self.caller_name = "CNVKit"
		self.output_folder = options['variants-copynumber', sample['PatientID'], self.caller_name]
		self.final_output = os.path.join(self.output_folder, 'reference_cnv.cnn')
		self.full_output = [
			sample['SampleID'] + '.cns',
			sample['SampleID'] + '.cnr',
			sample['SampleID'] + '.targetcoverage.cnn'
		]
		self.full_output = [os.path.join(self.output_folder, fn) for fn in self.full_output]
	
	def runCallerWorkflow(self, sample):
		call_status = self.runBatchCommand(sample)


	def runBatchCommand(self, sample):
		reference_cnn = os.path.join(self.output_folder, "reference_cnv.cnn")

		batch_command = """{cnvkit} batch {tumor} \
			--normal {normal} \
			--targets {targets} \
			--fasta {reference} \
			--output-reference {refcnn} \
			--output-dir {results} \
			--diagram \
			--scatter""".format(
				cnvkit = self.program,
				tumor = sample['TumorBAM'],
				normal = sample['NormalBAM'],
				targets = sample['ExomeTargets'],
				reference = self.reference,
				refcnn = reference_cnn,
				results = self.output_folder
			)

		label = "Batch Command"
		output_result = self.runCallerCommand(batch_command, label, reference_cnn)

		return output_result