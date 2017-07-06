import os
from .basicworkflow import Workflow

class MuSE(Workflow):

	def setCustomEnvironment(self, sample, options):
		self.caller_name = "MuSE"
		#Define relevant filenames
		self.call_output = '.'.join(self.abs_prefix.split('.')[:-1]) + '.MuSE.txt' #remove caller name form prefix
		self.sump_output = os.path.splitext(self.call_output)[0] + ".vcf"

		self.full_output = [self.call_output, self.sump_output]
		self.final_output= self.sump_output

	def runCallerWorkflow(self, sample):
		#Run the caller commands
		self.call_output_status = self.runMuseCall(sample)
		self.sump_output_status = self.runMuseSump()

	def runMuseCall(self, sample):

		call_command = "{program} call -O {prefix} -f {reference} {tumor} {normal}".format(
			program 	= self.program,
			reference 	= self.reference,
			prefix 		= self.call_output.replace('.MuSE.txt', ''),
			tumor 		= sample['TumorBAM'],
			normal 		= sample['NormalBAM'])
		label = "MuSE Call"
		output_result = self.runCallerCommand(call_command, label, self.call_output)

		return output_result

	def runMuseSump(self):
		
		sump_command = "{program} sump -I {call_output} -E -D {dbSNP} -O {output}".format(
			program 	= self.program,
			prefix 		= self.abs_prefix,
			dbSNP 		= self.dbSNP,
			call_output = self.call_output,
			output 		= self.sump_output)
		label = "MuSE Sump"
		output_result = self.runCallerCommand(sump_command, label, self.sump_output)

		return output_result