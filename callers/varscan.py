import os
from .basicworkflow import Workflow

class Varscan(Workflow):

	def setCustomEnvironment(self, sample, options):
		self.caller_name = "Varscan"
		self.raw_snps   = self.abs_prefix + '.snp.vcf'
		self.raw_indels = self.abs_prefix + '.indel.vcf'

		self.somatic_hc = self.abs_prefix + '.snp.Somatic.hc.vcf'
		self.somatic_lc = self.abs_prefix + '.snp.Somatic.lc.vcf'
		self.germline   = self.abs_prefix + '.snp.Germline.vcf'
		self.germline_hc= self.abs_prefix + '.snp.Germline.hc.vcf'
		self.loh        = self.abs_prefix + '.snp.LOH.vcf'
		self.loh_hc     = self.abs_prefix + '.snp.LOH.hc.vcf'

		self.full_output = [
			self.raw_snps, 
			self.raw_indels,
			self.somatic_hc,
			#self.somatic_lc,
			self.germline,
			self.germline_hc,
			self.loh,
			self.loh_hc]

		self.final_output = self.somatic_hc
	
	def runCallerWorkflow(self, sample):

		pileup_status = self.generateSinglePileup(sample)
		pileup_file = pileup_status['outputFiles']
		variant_discovery_status = self.runSingleVariantDiscovery(pileup_file)
		#self.runDoubleVariantDiscovery(sample)
		processing_status = self.postProcessing()

	def generateSinglePileup(self, sample):
		pileup = os.path.join(
			self.temp_folder,
			"{}_vs_{}.intermediate.mpileup".format(
				sample['NormalID'],
				sample['SampleID']
			)
		)
		# command = """samtools mpileup -f {reference} -q 1 -B {normal} {tumor} > {pileup}""".format(
		command = """samtools mpileup \
					-f {reference} \
					-q 1 \
					-B {normal} \
					{tumor}""".format(
			reference = self.reference,
			normal 	= sample['NormalBAM'],
			tumor 	= sample['TumorBAM'],
			pileup 	= pileup
		)
		#output_result = self.runCallerCommand(command, "GenerateSinglePileup", pileup, output_filename = pileup, show_output = True)
		#Need to use os.system due to how samtools generates mpileups and the size of the output.
		if not os.path.exists:
			os.system(command + " > {}".format(pileup))
		output_result = {'outputFiles': pileup}

		return output_result

	def runSingleVariantDiscovery(self, pileup):
		expected_output = [self.raw_snps, self.raw_indels]
		command = """java -jar {program} somatic {pileup} {output} \
			--mpileup 1 \
			--min-coverage 8 \
			--min-coverage-normal 8 \
			--min-coverage-tumor 6 \
			--min-var-freq 0.10 \
			--min-freq-for-hom 0.75 \
			--normal-purity 1.0 \
			--tumor-purity 1.00 \
			--p-value 0.99 \
			--somatic-p-value 0.05 \
			--strand-filter 0 \
			--output-vcf""".format(
				program = self.program,
				pileup  = pileup,
				output  = self.abs_prefix)
		print(command)
		output_result = self.runCallerCommand(command, "RunSingleVariantDiscovery", expected_output)

		return output_result
	
	def runDoubleVariantDiscovery(self, sample):

		normal_pileup = self.generatePileup(sample['NormalBAM'], sample['NormalID'])
		tumor_pileup = self.generatePileup(sample['TumorBAM'], sample['SampleID'])
		command = """java {memory} -jar {varscan} somatic {normal} {tumor} \
			--output-snp {snp} \
			--output-indel {indel} \
			--output-vcf 1""".format(
				varscan = self.program,
				memory  = self.max_memory_usage,
				normal  = normal_pileup,
				tumor   = tumor_pileup,
				snp     = self.raw_snps,
				indel   = self.raw_indels,
				mc      = self.min_coverage)
		output_files = [self.raw_snps, self.raw_indels]
		label = "Variant Discovery"
		status = self.runCallerCommand(command, label, output_files)

		return status

	def postProcessing(self):
		expected_output = [self.somatic_hc, self.germline, self.germline_hc, self.loh, self.loh_hc]

		process_command = """java -jar {varscan} processSomatic {vcf} \
			--min-tumor-freq 0.10 \
			--max-normal-freq 0.05 \
			--p-value 0.07""".format(
				memory  = self.max_memory_usage,
				varscan = self.program,
				vcf     = self.raw_snps
			)
		label = "Varscan Postprocessing"
		output_result = self.runCallerCommand(process_command, label, expected_output)
		return output_result

	def renameOutputFiles(self):
		pass


class VarscanCopynumber(Workflow):

	def setCustomEnvironment(self, sample, options):
		self.caller_name = "Varscan"

		self.output_folder = options['variants-copynumber', sample['PatientID'], self.caller_name]

		self.base_prefix = "{normal}_vs_{tumor}.{prefix}".format(
			tumor=sample['SampleID'],
			normal=sample['NormalID'],
			prefix=self.caller_name.lower()
		)

		self.abs_prefix = os.path.join(
			self.output_folder,
			self.base_prefix
		)

		self.rscript_filename   = os.path.join(self.output_folder, "{0}.varscan_CBS.r".format(sample['PatientID']))
		self.copynumber_output  = self.abs_prefix + '.copynumber'           # [prefix].copynumber
		self.called_copynumbers = self.copynumber_output + '.called'    # [prefix].copynumber.called
		self.called_homdels     = self.called_copynumbers + '.homdel'   # [prefix].copynumber.called.homdel
		self.copynumber_segments =self.called_copynumbers + '.segments' # [prefix].copynumber.called.segments
		self.full_output = [
			self.copynumber_output,
			self.called_copynumbers,
			self.called_homdels,
			self.copynumber_segments
		]
		self.final_output = self.copynumber_segments

	def runCallerWorkflow(self, sample):
		normal_pileup   = self.generatePileup(sample['NormalBAM'], sample['NormalID'])
		tumor_pileup    = self.generatePileup(sample['TumorBAM'], sample['SampleID'])

		copynumber_ratio_status = self.callCopynumberRatios(normal_pileup, tumor_pileup)
		copynumber_caller_status  = self.copyCaller()

		segmentation_status = self.circularBinarySegmentation()

	def circularBinarySegmentation(self):
		# NEED to strip first row of table before r script
		"""
				circular_binary_segmentation(
			varscan_prefix + '.copynumber.called', 
			varscan_prefix + '.copynumber.called.segments',
			rscript_file)
		"""
		script_contents = """library(DNAcopy)
			cn <- read.table("{input_file}",header=TRUE)
			CNA.object <-CNA( genomdat = cn[,6], chrom = cn[,1], maploc = cn[,2], data.type = 'logratio')
			CNA.smoothed <- smooth.CNA(CNA.object)
			segs <- segment(CNA.smoothed, verbose=0, min.width=2)
			segs2 = segs$output
			write.table(segs2[,2:6], file="{output_file}", row.names=F, col.names=F, quote=F, sep="\\t")""".format(
			input_file = self.called_copynumbers,
			output_file = self.copynumber_segments)

		with open(self.rscript_filename, 'w') as file1:
			for line in script_contents.splitlines():
				file1.write(line + '\n')

		command = """Rscript {script}""".format(script = self.rscript_filename)
		label = "Circulr Binary Segmentation"
		output_result = self.runCallerCommand(command, label, self.copynumber_segments)
		return output_result

	def callCopynumberRatios(self, normal_pileup, tumor_pileup):
		# -------------------------- Varscan (copynumber - caller) -----------------------------
		copynumber_command = """java {memory} -jar {varscan} copynumber {normal} {tumor} {prefix} \
			--min-base_qual {mbq} \
			--min-map-qual {mmq}""".format(
				mmq 	= self.min_mapping_quality,
				mbq 	= self.min_base_quality,
				normal 	= normal_pileup,
				tumor 	= tumor_pileup,
				varscan = self.program,
				memory 	= self.max_memory_usage,
				prefix 	= self.abs_prefix
		)
		label = "Varscan Copynumber Ratios"
		output_result = self.runCallerCommand(copynumber_command, label, self.copynumber_output)

		return output_result

	def copyCaller(self):
		command = """java {memory} -jar {varscan} copyCaller {copynumber} \
		--output-file {called} \
		--output-homdel-file {homdels} \ 
		--min-coverage {mq}""".format(
			memory 		= self.max_memory_usage,
			varscan 	= self.program,
			copynumber 	= self.copynumber_output,
			called 		= self.called_copynumbers,
			homdels 	= self.called_homdels,
			mq 			= self.min_coverage)
		label = "Varscan CopyCaller"
		output_result = self.runCallerCommand(command, label, self.called_copynumbers)

		return output_result