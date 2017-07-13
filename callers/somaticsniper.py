import os
from .basicworkflow import Workflow


class SomaticSniper(Workflow):

	def calculateConfidence(self, vcf_file):
		command = """perl {script} \
			--snp-file {snpfile} \
			--min-mapping-quality {mmq} \
			--min-somatic-score {ss} \
			--lq-output {lq} \
			--out-file {hq}""".format(
				script 	= self.hc_script,
				mmq 	= self.min_mapping_quality,
				ss 		= self.min_somatic_quality,
				snpfile	= vcf_file,
				hq 		= self.hq_variants,
				lq 		= self.lq_variants
		)

		label = "SomaticSniper: Filter lq Variants"
		output_result = self.runCallerCommand(command, label, [self.hq_variants, self.lq_variants])
		return output_result



	def generatePileup(self, bam_file, bam_name):

		raw_pileup_file = os.path.join(self.temp_folder, "{0}.somaticsniper.pileup".format(bam_name))
		filtered_output_file = raw_pileup_file + '.pileup'

		#samtools_command = "{samtools} pileup -cvi -f {reference} {bam} > {pileup}".format(
		samtools_command = "{samtools} pileup -cvi -f {reference} {bam}".format(
			samtools = self.samtools_program,
			reference = self.reference,
			bam = bam_file,
			pileup = raw_pileup_file)

		#filter_command = """samtools.pl varFilter -Q {basequality} {inputpileup} > {outputpileup}""".format(
		filter_command = """samtools.pl varFilter -Q {basequality} {inputpileup}""".format(
			basequality = self.min_base_quality,
			inputpileup = raw_pileup_file,
			outputpileup = filtered_output_file)
		# samtools.pl varFilter raw.pileup | awk '$6>=20' > final.pileup
		label = "SomaticSniper: Generate Pileup File"
		self.runCallerCommand(samtools_command, label, raw_pileup_file, output_filename = raw_pileup_file)

		label = "SomaticSniper: Filter Pileup File"
		self.runCallerCommand(filter_command, label, filtered_output_file, output_filename = filtered_output_file)

		return filtered_output_file

	def setCustomEnvironment(self, sample, options):
		# Define the relevant filenames
		self.caller_name = "SomaticSniper"
		somaticsniper_folder= self.program # for readability
		script_folder 		= os.path.join(somaticsniper_folder, 'src', 	'scripts')
		self.program 		= os.path.join(somaticsniper_folder, 'build', 	'bin', 	'bam-somaticsniper')
		self.readcount_program 	= options['Programs']['bam-readcount']
		self.samtools_program 	= options['Programs']['samtools-0.1.6']


		self.snpfilter_script= os.path.join(script_folder, 'snpfilter.pl')
		self.readcount_script= os.path.join(script_folder, 'prepare_for_readcount.pl')
		self.hc_script 		=  os.path.join(script_folder, 'highconfidence.pl')
		self.fpfilter  		=  os.path.join(script_folder, 'fpfilter.pl')

		self.raw_variants 	= self.abs_prefix + '.vcf'
		self.hq_variants 	= self.abs_prefix + '.hq.vcf'
		self.lq_variants 	= self.abs_prefix + '.lq.vcf'

	def readcounts(self, loh_file, tumor_bam):
		""" Expected Output:
		"""
		prepare_readcount_output 	= loh_file + '.pos'
		readcount_output 			= loh_file + '.readcounts.rc'
		# Prepare readcounts
		pr_command = """perl {script} \
			--snp-file {inputfile} \
			--out-file {output}""".format(
				script 		= self.readcount_script,
				inputfile 	= loh_file,
				output 		= prepare_readcount_output
			)
		label = "SomaticSniper: Prepare Readcounts"
		preparation_result = self.runCallerCommand(pr_command, label, prepare_readcount_output)
		# Readcounts

		readcount_command = """{program} \
			-b {mbq} \
			-q 1 \
			-f {reference} \
			-l {proutput} \
			{tumor} > {output}""".format(
				program 	= self.readcount_program, # different from readcount script
				mbq 		= self.min_base_quality,
				reference 	= self.reference,
				proutput 	= prepare_readcount_output,
				tumor 		= tumor_bam,
				output 		= readcount_output
			)
		label = "SomaticSniper: Generate Readcounts"
		# Note: The output for this command will be saved in the console log file.
		readcount_result = self.runCallerCommand(readcount_command, label, expected_output = readcount_output, output_filename = readcount_output)
		#self._captureReadcountOutput(readcount_output)

		return readcount_result

	def removeFalsePositives(self, loh_file, readcounts):
		false_positive_output = loh_file + '.fp_pass'
		fp_command = """perl {script} \
			--snp-file {snpfile} \
			--readcount-file {readcounts}""".format(
			script 		= self.fpfilter,
			snpfile 	= loh_file,
			readcounts 	= readcounts
		)

		label = "SomaticSniper: Remove False Positives"
		output_status = self.runCallerCommand(fp_command, label, false_positive_output)
		return output_status

	def removeLOH(self, vcf_file, pileup_file, output_file):
		# -------------------------- Filter and remove LOH --------------------------------
		command = """perl {snpfilter} \
			--snp-file {vcf} \
			--indel-file {pileup} \
			--out-file {output}""".format(
				snpfilter = self.snpfilter_script,
				vcf 	= vcf_file,
				pileup 	= pileup_file,
				output 	= output_file
			)
		label = "SomaticSniper: removeLOH"
		self.runCallerCommand(command, label, output_file)
		return output_file

	def runCallerWorkflow(self, sample):
		# Will print "Couldn't find single-end mapping quality. Check to see if the SM tag is in BAM."
		# This doesn't invalidate results, but try not to use single-end mapping quality in output
		if False:
			self.variant_discovery_status = self.runVariantDiscovery(sample)
		else:
			self._overwriteExistingFiles()
			self.variant_discovery_status = self._downloadGdcFile(sample, self.raw_variants)
		# Generate pileup files
		normal_pileup_file = self.generatePileup(sample['NormalBAM'], 	'normal') #returns pileup filename
		tumor_pileup_file  = self.generatePileup(sample['TumorBAM'], 	'tumor')

		# Filter LOH
		_intermediate_loh_filtered_output = self.abs_prefix + ".SNPfilter.intermediate"
		loh_filtered_output = self.abs_prefix + ".SNPfilter.final"
		_intermediate_file  = self.removeLOH(self.raw_variants, normal_pileup_file, _intermediate_loh_filtered_output)
		loh_filtered_output = self.removeLOH(_intermediate_loh_filtered_output, tumor_pileup_file, loh_filtered_output)
		# loh_filtered_output = _intermediate_file

		readcount_status 			= self.readcounts(loh_filtered_output, sample['TumorBAM'])
		readcounts 					= readcount_status['outputFiles']
		false_positive_status 		= self.removeFalsePositives(loh_filtered_output, readcounts)
		false_positive_output 		= false_positive_status['outputFiles']
		high_confidence_status	 	= self.calculateConfidence(false_positive_output)

		self.full_output = [self.raw_variants, self.hq_variants, self.lq_variants]
		self.final_output = self.hq_variants

	def runVariantDiscovery(self, sample):
		"""	Detects raw variants with somaticniper.
			Parameters
			----------
				sample: dict<>
					dictionary mapping sample-specific files.
		"""
		# -------------------------------- Variant Discovery Command -----------------------------
		somaticsniper_command = """{program} \
			-q 1 \
			-Q 15 \
			-s 0.01 \
			-T 0.85 \
			-N 2 \
			-r 0.001 \
			-G \
			-L \
			-n NORMAL \
			-t TUMOR \
			-F vcf \
			-f {reference} {tumor} {normal} {outputfile}""".format(
				program     = self.program,
				reference   = self.reference,
				tumor       = sample['TumorBAM'],
				normal      = sample['NormalBAM'],
				outputfile  = self.raw_variants
			)
		label = "SomaticSniper Variant Discovery"
		output_result = self.runCallerCommand(somaticsniper_command, label, self.raw_variants)
		return output_result