import vcf
import shutil
from pprint import pprint
import os

from settings import Settings

class TruthsetPipeline:
	def __init__(self, samples, options, **kwargs):
		truthset_type = kwargs['truthset_type']
		for sample in samples:
			self.runWorkflow(
				sample = sample,
				workflow_options = options,
				truthset_type = truthset_type
		)

	def runWorkflow(self, sample, workflow_options, truthset_type):

		truthset_indel = Truthset(sample, workflow_options, 'indel', truthset_type)
		truthset_snp   = Truthset(sample, workflow_options, 'snp',   truthset_type)

		result = {
			'filename:indel': truthset_indel.filename,
			'filename:snp':   truthset_snp.filename
		}
		return result


class Truthset:
	def __init__(self, sample, truthset_options, callset_type, truthset_type, **kwargs):
		########################### Define Common Attributes ##################
		self.debug = True
		if self.debug:

			print("\ttruthset_type = {}".format(truthset_type))

		# Parameters to use when generating the truthsets.
		self.truthset_type = truthset_type

		self.indel_intersection = kwargs.get('indel_intersection', 2)
		self.snp_intersection   = kwargs.get(  'snp_intersection', 5)

		self.min_tumor_vaf  = 0.08
		self.max_normal_vaf = 0.03

		self.gatk_program   = truthset_options['Programs']['GATK']
		self.picard_program = truthset_options['Programs']['Picard']
		self.reference      = truthset_options['Reference Files']['reference genome']
		self.output_folder = truthset_options.getPipelineFolder('truthset')

		######################## Generate the Truthset ########################
		self.filename = self.runWorkflow(
			sample = sample,
			callset_type = callset_type,
			truthset_type = truthset_type,
			**kwargs
		)
	
	def runWorkflow(self, sample, callset_type, truthset_type, **kwargs):
		"""
			Parameters
			----------
				callset_type: {'indel', 'snp'}
		"""
		patientId = sample['PatientID']
		############################## Define Filenames ##############################
		raw_callset_filename = os.path.join(
			self.output_folder,
			"{}.{}.raw.vcf".format(patientId, callset_type)
		)
		final_truthset_filename = os.path.join(
			self.output_folder,
			"{}.{}.{}.final.truthset.vcf".format(patientId, truthset_type, callset_type)
		)

		############################## Prepare the Truthset ###########################
		raw_callset = getSampleCallset(patientId, 'original-fixed-split', callset_type)

		if len(raw_callset) == 1:
			shutil.copy2(raw_callset.pop(), raw_callset_filename)
		elif len(raw_callset) > 1:
			GATK_MERGE_CALLSET(
				callset = raw_callset,
				filename = raw_callset_filename
			)
		else:
			message = "The callset is empty!"
			raise ValueError(message)

		# Generate two truthsets, one for snps and one for indels.
		final_truthset_filename = self._generateTruthset(
			input_vcf = raw_callset_filename,
			output_vcf = final_truthset_filename,
			training_type = truthset_type,
			callset_type = callset_type,
			**kwargs
		)

		return final_truthset_filename

	def _generateTruthset(self, input_vcf, output_vcf, truthset_type, **kwargs):
		""" Generates a truthset
			Parameters
			----------
				sample
				options
				training_type: {'RNA-seq', 'VAF', 'Intersection'}
		 """
		if self.debug:
			print("Generating Truthset...")
			print("\tInput: ", input_vcf)
			print("\tOutput:", output_vcf)

		with open(input_vcf, 'r') as input_file:
			reader = vcf.Reader(input_file)
			with open(output_vcf, 'w') as output_file:
				writer = vcf.Writer(output_file, reader)

				for index, record in enumerate(reader):
					validation_status = self._getRecordValidationStatus(
						record = record, 
						truthset_type = truthset_type
					)
					is_valid = validation_status['status']
					
					if is_valid:
						writer.write_record(record)

		return output_vcf

	def _getRecordValidationStatus(self, record, truthset_type):
		""" Determines if a given record is within the truthset.
			Returns
			-------
				dict<>
					* 'chrom': 
					* 'position': 
					* 'validation method': 
					* 'validation status': 
		"""
		if truthset_type == 'intersection':
			record_status = self._fromIntersection(record)
		elif truthset_type == 'rna':
			# Assume the record is from the vcf from the RNA-seq pipeline.
			# Assume all RNA-seq variants are true variants.
			record_status = self._fromRna(record)
		elif truthset_type == 'vaf':
			record_status = self._fromVaf(record)
		else:
			message = "The training type is not supported: '{}'".format(truthset_type)
			raise ValueError(message)

		result = {
			'chrom': record.CHROM,
			'position': record.POS,
			'method': truthset_type,
			'status': record_status
		}

		return result

	def _fromIntersection(self, record):

		_separator = '-'

		# GATK includes callers that filtered the variant, designated by a "filterIn{caller}" label.
		callset = [i for i in record.INFO['set'].split(_separator) if 'filter' not in i.lower()]

		# MuSE and Somaticsniper are snp-only, and should not be counted in the intersection of indel callsets
		if record.is_indel:
			num_callers_in_intersection = self.indel_intersection
		else:
			num_callers_in_intersection = self.snp_intersection

		_is_intersection = len(callset) == 1 and callset[0] == 'Intersection'

		_is_in_n_sets = len(callset) >= num_callers_in_intersection
		
		validation_status = _is_intersection or _is_in_n_sets

		row = {
			'chrom': record.CHROM,
			'position': record.POS,
			'validation method': 'Intersection',
			'validation status': validation_status
		}
		return row

	@staticmethod
	def _fromRna(record):
		is_valid = not record.FILTER # True if FILTER is empty or is None
		row = {
			'chrom': record.CHROM,
			'position': record.POS,
			'validation method': 'RNA-seq',
			'validation status': is_valid
		}
		return row

	def _fromVaf(self, record):
		sample_vaf = self._getSampleVAF(record)
		normal_vaf = sample_vaf['NORMAL']['vaf']
		tumor_vaf = sample_vaf['TUMOR']['vaf']

		# Filter out variants that were rejected by a caller
		# This is only needed for variants in the callsets of a single caller
		filter_out = False
		# Validate variants according to VAF status.
		# This is only for testing purposes, a better version should be used later.
		high_tumor_vaf  = tumor_vaf  >= self.min_tumor_vaf
		high_normal_vaf = normal_vaf >= self.max_normal_vaf
		_somatic_vaf = (high_tumor_vaf and not high_normal_vaf)
		_somatic_vaf = _somatic_vaf or (not high_tumor_vaf and normal_vaf == 0.0)
		if _somatic_vaf and not filter_out:
			validation_status = 1  # Somatic
		elif (not high_tumor_vaf and (normal_vaf > 0.0 and not high_normal_vaf)) or filter_out:
			validation_status = 0  # Non-Somatic
		else: validation_status = 2  # 'Unknown'

		validation_status = int(validation_status == 1)

		row = {
			'chrom': record.CHROM,
			'position': record.POS,
			'validation method': 'VAF',
			'validation status': validation_status
		}
		return row

	@staticmethod
	def _getSampleVAF(record):
		""" Use DP4 instead of DP
			Parameters
			----------
				record: from pyVCF reader
			Returns
			-------
				result: dict<>
				* {NORMAL, SAMPLE}:
					* 'alleles': THe number of alternate alleles.
					* 'reads': The number of reads used to calculate the vaf.
					* 'vaf': The variant allele frequency.
		"""
		response = dict()

		for sample in record.samples:
			sample_fields = sample.data._asdict()

			if 'DP4' in sample_fields:
				caller = 'somaticsniper'
			elif 'ALT_F1R2' in sample_fields:
				caller = 'mutect'
			elif 'FREQ' in sample_fields:
				caller = 'varscan'
			elif 'DP' in sample_fields:
				caller = 'muse'
			else:
				caller = 'strelka'
			caller = caller.lower()

			if caller == 'somaticsniper':

				sample_reads = sum(sample_fields['DP4'])
				sample_alleles = sum(sample_fields['DP4'][2:])
				sample_vaf = sample_alleles / sample_reads

			elif caller == 'mutect':
				sample_reads = sample_fields['ALT_F1R2'] + sample_fields['ALT_F2R1']
				sample_reads += sample_fields['REF_F1R2'] + sample_fields['REF_F2R1']
				sample_alleles = sample_fields['ALT_F1R2'] + sample_fields['ALT_F2R1']
				sample_vaf = sample_fields['AF']

			elif caller == 'muse':
				sample_reads = sample_fields['DP']
				sample_alleles = sample_fields['AD'][1]
				sample_vaf = sample_alleles / sample_reads

			elif caller == 'strelka':
				sample_ref = record.REF
				alleles = [i for i in ['A', 'C', 'G', 'T'] if i != sample_ref]
				sample_reads = sum([sample_fields[i + 'U'][1] for i in (alleles + [sample_ref])])
				sample_alleles = sum([sample_fields[i + 'U'][1] for i in alleles])
				if sample_reads == 0:
					sample_vaf = 0
				else:
					sample_vaf = sample_alleles / sample_reads

			elif caller == 'varscan':
				sample_reads = sample_fields['DP']
				sample_alleles = sample_fields['AD']
				sample_vaf = float(sample_fields['FREQ'].strip('%')) / 100
			else:
				message = "WARNING (getSampleVAF): {0} is not a caller!".format(caller)
				raise ValueError(message)

			result = {
				'reads': sample_reads,
				'alleles': sample_alleles,
				'vaf': sample_vaf
			}
			response[sample.sample] = result

		return response

if __name__ == "__main__":
	# Generate a truthset
	base_folder = ""
	default_config_filename = os.path.join(base_folder, "0_config_files", "pipeline_project_options.txt")
	truthset_options = Settings(default_config_filename)
	pipeline_type = 'rna'

	truthset_samples = []

	TruthsetPipeline(
		samples = truthset_samples,
		options = truthset_options,
		truthset_type = pipeline_type
	)