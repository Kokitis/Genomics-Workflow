import os
import vcf
from pprint import pprint
import configparser
PIPELINE_DIRECTORY = "/home/upmc/Documents/Variant_Discovery_Pipeline"
OPTIONS_FILENAME = os.path.join(PIPELINE_DIRECTORY, "0_config_files", "pipeline_project_options.txt")
OPTIONS = configparser.ConfigParser()
OPTIONS.read(OPTIONS_FILENAME)
def compareCallers(left, right):
	""" Uses GATK CombineVariants to merge and compare the output of two callers.
		Parameters
		----------
			left: string
			right: string
	"""
	gatk_program = OPTIONS['Programs']['GATK']
	reference = OPTIONS['Reference Files']['reference genome']
	output_file = "left-right-output.vcf"

	order = "left,right" #ordered by VAF confidence
	
	prioritize_command = """java -jar {gatk} \
		-T CombineVariants \
		-R {reference} \
		--variant:left {left} \
		--variant:right {right} \
		-o {output} \
		-genotypeMergeOptions PRIORITIZE \
		-priority {rod}""".format(
			gatk = gatk_program,
			reference = reference,
			left = left,
			right = right,
			output = output_file,
			rod = order)
	uniqueify_command = """java -jar {gatk} \
		-T CombineVariants \
		-R {reference} \
		--variant {left} \
		--variant {right} \
		-o {output} \
		-genotypeMergeOptions UNIQUIFY""".format(
			gatk = gatk_program,
			reference = reference,
			left = left,
			right = right,
			output = output_file,
			rod = order)
	os.system(uniqueify_command)

def countVariants(filename):
	""" Counts the number of variants detected per Chromosome."""
	chromosomes = dict()
	with open(filename, 'r') as file1:
		reader = vcf.Reader(file1)

		for record in reader:
			if record.CHROM not in chromosomes:
				chromosomes[record.CHROM] = 1
			else:
				chromosomes[record.CHROM] += 1
	return chromosomes

	pprint(chromosomes)
def compareOutput(left, right):
	print("Left File: ", left)
	print("Right File: ", right)
	lresults = countVariants(left)
	rresults = countVariants(right)

	for key, l in sorted(lresults.items()):
		print(key, '\t', l,'\t', rresults.get(key, 'N/A'))
if __name__ == "__main__":
	gdc_folder = "/home/upmc/Documents/TCGA-ESCA/TCGA-2H-A9GF/somatic_variants/GDC"
	gdc_filename = os.path.join(gdc_folder, "TCGA-2H-A9GF-01A-11D-A37C-09_TCGA-2H-A9GF-11A-11D-A37F-09_somaticsniper.vcf")
	#gdc_result = countVariants(gdc_filename)

	folder = "/home/upmc/Documents/Variant_Discovery_Pipeline/3_called_variants/TCGA-2H-A9GF-CHR1/SomaticSniper"
	filename = os.path.join(folder, "TCGA-2H-A9GF-CHR1-11A_vs_TCGA-2H-A9GF-CHR1-01A.somaticsniper.vcf")
	compareOutput(gdc_filename, filename)