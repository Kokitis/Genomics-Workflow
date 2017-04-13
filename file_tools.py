import os
import vcf

PIPELINE_DIRECTORY = "/home/upmc/Documents/Variant_Discovery_Pipeline"
OPTIONS_FILENAME = os.path.join(PIPELINE_DIRECTORY, "0_config_files", "pipeline_project_options.txt")

 def compareCallers(left, right):
 	""" Uses GATK CombineVariants to merge and compare the output of two callers.
 		Parameters
 		----------
 			left: string
 			right: string
 	"""
 	gatk_program = options['Programs']['GATK']
 	reference = options['Reference Files']['reference genome']
 	output_file = ""

 	order = "mutect,varscan,strelka,muse,somaticsniper" #ordered by VAF confidence
		
	command = """java -jar "{gatk}" \
		-T CombineVariants \
		-R "{reference}" \
		--variant:muse "{muse}" \
		--variant:mutect "{mutect}" \
		--variant:varscan "{varscan}" \
		--variant:somaticsniper "{ss}" \
		--variant:strelka "{strelka}" \
		-o "{output}" \
		-genotypeMergeOptions PRIORITIZE \
		-priority {rod}""".format(
			gatk = gatk_program,
			reference = reference,
			left = left,
			right = right,
			output = output_file,
			rod = order)
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

	pprint(chromosomes)

if __name__ == "__main__":
	filename = ""
	countVariants(filename)