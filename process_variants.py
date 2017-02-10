import os
import collections
import vcf
import csv
from pprint import pprint

def generateVennDiagram(filename = None, script_file = None, plot_file = None):
	"""
		TCGA-2H-A9GF
		area1 = 985,  n=4:, n=5;
		area2 = 618,
		area3 = 1140,
		area4 = 849,
		area5 = 1367,
		n12 = 496,
		n13 = 747,
		n14 = 748,
		n15 = 811,
		n23 = 442,
		n24 = 495,
		n25 = 432,
		n34 = 658,
		n35 = 734,
		n45 = 675,
		n123 = 434,
		n124 = 470,
		n125 = 420,
		n134 = 629,
		n135 = 686,
		n145 = 637,
		n234 = 426,
		n235 = 406,
		n245 = 407,
		n345 = 602,
		n1234 = 420,
		n1235 = 400,
		n1245 = 401,
		n1345 = 580,
		n2345 = 391,
		n12345 = 386,
	"""
	if filename is None:
		filename = "F:\\Combined_Variants\\output\\muse_mutect_somaticsniper_strelka_varscan_Uniquify.vcf"
	if script_file is None:
		script_file = "F:\\Combined_Variants\\output\\rscript.R"
	if plot_file is None:
		plot_file = '"F:/Combined_Variants/output/venn_diagram.tiff"'
	all_caller_names = ['muse', 'mutect', 'somaticsniper', 'strelka', 'varscan']
	combs = list()
	for i in range(1, 6):
		combs += ['_'.join(i) for i in list(itertools.combinations(all_caller_names, i))]
	callers = {k:0 for k in combs}

	with open(filename, 'r') as vcf_file:
		vcf_reader = vcf.Reader(vcf_file)
		for record in vcf_reader:
			category = record.INFO['set']
			category = '_'.join(sorted([i for i in category.split('-') if 'filter' not in i]))
			if category not in callers: callers[category] = 0
			callers[category] += 1

	caller_sets = {'Intersection': callers['Intersection']}
	for key in combs:
		caller_sets[key] = callers['Intersection']
		superset = set(key.split('_'))
		for subset, value in callers.items():
			subset = set(subset.split('_'))
			if superset.issubset(subset):
				caller_sets[key] += value
	r_script = """
		library(VennDiagram)
		venn.plot <- draw.quintuple.venn(
		area1 = {muse},
		area2 = {mutect},
		area3 = {somaticsniper},
		area4 = {strelka},
		area5 = {varscan},
		n12 = {muse_mutect},
		n13 = {muse_somaticsniper},
		n14 = {muse_strelka},
		n15 = {muse_varscan},
		n23 = {mutect_somaticsniper},
		n24 = {mutect_strelka},
		n25 = {mutect_varscan},
		n34 = {somaticsniper_strelka},
		n35 = {somaticsniper_varscan},
		n45 = {strelka_varscan},
		n123 = {muse_mutect_somaticsniper},
		n124 = {muse_mutect_strelka},
		n125 = {muse_mutect_varscan},
		n134 = {muse_somaticsniper_strelka},
		n135 = {muse_somaticsniper_varscan},
		n145 = {muse_strelka_varscan},
		n234 = {mutect_somaticsniper_strelka},
		n235 = {mutect_somaticsniper_varscan},
		n245 = {mutect_strelka_varscan},
		n345 = {somaticsniper_strelka_varscan},
		n1234 = {muse_mutect_somaticsniper_strelka},
		n1235 = {muse_mutect_somaticsniper_varscan},
		n1245 = {muse_mutect_strelka_varscan},
		n1345 = {muse_somaticsniper_strelka_varscan},
		n2345 = {mutect_somaticsniper_strelka_varscan},
		n12345 = {Intersection},
		category = c("MuSE", "MuTect2", "SomaticSniper", "Strelka", "Varscan"),
		fill = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"),
		cat.col = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"),
		cat.cex = 2,
		margin = 0.05,
		cex = c(1.5, 1.5, 1.5, 1.5, 1.5, 1, 0.8, 1, 0.8, 1, 0.8, 1, 0.8, 1, 0.8,
		1, 0.55, 1, 0.55, 1, 0.55, 1, 0.55, 1, 0.55, 1, 1, 1, 1, 1, 1.5),
		ind = TRUE
		);

		# Writing to file
		tiff(filename = {filename}, compression = "lzw");
		grid.draw(venn.plot);
		dev.off();
	""".format(**caller_sets, filename = plot_file)
	if script_file:
		with open(script_file, 'w') as outfile:
			outfile.write(r_script)


class HarmonizeVCFs:
	def __init__(self, variants, output_folder, options, tag = None,):
		"""
			Parameters
			----------
				variants: dict<>
					A dictionary mapping a caller to its output file.
		Output
		------
			output_folder/harmonized_vcfs/
				Each vcf file will be modified to include any field present in at least one of the raw vcf files.
				These new fields will be empty if no equivilent value is available in the raw vcf file.
			output_folder/merged_vcfs/
				This will contain the output of GATK CombineVariants. Callers are saved in the 
			
		"""
		if tag: self.tag = '.' + tag
		else: self.tag = ""
		self.callers = None #_getMergedRecord has hardcoded values
		self._define_caller_formats()
		
		harmonized_folder = os.path.join(output_folder, "harmonized_vcfs")
		merged_folder = os.path.join(output_folder, "merged_vcfs")
		print("Output Folder: ", output_folder)
		print("Harmonizing the vcfs...")
		variants = self._harmonizeVCFs(variants, output_folder = harmonized_folder)
		print("Merging vcfs...")
		merged_vcf = self.merge_vcfs(variants, output_folder = merged_folder) #GATK
		merged_vcf = self.validate_variants(merged_vcf)
		print("Formatting as a table...")
		merged_table = self.toTable(merged_vcf)
	#----------------------------------- Basic ----------------------------------------
	def _define_caller_formats(self):
		muse_format = ["GT", "DP", "AD", "BQ", "SS", 'READS', 'VAF']
		mutect_format = "GT:AD:AF:ALT_F1R2:ALT_F2R1:FOXOG:QSS:REF_F1R2:REF_F2R1:READS:VAF".split(':')
		somaticsniper_format = ["GT", "DP", "FDP", "SDP", "SUBDP", "AU", "CU", "GU", "TU","READS", "VAF"]
		strelka_format = "GT:IGT:DP:DP4:BCOUNT:GQ:JGQ:VAQ:BQ:MQ:AMQ:SS:SSC:READS:VAF".split(":")
		varscan_format = ['GT', 'GQ', 'DP', 'RD', 'AD', 'FREQ', 'DP4', 'READS', 'VAF'] 
		all_formats = muse_format + mutect_format + somaticsniper_format + strelka_format + varscan_format
		all_formats = ['GT'] + sorted([ i for i in set(all_formats) if i != 'GT']) #GT has to be first
		all_formats.remove('AF')
		all_formats.remove('FREQ')
		self.muse_tuple = collections.namedtuple('museData', muse_format)
		self.mutect_tuple = collections.namedtuple('mutectData', mutect_format)
		self.strelka_tuple = collections.namedtuple('strelkaData', strelka_format)
		self.ss_tuple = collections.namedtuple('SSData', somaticsniper_format)
		self.varscan_tuple = collections.namedtuple('varscanData', varscan_format)
		
		self.harmonized_tuple = collections.namedtuple('Format', all_formats)
	
	def _get_merged_data(self, fields):
		""" """
		harm_fields = list(self.harmonized_tuple._fields)
		for f in list(fields.keys()):
			if f not in harm_fields: fields.pop(f)
		_fields = {f:fields.get(f, '.') for f in harm_fields}
		_data = self.harmonized_tuple(**_fields)
		return _data
	#-------------------------------- Harmonize ---------------------------------------
	@staticmethod
	def getSampleVAF(sample, caller, sample_ref = None):
		""" Use DP4 instead of DP
			Parameters
			----------
				sample:
				caller: {'muse', 'mutect', 'somaticsniper', 'strelka', 'varscan'}
				sample_ref: {'A', 'C', 'G', 'T'}
					Only required for strelka.
			Returns
			-------
				result: dict<>
					* 'alleles': THe number of alternate alleles.
					* 'reads': The number of reads used to calculate the vaf.
					* 'vaf': The variant allele frequency.
		"""
		filter_out = False
		if caller == 'somaticsniper':
			sample_reads = sum(sample['DP4'])
			sample_alleles = sum(sample['DP4'][2:])
			sample_vaf = sample_alleles / sample_reads

		elif caller == 'mutect':
			#print(sample)
			sample_reads = sample['ALT_F1R2'] + sample['ALT_F2R1'] + sample['REF_F1R2'] + sample['REF_F2R1']
			sample_alleles = sample['ALT_F1R2'] + sample['ALT_F2R1']
			sample_vaf = sample['AF']

		elif caller == 'muse':
			sample_reads = sample['DP']
			sample_alleles = sample['AD'][1]
			sample_vaf = sample_alleles / sample_reads

		elif caller == 'strelka':
			alleles = [i for i in ['A', 'C', 'G', 'T'] if i != sample_ref]
			sample_reads = sum([sample[i+'U'][1] for i in (alleles + [sample_ref])])
			sample_alleles = sum([sample[i+'U'][1] for i in alleles])
			if sample_reads == 0: sample_vaf = 0
			else:
				sample_vaf = sample_alleles / sample_reads

		elif caller == 'varscan':
			sample_reads = sample['DP']
			sample_alleles = sample['AD']
			sample_vaf = float(sample['FREQ'].strip('%')) / 100

		result = {
			'ALLELES': sample_alleles,
			'READS': sample_reads,
			'VAF': float("{0:.5f}".format(sample_vaf)),
		}
		return result

	def _harmonizeRecord(self, record, caller):
		new_samples = list()
		for s in record.samples:
			sample_data = s.data._asdict()
			sample_vaf = self.getSampleVAF(s, caller, sample_ref = record.REF)
			sample_data = dict(list(sample_data.items()) + list(sample_vaf.items()))
			s.data = self._get_merged_data(sample_data)
			new_samples.append(s)


		new_record = vcf.model._Record(
			CHROM = record.CHROM,
			POS = record.POS,
			ID = record.ID,
			REF = record.REF,
			ALT = record.ALT,
			QUAL = record.QUAL,
			FILTER = record.FILTER,
			INFO = record.INFO,
			FORMAT = ':'.join(new_samples[0].data._fields),
			sample_indexes = record._sample_indexes,
			samples = new_samples)
		return new_record
	
	@staticmethod
	def _getHeaders(formats, infos, caller):
		new_formats = formats
		new_infos = infos
		_format_tuple = collections.namedtuple('Format', ['id', 'num', 'type', 'desc'])
		#_Format = collections.namedtuple('Format', ['id', 'num', 'type', 'desc'])
		##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
		##FORMAT=<ID=AF,Number=1,Type=Float,Description="Allele fraction of the event in the tumor">
		##FORMAT=<ID=FREQ,Number=1,Type=String,Description="Variant allele frequency">
		if caller == 'strelka':
			new_formats['GT'] = _format_tuple('GT', '1', 'String', 'Genotype')
		elif caller == 'varscan':
			new_formats['DP4'] = _format_tuple('DP4', '4', 'Integer', 'ref/fwd, ref/rev, var/fwd, var/rev')

		new_formats['READS'] = _format_tuple('READS', '1', 'Integer', 'Total Quality Reads')
		new_formats['VAF'] = _format_tuple('VAF', '1', 'Float', 'Variant Allele Frequency')
		new_formats['ALLELES'] = _format_tuple('ALLELES', '1', 'String', 'The number of alternate alleles')
		#new_formats['VS'] = _format_tuple('VS', '1', 'Integer', 'The Validation Status of the mutation.')

		new_infos['VS'] = _format_tuple('VS', '1', 'Integer', 'Validation Status')

		headers = {
			'formats': new_formats,
			'infos':   new_infos
		}
		return headers
	
	def _harmonizeVCFs(self, sample, options, variants, output_folder):
		""" Harmonizes all vcf for a given patient.
			Returns
			-------
				variants: dict<caller: path>
					A dictionary pointing to the harmonized vcf for each caller.
		"""
		for caller, filename in variants.items():
			hfile = self._harmonizeVCF(sample, options, filename, output_folder, caller)
			variants[caller] = hfile
		return variants
	
	def _harmonizeVCF(self, filename, output_folder, caller = None):
		""" Harmonizes the output VCFs from each caller.
			Parameters
			----------
				filename: string [PATH]
					The location of the output vcf. The name of the caller should be in the file's name.
					Format: {NormalID}_vs_{TumorID}.{CallerName}.vcf
				output_folder: string [PATH]
					The folder to save the output to.
				caller: {}; defualt None
			Returns
			-------
				output_filename: string
					Path the the harmonized VCF.
					Format: {NormalID}_vs_{TumorID}.{CallerName}.{TAG}.harmonized.vcf
			Class Attributes
			----------------
				self.tag: string
					A tag to append to include when naming the file.
			
		"""        
		if caller: pass
		elif 'muse' in filename: caller = 'muse'
		elif 'mutect' in filename: caller = 'mutect'
		elif 'somaticsniper' in filename: caller = 'somaticsniper'
		elif 'strelka' in filename: caller = 'strelka'
		else: caller = 'varscan'

		source_folder, basename = os.path.split(filename)
		output_filename = "{normal}_vs_{tumor}.{caller}.harmonized.vcf".format(sample['NormalID'], sample['SampleID'], caller)
		output_filename = os.path.join(output_folder, output_filename)

		with open(filename, 'r') as input_file:
			vcf_reader = vcf.Reader(input_file)
			headers = self._getHeaders(vcf_reader.formats, vcf_reader.infos, caller)
			vcf_reader.infos   = headers['infos']
			vcf_reader.formats = headers['formats']
			#vcf_reader.formats = self._getFormats(vcf_reader.formats, caller)
			with open(output_filename, 'w') as output_file:
				vcf_writer = vcf.Writer(output_file, vcf_reader)
				for index, record in enumerate(vcf_reader):
					new_record = self._harmonizeRecord(record, caller)
					vcf_writer.write_record(new_record)
		return output_filename
	#---------------------------------- Merge -----------------------------------------
	#Validate
	def _filterOut(self, record, caller = None):
		filter_out = not record.is_snp
		if record.FILTER is None:
			return filter_out
		filters = [i for i in record.FILTER if i not in ['PASS', 'Tier1', 'Tier2', 'Tier3', 'Tier4']]
		#Caller-specific filters
		if caller == 'muse':
			if len(record.FILTER) != 0 and record.FILTER[0] not in filters:
				filter_out = True
		elif caller == 'mutect':
			if len(filters) != 0: filter_out = True
		else:
			#General caller filter
			if len(filters) != 0: filter_out = True
		if record.INFO['set'] == 'FilteredInAll': filter_out = True

		return filter_out

	def _getValidationStatus(self, record, caller = None):
		"""
			Should use better Criteria.
			
			Use the intersection? (n = 4 or 5)
			Compare against dbSNP, COSMIC?
			Lower score if filtered by one caller but passed y others?
		"""

		normal = [i for i in record.samples if i.sample == 'NORMAL'][0]
		tumor  = [i for i in record.samples if i.sample == 'TUMOR'][0]
		tumor_vaf = tumor.data.VAF
		normal_vaf= normal.data.VAF

		#Filter out variants that were rejected by a caller
		#This is only needed for variants in the callsets of a single caller
		filter_out = self._filterOut(record, caller)
		
		#Validate variants according to VAF status.
		#This is only for testing purposes, a better version should be used later.
		_somatic_vaf = (tumor_vaf >= 0.08 and normal_vaf < 0.03) or (tumor_vaf < 0.08 and normal_vaf == 0.0)
		if _somatic_vaf and not filter_out:
			validation = 'Somatic'
		elif (tumor_vaf < 0.08 and (normal_vaf > 0.0 and normal_vaf < 0.03)) or filter_out:
			validation = 'Unknown'
		else: validation = 'Non-Somatic'
			
		if validation == 'Non-Somatic': validation = 0
		elif validation == 'Somatic': validation = 1
		else: validation = 2
		
		#Apply any additional filtering techniques.
		#if validation != 0:
		#    validation = int(len(record.INFO['set'].split('_')) == 5)

		#print(validation, normal_vaf, tumor_vaf, filter_out)

		return validation
	
	#Merge
	def merge_vcfs(self,sample, options, variants, output_folder):
		""" Uses GATK CombineVariants to merge the calls from each caller into a single file.
			Parameters
			----------
				variants: dict<caller, path>
					A dictionary linkng each caller to its harmonized output.
					Format: {NormalID}_vs_{TumorID}.{CallerName}.{TAG}.harmonized.vcf
			Returns
			-------
				Output_filename: string
					Format: {NormalID}_vs_{TumorID}.{CallerName}.{TAG}.merged.vcf
		"""
		if os.name == 'nt':
			gatk = "C:\\Users\\Deitrickc\\Downloads\\Genomic Programs\\GenomeAnalysisTK-3.7\\GenomeAnalysisTK.jar"
			reference = "G:\\Pipeline Files\\Reference\\GRCh38.d1.vd1.fa"
		else:
			gatk = options['Programs']['GATK']
			reference = options['Reference Files']['reference genome']
		
		#Check GATK
		if not os.path.isfile(gatk):
			print("ERROR: Could not locate GATK at ", gatk)
		if not os.path.isfile(reference):
			print("ERROR: Could not locate reference at ", reference)

		output_file = "{normal}_vs_{tumor}.merged.vcf".format(sample['NormalID'], sample['SampleID'])
		output_file = os.path.join(output_folder, output_file)
		
		order = "mutect,varscan,strelka,muse,somaticsniper" #ordered by VAF confidence
		order = ",".join([i for i in order.split(',') if i in variants])

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
			-priority {rod}"""
		command = command.format(
				gatk = gatk,
				reference = reference,
				muse = variants['muse'],
				mutect = variants['mutect'],
				ss = variants['somaticsniper'],
				strelka = variants['strelka'],
				varscan = variants['varscan'],
				variants = templates,
				rod = order,
				output = output_file)

		if not os.path.isfile(output_file):
			print(command)
			os.system(command)
		return output_file
		
	def _getMergedRecord(self, record):
		#Modify the INFO field to add the validation status and format the 'sets' field
		#infotuple = collections.namedtuple('INFO', sorted(record.INFO) + ['VS'])
		
		
		new_info = record.INFO
		if new_info['set'] == 'Intersection':
			caller_sets = sorted(['muse', 'mutect', 'somaticsniper', 'strelka', 'varscan'])
		else:
			caller_sets = sorted(i for i in new_info['set'].split('-') if 'filter' not in i)
		caller_sets = '_'.join(caller_sets)
		new_info['set'] = caller_sets
		new_info['VS'] = self._getValidationStatus(record)
		#new_info = infotuple(**new_info)

		new_record = vcf.model._Record(
			CHROM = record.CHROM,
			POS = record.POS,
			ID = record.ID,
			REF = record.REF,
			ALT = record.ALT,
			QUAL = record.QUAL,
			FILTER = record.FILTER,
			INFO = new_info,#record.INFO,
			FORMAT = record.FORMAT,
			sample_indexes = record._sample_indexes,
			samples = record.samples)
		
		return new_record
	 
	def validate_variants(self, merged_vcf):
		source_folder, basename = os.path.split(merged_vcf)
		basename = basename.split('.')[0] + self.tag + '.validated.vcf'
		output_filename = os.path.join(source_folder, basename)
		with open(merged_vcf, 'r') as file1:
			vcf_reader = vcf.Reader(file1)
			with open(output_filename, "w") as output:
				vcf_writer = vcf.Writer(output, vcf_reader)
				
				for record in vcf_reader:
					#validation_status = self._getValidationStatus(record)
					new_record = self._getMergedRecord(record)
					vcf_writer.write_record(new_record)
		return output_filename
			
	def toTable(self, merged_vcf):
		table = list()
		stats = collections.defaultdict(list)
		with open(merged_vcf, 'r') as file1:
			vcf_reader = vcf.Reader(file1)
			for record in vcf_reader:
				normal = [s for s in record.samples if s.sample == 'NORMAL'][0].data._asdict()
				tumor  = [s for s in record.samples if s.sample == 'TUMOR'][0].data._asdict()
				variant_allele = record.ALT[0]
				if variant_allele is None or len(variant_allele) > 1: continue
				row = {
					"Patient": "TCGA-2H-A9GF",
					"Validation_Status": record.INFO['VS'],
					"Mutation_Set": record.INFO['set'],
					"Reference_Allele": record.REF,
					"Variant_Allele": variant_allele,
					"Depth_Tumor": tumor['READS'],
					"Depth_Normal": normal['READS'],
					"Vaf_Tumor": tumor['VAF']*100,
					"Vaf_Normal": normal['VAF']*100
				}
				
				if row['Mutation_Set'] == 'FilteredInAll': continue
				elif row['Validation_Status'] == 2: continue
				
				stats[row['Mutation_Set']].append(row['Validation_Status'])
				table.append(row)
		print("Validation per set:")

		for key, series in sorted(stats.items()):
			
			_total_validated = series.count(1)
			_total_detected = len(series)
			_ratio = _total_validated / _total_detected
			print("{0:<45}\t{1}\t{2:.1%}".format(key, _total_detected, _ratio))
				
		table_filename = merged_vcf + '.table'
		header = sorted(table[0].keys())
		with open(table_filename, 'w') as outputfile:
			writer = csv.DictWriter(outputfile, delimiter = '\t', fieldnames = header)
			writer.writeheader()
			writer.writerows(table)
		return table_filename
				
	def expandTable(table_file, callers, caller_sets):
		""" Adds and initializes columns to a dataframe based on the passed variable list.
		Parameters
		----------
			df: pandas.DataFrame
				The dataframe to modify
			ggname: string
				The name of an existing column in the df.
			filters: list<str>
		
		"""
		#The subtype is just the reference and alternate alleles concatenated
		for line in table:
			row = line
			row['subtype'] = row['Reference_Allele'] + row['Variant_Allele']
			for caller in callers:
				row['caller'] = int(caller in row['Mutation_Set'])
			for caller_set in caller_sets:
				row[caller_set] = int(row['Mutation_Set'] == caller_set)

		return df

def GetVariantList(sample, options):

	vcf_folder = os.path.join(options['Pipeline Options'], sample['PatientID'])

	patientID = sample['PatientID']
	normalID  = sample['NormalID']
	tumorID   = sample['SampleID']

	variants   = {
		#os.path.join(options['output']['Bambino'].format(patient = patientID),"{0}_vs_{1}.bambino.vcf".format(normalID, tumorID)),
		#os.path.join(options['output']['Haplotypecaller'].format(patient = patientID), "{0}_vs_{1}.raw.snps.indels.vcf".join(normalID, tumorID)),
		'muse':          os.path.join(vcf_folder, "{0}_vs_{1}.MuSE.vcf".format(normalID, tumorID)),
		'mutect2':       os.path.join(vcf_folder, "{0}_vs_{1}.mutect2.vcf".format(normalID, tumorID)),
		'somaticsniper': os.path.join(vcf_folder, "{0}_vs_{1}.somaticsniper.hq.vcf".format(normalID, tumorID)),
		'strelka-indel': os.path.join(vcf_folder,"results", "{0}_vs_{1}.passed.somatic.indels.strelka.vcf".format(normalID, tumorID)),
		'strelka-snv':   os.path.join(vcf_folder,"results", "{0}_vs_{1}.passed.somatic.snvs.strelka.vcf".format(normalID, tumorID)),
		'varscan-indel':       os.path.join(vcf_folder, "{0}_vs_{1}.snp.Somatic.hc".format(normalID, tumorID)),
		'varscan-snv': os.path.join()
	}

	#variants = {k:v for k, v in variants.items() if os.path.isfile(v)}

	return variants

def VariantEffectPredictor(sample, options, variants, output_folder):
	""" VariantEffectPredictor uses terminology as defined by the Sequence Ontology Project
	"""
	#perl variant_effect_predictor.pl --cache -i input.txt -o output.txt	

	reference = options['Reference Files']['reference genome']
	program = options['Programs']['varianteffectpredictor']

	pprint(variants)
	for index, source in enumerate(variants):
		_, fn = os.path.split(source)
		destination = os.path.join(output_folder, fn + '.VEP_annotated.vcf')

		command = """perl {vep} \
			--input_file {inputfile} \
			--output_file {output} \
			--fasta {reference} \
			--species homo_sapiens \
			--assembly GRCh38 \
			--format vcf \
			--cache \
			--html \
			--symbol \
			--biotype \
			--total_length \
			--numbers \
			--fields Consequence,Codons,Aminoacids,Gene,SYMBOL,Feature,EXON,PolyPhen,SIFT,Protein_position,BIOTYPE \
			--vcf""".format(
				vep = program,
				inputfile = source,
				output = destination,
				reference = reference)
		#print(command)
		Terminal(command, show_output = True)

	vep_log = {

	}
	return vep_log

def vcftomaf(vcf_file = None):
	vcftomaf_script = "/home/upmc/Programs/vcf2maf-1.6.12/vcf2maf.pl"
	output_folder = "/home/upmc/Documents/Variant_Discovery_Pipeline/8_combined_variants"
	output_file = os.path.join(output_folder, "maf_test.maf")

	if vcf_file is None:
		vcf_file = os.path.join(output_folder, "TCGA-2H-A9GF-01A-11D-A37C-09_TCGA-2H-A9GF-11A-11D-A37F-09_mutect_annotated.vcf")

	#perl vcf2maf.pl --input-vcf tests/test.vcf --output-maf tests/test.vep.maf --tumor-id WD1309 --normal-id NB1308

	command = """perl {script} \
				--input-vcf {vcf_file} \
				--output-maf {output} \
				--tumor-id {tumor} \
				--normal-id {normal}""".format(
					script = vcftomaf_script,
					vcf_file = vcf_file,
					output = output_file)

def VariantAnnotator(sample, options):
	"""
		 java -jar GenomeAnalysisTK.jar \
		-R reference.fasta \
		-T VariantAnnotator \
		-I input.bam \
		-o output.vcf \
		-A Coverage \
		-V input.vcf \
		-L input.vcf \
   --dbsnp dbsnp.vcf
	"""

	va_command =  """java {memory} -jar {GATK} \
		-R {reference} \
		-T VariantAnnotator \
		-I {bam} \
		-o output.vcf \
		-A Coverage \
		-V input.vcf \
		-L input.vcf \
		--dbsnp dbsnp.vcf""".format(
		memory = options['Parameters']['JAVA_MAX_MEMORY_USAGE'],
		GATK = gatk_program,
		reference = reference,
		bam = "")

def VariantRecalibrator(sample, options, variant_type):
	"""	Recalibrates variant quality scores.
		Example Usage
		-------------
		java -jar GenomeAnalysisTK.jar \ 
		-T VariantRecalibrator \ 
		-R reference.fa \ 
		-input raw_variants.vcf \ 
		-resource:hapmap,known=false,training=true,truth=true,prior=15.0 hapmap.vcf \ 
		-resource:omni,known=false,training=true,truth=true,prior=12.0 omni.vcf \ 
		-resource:1000G,known=false,training=true,truth=false,prior=10.0 1000G.vcf \ 
		-resource:dbsnp,known=true,training=false,truth=false,prior=2.0 dbsnp.vcf \ 
		-an DP \ 
		-an QD \ 
		-an FS \ 
		-an SOR \ 
		-an MQ \
		-an MQRankSum \ 
		-an ReadPosRankSum \ 
		-an InbreedingCoeff \
		-mode SNP \ 
		-tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 \ 
		-recalFile recalibrate_SNP.recal \ 
		-tranchesFile recalibrate_SNP.tranches \ 
		-rscriptFile recalibrate_SNP_plots.R 

		Reference
		---------
			https://software.broadinstitute.org/gatk/guide/article?id=2805
	"""
	#-------------------------- Generate the sample-specific options ---------------------------
	prefix = "{0}_vs_{1}".format(sample['NormalID'], sample['SampleID'])
	reference = options['Reference Files']['reference genome']

	#------------------------ Build the recalibration model ----------------------

	snp_build_command = """java {memory} -jar {GATK} \
		-T VariantRecalibrator \
		-R {reference} \
		-input {raw_variants} \
		-an QD \
		-an FS \
		-an SOR \
		-an MQ \
		-mode SNP \
		-resource:hapmap,known=false,training=true,truth=true,prior=15.0 {hapmap} \
		-resource:omni,known=false,training=true,truth=true,prior=12.0 {OMNI} \
		-resource:1000G,known=false,training=true,truth=false,prior=10.0 {OTG} \
		-resource:dbsnp,known=true,training=false,truth=false,prior=2.0 {dbSNP} \
		-tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 \
		-recalFile {prefix}.recalibrate_SNP.recal \
		-tranchesFile {prefix}.recalibrate_SNP.tranches \
		-rscriptFile {prefix}.recalibrate_SNP_plots.R""".format(
			memory       = options['Parameters']['JAVA_MAX_MEMORY_USAGE'],
			GATK         = options['Programs']['GATK'],
			reference    = reference,
			prefix       = prefix,
			raw_variants = "",
			mode         = variant_type.upper(),

			hapmap       = options['Reference Files']['hapmap'],
			OMNI         = options['Reference Files']['omni'],
			OTG          = options['Reference Files']['1000G'],
			dbSNP        = options['Reference Files']['dbSNP'])

	#------------------------- Recalibrate the SNP quality scores ------------------------
	#ts_filter_level = 99.9 recommended by GATK best practices
	#Retrieves 99.9% of true sites, includes an amount of false positives
	snp_apply_command = """java {memory} -jar {gatk} \
		-T ApplyRecalibration \
		-R {reference} \
		-input {raw_variants} \
		-mode SNP \
		--ts_filter_level 99.9 \
		-recalFile {prefix}.recalibrate_SNP.recal \
		-tranchesFile {prefix}.recalibrate_SNP_plots.R \
		-o {prefix}.recalibrated_snps_raw_indels.vcf""".format(
			memory = options['Parameters']['JAVA_MAX_MEMORY_USAGE'],
			gatk = gatk_program,
			reference = reference,
			raw_variants = "",
			prefix = prefix)
	#----------------------------- Recalibrate Indels ----------------------------------
	indel_build_command = """java {memory} -jar {gatk} \
		-T VariantRecalibrator \
		-R {reference} \
		-input {prefix}.recalibrated_snps_raw_indels.vcf \
		-an QD \
		-an FS \
		-an SOR \
		-an MQ \
		-mode INDEL \
		--maxGaussians 4 \
		-resource:dbsnp,known=true,training=false,truth=false,prior=2.0 {dbSNP} \
		-resource:mills,known=false,training=true,truth=true,prior=12.0 {mills} \
		-tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 \
		-recalFile {prefix}.recalibrate_INDEL.recal \
		-tranchesFile {prefix}.recalibrate_INDEL.tranches \
		-rscriptFile {prefix}.recalibrate_INDEL_plots.R""".format(
			memory       = options['Parameters']['JAVA_MAX_MEMORY_USAGE'],
			GATK         = options['Programs']['GATK'],
			reference    = reference,
			prefix       = prefix,
			raw_variants = "",

			mills        = "",
			dbSNP        = options['Reference Files']['dbSNP'])
	
	indel_apply_command = """java {memory} -jar {gatk} \
		-T ApplyRecalibration \
		-R {reference} \
		-input {raw_variants} \
		-mode INDEL \
		--ts_filter_level 99.9 \
		-recalFile {prefix}.recalibrate_INDEL.recal \
		-tranchesFile {prefix}.recalibrate_INDEL_plots.R \
		-o {prefix}.recalibrated_variants.vcf""".format(
			memory = options['Parameters']['JAVA_MAX_MEMORY_USAGE'],
			gatk = gatk_program,
			reference = reference,
			raw_variants = "",
			prefix = prefix)

def ProcessPipelineMethod(sample, options):
	"""
		1. Copy raw VCF Files to the processing directory
			/Processing Directory
				/PatientID
					/raw_vcfs
					/harmonized_vcfs
					/processed_vcfs
		2. Harmonize VCFs
		3. Merge VCFs
		4. Generate Truthset
	"""
	#-------------------------- make all of the relevant folders ----------------------------
	pipeline_folder = os.path.join(os.getcwd(), '8_combined_variants') #Holds folders for each patient
	patient_folder = os.path.join(pipeline_folder, sample['PatientID']) #holds folders for an individual patient
	raw_variant_folder = os.path.join(patient_folder, 'raw_variants')
	annotated_variants_folder = os.path.join(patient_folder, 'annotated_variants')

	for _dir in [pipeline_folder, patient_folder, raw_variant_folder, annotated_variants_folder]:
		if not os.path.isdir(_dir):
			os.mkdir(_dir)
	
	#-------------------------- Get a list of all the variants used in this analysis ---------------------
	sample_variants = GetVariantList(sample, options)
	
	
	#------------------------------ Annotate Variants ----------------------------------------
	raw_variants = list()

	for key, source in sample_variants.items():
		folder, file_name = os.path.split(source)
		destination = os.path.join(raw_variant_folder, file_name)
		shutil.copy2(source, destination)
		raw_variants.append(destination)


	vep_log = VariantEffectPredictor(sample, options, raw_variants, annotated_variants_folder)

def generate_truthset(variants, merged_variants, training_type):
	""" Generates a truthset based on the passed variants/
		Option 1: intersection of all 5 callers.
		Option 2: Dream-SEQ data
		Option 3: RNA-seq?
	"""

	if training_type == 'DREAM': pass
	elif training_type == 'intersection':
		#Use the intersection of the callers
		#Open the merged/combined vcf file and mark all variants in the intersection as 'true'
		pass

def ProcessVariants(sample, options):
	"""
		1. Copy raw VCF Files to the processing directory
			/Processing Directory
				/PatientID
					/raw_vcfs
					/harmonized_vcfs
					/processed_vcfs
		2. Harmonize VCFs
		3. Merge VCFs
		4. Generate Truthset
	"""
	output_folder = ""
	raw_variants = GetVariantList(sample, options)

	#Merge the vcf files with gatk combineVariants.
	#Output is in output_folder/merged_vcfs/
	HarmonizeVCFs(raw_variants, output_folder)

if __name__ == "__main__":
	sample = {}
	options = configparser.ConfigParser()
	options.read(filename)
	ProcessVariants()





