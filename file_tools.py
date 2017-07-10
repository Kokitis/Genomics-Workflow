import os
import csv
import vcf
from pprint import pprint
import configparser
PIPELINE_DIRECTORY = "/home/upmc/Documents/Variant_Discovery_Pipeline"
OPTIONS_FILENAME = os.path.join(PIPELINE_DIRECTORY, "0_config_files", "pipeline_project_options.txt")
OPTIONS = configparser.ConfigParser()
OPTIONS.read(OPTIONS_FILENAME)

def readTSV(filename, returnFieldnames = False):
    with open(filename, 'r') as file1:
        reader = csv.DictReader(file1, delimiter = "\t")
        fieldnames = reader.fieldnames
        reader = list(reader)

    if returnFieldnames:
        return reader, fieldnames
    else:
        return reader

def writeTSV(table, filename, fieldnames = None):
    if fieldnames is None:
        fieldnames = sorted(table[0].keys())

    with open(filename, 'w', newline = "") as tsv_file:
        writer = csv.DictWriter(tsv_file, delimiter = "\t", fieldnames = fieldnames, restval = "")
        writer.writeheader()
        writer.writerows(table)

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

def compareOutput(left, right):
    print("Left File: ", left)
    print("Right File: ", right)
    lresults = countVariants(left)
    rresults = countVariants(right)

    for key, l in sorted(lresults.items()):
        print(key, '\t', l,'\t', rresults.get(key, 'N/A'))

def sortRefSeq(filename):
    def parseChrom(chrom):
        if len(chrom) == 4:
            number = chrom[-1]
            if number in "1234567890":
                number = "0" + number

        elif len(chrom) == 5:
            number = chrom[-2:]
        else:
            number = chrom[3:]

        return number
    print("Reading ", filename)
    table, fieldnames = readTSV(filename, True)

    print("Sorting the table...")

    table = sorted(table, key = lambda s: (parseChrom(s['chrom']), int(s['txStart'])))

    #Remove non-standard contigs, as they conflict with the reference genome.
    reduced_table = [row for row in table if (len(row['chrom']) < 6) and row['chrom'] != 'chrM']

    #the contig "chrM" needs to go after contigs "chrX" and "chrY"
    m_table = [row for row in table if row['chrom'] == 'chrM']
    final_table = reduced_table + m_table


    print("Writing the table...")
    writeTSV(final_table, filename + '.sorted.tsv', fieldnames = fieldnames)


if __name__ == "__main__" and False:
    gdc_folder = "/home/upmc/Documents/TCGA-ESCA/TCGA-2H-A9GF/somatic_variants/GDC"
    gdc_filename = os.path.join(gdc_folder, "TCGA-2H-A9GF-01A-11D-A37C-09_TCGA-2H-A9GF-11A-11D-A37F-09_varscan.vcf")
    #gdc_result = countVariants(gdc_filename)

    folder = "/home/upmc/Documents/Variant_Discovery_Pipeline/3_called_variants/TCGA-2H-A9GF-CHR1/SomaticSniper"
    filename = os.path.join(folder, "TCGA-2H-A9GF-CHR1-11A_vs_TCGA-2H-A9GF-CHR1-01A.somaticsniper.vcf")
    compareOutput(gdc_filename, filename)
elif False:
    filename = "/home/upmc/Documents/Reference/GRCh38-hg38-NCBIRefSeq-UCSCRefSeq-allFieldsFromTable-WholeGene.txt"
    sortRefSeq(filename)